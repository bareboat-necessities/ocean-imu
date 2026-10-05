#!/usr/bin/env python3
"""Conventional heave baselines for the OU-III comparison.

Two comparators, both driven by the same measurement-only vertical acceleration
(Mahony proxy, paper Eq. 73) and optionally the same measurement-only period
estimate (paper Eq. 80) that OU-III uses. Reference motion is used only to
score and to select frozen parameters on a development split.

  FDDI  frequency-domain double integration (offline, zero-phase, non-causal):
        Y(f) = -W(f) A(f) / (2 pi f)^2, raised-cosine high-pass taper W.
        A record-level reference for linear processing, not a real-time method.
  HPDI  causal high-pass double integrator, src/hpdi/HeaveHPDI.h, replayed via
        the hpdi_replay binary. H(s) = [1 - R(s)/D(s)] / s^2 with Butterworth
        D of order n and m DC zeros; m == n is the classic cascaded high-pass.

Subcommands
  synth      write synthetic vertical-only records + manifest (smoke tests only)
  selfcheck  verify hpdi_replay against scipy's bilinear transform of H(s)
  tune       grid-search both baselines on the dev split, freeze parameters
  evaluate   score frozen parameters on the eval split; optional paired
             bootstrap against another method's per-record results

Manifest CSV columns: path, family, hs, split (dev|eval), seed (optional).
Record CSV columns (names configurable): t, a_up, z_ref, Tz_hat (optional).
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import json
import math
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy import fft as sfft

# --------------------------------------------------------------------------- data


@dataclass
class Entry:
    path: Path
    family: str
    hs: float
    split: str
    seed: str = ""


@dataclass
class Record:
    t: np.ndarray
    a: np.ndarray
    z: np.ndarray | None
    tz: np.ndarray | None


@dataclass
class Columns:
    t: str = "t"
    a: str = "a_up"
    z: str = "z_ref"
    tz: str = "Tz_hat"
    flip_z: bool = False


def read_manifest(path: Path) -> list[Entry]:
    out = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            p = Path(row["path"])
            if not p.is_absolute():
                p = path.parent / p
            out.append(Entry(p, row["family"], float(row["hs"]), row["split"].strip().lower(),
                             row.get("seed", "") or ""))
    if not out:
        sys.exit(f"empty manifest {path}")
    return out


def load_record(path: Path, cols: Columns) -> Record:
    with open(path) as fh:
        header = [h.strip().strip('"') for h in fh.readline().strip().split(",")]
    idx = {name: i for i, name in enumerate(header)}
    for need in (cols.t, cols.a):
        if need not in idx:
            sys.exit(f"{path}: column '{need}' not found in {header}")
    use = [idx[cols.t], idx[cols.a]]
    has_z, has_tz = cols.z in idx, cols.tz in idx
    if has_z:
        use.append(idx[cols.z])
    if has_tz:
        use.append(idx[cols.tz])
    data = np.loadtxt(path, delimiter=",", skiprows=1, usecols=use, ndmin=2)
    t, a = data[:, 0], data[:, 1]
    z = data[:, 2] if has_z else None
    if z is not None and cols.flip_z:
        z = -z
    tz = data[:, -1] if has_tz else None
    if np.any(np.diff(t) <= 0):
        sys.exit(f"{path}: time not strictly increasing")
    return Record(t, a, z, tz)


def score_mask(t: np.ndarray, score_last: float) -> np.ndarray:
    return t >= t[-1] - score_last


def pct_hs(rms: float, hs: float) -> float:
    return 100.0 * rms / hs


# --------------------------------------------------------------------------- FDDI


class FDDI:
    """Zero-phase FFT double integration with reflect padding.

    The record is mean-removed, reflect-padded by pad_s at both ends, and the
    padding (only) is cosine-tapered to zero, so scored samples near the record
    end are not attenuated by an analysis window.
    """

    def __init__(self, t: np.ndarray, a: np.ndarray, pad_s: float = 120.0):
        dt = np.diff(t)
        self.fs = 1.0 / float(np.median(dt))
        self.uniform = bool(np.max(np.abs(dt * self.fs - 1.0)) < 1e-3)
        self.t = t
        if self.uniform:
            tu, au = t, a
        else:  # resample to a uniform grid, map the result back afterwards
            tu = np.arange(t[0], t[-1], 1.0 / self.fs)
            au = np.interp(tu, t, a)
        self.tu = tu
        good = np.isfinite(au)
        if not good.all():
            au = np.interp(tu, tu[good], au[good])
        x = au - au.mean()
        n = len(x)
        npad = min(int(round(pad_s * self.fs)), n - 1)
        xp = np.pad(x, npad, mode="reflect")
        if npad > 0:
            ramp = 0.5 * (1.0 - np.cos(np.pi * np.arange(npad) / npad))
            xp[:npad] *= ramp
            xp[-npad:] *= ramp[::-1]
        self.n, self.npad = n, npad
        self.nfft = sfft.next_fast_len(len(xp), real=True)
        self.A = sfft.rfft(xp, self.nfft)
        self.f = sfft.rfftfreq(self.nfft, 1.0 / self.fs)
        with np.errstate(divide="ignore"):
            self.inv_w2 = np.where(self.f > 0, -1.0 / (2 * np.pi * self.f) ** 2, 0.0)

    @staticmethod
    def weight(f: np.ndarray, f1: float, f2: float, fh1: float | None = None, fh2: float | None = None):
        w = np.clip((f - f1) / max(f2 - f1, 1e-12), 0.0, 1.0)
        w = 0.5 * (1.0 - np.cos(np.pi * w))
        if fh1 is not None and fh2 is not None:
            v = np.clip((fh2 - f) / max(fh2 - fh1, 1e-12), 0.0, 1.0)
            w *= 0.5 * (1.0 - np.cos(np.pi * v))
        return w

    def displacement(self, f1: float, taper: float, fh1: float | None = None, fh2: float | None = None):
        W = self.weight(self.f, f1, f1 * taper, fh1, fh2)
        y = sfft.irfft(self.A * W * self.inv_w2, self.nfft)[self.npad:self.npad + self.n]
        return y if self.uniform else np.interp(self.t, self.tu, y)


def fddi_cutoff(cfg: dict, rec: Record, mask: np.ndarray) -> float:
    if cfg["mode"] == "fixed":
        return cfg["value"]
    tz = rec.tz[mask] if rec.tz is not None else np.array([])
    tz = tz[np.isfinite(tz) & (tz > 0)]
    if tz.size == 0:
        raise ValueError("period mode needs a finite Tz_hat column")
    return cfg["value"] / float(np.median(tz))


def fddi_grid(has_tz: bool) -> list[dict]:
    grid = []
    for taper in (1.25, 1.5, 2.0):
        for f1 in np.geomspace(0.005, 0.2, 24):
            grid.append({"mode": "fixed", "value": float(f1), "taper": taper})
        if has_tz:
            for r in np.geomspace(0.05, 1.0, 24):
                grid.append({"mode": "period", "value": float(r), "taper": taper})
    return grid


def fddi_scores(entries: list[Entry], cols: Columns, grid: list[dict], score_last: float, pad_s: float):
    """Returns array [config, record] of % H_s."""
    out = np.full((len(grid), len(entries)), np.nan)
    for j, e in enumerate(entries):
        rec = load_record(e.path, cols)
        if rec.z is None:
            sys.exit(f"{e.path}: reference column '{cols.z}' missing")
        mask = score_mask(rec.t, score_last)
        eng = FDDI(rec.t, rec.a, pad_s)
        for i, cfg in enumerate(grid):
            try:
                f1 = fddi_cutoff(cfg, rec, mask)
            except ValueError:
                continue
            y = eng.displacement(f1, cfg["taper"])
            err = y[mask] - rec.z[mask]
            out[i, j] = pct_hs(float(np.sqrt(np.nanmean(err ** 2))), e.hs)
    return out


# --------------------------------------------------------------------------- HPDI


def hpdi_grid(has_tz: bool) -> list[dict]:
    grid = []
    for n in range(3, 7):
        for m in range(3, n + 1):
            for fc in np.geomspace(0.005, 0.15, 16):
                grid.append({"n": n, "m": m, "mode": "fixed", "value": float(fc)})
            if has_tz:
                for r in np.geomspace(0.03, 0.6, 16):
                    grid.append({"n": n, "m": m, "mode": "period", "value": float(r)})
    return grid


def replay_args(binary: Path, cols: Columns, opts: argparse.Namespace) -> list[str]:
    args = [str(binary), "--t-col", cols.t, "--a-col", cols.a, "--z-col", cols.z, "--tz-col", cols.tz,
            "--smoothing-periods", str(opts.smoothing_periods), "--refresh", str(opts.refresh),
            "--fallback-period", str(opts.fallback_period)]
    if cols.flip_z:
        args.append("--flip-z")
    if opts.float32:
        args.append("--float")
    return args


def hpdi_scores(entries: list[Entry], cols: Columns, grid: list[dict], opts: argparse.Namespace):
    binary = Path(opts.hpdi_bin)
    if not binary.exists():
        sys.exit(f"hpdi_replay binary not found: {binary} (build tools/hpdi_replay.cpp)")
    with tempfile.TemporaryDirectory() as td:
        gpath = Path(td) / "grid.txt"
        gpath.write_text("".join(f"{g['n']} {g['m']} {g['mode']} {g['value']!r}\n" for g in grid))

        def one(e: Entry) -> np.ndarray:
            cmd = replay_args(binary, cols, opts) + ["--in", str(e.path), "--grid", str(gpath),
                                                     "--score-last", str(opts.score_last)]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                raise RuntimeError(f"hpdi_replay failed on {e.path}: {res.stderr.strip()}")
            rows = list(csv.DictReader(res.stdout.splitlines()))
            if len(rows) != len(grid):
                raise RuntimeError(f"{e.path}: expected {len(grid)} rows, got {len(rows)}")
            return np.array([pct_hs(float(r["rms_m"]), e.hs) for r in rows])

        with cf.ThreadPoolExecutor(max_workers=opts.jobs) as ex:
            cols_out = list(ex.map(one, entries))
    return np.stack(cols_out, axis=1)


# --------------------------------------------------------------------------- selection


def select(scores: np.ndarray, grid: list[dict], keep=lambda g: True) -> tuple[dict, float]:
    obj = np.nanmean(scores, axis=1)
    cand = [i for i, g in enumerate(grid) if keep(g) and np.isfinite(obj[i])]
    if not cand:
        raise ValueError("no admissible configuration")
    best = min(cand, key=lambda i: obj[i])
    return dict(grid[best]), float(obj[best])


def edge_warning(name: str, best: dict, grid: list[dict]) -> None:
    same = [g["value"] for g in grid if g["mode"] == best["mode"]
            and all(g.get(k) == best.get(k) for k in ("n", "m", "taper") if k in best)]
    if same and best["value"] in (min(same), max(same)):
        print(f"  warning: {name} optimum at grid edge ({best['mode']} {best['value']:.4g}); widen the grid")


def cmd_tune(opts: argparse.Namespace) -> None:
    cols = columns_from(opts)
    entries = [e for e in read_manifest(Path(opts.manifest)) if e.split == "dev"]
    if not entries:
        sys.exit("no dev records in manifest")
    has_tz = load_record(entries[0].path, cols).tz is not None
    print(f"tuning on {len(entries)} dev records; period-scaled modes {'on' if has_tz else 'off (no T_z column)'}")

    fgrid = fddi_grid(has_tz)
    fs = fddi_scores(entries, cols, fgrid, opts.score_last, opts.pad)
    fbest, fobj = select(fs, fgrid)
    edge_warning("FDDI", fbest, fgrid)

    hgrid = hpdi_grid(has_tz)
    hs = hpdi_scores(entries, cols, hgrid, opts)
    hbest, hobj = select(hs, hgrid)
    cbest, cobj = select(hs, hgrid, lambda g: g["m"] == g["n"])
    pbest, pobj = select(hs, hgrid, lambda g: g["m"] < g["n"])
    for name, b in (("HPDI", hbest), ("HPDI-classic", cbest), ("HPDI-compensated", pbest)):
        edge_warning(name, b, hgrid)

    frozen = {
        "protocol": {"score_last_s": opts.score_last, "dev_records": [str(e.path) for e in entries],
                     "objective": "mean normalized vertical RMS, % H_s, over dev records",
                     "fddi_pad_s": opts.pad, "hpdi_smoothing_periods": opts.smoothing_periods,
                     "hpdi_refresh_s": opts.refresh, "hpdi_fallback_period_s": opts.fallback_period,
                     "hpdi_float32": bool(opts.float32)},
        "FDDI": {**fbest, "dev_pct_hs": fobj},
        "HPDI-classic": {**cbest, "dev_pct_hs": cobj},
        "HPDI-compensated": {**pbest, "dev_pct_hs": pobj},
        "HPDI-best": {**hbest, "dev_pct_hs": hobj},
    }
    Path(opts.out).write_text(json.dumps(frozen, indent=2))
    for k in ("FDDI", "HPDI-classic", "HPDI-compensated"):
        print(f"  {k:18s} dev {frozen[k]['dev_pct_hs']:.3f} % H_s  {describe(frozen[k])}")
    print(f"frozen parameters -> {opts.out}")


def describe(cfg: dict) -> str:
    unit = "Hz" if cfg["mode"] == "fixed" else "x 1/T_z"
    head = f"n={cfg['n']} m={cfg['m']} " if "n" in cfg else f"taper={cfg['taper']} "
    return f"{head}{cfg['mode']} {cfg['value']:.4g} {unit}"


# --------------------------------------------------------------------------- evaluate


def paired_compare(rows: list[dict], other_csv: Path, other_method: str, family: str, reps: int) -> None:
    """Seed-level paired contrast, primary-endpoint style: per seed, mean % H_s over
    the family's records; bootstrap of seed differences (PCG64, seed 20260317) and
    exact sign-flip test. Differences are baseline minus the other method."""
    other: dict[tuple[str, float], float] = {}
    with open(other_csv, newline="") as fh:
        for r in csv.DictReader(fh):
            if r["method"] == other_method and r["family"] == family:
                other[(r["seed"], round(float(r["hs"]), 6))] = float(r["pct_hs"])
    methods = sorted({r["method"] for r in rows})
    print(f"\npaired vs {other_method}, family {family} (baseline - {other_method}, pp of H_s)")
    for meth in methods:
        per_seed: dict[str, list[float]] = {}
        for r in rows:
            if r["method"] != meth or r["family"] != family:
                continue
            key = (r["seed"], round(float(r["hs"]), 6))
            if key in other:
                per_seed.setdefault(r["seed"], []).append(float(r["pct_hs"]) - other[key])
        counts = {len(v) for v in per_seed.values()}
        if len(per_seed) < 2 or len(counts) != 1:
            print(f"  {meth}: need >= 2 seeds with matching record sets (got {len(per_seed)})")
            continue
        d = np.array([np.mean(v) for v in per_seed.values()])
        rng = np.random.Generator(np.random.PCG64(20260317))
        boot = rng.choice(d, size=(reps, d.size), replace=True).mean(axis=1)
        lo, hi = np.percentile(boot, [2.5, 97.5])
        n = d.size
        if n <= 20:
            signs = ((np.arange(2 ** n)[:, None] >> np.arange(n)) & 1) * 2 - 1
            stats = np.abs((signs * d).mean(axis=1))
            p = float(np.mean(stats >= abs(d.mean()) - 1e-12))
            ptxt = f"exact sign-flip p={p:.4f}"
        else:
            ptxt = "sign-flip skipped (n > 20)"
        print(f"  {meth:18s} {d.mean():+.3f} [{lo:+.3f}, {hi:+.3f}]  n_seeds={n}  {ptxt}")


def cmd_evaluate(opts: argparse.Namespace) -> None:
    cols = columns_from(opts)
    frozen = json.loads(Path(opts.params).read_text())
    entries = [e for e in read_manifest(Path(opts.manifest)) if e.split == opts.split]
    if not entries:
        sys.exit(f"no '{opts.split}' records in manifest")
    rows: list[dict] = []

    fcfg = frozen["FDDI"]
    for e in entries:
        rec = load_record(e.path, cols)
        mask = score_mask(rec.t, opts.score_last)
        y = FDDI(rec.t, rec.a, opts.pad).displacement(fddi_cutoff(fcfg, rec, mask), fcfg["taper"])
        rms = float(np.sqrt(np.nanmean((y[mask] - rec.z[mask]) ** 2)))
        rows.append(result_row("FDDI", fcfg, e, rms))

    hnames = [k for k in ("HPDI-classic", "HPDI-compensated") if k in frozen]
    hgrid = [frozen[k] for k in hnames]
    hs = hpdi_scores(entries, cols, [{k: g[k] for k in ("n", "m", "mode", "value")} for g in hgrid], opts)
    for i, name in enumerate(hnames):
        for j, e in enumerate(entries):
            rows.append(result_row(name, hgrid[i], e, hs[i, j] * e.hs / 100.0))

    with open(opts.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"per-record results -> {opts.out}")
    summarize(rows)
    if opts.compare_csv:
        for fam in sorted({r["family"] for r in rows}):
            paired_compare(rows, Path(opts.compare_csv), opts.compare_method, fam, opts.bootstrap)


def result_row(method: str, cfg: dict, e: Entry, rms: float) -> dict:
    return {"method": method, "config": describe(cfg), "family": e.family, "hs": e.hs, "seed": e.seed,
            "path": str(e.path), "rms_m": f"{rms:.6g}", "pct_hs": f"{pct_hs(rms, e.hs):.4f}"}


def summarize(rows: list[dict]) -> None:
    keys = sorted({(r["method"], r["family"], float(r["hs"])) for r in rows})
    print(f"\n{'method':18s} {'family':10s} {'Hs':>6s} {'mean %Hs':>9s} {'sd':>6s} {'n':>3s}")
    for meth, fam, hs in keys:
        v = np.array([float(r["pct_hs"]) for r in rows
                      if r["method"] == meth and r["family"] == fam and float(r["hs"]) == hs])
        sd = v.std(ddof=1) if v.size > 1 else float("nan")
        print(f"{meth:18s} {fam:10s} {hs:6.2f} {v.mean():9.3f} {sd:6.3f} {v.size:3d}")


# --------------------------------------------------------------------------- selfcheck


def hpdi_tf(n: int, m: int, fc: float):
    w = 2 * np.pi * fc
    c = [1.0]
    g = np.pi / (2 * n)
    for k in range(1, n + 1):
        c.append(c[-1] * math.cos((k - 1) * g) / math.sin(k * g))
    den = [c[k] * w ** (n - k) for k in range(n, -1, -1)]                       # descending
    num = [c[p + 2] * w ** (n - p - 2) if p + 2 >= m else 0.0 for p in range(n - 2, -1, -1)]
    return np.array(num), np.array(den)


def cmd_selfcheck(opts: argparse.Namespace) -> None:
    from scipy import signal

    fs, dur = 200.0, 900.0
    t = np.arange(int(fs * dur)) / fs
    rng = np.random.default_rng(7)
    a = (-(2 * np.pi / 7) ** 2 * np.sin(2 * np.pi / 7 * t) - 0.4 * (2 * np.pi / 4) ** 2 * np.sin(2 * np.pi / 4 * t + 1)
         + 0.03 + 0.01 * rng.standard_normal(t.size))
    ok = True
    with tempfile.TemporaryDirectory() as td:
        rec = Path(td) / "rec.csv"
        np.savetxt(rec, np.column_stack([t, a]), delimiter=",", header="t,a_up", comments="", fmt="%.9g")
        for n, m, fc in ((4, 3, 0.04), (5, 5, 0.03), (6, 3, 0.02), (3, 2, 0.05)):
            out = Path(td) / "out.csv"
            subprocess.run([opts.hpdi_bin, "--in", str(rec), "--tz-col", "none", "--dump", str(out),
                            "--n", str(n), "--m", str(m), "--mode", "fixed", "--value", str(fc)], check=True)
            y = np.loadtxt(out, delimiter=",", skiprows=1, usecols=1)
            num, den = hpdi_tf(n, m, fc)
            z, p, k = signal.tf2zpk(num, den)
            sos = signal.zpk2sos(*signal.bilinear_zpk(z, p, k, fs))
            ref = signal.sosfilt(sos, a)
            sel = t > 500  # past the slowest start-up transient (n = 6)
            rel = np.sqrt(np.mean((y[sel] - ref[sel]) ** 2)) / np.sqrt(np.mean(ref[sel] ** 2))
            good = rel < 1e-6
            ok &= good
            print(f"  n={n} m={m} fc={fc}: relative RMS difference {rel:.2e} {'ok' if good else 'FAIL'}")
    if not ok:
        sys.exit(1)
    print("selfcheck passed: hpdi_replay == bilinear(H(s))")


# --------------------------------------------------------------------------- synth


def jonswap(f: np.ndarray, hs: float, tp: float, gamma: float = 3.3) -> np.ndarray:
    fp = 1.0 / tp
    sig = np.where(f <= fp, 0.07, 0.09)
    s = f ** -5 * np.exp(-1.25 * (fp / f) ** 4) * gamma ** np.exp(-((f - fp) ** 2) / (2 * sig ** 2 * fp ** 2))
    df = np.gradient(f)
    return s * (hs ** 2 / 16.0) / np.sum(s * df)


def cmd_synth(opts: argparse.Namespace) -> None:
    """Vertical-only JONSWAP records with constant + random-walk bias and white
    noise. No vessel response, tilt or Mahony proxy: plumbing tests only."""
    outdir = Path(opts.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fs, dur = opts.fs, opts.duration
    t = np.arange(int(fs * dur)) / fs
    seas = [(0.27, 3.0), (1.5, 5.7), (4.0, 8.5), (8.5, 11.4)]
    rows = []
    for split, seeds in (("dev", range(opts.seeds)), ("eval", range(100, 100 + opts.seeds))):
        for seed in seeds:
            for hs, tp in seas:
                rng = np.random.default_rng(seed * 1000 + int(hs * 100))
                f = np.linspace(0.03, 1.0, 900)
                S = jonswap(f, hs, tp)
                amp = np.sqrt(2 * S * np.gradient(f))
                ph = rng.uniform(0, 2 * np.pi, f.size)
                W = 2 * np.pi * f
                z = np.zeros_like(t)
                acc = np.zeros_like(t)
                for A, w, p in zip(amp, W, ph):
                    s = np.sin(w * t + p)
                    z += A * s
                    acc -= A * w * w * s
                bias = rng.normal(0, 0.05) + np.cumsum(rng.normal(0, opts.bias_rw / np.sqrt(fs), t.size))
                a_meas = acc + bias + rng.normal(0, opts.noise, t.size)
                m0, m2 = np.sum(S * np.gradient(f)), np.sum(S * W ** 2 * np.gradient(f))
                tz = np.full_like(t, 2 * np.pi * np.sqrt(m0 / m2))
                name = outdir / f"jonswap_hs{hs:g}_{split}_s{seed}.csv"
                np.savetxt(name, np.column_stack([t, a_meas, z, tz]), delimiter=",",
                           header="t,a_up,z_ref,Tz_hat", comments="", fmt="%.8g")
                rows.append({"path": name.name, "family": "JONSWAP", "hs": hs, "split": split, "seed": seed})
    with open(outdir / "manifest.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["path", "family", "hs", "split", "seed"], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} records + manifest.csv to {outdir}")


# --------------------------------------------------------------------------- CLI


def columns_from(opts: argparse.Namespace) -> Columns:
    return Columns(opts.t_col, opts.a_col, opts.z_col, opts.tz_col, opts.flip_z)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--manifest", required=True)
        p.add_argument("--hpdi-bin", default="build/hpdi_replay")
        p.add_argument("--t-col", default="t")
        p.add_argument("--a-col", default="a_up")
        p.add_argument("--z-col", default="z_ref")
        p.add_argument("--tz-col", default="Tz_hat")
        p.add_argument("--flip-z", action="store_true", help="reference is NED-down")
        p.add_argument("--score-last", type=float, default=900.0)
        p.add_argument("--pad", type=float, default=120.0, help="FDDI reflect padding, s")
        p.add_argument("--smoothing-periods", type=float, default=0.75)
        p.add_argument("--refresh", type=float, default=0.1)
        p.add_argument("--fallback-period", type=float, default=6.0)
        p.add_argument("--float32", action="store_true", help="replay HPDI in float, as on ESP32-S3")
        p.add_argument("--jobs", type=int, default=os.cpu_count() or 2)

    p = sub.add_parser("tune")
    common(p)
    p.add_argument("--out", default="hpdi_fddi_frozen.json")
    p.set_defaults(func=cmd_tune)

    p = sub.add_parser("evaluate")
    common(p)
    p.add_argument("--params", default="hpdi_fddi_frozen.json")
    p.add_argument("--split", default="eval")
    p.add_argument("--out", default="hpdi_fddi_results.csv")
    p.add_argument("--compare-csv", help="long CSV with method,family,hs,seed,pct_hs")
    p.add_argument("--compare-method", default="OU-III")
    p.add_argument("--bootstrap", type=int, default=10000)
    p.set_defaults(func=cmd_evaluate)

    p = sub.add_parser("selfcheck")
    p.add_argument("--hpdi-bin", default="build/hpdi_replay")
    p.set_defaults(func=cmd_selfcheck)

    p = sub.add_parser("synth")
    p.add_argument("--outdir", default="synthetic_heave")
    p.add_argument("--seeds", type=int, default=2)
    p.add_argument("--fs", type=float, default=200.0)
    p.add_argument("--duration", type=float, default=1200.0)
    p.add_argument("--noise", type=float, default=0.0148)
    p.add_argument("--bias-rw", type=float, default=5e-4)
    p.set_defaults(func=cmd_synth)

    opts = ap.parse_args()
    opts.func(opts)


if __name__ == "__main__":
    main()
