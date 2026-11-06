#!/usr/bin/env python3
"""Validation runner of the standalone shared-modes QMEKF (EXPERIMENTAL).

Companion of docs/kalman-carma-shared-modes.md and
doc/kalman_carma_shared_modes/.  Drives
tests/kalman_carma_shared_modes/kalman_carma_shared_modes-sim (dataset replay)
and kalman_carma_shared_modes-bench (paired Monte Carlo, carma3_bench's
engine) and writes reports/results/kalman_carma_shared_modes/.

Stages (resumable; each appends to its own CSV):

    smoke     one paired row (JONSWAP, seed 1) with the reduction checks
    dataset   replay of the pinned vessel-RAO records (oceanography-waves-lib),
              CSM, CSM_OFF and the CARMA3 baseline A41, default seeds
    run       final paired Monte Carlo: 18 scenarios x 50 seeds (1..50),
              filters A41 (baseline), SHM:est (research Option C), CSM
              (deployable default), CSM:off (its fallback model)
    robust    coupling-mismatch injections: reversed, rotated, phase-shifted,
              doubled and true external priors, with and without the
              confidence, on three seas x 30 seeds; frozen keel-frame coupling
              under heading changes
    traces    1 Hz confidence traces and 10 Hz error time series, seed 1
    spectra   Welch spectra of truth and estimates, seeds 1..10, three seas
    timing    host update time, one job on an otherwise idle machine
    report    statistics, figures (SVG + PGF), LaTeX tables, manifest.json

Seeds: the only fitted hyper-parameters (estimator horizon 240 s, envelope band
0.15 f_op) were chosen on screening seeds 1001-1020 by the coupled-prior
study; the confidence thresholds are fixed a priori from the statistics of the
mismatch test (E[D^2] = 2 for a consistent prior).  Seeds 1..50 never entered
a choice.  Uncertainty: paired bootstrap over independent realizations (seeds),
never over the autocorrelated 200 Hz samples.

Run: python3 tools/carma_shared_modes_validation.py {smoke,dataset,run,robust,traces,spectra,timing,report,all} [--jobs 4]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TDIR = ROOT / "tests" / "kalman_carma_shared_modes"
SIM = TDIR / "kalman_carma_shared_modes-sim"        # dataset replay
BENCH = TDIR / "kalman_carma_shared_modes-bench"    # paired Monte Carlo
OUT = ROOT / "reports" / "results" / "kalman_carma_shared_modes"
RAW = OUT / "raw_runs.csv"
ROB = OUT / "robustness_runs.csv"
TIM = OUT / "timing_runs.csv"
SMOKE = OUT / "smoke_runs.csv"
SCREEN = OUT / "screen_runs.csv"
DATASET = OUT / "dataset_summary.csv"
TRACE_DIR = OUT / "traces"
SPEC_DIR = OUT / "spectra"
FIG_DIR = OUT / "figures"

sys.path.insert(0, str(Path(__file__).resolve().parent))
import carma3_coupled_benchmark as research  # noqa: E402  (scenario table of the coupled-prior study)

J29 = research.J29
SCENARIOS = dict(research.SCENARIOS)
FINAL_FILTERS = "A41,A41@C,SHM:est,CSM,CSM:off"
BASE = "A41"
ROB_SEAS = [J29, "mixed_swell_wind", "pm_fp0.25_hs1"]
ROB_FILTERS = ("A41,CSM:off,CSM/ext-true,CSM:noconf/ext-true,CSM/ext-flip,CSM:noconf/ext-flip,"
               "CSM/ext-rot90,CSM:noconf/ext-rot90,CSM/ext-ph45,CSM:noconf/ext-ph45,CSM/ext-amp2,CSM:noconf/ext-amp2")
ROB_TURNS = f"{J29}+turns"
ROB_TURN_FILTERS = "A41,CSM,CSM:noconf,CSM/fz400/keel,CSM:noconf/fz400/keel"
TRACE_SEAS = [J29, "dir_change", f"{J29}+turns", f"{J29}+on400w60+off900w60", f"{J29}+maneuver", "mixed_swell_wind",
              "crossing", "calm"]
TRACE_FILTERS = "A41,SHM:est,CSM"
TRACE_EXTRA = {J29: "A41,CSM/ext-flip,CSM:noconf/ext-flip"}
SPEC_SEAS = [J29, "mixed_swell_wind", "crossing"]
SPEC_FILTERS = ["A41", "SHM:est", "CSM"]
TIMING_FILTERS = ["A41", "SHM:est", "CSM", "CSM:off"]
LABEL = {"A41": "CARMA3 baseline (A41)", "A41@C": "CARMA3 baseline, Evidence excitation (A41@C)", "SHM:est": "Option C research (SHM:est)", "CSM": "Shared modes (CSM)",
         "CSM:off": "CSM coupling off", "CSM:noconf": "CSM without confidence"}
COLORS = {"A41": "#555555", "A41@C": "#8c564b", "SHM:est": "#ff9896", "CSM": "#d62728", "CSM:off": "#1f77b4", "truth": "#000000"}


# ---------------------------------------------------------------------------
# running
# ---------------------------------------------------------------------------

def require(exe=None):
    for x in ([exe] if exe else [SIM, BENCH]):
        if not x.exists():
            sys.exit(f"missing {x.relative_to(ROOT)}: make -C tests/kalman_carma_shared_modes build")


def read_rows(path):
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f))


def cmd(sc, seed, filters, score_start=300.0):
    return ["--level", "full", "--scenario", sc, "--seed", str(seed), "--filters", filters,
            "--score-start", f"{score_start:g}"]


def run_rows(cmds, out_path, jobs, extra=("--extended", "--coupling"), env=None):
    require()
    extra = list(extra)
    hdr = subprocess.run([str(BENCH), "--header"] + extra, capture_output=True, text=True, check=True).stdout.strip()
    new_file = not out_path.exists()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    failures = 0
    with out_path.open("a") as out:
        if new_file:
            out.write(hdr + "\n")
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futs = {ex.submit(subprocess.run, [str(BENCH)] + extra + c, capture_output=True, text=True,
                              env=env): c for c in cmds}
            for k, fut in enumerate(as_completed(futs), 1):
                c, p = futs[fut], fut.result()
                if p.returncode != 0:
                    failures += 1
                    print(f"ERROR {' '.join(c)}: rc={p.returncode} {p.stderr[-400:]}", flush=True)
                    continue
                out.write(p.stdout)
                out.flush()
                if k % 20 == 0 or k == len(cmds):
                    print(f"  {k}/{len(cmds)}", flush=True)
    return failures


def pending(rows, work):
    have = {(r["scenario"], r["seed"]) for r in rows}
    return [w for w in work if (w[0], str(w[1])) not in have]


def stage_smoke(args):
    if SMOKE.exists():
        SMOKE.unlink()
    run_rows([cmd(f"{J29}+T600", 1, "A41,A41@C,SHM:est,CSM:noconf@A,CSM,CSM:off", 200.0)], SMOKE, 1)
    rows = {r["filter"]: r for r in read_rows(SMOKE)}

    def rel(a, b, k):
        x, y = float(rows[a][k]), float(rows[b][k])
        return abs(x - y) / max(abs(y), 1e-12)

    # With the excitation paired (Legacy, @A), CSM:noconf is the research
    # filter SHM:est; the two differ only by float summation order.
    for k in ("rms_h", "rms_x", "rms_y", "rms_z", "nees_x", "nees_y", "nees_z"):
        if k in rows["SHM:est"]:
            assert rel("CSM:noconf@A", "SHM:est", k) < 1e-3, ("CSM:noconf@A must reproduce SHM:est", k)
    # Coupling off is the CARMA(4,1) law of the baseline at the same excitation.
    for k in ("rms_h", "rms_z"):
        if k in rows["A41@C"]:
            assert rel("CSM:off", "A41@C", k) < 0.03, ("CSM:off must track A41@C", k)
    for f in ("CSM", "CSM:off"):
        assert rows[f]["failed"] == "0", f
    print("smoke ok:", {f: round(1000 * float(r["rms_h"]), 2) for f, r in rows.items()})


# Screening of the amplitude commit deadband (the only hyper-parameter this
# filter adds that is chosen from data), on the screening seeds 1001-1020 only.
# Criterion fixed before the run: the tolerance with the largest mean paired
# 3-D RMS improvement over CSM without deadband across the wave seas, subject
# to no sea's heave RMS worse by more than 1 %; ties go to the smaller value.
SCREEN_SEAS = list(research.SCREEN_SCENARIOS) + ["small_amplitude", "jonswap_fp0.45_hs0.4_g3.3"]
SCREEN_TOLS = [0.1, 0.2, 0.3]


def stage_screen(args):
    filters = "A41,CSM," + ",".join(f"CSM/atol{t:g}" for t in SCREEN_TOLS)
    work = [(sc, s) for sc in SCREEN_SEAS for s in range(1001, 1021)]
    todo = pending(read_rows(SCREEN), work)
    print(f"screen: {len(todo)} of {len(work)} rows to run", flush=True)
    run_rows([cmd(sc, s, filters) for sc, s in todo], SCREEN, args.jobs)
    df = load(SCREEN)
    if df.empty:
        return
    df["rms3"] = np.sqrt(df.rms_x ** 2 + df.rms_y ** 2 + df.rms_z ** 2)
    best, best_gain = 0.0, 0.0
    for t in SCREEN_TOLS:
        f = f"CSM/atol{t:g}"
        gains, worst_z = [], 1.0
        for sc in SCREEN_SEAS:
            if sc in ("calm", "stationary"):
                continue
            g = paired(df, sc, f, "rms3", base="CSM")
            z = paired(df, sc, f, "rms_z", base="CSM")
            if g and z:
                gains.append(g["change"])
                worst_z = min(worst_z, z["change"])
        mean_gain = float(np.mean(gains)) if gains else float("nan")
        ok = worst_z >= -0.01
        print(f"  atol {t:g}: mean 3-D change vs CSM {100 * mean_gain:+.2f} %, worst heave {100 * worst_z:+.2f} %"
              f" -> {'admissible' if ok else 'rejected'}")
        if ok and mean_gain > best_gain + 1e-4:
            best, best_gain = t, mean_gain
    print(f"screen: selected amplitude tolerance {best:g} (mean 3-D gain {100 * best_gain:+.2f} %)")


def stage_run(args):
    work = [(sc, s) for sc in SCENARIOS for s in range(1, args.seeds + 1)]
    todo = pending(read_rows(RAW), work)
    print(f"run: {len(todo)} of {len(work)} rows to run", flush=True)
    run_rows([cmd(sc, s, FINAL_FILTERS) for sc, s in todo], RAW, args.jobs)


def stage_robust(args):
    work = [(sc, s, ROB_FILTERS) for sc in ROB_SEAS for s in range(1, args.robust_seeds + 1)]
    work += [(ROB_TURNS, s, ROB_TURN_FILTERS) for s in range(1, args.robust_seeds + 1)]
    have = {(r["scenario"], r["seed"]) for r in read_rows(ROB)}
    todo = [w for w in work if (w[0], str(w[1])) not in have]
    print(f"robust: {len(todo)} of {len(work)} rows to run", flush=True)
    run_rows([cmd(sc, s, f) for sc, s, f in todo], ROB, args.jobs)


def stage_traces(args):
    TRACE_DIR.mkdir(parents=True, exist_ok=True)
    cmds = []
    for sc in TRACE_SEAS:
        tag = re.sub(r"[^A-Za-z0-9_.+-]", "_", sc)
        for fl, sfx in ((TRACE_FILTERS, ""), (TRACE_EXTRA.get(sc), "_ext")):
            if not fl:
                continue
            p = TRACE_DIR / f"{tag}{sfx}_ts.csv"
            env = dict(os.environ, CSM_TRACE=str(TRACE_DIR / f"{tag}{sfx}_conf_"))
            cmds.append((cmd(sc, 1, fl) + ["--timeseries", str(p), "--decimate", "20"], env))
    require()
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(subprocess.run, [str(BENCH), "--extended", "--coupling"] + c, capture_output=True, text=True,
                          env=e) for c, e in cmds]
        for f in futs:
            p = f.result()
            if p.returncode != 0:
                print("ERROR trace:", p.stderr[-300:])
    print(f"traces: {len(cmds)} runs -> {TRACE_DIR.relative_to(ROOT)}")


def welch_cross(x, y, fs, nper=4096):
    from scipy.signal import csd
    f, S = csd(x, y, fs=fs, nperseg=nper, detrend="linear")
    return f, S


def stage_spectra(args):
    SPEC_DIR.mkdir(parents=True, exist_ok=True)
    import tempfile
    for sc in SPEC_SEAS:
        acc = {}
        nseeds = 0
        with tempfile.TemporaryDirectory() as td:
            cmds = []
            for s in range(1, args.spectra_seeds + 1):
                p = Path(td) / f"s{s}.csv"
                cmds.append((p, cmd(sc, s, ",".join(SPEC_FILTERS)) + ["--timeseries", str(p), "--decimate", "10"]))
            with ThreadPoolExecutor(max_workers=args.jobs) as ex:
                list(ex.map(lambda pc: subprocess.run([str(BENCH)] + pc[1], capture_output=True, text=True), cmds))
            for p, _ in cmds:
                if not p.exists():
                    continue
                import pandas as pd
                d = pd.read_csv(p)
                d = d[d.t >= 300.0]
                fs = 1.0 / float(np.median(np.diff(d.t.values)))
                series = {"truth": (d.x_true.values, d.y_true.values, d.z_true.values)}
                for fl in SPEC_FILTERS:
                    series[fl] = (d[f"x_{fl}"].values, d[f"y_{fl}"].values, d[f"z_{fl}"].values)
                for name, (x, y, z) in series.items():
                    comps = {"xx": (x, x), "yy": (y, y), "zz": (z, z), "xz": (x, z), "yz": (y, z), "xy": (x, y)}
                    if name != "truth":
                        tx, ty, tz = series["truth"]
                        comps.update({"ex": (x - tx, x - tx), "ey": (y - ty, y - ty), "ez": (z - tz, z - tz)})
                    for k, (u, v) in comps.items():
                        f, S = welch_cross(u, v, fs)
                        acc.setdefault((name, k), 0)
                        acc[(name, k)] = acc[(name, k)] + S
                nseeds += 1
        if not nseeds:
            continue
        rows = {"f": f.tolist()}
        for (name, k), S in acc.items():
            rows[f"{name}|{k}|re"] = (np.real(S) / nseeds).tolist()
            rows[f"{name}|{k}|im"] = (np.imag(S) / nseeds).tolist()
        tag = re.sub(r"[^A-Za-z0-9_.+-]", "_", sc)
        with (SPEC_DIR / f"{tag}.csv").open("w", newline="") as fh:
            w = csv.writer(fh)
            keys = list(rows)
            w.writerow(keys)
            for i in range(len(rows["f"])):
                w.writerow([f"{rows[k][i]:.6e}" for k in keys])
        print(f"spectra {sc}: {nseeds} seeds", flush=True)


def stage_timing(args):
    if TIM.exists():
        TIM.unlink()
    run_rows([cmd(J29, s, f) for s in range(1, 6) for f in TIMING_FILTERS], TIM, 1)


def stage_dataset(args):
    """Replay of the pinned vessel-RAO records with the default seeds: CSM and
    CSM_OFF (kalman_carma_shared_modes-sim) and the CARMA3 baseline A41
    (tests/kalman_carma/carma3-sim, the same records, seeds and adapter); the
    metrics are the replays' own VALIDATION_METRICS lines (last 900 s)."""
    require(SIM)
    carma3_sim = ROOT / "tests" / "kalman_carma" / "carma3-sim"
    if not carma3_sim.exists():
        sys.exit("missing tests/kalman_carma/carma3-sim: make -C tests/kalman_carma build")
    runs = [("A41", carma3_sim, {"W3D_CARMA3_VARIANT": "A41"}),
            ("CSM_OFF", SIM, {"W3D_CSM_VARIANT": "CSM_OFF"}),
            ("CSM", SIM, {"W3D_CSM_VARIANT": "CSM"})]
    rows = []
    for variant, exe, extra in runs:
        env = {k: v for k, v in os.environ.items() if not k.startswith("W3D_")}
        env.update(extra, W3D_WRITE_TIMESERIES="1")
        p = subprocess.run([str(exe)], cwd=exe.parent, capture_output=True, text=True, env=env)
        if p.returncode != 0:
            print(f"ERROR dataset {variant}: rc={p.returncode} {p.stderr[-300:]}")
        for line in p.stdout.splitlines():
            if not line.startswith("VALIDATION_METRICS family="):
                continue
            kv = dict(tok.split("=", 1) for tok in line.split()[1:] if "=" in tok)
            row = {"variant": variant, "record": kv.get("input", "").replace("wave_data_", "").replace(".csv", "")}
            for k in ("disp_x_rms_m", "disp_y_rms_m", "disp_z_rms_m", "disp_3d_rms_m", "disp_z_pct_hs",
                      "disp_3d_pct_refmax", "roll_rms_deg", "pitch_rms_deg", "yaw_rms_deg", "accel_3d_rms_mps2",
                      "accel_bias_3d_rms_mps2", "gyro_bias_3d_rms_radps"):
                try:
                    row[k] = float(kv[k])
                except (KeyError, ValueError):
                    row[k] = float("nan")
            rows.append(row)
    if not rows:
        return
    OUT.mkdir(parents=True, exist_ok=True)
    with DATASET.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"dataset: {len(rows)} record replays -> {DATASET.relative_to(ROOT)}")
    # one representative time series per filter for the tracking figures
    TRACE_DIR.mkdir(parents=True, exist_ok=True)
    import pandas as pd
    for variant, pat in (("A41", ROOT / "tests" / "kalman_carma" / "w3d_jonswap_H1.500*_fusion_CARMA3_A41.csv"),
                         ("CSM", TDIR / "w3d_jonswap_H1.500*_fusion_CARMA_SHARED_MODES_CSM.csv")):
        src = next(iter(sorted(pat.parent.glob(pat.name))), None)
        if src:
            d = pd.read_csv(src)
            step = max(1, int(round(0.1 / float(np.median(np.diff(d.time.values))))))
            keep = [c for c in d.columns if c.startswith(("time", "disp_", "vel_", "roll_", "pitch_", "yaw_", "acc_bias",
                                                           "gyro_bias", "tau_applied", "sigma_a_applied", "freq_tracker"))]
            d[keep].iloc[::step].to_csv(TRACE_DIR / f"dataset_jonswap_H1.5_{variant}.csv", index=False, float_format="%.6g")


def dataset_metrics(path, t0=None):
    import pandas as pd
    d = pd.read_csv(path)
    cols = set(d.columns)
    tcol = d.columns[0]
    t = d[tcol].values
    if t0 is None:
        t0 = t[-1] - 900.0 if t[-1] > 1000 else t[0] + 0.25 * (t[-1] - t[0])
    m = t >= t0
    out = {}
    for ax in "xyz":
        a, b = f"disp_ref_{ax}", f"disp_est_{ax}"
        if a in cols and b in cols:
            e = d[b].values[m] - d[a].values[m]
            out[f"rms_{ax}"] = float(np.sqrt(np.mean(e ** 2)))
        a, b = f"vel_ref_{ax}", f"vel_est_{ax}"
        if a in cols and b in cols:
            e = d[b].values[m] - d[a].values[m]
            out[f"rmsv_{ax}"] = float(np.sqrt(np.mean(e ** 2)))
    for c in ("roll", "pitch", "yaw"):
        a, b = f"{c}_ref", f"{c}_est"
        if a in cols and b in cols:
            e = (d[b].values[m] - d[a].values[m] + 180.0) % 360.0 - 180.0
            out[f"{c}_rms_deg"] = float(np.sqrt(np.mean(e ** 2)))
    return out


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------

def paired(df, sc, cand, col, base=BASE, n_boot=4000):
    a = df[(df.scenario == sc) & (df["filter"] == base)].set_index("seed")[col]
    b = df[(df.scenario == sc) & (df["filter"] == cand)].set_index("seed")[col]
    idx = a.index.intersection(b.index)
    a, b = a.loc[idx].values.astype(float), b.loc[idx].values.astype(float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 3:
        return None
    rng = np.random.default_rng(20261106)
    k = rng.integers(0, len(a), (n_boot, len(a)))
    bs = 1.0 - b[k].mean(axis=1) / a[k].mean(axis=1)
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return {"base": float(a.mean()), "cand": float(b.mean()), "change": float(1.0 - b.mean() / a.mean()),
            "ci": [float(lo), float(hi)], "n": int(len(a))}


def load(path):
    import pandas as pd
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def summarize(df, filters, scenarios, base=BASE):
    out = {}
    for sc in scenarios:
        sub = df[df.scenario == sc]
        if sub.empty:
            continue
        out[sc] = {}
        for f in filters:
            g = sub[sub["filter"] == f]
            if g.empty:
                continue
            r = {"n": int(len(g)), "failed": int(g.failed.sum())}
            for c in ("rms_x", "rms_y", "rms_z", "rms_h", "rmsv_x", "rmsv_y", "rmsv_z", "peak_x", "peak_y", "peak_z",
                      "lf_x", "lf_y", "lf_z", "su_rms_x", "su_rms_z", "roll_rms_deg", "pitch_rms_deg", "yaw_rms_deg",
                      "rmsa_x", "rmsa_y", "rmsa_z", "pct_z", "conv_z", "f60_z"):
                if c in g:
                    r[c] = float(g[c].mean())
            for c in ("nees_pv6", "nees_pos3", "nis", "att_nees", "nees9_native"):
                if c in g:
                    v = g[c].astype(float)
                    r[c + "_mean"] = float(v.mean())
                    r[c + "_median"] = float(v.median())
            for c in ("cov95_x", "cov95_y", "cov95_z", "nis_in95", "wi_ac1_z", "wi_cc_xz", "cpl_coh_xz", "cpl_ph_xz"):
                if c in g:
                    r[c] = float(g[c].astype(float).mean())
            if f != base:
                for c in ("rms_x", "rms_y", "rms_z", "rms_h", "rmsv_z", "nees_pv6"):
                    p = paired(df, sc, f, c, base)
                    if p:
                        r["chg_" + c] = p
            out[sc][f] = r
    return out


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

def git_rev():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    except OSError:
        return "unknown"


def manifest(stats):
    mk = (ROOT / "Makefile").read_text()
    ver = re.search(r"SIM_DATA_VERSION \?= (\S+)", mk)
    zipn = re.search(r"SIM_DATA_ZIP_NAME \?= (\S+)", mk)
    def sha(p):
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    inputs = {}
    for p in sorted(TDIR.glob("wave_data_*.csv"))[:40]:
        inputs[p.name] = sha(p)
    results = {p.name: sha(p) for p in sorted(OUT.glob("*.csv"))}
    return {
        "git_revision": git_rev(),
        "dataset": {"repository": "bareboat-necessities/oceanography-waves-lib",
                    "version": ver.group(1) if ver else None, "archive": zipn.group(1) if zipn else None,
                    "record_sha256": inputs},
        "monte_carlo": {"scenarios": SCENARIOS, "seeds": "1..50 (final), 1..30 (robustness), 1..10 (spectra)",
                        "screening_seeds": "1001..1020 (coupled-prior study; estimator horizon and band)",
                        "duration_s": "per scenario (tests/kalman_carma/CarmaBenchTruth.h; default 1200 s, scored from 300 s)",
                        "filters": FINAL_FILTERS, "robustness_filters": ROB_FILTERS, "turn_filters": ROB_TURN_FILTERS,
                        "sample_rate_hz": 200},
        "configuration": {"coupling": "estimated (240 s, band 0.15 f_op)",
                          "confidence": "src/kalman_carma_shared_modes/CarmaCouplingConfidence.h defaults"},
        "software": {"python": platform.python_version(), "numpy": np.__version__, "platform": platform.platform()},
        "result_sha256": results,
        "n_rows": stats.get("n_rows", {}),
    }


def fmt_chg(r):
    if not r:
        return "--"
    return f"{100 * r['change']:+.1f} [{100 * r['ci'][0]:+.1f}, {100 * r['ci'][1]:+.1f}]"


def stage_report(args):
    OUT.mkdir(parents=True, exist_ok=True)
    df = load(RAW)
    rob = load(ROB)
    stats = {"n_rows": {"final": int(len(df)), "robust": int(len(rob))}}
    md = ["# Shared-modes QMEKF validation (generated)\n",
          "Paired change against the CARMA3 baseline A41 (+ = better) with 95 % bootstrap intervals over seeds. "
          "NEES of [p, v] (6 = consistent in the MEAN; the median of a chi^2_6 variable is 5.35) and NIS.\n"]
    if not df.empty:
        F = FINAL_FILTERS.split(",")
        s = summarize(df, F, SCENARIOS)
        stats["final"] = s
        # The same rows paired against the same-excitation baseline A41@C:
        # CSM (Evidence) - A41@C is the effect of the shared modes alone.
        se = summarize(df, [f for f in F if f != BASE], SCENARIOS, base="A41@C")
        stats["final_vs_evidence"] = se
        md.append("## Final Monte Carlo against the same-excitation baseline A41@C\n")
        md.append("| scenario | CSM x / y / z (%) | CSM:off x / y / z (%) |")
        md.append("|---|---|---|")
        for sc, purpose in SCENARIOS.items():
            r = se.get(sc, {})
            cells = []
            for f in ("CSM", "CSM:off"):
                q = r.get(f, {})
                cells.append(" / ".join(f"{100 * q['chg_rms_' + a]['change']:+.1f}" if q.get("chg_rms_" + a) else "--"
                                        for a in "xyz"))
            md.append(f"| {purpose} | " + " | ".join(cells) + " |")
        md.append("")
        md.append("## Final Monte Carlo (50 seeds)\n")
        md.append("| scenario | A41 rms x/y/z (mm) | " + " | ".join(f"{f} x / y / z (%)" for f in F[1:]) + " | NEES6 mean (A41 / " + " / ".join(F[1:]) + ") |")
        md.append("|---" * (3 + len(F) - 1) + "|")
        for sc, purpose in SCENARIOS.items():
            r = s.get(sc)
            if not r or BASE not in r:
                continue
            b = r[BASE]
            cells = [f"{1000 * b['rms_x']:.1f} / {1000 * b['rms_y']:.1f} / {1000 * b['rms_z']:.1f}"]
            for f in F[1:]:
                q = r.get(f, {})
                cells.append(" / ".join(f"{100 * q['chg_rms_' + a]['change']:+.1f}" if q.get("chg_rms_" + a) else "--"
                                        for a in "xyz"))
            ne = " / ".join(f"{r[f]['nees_pv6_mean']:.1f}" if f in r and "nees_pv6_mean" in r[f] else "--" for f in F)
            md.append(f"| {purpose} | " + " | ".join(cells) + f" | {ne} |")
        md.append("")
    if not rob.empty:
        R = [f for f in ROB_FILTERS.split(",") if f != BASE]
        seas = ROB_SEAS + [ROB_TURNS]
        rs = summarize(rob, list(dict.fromkeys(R + ROB_TURN_FILTERS.split(","))), seas)
        stats["robust"] = rs
        md.append("## Coupling-mismatch robustness (30 seeds)\n")
        md.append("| injection | sea | rms h / z change vs A41 (%) | NEES6 mean / median | mean c |")
        md.append("|---|---|---|---|---|")
        for sc in seas:
            for f, q in rs.get(sc, {}).items():
                if f == BASE:
                    continue
                md.append(f"| {f} | {sc} | {fmt_chg(q.get('chg_rms_h'))} / {fmt_chg(q.get('chg_rms_z'))} | "
                          f"{q.get('nees_pv6_mean', float('nan')):.1f} / {q.get('nees_pv6_median', float('nan')):.1f} | -- |")
        md.append("")
    ds = load(DATASET)
    if not ds.empty:
        md.append("## Dataset replay (pinned vessel-RAO records, default seeds)\n")
        md.append(ds.to_markdown(index=False, floatfmt=".4f") if hasattr(ds, "to_markdown") else ds.to_string())
        md.append("")
    (OUT / "summary.md").write_text("\n".join(md) + "\n")
    (OUT / "stats.json").write_text(json.dumps(stats, indent=1))
    (OUT / "manifest.json").write_text(json.dumps(manifest(stats), indent=1))
    write_tex(stats)
    try:
        figures(stats, FIG_DIR, ("svg",))
    except Exception as e:  # noqa: BLE001 - figures must not hide the tables
        print("figures failed:", e)
    print(f"report -> {OUT.relative_to(ROOT)}")


# ---------------------------------------------------------------------------
# LaTeX tables (doc/kalman_carma_shared_modes)
# ---------------------------------------------------------------------------

TEX_DIR = ROOT / "doc" / "kalman_carma_shared_modes"


def tex(t):
    return t.replace("_", "\\_").replace("%", "\\%").replace("&", "\\&").replace("->", "$\\to$")


def write_tex(stats):
    TEX_DIR.mkdir(parents=True, exist_ok=True)
    s = stats.get("final")
    L = []
    if s:
        F = FINAL_FILTERS.split(",")
        L += ["\\begin{table*}[t]", "\\caption{Final paired Monte Carlo, 50 seeds per sea: baseline RMS (mm) and the "
              "change of each filter against it, $x$/$y$/$z$ in \\% (positive = better). Mean NEES of "
              "$[\\vct p,\\vct v]$ (6 = consistent).}", "\\label{tab:final}", "\\centering\\scriptsize",
              "\\begin{adjustbox}{max width=\\textwidth}", "\\begin{tabular}{l" + "c" * (len(F) + 1) + "}", "\\toprule",
              "Sea & A41 (mm) & " + " & ".join(tex(f) for f in F[1:]) + " & NEES " + "/".join(tex(f) for f in F) + " \\\\",
              "\\midrule"]
        for sc, purpose in SCENARIOS.items():
            r = s.get(sc)
            if not r or BASE not in r:
                continue
            b = r[BASE]
            cells = [f"{1000 * b['rms_x']:.1f}/{1000 * b['rms_y']:.1f}/{1000 * b['rms_z']:.1f}"]
            for f in F[1:]:
                q = r.get(f, {})
                cells.append("/".join(f"{100 * q['chg_rms_' + a]['change']:+.0f}" if q.get("chg_rms_" + a) else "--"
                                      for a in "xyz"))
            ne = "/".join(f"{r[f]['nees_pv6_mean']:.1f}" if f in r and "nees_pv6_mean" in r[f] else "--" for f in F)
            L.append(f"{tex(purpose)} & " + " & ".join(cells) + f" & {ne} \\\\")
        L += ["\\bottomrule", "\\end{tabular}", "\\end{adjustbox}", "\\end{table*}", ""]
    rs = stats.get("robust")
    if rs:
        L += ["\\begin{table*}[t]", "\\caption{Coupling-mismatch injections, 30 seeds: paired change of the horizontal / "
              "vertical RMS against A41 (\\%, 95\\,\\% interval) and the mean / median NEES of $[\\vct p,\\vct v]$, "
              "with and without the coupling confidence.}", "\\label{tab:robust}", "\\centering\\scriptsize",
              "\\begin{adjustbox}{max width=\\textwidth}", "\\begin{tabular}{llccc}", "\\toprule",
              "Injection & Sea & horizontal & heave & NEES mean/median \\\\", "\\midrule"]
        for sc in ROB_SEAS + [ROB_TURNS]:
            for f, q in rs.get(sc, {}).items():
                if f == BASE:
                    continue
                L.append(f"{tex(f)} & {tex(sc)} & {tex(fmt_chg(q.get('chg_rms_h')))} & {tex(fmt_chg(q.get('chg_rms_z')))} & "
                         f"{q.get('nees_pv6_mean', float('nan')):.1f}/{q.get('nees_pv6_median', float('nan')):.1f} \\\\")
        L += ["\\bottomrule", "\\end{tabular}", "\\end{adjustbox}", "\\end{table*}", ""]
    (TEX_DIR / "csm-results-generated.tex-part").write_text("\n".join(L) + "\n")


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def figures(stats, outdir, exts=("svg",)):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    import shutil
    from matplotlib.text import Text
    plt.rcParams.update({"font.size": 8, "axes.grid": True, "grid.alpha": 0.3, "svg.hashsalt": "csm",
                         "pgf.rcfonts": False})
    if not shutil.which(plt.rcParams["pgf.texsystem"]) and shutil.which("lualatex"):
        plt.rcParams["pgf.texsystem"] = "lualatex"
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    def tex_safe(fig):
        # PGF text is typeset by LaTeX: "^" and "%" are escaped by matplotlib,
        # a bare "_" outside math mode is not ("Missing $ inserted").
        for t in fig.findobj(Text):
            txt = t.get_text()
            if "_" in txt:
                parts = txt.split("$")
                for i in range(0, len(parts), 2):
                    parts[i] = parts[i].replace("\\_", "_").replace("_", "\\_")
                t.set_text("$".join(parts))

    def save(fig, name):
        fig.tight_layout()
        for e in sorted(exts, key=lambda x: x == "pgf"):   # raw text formats first, PGF last
            if e == "pgf":
                tex_safe(fig)
            fig.savefig(outdir / f"{name}.{e}", dpi=150)
        plt.close(fig)

    def replay(variant):
        """Dataset-replay time series (JONSWAP Hs 1.5 m): the committed
        decimated trace, else a replay CSV of this run (CI copies the
        replays into plots/kalman_carma_shared_modes) or of the test dir."""
        p = TRACE_DIR / f"dataset_jonswap_H1.5_{variant}.csv"
        if p.exists():
            return pd.read_csv(p)
        fam = "CARMA3_A41" if variant == "A41" else f"CARMA_SHARED_MODES_{variant}"
        pat = f"w3d_jonswap_H1.500*_fusion_{fam}.csv"
        src = next(iter(sorted(Path.cwd().glob(pat)) + sorted(TDIR.glob(pat))
                        + sorted((ROOT / "tests" / "kalman_carma").glob(pat))), None)
        if src is None:
            return None
        d = pd.read_csv(src)
        step = max(1, int(round(0.1 / float(np.median(np.diff(d.time.values))))))
        return d.iloc[::step].reset_index(drop=True)

    # A-D: dataset replay tracking (truth vs estimate), attitude errors, biases
    dc, da = replay("CSM"), replay("A41")
    if dc is not None:
        t = dc.iloc[:, 0].values
        w = (t >= t[-1] - 120.0)
        for kind, lab, unit, name in (("disp", "displacement", "m", "tracking_position"),
                                      ("vel", "velocity", "m/s", "tracking_velocity")):
            if f"{kind}_ref_x" not in dc:
                continue
            fig, axs = plt.subplots(3, 1, figsize=(7.0, 5.4), sharex=True)
            for i, ax in enumerate("xyz"):
                axs[i].plot(t[w], dc[f"{kind}_ref_{ax}"].values[w], color="k", lw=1.4, label="truth")
                if da is not None and len(da) == len(dc):
                    axs[i].plot(t[w], da[f"{kind}_est_{ax}"].values[w], color=COLORS["A41"], lw=0.9, label="A41")
                axs[i].plot(t[w], dc[f"{kind}_est_{ax}"].values[w], color=COLORS["CSM"], lw=0.9, ls="--", label="CSM")
                axs[i].set_ylabel(f"{['surge x', 'sway y', 'heave z'][i]} ({unit})")
            axs[0].legend(loc="upper right", ncol=3, fontsize=7)
            axs[2].set_xlabel("time (s)")
            fig.suptitle(f"Dataset replay JONSWAP Hs 1.5 m: {lab}, last 120 s", fontsize=9)
            save(fig, name)
        cols = [(f"{x}_ref", f"{x}_est") for x in ("roll", "pitch", "yaw")]
        if all(r in dc and e in dc for r, e in cols):
            fig, axs = plt.subplots(3, 1, figsize=(7.0, 5.0), sharex=True)
            for i, (r, e) in enumerate(cols):
                for d, f in ((da, "A41"), (dc, "CSM")):
                    if d is None:
                        continue
                    err = (d[e].values - d[r].values + 180.0) % 360.0 - 180.0
                    axs[i].plot(d.iloc[:, 0].values, err, color=COLORS[f], lw=0.8, label=f)
                axs[i].set_ylabel(f"{r.split('_')[0]} error (deg)")
            axs[0].legend(fontsize=7)
            axs[2].set_xlabel("time (s)")
            save(fig, "tracking_attitude")
        if "acc_bias_est_x" in dc:
            fig, axs = plt.subplots(2, 1, figsize=(7.0, 4.6), sharex=True)
            for i, ax in enumerate("xyz"):
                col = ["#1f77b4", "#2ca02c", "#d62728"][i]
                axs[0].plot(t, dc[f"acc_bias_{ax}"].values, color=col, lw=1.4, alpha=0.5, label=f"true {ax}")
                axs[0].plot(t, dc[f"acc_bias_est_{ax}"].values, color=col, lw=0.8, ls="--", label=f"CSM {ax}")
                axs[1].plot(t, 57.29578 * dc[f"gyro_bias_{ax}"].values, color=col, lw=1.4, alpha=0.5, label=f"true {ax}")
                axs[1].plot(t, 57.29578 * dc[f"gyro_bias_est_{ax}"].values, color=col, lw=0.8, ls="--", label=f"CSM {ax}")
            axs[0].set_ylabel("accelerometer bias (m/s$^2$)")
            axs[1].set_ylabel("gyroscope bias (deg/s)")
            axs[1].set_xlabel("time (s)")
            axs[0].legend(fontsize=6, ncol=3)
            save(fig, "bias_convergence")
    # E: spectra
    for sc in SPEC_SEAS:
        p = SPEC_DIR / f"{re.sub(r'[^A-Za-z0-9_.+-]', '_', sc)}.csv"
        if not p.exists():
            continue
        d = pd.read_csv(p)
        f = d.f.values
        m = (f > 0.03) & (f < 1.5)
        fig, axs = plt.subplots(2, 3, figsize=(10.0, 5.6))
        for j, k in enumerate(("xx", "yy", "zz")):
            for name in ["truth"] + SPEC_FILTERS:
                axs[0, j].loglog(f[m], d[f"{name}|{k}|re"].values[m], color=COLORS.get(name, "k"),
                                 lw=1.8 if name == "truth" else 1.0, label=name)
            for name in SPEC_FILTERS:
                key = f"{name}|e{k[0]}|re"
                if key in d:
                    axs[0, j].loglog(f[m], d[key].values[m], color=COLORS.get(name), lw=0.8, ls=":")
            axs[0, j].set_title(f"PSD {k[0]} (solid) and error PSD (dotted)", fontsize=8)
        for j, k in enumerate(("xz", "yz")):
            for name in ["truth"] + SPEC_FILTERS:
                S = d[f"{name}|{k}|re"].values + 1j * d[f"{name}|{k}|im"].values
                Sa = d[f"{name}|{k[0]}{k[0]}|re"].values
                Sb = d[f"{name}|{k[1]}{k[1]}|re"].values
                coh = np.abs(S) ** 2 / np.maximum(Sa * Sb, 1e-30)
                axs[1, j].semilogx(f[m], np.degrees(np.angle(S[m])), color=COLORS.get(name, "k"),
                                   lw=1.8 if name == "truth" else 1.0, label=name)
                if j == 0:
                    axs[1, 2].semilogx(f[m], coh[m], color=COLORS.get(name, "k"), lw=1.8 if name == "truth" else 1.0,
                                       label=name)
            axs[1, j].set_ylim(-180, 180)
            axs[1, j].set_title(f"{k[0]}-{k[1]} cross-spectral phase (deg)", fontsize=8)
        axs[1, 2].set_title("x-z magnitude-squared coherence", fontsize=8)
        axs[0, 0].legend(fontsize=6)
        for ax in axs[1]:
            ax.set_xlabel("frequency (Hz)")
        fig.suptitle(f"Spectra, {SCENARIOS.get(sc, sc)} (Welch, 10 seeds)", fontsize=9)
        save(fig, f"spectra_{re.sub(r'[^A-Za-z0-9]', '_', sc)}")
    # F / H: confidence and adaptation traces
    for sc in TRACE_SEAS + [J29 + "_ext"]:
        tag = re.sub(r"[^A-Za-z0-9_.+-]", "_", sc.replace("_ext", "")) + ("_ext" if sc.endswith("_ext") else "")
        files = sorted(TRACE_DIR.glob(f"{tag}_conf_*.csv"))
        ts = TRACE_DIR / f"{tag}_ts.csv"
        if not files:
            continue
        fig, axs = plt.subplots(4, 1, figsize=(7.2, 7.0), sharex=True)
        for fp in files:
            d = pd.read_csv(fp)
            lab = fp.stem.split("_conf_")[-1].replace("_", ":", 1).replace("_", "/")
            axs[0].plot(d.t, d.confidence, lw=1.0, label=lab)
            axs[1].plot(d.t, d.coh_xz, lw=0.9, label=lab)
            axs[2].plot(d.t, d.f_op, lw=0.9, label=lab)
            axs[2].plot(d.t, d.sigma_z, lw=0.9, ls="--")
            axs[3].semilogy(d.t, np.maximum(d.d2, 1e-3), lw=0.8, label=lab)
        axs[0].set_ylabel("confidence c")
        axs[1].set_ylabel("applied x-z coherence")
        axs[2].set_ylabel(r"$f_{op}$ (Hz), $\sigma_z$ (dashed)")
        axs[3].set_ylabel("mismatch $D^2$")
        axs[3].axhline(16.0, color="k", lw=0.6, ls=":")
        axs[3].axhline(6.0, color="k", lw=0.6, ls=":")
        axs[0].legend(fontsize=6)
        axs[3].set_xlabel("time (s)")
        fig.suptitle(f"Adaptation and coupling confidence: {SCENARIOS.get(sc.replace('_ext', ''), sc)}"
                     + (" (reversed external prior)" if sc.endswith("_ext") else ""), fontsize=9)
        save(fig, f"confidence_{re.sub(r'[^A-Za-z0-9]', '_', sc)}")
        if ts.exists():
            d = pd.read_csv(ts)
            fig, ax = plt.subplots(figsize=(7.2, 3.0))
            filt = [c[2:] for c in d.columns if c.startswith("x_") and c not in ("x_true", "x_total")]
            for fl in filt:
                e = np.sqrt((d[f"x_{fl}"] - d.x_true) ** 2 + (d[f"y_{fl}"] - d.y_true) ** 2)
                ax.plot(d.t, pd.Series(e).rolling(300, min_periods=10).apply(lambda v: np.sqrt(np.mean(v ** 2)), raw=True),
                        lw=0.9, label=fl)
            ax.set_ylabel("rolling 30 s horizontal RMS error (m)")
            ax.set_xlabel("time (s)")
            ax.legend(fontsize=6)
            save(fig, f"error_{re.sub(r'[^A-Za-z0-9]', '_', sc)}")
    # G: comparison across seas
    s = stats.get("final")
    if s:
        F = [f for f in FINAL_FILTERS.split(",") if f != BASE]
        seas = [sc for sc in SCENARIOS if sc in s and BASE in s[sc]]
        fig, axs = plt.subplots(1, 2, figsize=(10.0, 6.4), sharey=True)
        y = np.arange(len(seas))
        wbar = 0.8 / len(F)
        for j, col in enumerate(("rms_h", "rms_z")):
            for k, f in enumerate(F):
                v, lo, hi = [], [], []
                for sc in seas:
                    q = s[sc].get(f, {}).get("chg_" + col)
                    v.append(100 * q["change"] if q else np.nan)
                    lo.append(100 * (q["change"] - q["ci"][0]) if q else 0)
                    hi.append(100 * (q["ci"][1] - q["change"]) if q else 0)
                axs[j].barh(y + (k - len(F) / 2 + 0.5) * wbar, v, wbar, xerr=[lo, hi], color=COLORS.get(f), label=LABEL.get(f, f))
            axs[j].axvline(0, color="k", lw=0.6)
            axs[j].set_xlabel(("horizontal" if j == 0 else "heave") + " RMS change vs A41 (%), + = better")
        axs[0].set_yticks(y)
        axs[0].set_yticklabels([SCENARIOS[sc] for sc in seas], fontsize=7)
        axs[0].legend(fontsize=7, loc="lower left")
        save(fig, "comparison_rms")
        # I: consistency
        fig, axs = plt.subplots(1, 2, figsize=(10.0, 4.2))
        allf = FINAL_FILTERS.split(",")
        for k, f in enumerate(allf):
            ne = [s[sc][f].get("nees_pv6_mean", np.nan) if f in s[sc] else np.nan for sc in seas]
            nm = [s[sc][f].get("nees_pv6_median", np.nan) if f in s[sc] else np.nan for sc in seas]
            xs = np.arange(len(seas)) + (k - len(allf) / 2 + 0.5) * 0.18
            axs[0].semilogy(xs, ne, "o", color=COLORS.get(f), label=f"{f} mean", ms=4)
            axs[0].semilogy(xs, nm, "x", color=COLORS.get(f), ms=4)
            ni = [s[sc][f].get("nis_mean", np.nan) if f in s[sc] else np.nan for sc in seas]
            axs[1].plot(xs, ni, "o", color=COLORS.get(f), label=f, ms=4)
        axs[0].axhline(6.0, color="k", lw=0.8, ls="--")
        axs[0].set_ylabel("NEES [p, v] (o mean, x median; 6 = consistent mean)")
        axs[1].set_ylabel("mean accelerometer NIS (3 dof nominal)")
        for ax in axs:
            ax.set_xticks(np.arange(len(seas)))
            ax.set_xticklabels([SCENARIOS[sc] for sc in seas], rotation=70, fontsize=6)
        axs[0].legend(fontsize=6)
        save(fig, "consistency_nees_nis")
    # H: robustness bars
    rs = stats.get("robust")
    if rs:
        fig, axs = plt.subplots(1, len(ROB_SEAS) + 1, figsize=(12.0, 4.6), sharey=False)
        for j, sc in enumerate(ROB_SEAS + [ROB_TURNS]):
            r = rs.get(sc, {})
            fl = [f for f in r if f != BASE]
            v = [100 * r[f]["chg_rms_h"]["change"] if r[f].get("chg_rms_h") else np.nan for f in fl]
            cols = ["#d62728" if ":noconf" not in f else "#7f7f7f" for f in fl]
            axs[j].barh(np.arange(len(fl)), v, color=cols)
            axs[j].set_yticks(np.arange(len(fl)))
            axs[j].set_yticklabels(fl, fontsize=6)
            axs[j].axvline(0, color="k", lw=0.6)
            axs[j].set_title(SCENARIOS.get(sc, sc), fontsize=8)
            axs[j].set_xlabel("horizontal RMS change vs A41 (%)")
        save(fig, "robustness_rms")
    # J: cost
    size = ROOT / "reports" / "results" / "kalman_carma_shared_modes" / "esp32s3_size.csv"
    tim = load(TIM)
    if size.exists() or not tim.empty:
        fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.4))
        if not tim.empty:
            g = tim.groupby("filter")[["runtime_ns", "runtime_p99_ns"]].mean() / 1000.0
            g = g.reindex([f for f in TIMING_FILTERS if f in g.index])
            axs[0].bar(np.arange(len(g)) - 0.2, g.runtime_ns, 0.4, label="mean")
            axs[0].bar(np.arange(len(g)) + 0.2, g.runtime_p99_ns, 0.4, label="99th percentile")
            axs[0].set_xticks(np.arange(len(g)))
            axs[0].set_xticklabels(g.index, fontsize=7)
            axs[0].set_ylabel("host update time (us)")
            axs[0].legend(fontsize=7)
        if size.exists():
            z = pd.read_csv(size)
            x = np.arange(len(z))
            axs[1].bar(x - 0.2, z.code_bytes / 1000, 0.4, label="code (kB)")
            axs[1].bar(x + 0.2, (z.sizeof_qmekf + z.sizeof_wrapper + z.sizeof_top) / 1000, 0.4,
                       label="objects incl. heap (kB)")
            ax2 = axs[1].twinx()
            ax2.plot(x, z.max_stack_frame / 1000, "kD", label="largest frame (kB)")
            ax2.plot(x, z.update_chain_stack / 1000, "ks", mfc="none", label="static update chain (kB)")
            ax2.set_ylabel("stack (kB)")
            ax2.legend(fontsize=6, loc="upper right")
            axs[1].set_xticks(x)
            axs[1].set_xticklabels(z.variant, fontsize=7)
            axs[1].legend(fontsize=6, loc="upper left")
        save(fig, "cost_time_memory")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=["screen", "smoke", "dataset", "run", "robust", "traces", "spectra", "timing", "report",
                                      "figures", "all"])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--robust-seeds", type=int, default=30)
    ap.add_argument("--spectra-seeds", type=int, default=10)
    ap.add_argument("--outdir", default=None, help="figures stage: output directory")
    ap.add_argument("--ext", default="svg", help="figures stage: comma-separated formats (svg, pgf, png, pdf)")
    args = ap.parse_args(argv)
    if args.stage == "figures":
        stats = json.loads((OUT / "stats.json").read_text()) if (OUT / "stats.json").exists() else {}
        figures(stats, Path(args.outdir) if args.outdir else FIG_DIR, tuple(args.ext.split(",")))
        return 0
    stages = ["smoke", "dataset", "run", "robust", "traces", "spectra", "timing", "report"] if args.stage == "all" \
        else [args.stage]
    for st in stages:
        print(f"== {st}", flush=True)
        globals()[f"stage_{st}"](args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
