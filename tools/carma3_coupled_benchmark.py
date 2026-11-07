#!/usr/bin/env python3
"""Coupled three-axis CARMA benchmark driver (EXPERIMENTAL research).

Companion of docs/carma-coupled-three-axis.md and
doc/kalman_carma/kalman_carma_coupled.tex.  Every row runs
tests/kalman_carma/carma3_bench on one (scenario, seed): one truth and one
sensor stream fed identically to every filter of the row (paired design).

Seeds: the final and ablation runs use 1..50 / 1..30, which never enter a
choice.  Every subcommand is resumable.  The shared directional modes
(Experiment B) are the standalone filter of tools/carma_shared_modes_study.py.

    run        final comparison, 18 scenarios x 50 seeds
    ablate     co-spectrum only, modal coupling with oracle and reversed
               oracle priors, seven seas x 30 seeds
    spectra    time series (20 Hz) of seeds 1..10 on five seas, reduced to
               Welch auto- and cross-spectra of truth and estimates
    traces     time series of seed 1 for the figures
    timing     host update time, one job on an otherwise idle machine
    report     summary.md, paired_stats.json, figures (tools/carma3_coupled_report.py)

Run: python3 tools/carma3_coupled_benchmark.py {run,ablate,spectra,traces,timing,report} [--jobs 4]
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TDIR = ROOT / "tests" / "kalman_carma"
BENCH = TDIR / "carma3_bench"
OUT = ROOT / "reports" / "results" / "carma3_coupled"
RAW = OUT / "raw_runs.csv"
ABL_RAW = OUT / "ablation_runs.csv"
TIMING_RAW = OUT / "timing_runs.csv"
SPEC_DIR = OUT / "spectra"
TRACE_DIR = OUT / "traces"

J29 = "jonswap_fp0.29_hs0.8_g3.3"
# Final scenarios (study brief section 5.1) -> purpose.
SCENARIOS = {
    "calm": "calm water (afloat, no waves)",
    "stationary": "still IMU",
    "regular_f0.25_h1": "regular monochromatic wave",
    f"{J29}+sp10000": "unidirectional JONSWAP",
    J29: "directionally spread JONSWAP (cos^8)",
    "pm_fp0.25_hs1": "Pierson-Moskowitz",
    "swell": "long-period swell",
    "jonswap_fp0.45_hs0.4_g3.3": "short-period wind waves",
    "mixed_swell_wind": "mixed wind sea + swell",
    "bimodal": "two systems, different frequencies",
    "crossing": "two crossing systems (90 deg)",
    "dir_change": "changing wave direction (30 -> 120 deg)",
    "abrupt_change": "abrupt change of wave period",
    f"{J29}+on400w60+off900w60": "wave onset and cessation",
    f"{J29}+turns": "vessel heading changes (90 deg zigzag)",
    f"{J29}+maneuver": "vessel manoeuvring",
    "small_amplitude": "very low sea state (5 cm)",
    f"{J29}+magdist+scal": "magnetic disturbances, calibration residuals",
}
FINAL_FILTERS = "A41,MOD:est/nz,COR:est,COR:orc"
ABL_SCENARIOS = [J29, "pm_fp0.25_hs1", "mixed_swell_wind", "bimodal", "crossing", f"{J29}+turns", "calm"]
ABL_FILTERS = "A41,COR:est/xy,MOD:est,MOD:orc/nz,MOD:orc/nz/flip"
SPEC_SCENARIOS = [f"{J29}+sp10000", J29, "mixed_swell_wind", "bimodal", "crossing"]
SPEC_FILTERS = ["A41", "MOD:est/nz", "COR:est"]
TRACE_SCENARIOS = [J29, "dir_change", f"{J29}+on400w60+off900w60", f"{J29}+turns"]
TRACE_FILTERS = "A41,MOD:est/nz,COR:est"
# The committed traces keep only what the report's figures read
# (tools/carma3_coupled_report.py trace_plots): the J29 record at the bench's
# 10 Hz inside its plotted 690-770 s window, the others at 2 Hz (30 s rolling
# RMS and the applied x-z phase).
TRACE_FIG_FILTERS = ["A41", "COR:est", "MOD:est/nz"]
TRACE_J29_WINDOW = (690.0, 770.0)
TRACE_EVERY = 5
TIMING_FILTERS = ["A41", "MOD:est/nz", "COR:est"]


def require(path):
    if not path.exists():
        sys.exit(f"missing {path.relative_to(ROOT)}: build it first (make -C tests/kalman_carma build)")


def read_rows(path):
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f))


def header(extra):
    return subprocess.run([str(BENCH), "--header"] + extra, capture_output=True, text=True, check=True).stdout.strip()


def run_bench(cmds, out_path, jobs, extra=("--extended", "--coupling")):
    require(BENCH)
    extra = list(extra)
    hdr = header(extra)
    new_file = not out_path.exists()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("a") as out:
        if new_file:
            out.write(hdr + "\n")
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futs = {ex.submit(subprocess.run, [str(BENCH)] + extra + c, capture_output=True, text=True): c for c in cmds}
            done = 0
            for fut in as_completed(futs):
                c = futs[fut]
                p = fut.result()
                done += 1
                if p.returncode != 0:
                    print(f"ERROR {' '.join(c)}: rc={p.returncode} {p.stderr[-400:]}", flush=True)
                    continue
                if p.stderr.strip():
                    print(f"  note {' '.join(c[1:6])}: {p.stderr.strip()[-200:]}", flush=True)
                out.write(p.stdout)
                out.flush()
                if done % 20 == 0 or done == len(cmds):
                    print(f"  {done}/{len(cmds)}", flush=True)


def pending(rows, work):
    have = {(r["scenario"], r["seed"]) for r in rows}
    return [w for w in work if (w[0], str(w[1])) not in have]


def cmd(sc, seed, filters, score_start=300.0):
    return ["--level", "full", "--scenario", sc, "--seed", str(seed), "--filters", filters,
            "--score-start", f"{score_start:g}"]


def final(args):
    work = [(sc, s) for sc in SCENARIOS for s in range(1, args.seeds + 1)]
    todo = pending(read_rows(RAW), work)
    print(f"run: {len(todo)} of {len(work)} rows to run", flush=True)
    run_bench([cmd(sc, s, FINAL_FILTERS) for sc, s in todo], RAW, args.jobs)


def ablate(args):
    work = [(sc, s) for sc in ABL_SCENARIOS for s in range(1, args.ablation_seeds + 1)]
    todo = pending(read_rows(ABL_RAW), work)
    print(f"ablate: {len(todo)} of {len(work)} rows to run", flush=True)
    run_bench([cmd(sc, s, ABL_FILTERS) for sc, s in todo], ABL_RAW, args.jobs)


# ---------------------------------------------------------------------------
# Spectra: Welch auto/cross spectra of the truth and of every estimate.
# ---------------------------------------------------------------------------

def welch_csd(x, y, fs, nper):
    """One-sided cross-spectral density (x conj(y) convention: S_xy = E[X Y*])."""
    from scipy.signal import csd
    f, s = csd(x, y, fs=fs, nperseg=nper, noverlap=nper // 2, window="hann", detrend="constant")
    # scipy returns E[conj(X) Y]; conjugate to S_xy = E[X conj(Y)] (truth/model convention).
    return f, np.conj(s)


def spectra_one(ts_path, filters, t0=300.0):
    import pandas as pd
    df = pd.read_csv(ts_path)
    df = df[df["t"] >= t0]
    fs = 1.0 / np.median(np.diff(df["t"].values))
    nper = int(round(102.4 * fs))
    out = {}
    sets = {"truth": ("x_true", "y_true", "z_true")}
    for f in filters:
        sets[f] = (f"x_{f}", f"y_{f}", f"z_{f}")
    for name, cols in sets.items():
        X = [df[c].values - np.mean(df[c].values) for c in cols]
        S = {}
        for i in range(3):
            for j in range(i, 3):
                fr, s = welch_csd(X[i], X[j], fs, nper)
                S[(i, j)] = s
        out[name] = (fr, S)
        if name != "truth":
            E = [df[c].values - df[t].values for c, t in zip(cols, sets["truth"])]
            for i in range(3):
                fr, s = welch_csd(E[i], E[i], fs, nper)
                S[("err", i)] = s
    return out


def spectra(args):
    import shutil
    import tempfile
    SPEC_DIR.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix="carma3_coupled_ts_"))
    try:
        for sc in SPEC_SCENARIOS:
            dest = SPEC_DIR / f"{sc.replace('+', '_')}.csv"
            if dest.exists():
                continue
            cmds, paths = [], []
            for s in range(1, args.spectra_seeds + 1):
                p = tmp / f"{sc}_{s}.csv"
                paths.append(p)
                cmds.append(cmd(sc, s, ",".join(SPEC_FILTERS)) + ["--timeseries", str(p), "--decimate", "10"])
            with ThreadPoolExecutor(max_workers=args.jobs) as ex:
                list(ex.map(lambda c: subprocess.run([str(BENCH)] + c, capture_output=True, text=True, check=True), cmds))
            acc = {}
            fr = None
            for p in paths:
                o = spectra_one(p, SPEC_FILTERS)
                for name, (f, S) in o.items():
                    fr = f
                    for k, v in S.items():
                        acc.setdefault((name, k), []).append(v)
                p.unlink()
            with dest.open("w", newline="") as fh:
                w = csv.writer(fh)
                keys = sorted(acc.keys(), key=str)
                cols = ["f_hz"]
                for name, k in keys:
                    tag = f"{name}|{'err' + str(k[1]) if k[0] == 'err' else str(k[0]) + str(k[1])}"
                    cols += [tag + "|re", tag + "|im"]
                w.writerow(cols)
                means = {k: np.mean(np.array(v), axis=0) for k, v in acc.items()}
                for idx in range(len(fr)):
                    if fr[idx] > 2.5:
                        break
                    row = [f"{fr[idx]:.6g}"]
                    for k in keys:
                        row += [f"{means[k][idx].real:.6e}", f"{means[k][idx].imag:.6e}"]
                    w.writerow(row)
            print(f"  spectra {sc}: {len(paths)} seeds", flush=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def trim_trace(full: Path, out: Path, scenario: str) -> None:
    """Reduce a full bench trace to the columns, rate and window the figures read."""
    import pandas as pd
    t = pd.read_csv(full)
    keep = ["t", "x_true", "y_true", "z_true"]
    keep += [f"{c}_{n}" for n in TRACE_FIG_FILTERS for c in "xyz"]
    keep += [f"phxz_{n}" for n in TRACE_FIG_FILTERS] + ["dyaw_A41"]
    t = t[[k for k in keep if k in t.columns]]
    if scenario == J29:
        t = t[(t.t >= TRACE_J29_WINDOW[0]) & (t.t <= TRACE_J29_WINDOW[1])]
    else:
        t = t.iloc[::TRACE_EVERY]
    t.to_csv(out, index=False, float_format="%.6g")


def traces(args):
    TRACE_DIR.mkdir(parents=True, exist_ok=True)
    cmds, jobs = [], []
    for sc in TRACE_SCENARIOS:
        p = TRACE_DIR / f"{sc.replace('+', '_')}_seed1.csv"
        if p.exists():
            continue
        full = p.with_suffix(".full.csv")
        cmds.append(cmd(sc, 1, TRACE_FILTERS) + ["--timeseries", str(full), "--decimate", "20", "--extended", "--coupling"])
        jobs.append((full, p, sc))
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        list(ex.map(lambda c: subprocess.run([str(BENCH)] + c, capture_output=True, text=True, check=True), cmds))
    for full, p, sc in jobs:
        trim_trace(full, p, sc)
        full.unlink()
    print(f"traces: {len(cmds)} written", flush=True)


def timing(args):
    rows = read_rows(TIMING_RAW)
    have = {(r["scenario"], r["seed"], r["filter"]) for r in rows}
    cmds = []
    for sc in (J29, f"{J29}+jit1", f"{J29}+jit"):
        for s in (1, 2, 3):
            for f in TIMING_FILTERS:
                if (sc, str(s), f) not in have:
                    cmds.append(cmd(sc, s, f))
    print(f"timing: {len(cmds)} single-filter runs, one job", flush=True)
    run_bench(cmds, TIMING_RAW, 1)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=["run", "ablate", "spectra", "traces", "timing", "report"])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--ablation-seeds", type=int, default=30)
    ap.add_argument("--spectra-seeds", type=int, default=10)
    args = ap.parse_args(argv)
    if args.stage == "report":
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import carma3_coupled_report as rep
        return rep.main([])
    {"run": final, "ablate": ablate, "spectra": spectra, "traces": traces, "timing": timing}[args.stage](args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
