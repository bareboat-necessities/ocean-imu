#!/usr/bin/env python3
"""Three-axis CARMA benchmark driver (EXPERIMENTAL; docs/carma-three-axis.md).

Runs tests/kalman_carma/carma3_bench over a fixed grid of (level, scenario,
seed) with every filter of a row fed the same truth and sensor stream, the
pinned vessel-RAO dataset through the shared W3dSimulationRunner, and writes
the report.  Every subcommand is resumable: rows already in the raw CSV are
skipped.

    run          synthetic grid (complete INS and level 1)
    dataset      pinned 28 ft vessel-RAO records, 10 seeds: kalman_ou_iii-sim
                 and carma3-sim
    traces       time series for the figures
    truth        realized motion levels of every scenario (seed 1)
    timing       host update time on an otherwise idle machine (one job):
                 fixed, +-1 % (1 us clock) and +-25 % jittered intervals
    report       tables (summary.md), paired statistics (paired_stats.json)
                 and figures

Paired statistics (tools/carma3_report.py): improvement = 1 - mean RMS
(candidate) / mean RMS (OU-III) with a 95 % bootstrap interval over seeds
(10 000 resamples, fixed seed), win rate and Wilcoxon signed-rank p.

Run: python3 tools/carma3_benchmark.py {run,dataset,traces,truth,report} [--jobs 4]
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(__file__).resolve().parents[1]
TDIR = ROOT / "tests" / "kalman_carma"
BENCH = TDIR / "carma3_bench"
OUT = ROOT / "reports" / "results" / "carma3_benchmark"
RAW = OUT / "raw_runs.csv"
DATASET_RAW = OUT / "dataset_runs.csv"
TRUTH_RAW = OUT / "truth_motion.csv"
TIMING_RAW = OUT / "timing_runs.csv"
TRACE_DIR = OUT / "traces"
OU3_SIM = ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim"
CARMA3_SIM = TDIR / "carma3-sim"

FULL_FILTERS = "OU3,A41,A52,M45,M54,A52S,A41S,H41"
L1_FILTERS = "OU3,A41,A52,M45,M54,A52S,A41S,H41,A30"
J29 = "jonswap_fp0.29_hs0.8_g3.3"

# Principal seas: 100 paired seeds at both levels.
PRIMARY = [J29, "pm_fp0.25_hs1", "mixed_swell_wind", "changing_fp"]
# Secondary: 30 seeds, complete INS.  Headings are the wave propagation
# direction relative to the bow (0 following, 90 beam, 180 head).
SECONDARY = [
    f"{J29}+rel0", f"{J29}+rel90", f"{J29}+rel180", f"{J29}+rel135",
    f"{J29}+rel180+U2.5", f"{J29}+rel0+U2.5", f"{J29}+rel90+U2.5",
    f"{J29}+anchor", f"{J29}+drift", f"{J29}+maneuver",
    "bimodal", "broadband", "abrupt_change",
    "jonswap_fp0.20_hs1_g3.3", "jonswap_fp0.45_hs0.4_g3.3", "pm_fp0.15_hs2.5", "jonswap_fp0.10_hs4_g3.3",
    "bias_drift", "tilted", "large_rotation", "small_amplitude",
    f"{J29}+maggap", f"{J29}+jit", f"{J29}+hlegacy", "mixed_swell_wind+hlegacy",
    "stationary", "stationary+calbias",
]
LONG = [f"{J29}+T3600", "mixed_swell_wind+T3600"]
L1_SECONDARY = ["bimodal", f"{J29}+rel90", f"{J29}+rel180", f"{J29}+anchor", f"{J29}+maneuver"]


def require(path):
    if not path.exists():
        sys.exit(f"missing {path.relative_to(ROOT)}: build it first (make -C {path.parent.relative_to(ROOT)} build)")


def read_rows(path):
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f))


def header():
    return subprocess.run([str(BENCH), "--header"], capture_output=True, text=True, check=True).stdout.strip()


def run_bench(cmds, out_path, jobs):
    """Run carma3_bench for every argument list, appending rows to out_path."""
    require(BENCH)
    hdr = header()
    new_file = not out_path.exists()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("a") as out:
        if new_file:
            out.write(hdr + "\n")
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futs = {ex.submit(subprocess.run, [str(BENCH)] + c, capture_output=True, text=True): c for c in cmds}
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


def grid(args):
    work = []
    for sc in PRIMARY:
        for seed in range(1, args.primary_seeds + 1):
            work.append(("full", sc, seed, FULL_FILTERS))
            work.append(("l1", sc, seed, L1_FILTERS))
    for sc in SECONDARY:
        for seed in range(1, args.secondary_seeds + 1):
            work.append(("full", sc, seed, FULL_FILTERS + (",A30" if sc == f"{J29}+rel90" else "")))
    for sc in L1_SECONDARY:
        for seed in range(1, args.secondary_seeds + 1):
            work.append(("l1", sc, seed, L1_FILTERS))
    for sc in LONG:
        for seed in range(1, args.long_seeds + 1):
            work.append(("full", sc, seed, "OU3,A41,A52,M45,M54,A52S"))
    # Coupled x-y prior with the ORACLE horizontal acceleration correlation of
    # the default sea (tools/carma3_analysis.py coupled_xy_r): a diagnostic
    # upper bound, not a candidate (the estimator does not measure r).
    r = oracle_xy_correlation()
    for seed in range(1, args.secondary_seeds + 1):
        work.append(("full", J29, seed, f"A41,A41C:{r:.3f}"))
    return work


def oracle_xy_correlation():
    import carma3_analysis as ca
    f = np.linspace(1e-3, 5.0, 5000)
    S = ca.cross_spectrum(J29, f, 0.0)
    w4 = (2 * np.pi * f) ** 4
    sxx = np.trapezoid(w4 * np.real(S[:, 0, 0]), f)
    syy = np.trapezoid(w4 * np.real(S[:, 1, 1]), f)
    sxy = np.trapezoid(w4 * np.real(S[:, 0, 1]), f)
    return float(sxy / np.sqrt(sxx * syy))


def cmd_run(args):
    done = {(r["level"], r["scenario"], r["seed"], r["filter"]) for r in read_rows(RAW)}
    cmds = []
    for level, sc, seed, filters in grid(args):
        todo = [f for f in filters.split(",") if (level, sc, str(seed), f) not in done]
        if not todo:
            continue
        cmds.append(["--level", level, "--scenario", sc, "--seed", str(seed), "--filters", ",".join(todo)])
    # Longest first keeps the pool busy at the end.
    cmds.sort(key=lambda c: ("T3600" not in c[3], c[1] != "full"))
    print(f"{len(cmds)} bench invocations to do", flush=True)
    run_bench(cmds, RAW, args.jobs)


# ---------------------------------------------------------------------------
# Pinned dataset
# ---------------------------------------------------------------------------

DATASET_FIELDS = ["family", "seed", "input", "disp_z_rms_m", "disp_z_pct_hs", "disp_x_rms_m", "disp_y_rms_m",
                  "disp_x_pct_hs", "disp_y_pct_hs", "roll_rms_deg", "pitch_rms_deg", "yaw_rms_deg",
                  "accel_z_rms_mps2"]


def parse_metrics(text):
    rows = []
    for line in text.splitlines():
        if line.startswith("VALIDATION_METRICS "):
            rows.append(dict(kv.split("=", 1) for kv in line.split()[1:] if "=" in kv))
    return rows


def dataset_invocation(cmd, extra, seed):
    files = sorted(str(p) for p in TDIR.glob("wave_data_*.csv"))
    if not files:
        sys.exit("no dataset in tests/kalman_carma: make ensure-sim-data")
    env = {k: v for k, v in os.environ.items() if not k.startswith("W3D_")}
    env.update(W3D_VALIDATION_WINDOW_SEC="900", W3D_WRITE_TIMESERIES="0", W3D_COLLECT_ALL_GATES="1", **extra)
    if seed is not None:
        env["W3D_SEED"] = str(seed)
    argv = cmd + sum((["--input", f] for f in files), [])
    p = subprocess.run(argv, capture_output=True, text=True, env=env, cwd=str(TDIR))
    return p, parse_metrics(p.stdout)


def cmd_dataset(args):
    for s in (OU3_SIM, CARMA3_SIM):
        require(s)
    done = {(r["family"], r["seed"]) for r in read_rows(DATASET_RAW)}
    work = []
    for seed in range(1, args.seeds + 1):
        work.append(("OU3", seed, [str(OU3_SIM)], {}))
        for v in args.variants.split(","):
            work.append((v, seed, [str(CARMA3_SIM)], {"W3D_CARMA3_VARIANT": v}))
    work = [w for w in work if (w[0], str(w[1])) not in done]
    print(f"{len(work)} dataset invocations to do", flush=True)
    new_file = not DATASET_RAW.exists()
    DATASET_RAW.parent.mkdir(parents=True, exist_ok=True)
    with DATASET_RAW.open("a") as out:
        wr = csv.DictWriter(out, fieldnames=DATASET_FIELDS, extrasaction="ignore")
        if new_file:
            wr.writeheader()
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(dataset_invocation, cmd, extra, seed): (fam, seed) for fam, seed, cmd, extra in work}
            for fut in as_completed(futs):
                fam, seed = futs[fut]
                p, metrics = fut.result()
                if p.returncode not in (0, 1) or not metrics:
                    print(f"ERROR {fam} seed {seed}: rc={p.returncode} {p.stderr[-300:]}", flush=True)
                    continue
                for m in metrics:
                    m["family"] = fam
                    m["seed"] = str(seed)
                    m["input"] = Path(m.get("input", "")).name
                    wr.writerow(m)
                out.flush()
                print(f"  {fam} seed {seed}: {len(metrics)} records", flush=True)


# ---------------------------------------------------------------------------
# Traces and truth
# ---------------------------------------------------------------------------

TRACES = [
    ("full", J29, 1, "OU3,A41,A52"),
    ("l1", J29, 1, "OU3,A41,A52"),
    ("full", "mixed_swell_wind", 1, "OU3,A41,A52"),
    ("full", f"{J29}+maneuver", 1, "OU3,A41,A52"),
    ("full", "stationary", 1, "OU3,A41,A52"),
    ("full", f"{J29}+rel90", 1, "OU3,A41,A52"),
]


def trace_path(level, sc, seed):
    return TRACE_DIR / f"{level}_{sc.replace('+', '_')}_s{seed}.csv"


def cmd_traces(args):
    require(BENCH)
    TRACE_DIR.mkdir(parents=True, exist_ok=True)
    for level, sc, seed, filters in TRACES:
        path = trace_path(level, sc, seed)
        subprocess.run([str(BENCH), "--level", level, "--scenario", sc, "--seed", str(seed), "--filters", filters,
                        "--timeseries", str(path), "--decimate", "4"], check=True, capture_output=True)
        print("wrote", path.relative_to(ROOT))


def cmd_truth(args):
    require(BENCH)
    scenarios = list(dict.fromkeys(PRIMARY + SECONDARY))
    rows = []
    for sc in scenarios:
        out = subprocess.run([str(BENCH), "--level", "truth", "--scenario", sc, "--seed", "1"],
                             capture_output=True, text=True, check=True).stdout
        kv = dict(t.split("=", 1) for t in out.split())
        rows.append(kv)
    with TRUTH_RAW.open("w") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    print("wrote", TRUTH_RAW.relative_to(ROOT))


def cmd_timing(args):
    """Serial runs on an idle host: fixed, +-1 % on a 1 us clock and +-25 %
    jittered sample interval."""
    if TIMING_RAW.exists():
        TIMING_RAW.unlink()
    cmds = [["--level", "full", "--scenario", sc, "--seed", str(seed), "--filters", "OU3,A41,A52,M45,M54,A52S"]
            for sc in (J29, f"{J29}+jit1", f"{J29}+jit") for seed in range(1, args.seeds + 1)]
    run_bench(cmds, TIMING_RAW, 1)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--primary-seeds", type=int, default=100)
    p.add_argument("--secondary-seeds", type=int, default=30)
    p.add_argument("--long-seeds", type=int, default=20)
    p.set_defaults(fn=cmd_run)
    p = sub.add_parser("dataset")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--seeds", type=int, default=10)
    p.add_argument("--variants", default="A41,A52,M45,M54,A52S,H41")
    p.set_defaults(fn=cmd_dataset)
    p = sub.add_parser("traces")
    p.set_defaults(fn=cmd_traces)
    p = sub.add_parser("truth")
    p.set_defaults(fn=cmd_truth)
    p = sub.add_parser("timing")
    p.add_argument("--seeds", type=int, default=3)
    p.set_defaults(fn=cmd_timing)
    p = sub.add_parser("report")
    p.set_defaults(fn=lambda a: __import__("carma3_report").report())
    args = ap.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
