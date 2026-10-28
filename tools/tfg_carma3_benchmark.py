#!/usr/bin/env python3
"""Two-frame Lie-group CARMA filter benchmark (EXPERIMENTAL; docs/tfg-carma3-research.md).

Runs tests/kalman_carma/carma3_bench --extended over a fixed grid of
(scenario, seed).  Every filter of a row is fed the same truth and the same
sensor-noise realization, so every comparison is paired.  The primary pair is
CARMA3 (quaternion MEKF, A41) against TFG-CARMA3 (two-frame Lie group, T41),
both behind the same wrapper code; OU-III (OU3) is the shipping baseline;
T41R0 / T41RA / T41NW are the reset and covariance-write ablations.  The
CARMA(5,2) pair (A52, T52) is a separate follow-up grid.  Every subcommand is
resumable: rows already in the raw CSV are skipped.

    run       synthetic grid -> reports/results/tfg_carma3/raw_runs.csv
    followup  CARMA(5,2) grid -> followup_runs.csv
    mechanism post-hoc test of the mixed-sea result -> mechanism_runs.csv
    truth     realized motion levels of every condition -> truth_levels.csv
    dataset   pinned 28 ft vessel-RAO records through carma3-sim (A41, T41,
              T41R0) and kalman_ou_iii-sim (OU3), 10 seeds -> dataset_runs.csv
    traces    time series of selected runs for the figures
    timing    host update time on an otherwise idle machine (one job)
    report    tables, paired statistics and figures (tools/tfg_carma3_report.py)

Run: python3 tools/tfg_carma3_benchmark.py {run,followup,mechanism,truth,dataset,traces,timing,report} [--jobs 4]
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(__file__).resolve().parents[1]
TDIR = ROOT / "tests" / "kalman_carma"
BENCH = TDIR / "carma3_bench"
CARMA3_SIM = TDIR / "carma3-sim"
OU3_SIM = ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim"
OUT = ROOT / "reports" / "results" / "tfg_carma3"
RAW = OUT / "raw_runs.csv"
FOLLOWUP_RAW = OUT / "followup_runs.csv"
MECHANISM_RAW = OUT / "mechanism_runs.csv"
TRUTH_RAW = OUT / "truth_levels.csv"
DATASET_RAW = OUT / "dataset_runs.csv"
TIMING_RAW = OUT / "timing_runs.csv"
TRACE_DIR = OUT / "traces"

FILTERS = "OU3,A41,T41,T41R0,T41RA,T41NW"
LONG_FILTERS = "OU3,A41,T41,T41R0"
FOLLOWUP_FILTERS = "A41,T41,A52,T52,T52R0"
J29 = "jonswap_fp0.29_hs0.8_g3.3"

# Principal seas: 100 paired seeds.
PRIMARY = [J29, "pm_fp0.25_hs1", "mixed_swell_wind", "changing_fp"]
# Scenario matrix: 30 paired seeds each.  Headings are the wave propagation
# direction relative to the bow (0 following, 90 beam, 180 head).
SECONDARY = [
    "broadband", "bimodal", "abrupt_change", "jonswap_fp0.20_hs1_g3.3", "pm_fp0.15_hs2.5",
    f"{J29}+rel0", f"{J29}+rel45", f"{J29}+rel90", f"{J29}+rel135", f"{J29}+rel180",
    f"{J29}+anchor", f"{J29}+drift", f"{J29}+rel180+U2.5", f"{J29}+rel0+U2.5", f"{J29}+maneuver",
    "stationary", "stationary+calbias", "small_amplitude", "jonswap_fp0.40_hs0.02_g3.3",
    "large_rotation", "tilted", f"{J29}+turns",
    f"{J29}+maggap", f"{J29}+T1800+maggap900",
    f"{J29}+gbdrift", "bias_drift", f"{J29}+temp",
    f"{J29}+jit", f"{J29}+jit1",
]
# 3600 s runs; the last loses the magnetometer from 1200 s to the end (2400 s
# without heading information).
LONG = [f"{J29}+T3600", "mixed_swell_wind+T3600", f"{J29}+T3600+maggap2400"]

# Post hoc (designed after the principal-sea result; docs/tfg-carma3-research.md
# section 8): the mixed sea, where the pair separates, with the magnetometer
# removed for 300 s from T/3, so the heading information each filter gains
# from the accelerometer alone is measured directly; and JONSWAP with a
# calibrated accelerometer (+-1 mg turn-on bias instead of up to 8 mg), the
# test of whether the common attitude overconfidence is the bias prior.
MECHANISM = ["mixed_swell_wind+maggap", f"{J29}+calbias"]
MECHANISM_FILTERS = "OU3,A41,T41,T41R0,T41RA"

TRACES = [
    (J29, 1), ("mixed_swell_wind", 1), ("stationary", 1), (f"{J29}+maggap", 1), (f"{J29}+T3600+maggap2400", 1),
    (f"{J29}+T1800+maggap900", 1), (f"{J29}+turns", 1), ("large_rotation", 1), ("changing_fp", 1),
    ("mixed_swell_wind+maggap", 1),
]


def require(path):
    if not path.exists():
        sys.exit(f"missing {path.relative_to(ROOT)}: build it first (make -C {path.parent.relative_to(ROOT)} build)")


def read_rows(path):
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f))


def header():
    out = subprocess.run([str(BENCH), "--header", "--extended"], capture_output=True, text=True, check=True).stdout
    return out.splitlines()[0].strip()


def run_bench(cmds, out_path, jobs):
    """Run carma3_bench for every argument list, appending rows to out_path."""
    require(BENCH)
    hdr = header()
    new_file = not out_path.exists()
    if not new_file:
        with out_path.open() as f:
            if f.readline().strip() != hdr:
                sys.exit(f"{out_path.relative_to(ROOT)}: header differs from the current carma3_bench; "
                         "move it aside before resuming")
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


def pending(work, raw):
    done = {(r["scenario"], r["seed"]) for r in read_rows(raw)}
    return [w for w in work if (w[0], str(w[1])) not in done]


def cmds_for(work):
    return [["--extended", "--level", "full", "--scenario", sc, "--seed", str(seed), "--filters", flt]
            for sc, seed, flt in work]


def cmd_run(args):
    work = []
    for sc in PRIMARY:
        work += [(sc, s, FILTERS) for s in range(1, args.primary_seeds + 1)]
    for sc in SECONDARY:
        work += [(sc, s, FILTERS) for s in range(1, args.secondary_seeds + 1)]
    for sc in LONG:
        work += [(sc, s, LONG_FILTERS) for s in range(1, args.long_seeds + 1)]
    work = pending(work, RAW)
    print(f"{len(work)} runs to do", flush=True)
    run_bench(cmds_for(work), RAW, args.jobs)


def cmd_followup(args):
    work = [(sc, s, FOLLOWUP_FILTERS) for sc in PRIMARY for s in range(1, args.secondary_seeds + 1)]
    work += [(sc, s, FOLLOWUP_FILTERS) for sc in ("stationary", f"{J29}+maggap", "bimodal")
             for s in range(1, args.secondary_seeds + 1)]
    work = pending(work, FOLLOWUP_RAW)
    print(f"{len(work)} follow-up runs to do", flush=True)
    run_bench(cmds_for(work), FOLLOWUP_RAW, args.jobs)


def cmd_mechanism(args):
    work = [(sc, s, MECHANISM_FILTERS) for sc in MECHANISM for s in range(1, args.secondary_seeds + 1)]
    work = pending(work, MECHANISM_RAW)
    print(f"{len(work)} mechanism runs to do", flush=True)
    run_bench(cmds_for(work), MECHANISM_RAW, args.jobs)


def cmd_truth(args):
    """Realized wave motion of every condition (no filter), three seeds each."""
    require(BENCH)
    rows = []
    for sc in dict.fromkeys(PRIMARY + SECONDARY + MECHANISM):
        for seed in (1, 2, 3):
            out = subprocess.run([str(BENCH), "--level", "truth", "--scenario", sc, "--seed", str(seed)],
                                 capture_output=True, text=True, check=True).stdout
            rows.append(dict(t.split("=", 1) for t in out.split()))
    TRUTH_RAW.parent.mkdir(parents=True, exist_ok=True)
    with TRUTH_RAW.open("w") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    print(f"truth: {len(rows)} rows", flush=True)


def cmd_dataset(args):
    import carma3_benchmark as cb
    for s in (OU3_SIM, CARMA3_SIM):
        require(s)
    done = {(r["family"], r["seed"]) for r in read_rows(DATASET_RAW)}
    work = []
    for seed in range(1, args.seeds + 1):
        work.append(("OU3", seed, [str(OU3_SIM)], {}))
        for v in ("A41", "T41", "T41R0"):
            work.append((v, seed, [str(CARMA3_SIM)], {"W3D_CARMA3_VARIANT": v}))
    work = [w for w in work if (w[0], str(w[1])) not in done]
    print(f"{len(work)} dataset invocations to do", flush=True)
    new_file = not DATASET_RAW.exists()
    DATASET_RAW.parent.mkdir(parents=True, exist_ok=True)
    with DATASET_RAW.open("a") as out:
        wr = csv.DictWriter(out, fieldnames=cb.DATASET_FIELDS, extrasaction="ignore")
        if new_file:
            wr.writeheader()
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(cb.dataset_invocation, cmd, extra, seed): (fam, seed) for fam, seed, cmd, extra in work}
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


def cmd_traces(args):
    require(BENCH)
    TRACE_DIR.mkdir(parents=True, exist_ok=True)
    cmds = []
    for sc, seed in TRACES:
        name = sc.replace("+", "_")
        ts = TRACE_DIR / f"{name}_s{seed}.csv"
        if ts.exists():
            continue
        cmds.append([str(BENCH), "--header", "--extended", "--level", "full", "--scenario", sc, "--seed", str(seed),
                     "--filters", "OU3,A41,T41,T41R0", "--timeseries", str(ts), "--decimate", "10"])
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for p in ex.map(lambda c: subprocess.run(c, capture_output=True, text=True), cmds):
            if p.returncode != 0:
                print(f"ERROR trace: {p.stderr[-300:]}", flush=True)
    print(f"{len(cmds)} traces written", flush=True)


def cmd_timing(args):
    """Host update time, one process at a time on an otherwise idle machine."""
    require(BENCH)
    rows = []
    for mod in ("", "+jit1", "+jit"):
        for seed in (1, 2, 3):
            p = subprocess.run([str(BENCH), "--level", "full", "--scenario", f"{J29}{mod}", "--seed", str(seed),
                                "--filters", "OU3,A41,T41,T41R0,T41RA"], capture_output=True, text=True, check=True)
            hdr = header().split(",")
            for line in p.stdout.strip().splitlines():
                r = dict(zip(hdr, line.split(",")))
                rows.append({"clock": mod.lstrip("+") or "fixed", "seed": seed, "filter": r["filter"],
                             "runtime_ns": r["runtime_ns"], "memory_bytes": r["memory_bytes"]})
    TIMING_RAW.parent.mkdir(parents=True, exist_ok=True)
    with TIMING_RAW.open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"timing: {len(rows)} rows", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "followup", "mechanism", "truth", "dataset", "traces", "timing",
                                        "report"])
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--primary-seeds", type=int, default=100)
    ap.add_argument("--secondary-seeds", type=int, default=30)
    ap.add_argument("--long-seeds", type=int, default=20)
    ap.add_argument("--seeds", type=int, default=10, help="dataset seeds")
    args = ap.parse_args()
    if args.command == "report":
        import tfg_carma3_report
        tfg_carma3_report.main()
        return
    {"run": cmd_run, "followup": cmd_followup, "mechanism": cmd_mechanism, "truth": cmd_truth,
     "dataset": cmd_dataset, "traces": cmd_traces, "timing": cmd_timing}[args.command](args)


if __name__ == "__main__":
    main()
