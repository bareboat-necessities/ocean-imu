#!/usr/bin/env python3
"""Anisotropy ablation of the OU-III stationary prior against its regularizer.

OU-III scales the stationary acceleration standard deviation horizontally by
S_factor = 1.87 and leaves the integral pseudo-measurement noise isotropic
(R_S_xy_factor = 1).  The similarity theorem the r_S schedule is derived from
gives sigma_S ~ sigma_aw tau^3 componentwise, so the axis-consistent choice
would be R_S_xy_factor = S_factor.  The deployed pair is not that point.

This driver sweeps both constants together over the four stationary JONSWAP
records and several IMU/initialization seeds, and reports paired per-seed
differences against the deployed pair.  Arms named "consistent" satisfy
rho_xy = S_sigma; the others isolate one knob at a time.

    python3 tools/ou_anisotropy_ablation.py --seeds 1,2,3,4,5

Results and reading: docs/ou-iii-anisotropy-consistency.md.  Note that
setRSXYFactor() clamped rho_xy to 1 before this study, so rho_xy > 1 arms run
against an older build silently reproduce the deployed arm.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import os
import re
import statistics as st
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SIM = REPO_ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim"

RECORDS = (
    "wave_data_jonswap_H0.270_L14.047_A30.00_P60.00.csv",
    "wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv",
    "wave_data_jonswap_H4.000_L112.766_A30.00_P30.00.csv",
    "wave_data_jonswap_H8.500_L202.839_A-30.00_P72.00.csv",
)

# label -> (S_factor, R_S_xy_factor).  A is the deployed pair and the baseline
# every other arm is paired against.
ARMS: dict[str, tuple[str, str]] = {
    "A_S1.87_rho1.00_deployed":   ("1.87", "1.00"),
    "B_S1.87_rho1.87_consistent": ("1.87", "1.87"),
    "C_S1.00_rho1.00_consistent": ("1.00", "1.00"),
    "D_S1.00_rho1.87_control":    ("1.00", "1.87"),
    "E_S0.80_rho1.00_inconsist":  ("0.80", "1.00"),
    "F_S0.80_rho0.80_consistent": ("0.80", "0.80"),
    "G_S1.87_rho0.53_antilaw":    ("1.87", "0.5348"),
}
BASE = "A_S1.87_rho1.00_deployed"

METRICS = ("disp_x_rms_m", "disp_y_rms_m", "disp_z_rms_m", "disp_3d_rms_m",
           "disp_z_pct_hs", "tau_applied_s", "sigma_applied_mps2")
REPORTED = METRICS[:4]


def run_one(job: tuple[str, str, int], window_s: int) -> dict:
    arm, record, seed = job
    s_factor, rho = ARMS[arm]
    env = os.environ.copy()
    env.update({
        "OU_III_S_FACTOR": s_factor,
        "OU_III_R_S_XY_FACTOR": rho,
        "W3D_SEED": str(seed),
        "W3D_WRITE_TIMESERIES": "0",          # keeps parallel arms from colliding
        "W3D_VALIDATION_WINDOW_SEC": str(window_s),
    })
    # A non-zero return code only means a deployed quality gate tripped, which
    # the off-default seeds are not gauged for.  The metrics line is emitted
    # either way, and the gate verdict is recorded alongside it.
    proc = subprocess.run([str(SIM), "--input", record], cwd=SIM.parent, env=env,
                          capture_output=True, text=True, check=False)
    try:
        line = next(l for l in proc.stdout.splitlines()
                    if l.startswith("VALIDATION_METRICS"))
    except StopIteration:
        raise SystemExit(f"no metrics from {arm} {record} seed={seed}:\n{proc.stdout[-2000:]}")
    row = {"arm": arm, "record": record.split("_")[3], "seed": seed,
           "gate": "PASS" if "QUALITY_GATE: PASS=1" in proc.stdout else "FAIL"}
    for m in METRICS:
        hit = re.search(rf"\b{m}=([-\d.eE+]+)", line)
        row[m] = float(hit.group(1)) if hit else float("nan")
    print(f"  {arm:29} {row['record']:7} seed={seed} "
          f"3D={row['disp_3d_rms_m']:.4f} gate={row['gate']}", flush=True)
    return row


def report(rows: list[dict]) -> None:
    by: dict[tuple[str, str], dict] = defaultdict(dict)
    for r in rows:
        by[(r["record"], str(r["seed"]))][r["arm"]] = {
            m: float(r[m]) for m in REPORTED}
    records = sorted({r["record"] for r in rows})
    arms = [a for a in ARMS if a != BASE and any(r["arm"] == a for r in rows)]

    def deltas(arm: str, metric: str, record: str | None = None) -> list[float]:
        return [100.0 * (v[arm][metric] - v[BASE][metric]) / v[BASE][metric]
                for k, v in by.items()
                if arm in v and BASE in v and (record is None or k[0] == record)]

    print("\nPaired per-seed change against the deployed pair, percent "
          "(negative is better).")
    print("Bracketed range is over seeds; * marks one sign on every seed.\n")
    head = f"{'record':8} {'arm':29} " + " ".join(
        f"{m.split('_')[1]:>22}" for m in REPORTED)
    print(head)
    print("-" * len(head))
    for record in records:
        for arm in arms:
            cells = []
            for m in REPORTED:
                d = deltas(arm, m, record)
                if not d:
                    cells.append(" " * 22)
                    continue
                star = "*" if all(x > 0 for x in d) or all(x < 0 for x in d) else " "
                cells.append(
                    f"{st.fmean(d):+7.2f} [{min(d):+6.1f},{max(d):+6.1f}]{star}".rjust(22))
            print(f"{record:8} {arm:29} " + " ".join(cells))
        print()

    print("Pooled over records:\n")
    for arm in arms:
        cells = [f"{m.split('_')[1]}={st.fmean(deltas(arm, m)):+7.2f}%"
                 for m in REPORTED]
        print(f"{arm:29} " + "  ".join(cells))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", default="1,2,3,4,5",
                    help="comma-separated W3D_SEED values (default 1,2,3,4,5)")
    ap.add_argument("--arms", default="",
                    help="comma-separated arm labels or prefixes; default all")
    ap.add_argument("--jobs", type=int, default=4, help="concurrent simulator runs")
    ap.add_argument("--window-sec", type=int, default=900,
                    help="trailing scoring window (default 900)")
    ap.add_argument("--output", type=Path,
                    default=REPO_ROOT / "reports" / "results" / "ou_anisotropy" /
                    "ou_anisotropy_raw.csv")
    args = ap.parse_args()

    if not SIM.is_file():
        raise SystemExit(f"{SIM} not built; run make -C tests/kalman_ou_iii build")
    missing = [r for r in RECORDS if not (SIM.parent / r).is_file()]
    if missing:
        raise SystemExit(f"missing records {missing}; run make fetch-sim-data")

    seeds = [int(s) for s in args.seeds.split(",") if s.strip()]
    selected = list(ARMS)
    if args.arms:
        wanted = [a.strip() for a in args.arms.split(",") if a.strip()]
        selected = [a for a in ARMS
                    if any(a == w or a.startswith(w) for w in wanted)]
        if BASE not in selected:
            selected.insert(0, BASE)     # every arm is reported against it

    jobs = [(a, r, s) for a in selected for r in RECORDS for s in seeds]
    print(f"{len(jobs)} runs: {len(selected)} arms x {len(RECORDS)} records "
          f"x {len(seeds)} seeds", flush=True)
    with cf.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows = list(pool.map(lambda j: run_one(j, args.window_sec), jobs))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["arm", "record", "seed", "gate", *METRICS])
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {args.output}")
    report(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
