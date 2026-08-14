#!/usr/bin/env python3
"""Sweep of the MEKF variances the OU families take through their constructors.

Kalman3D_Wave_OU_III and Kalman3D_Wave_OU_II are each handed seven variances at
construction -- the three sensor sigmas (sigma_a, sigma_g, sigma_m), Pq0, Pb0,
b0, and the family's pseudo-measurement regularizer -- and none had ever been
swept.  Four of the seven were not even reachable: both wrappers' begin()
called the three-argument initialize(), leaving initialize_ext() dead code.
They are Config fields now, and both simulators expose all seven:

    SF_SIGMA_A_SCALE      scale on the harness sigma_a (2.8x injected accel white)
    SF_SIGMA_G_SCALE      scale on the harness sigma_g (2.0x injected gyro white)
    SF_SIGMA_M_SCALE      scale on the harness sigma_m (1.2x injected mag white)
    SF_PQ0                initial attitude-error variance, rad^2
    SF_PB0                initial gyro-bias variance, (rad/s)^2
    SF_GYRO_BIAS_RW_VAR   gyro-bias random-walk variance density, (rad/s)^2/s
    SF_RS_NOISE_VAR       OU-III: initial integral pseudo-measurement variance
    SF_RP0_NOISE_VAR      OU-II: initial position pseudo-measurement variance
    SF_RV0_NOISE_VAR      OU-II: initial velocity pseudo-measurement variance

The default run is a coordinate sweep: one axis moved at a time about the
deployed point, every arm paired against the deployed arm on the same
(record, seed).  Joint arms for the refine rounds are given directly:

    python3 tools/ou_qmekf_variance_study.py --family ou_ii --seeds 1,2,3,4,5
    python3 tools/ou_qmekf_variance_study.py --arms \
        'g0.1:SF_SIGMA_G_SCALE=0.1' \
        'g0.1_m2:SF_SIGMA_G_SCALE=0.1,SF_SIGMA_M_SCALE=2'

Results and reading: docs/ou-iii-qmekf-variances.md and
docs/ou-ii-qmekf-variances.md.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import math
import os
import re
import statistics as st
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# Both OU families take the same seven variances and differ only in the last
# ones: OU-III regularizes the third integral with a single r_S, OU-II
# regularizes position and velocity with r_p0 and r_v0.  Everything else --
# the harness sigmas, the startup covariances, the gyro-bias random walk -- is
# shared, including the units mismatch on sigma_g.
FAMILIES = {
    "ou_iii": {
        "sim": REPO_ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim",
        "tail_axes": {
            "rs0": ("SF_RS_NOISE_VAR", (0.015, 0.15, 15.0, 150.0)),
        },
    },
    "ou_ii": {
        "sim": REPO_ROOT / "tests" / "kalman_ou_ii" / "kalman_ou_ii-sim",
        "tail_axes": {
            "rp0": ("SF_RP0_NOISE_VAR", (0.015, 0.15, 15.0, 150.0)),
            "rv0": ("SF_RV0_NOISE_VAR", (0.003, 0.03, 3.0, 30.0)),
        },
    },
}

RECORDS_BY_FAMILY = {
    "jonswap": (
        "wave_data_jonswap_H0.270_L14.047_A30.00_P60.00.csv",
        "wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv",
        "wave_data_jonswap_H4.000_L112.766_A30.00_P30.00.csv",
        "wave_data_jonswap_H8.500_L202.839_A-30.00_P72.00.csv",
    ),
    "pmstokes": (
        "wave_data_pmstokes_H0.270_L14.047_A30.00_P60.00.csv",
        "wave_data_pmstokes_H1.500_L50.710_A-30.00_P120.00.csv",
        "wave_data_pmstokes_H4.000_L112.766_A30.00_P30.00.csv",
        "wave_data_pmstokes_H8.500_L202.839_A-30.00_P72.00.csv",
    ),
}

BASE = "deployed"

# Coordinate sweep about the deployed point.  Every axis carries the deployed
# value implicitly through the shared baseline arm, so it is not repeated here.
# The family's own regularizer axes are appended by make_axes().
COMMON_AXES: dict[str, tuple[str, tuple[float, ...]]] = {
    "sga": ("SF_SIGMA_A_SCALE",    (0.36, 0.5, 0.71, 1.41, 2.0, 3.0)),
    "sgg": ("SF_SIGMA_G_SCALE",    (0.02, 0.035, 0.056, 0.1, 0.18, 0.32, 0.56, 2.0)),
    "sgm": ("SF_SIGMA_M_SCALE",    (0.25, 0.5, 2.0, 4.0, 8.0)),
    "pq0": ("SF_PQ0",              (5e-5, 5e-3, 5e-2)),
    "pb0": ("SF_PB0",              (1e-8, 1e-7, 1e-5, 1e-4)),
    "brw": ("SF_GYRO_BIAS_RW_VAR", (1e-13, 1e-12, 1e-10, 1e-9, 1e-8)),
}


def make_axes(family: str) -> dict[str, tuple[str, tuple[float, ...]]]:
    return {**COMMON_AXES, **FAMILIES[family]["tail_axes"]}

METRICS = ("roll_rms_deg", "pitch_rms_deg", "yaw_rms_deg",
           "disp_z_rms_m", "disp_3d_rms_m", "disp_z_pct_hs",
           "accel_bias_3d_rms_mps2", "gyro_bias_3d_rms_radps",
           "tau_applied_s", "sigma_applied_mps2")

# Attitude is what this study is for; displacement and the two bias channels
# are carried as guards, because a variance that buys tilt by unloading the
# error into the accelerometer bias has not bought anything.
REPORTED = ("roll_rms_deg", "pitch_rms_deg", "yaw_rms_deg",
            "disp_z_rms_m", "disp_3d_rms_m", "accel_bias_3d_rms_mps2")


def parse_arms(specs: list[str]) -> dict[str, dict[str, str]]:
    """'label:K=V,K=V' -> {label: {K: V}}.  A bare 'K=V,...' labels itself."""
    arms: dict[str, dict[str, str]] = {}
    for spec in specs:
        label, sep, body = spec.partition(":")
        if not sep:
            label, body = spec, spec
        env = {}
        for item in body.split(","):
            item = item.strip()
            if not item:
                continue
            key, eq, value = item.partition("=")
            if not eq:
                raise SystemExit(f"bad arm term {item!r} in {spec!r}; want KEY=VALUE")
            env[key.strip()] = value.strip()
        if not env:
            raise SystemExit(f"arm {spec!r} sets nothing")
        arms[label] = env
    return arms


def axis_arms(selected: list[str], axes: dict) -> dict[str, dict[str, str]]:
    arms: dict[str, dict[str, str]] = {}
    for key in selected:
        if key not in axes:
            raise SystemExit(f"unknown axis {key!r}; pick from {sorted(axes)}")
        env_name, values = axes[key]
        for v in values:
            arms[f"{key}_{v:g}"] = {env_name: f"{v:g}"}
    return arms


def run_one(job: tuple[str, str, int], arms: dict[str, dict[str, str]],
            base_env: dict[str, str], window_s: int, sim: Path) -> dict:
    arm, record, seed = job
    env = os.environ.copy()
    env.update(base_env)
    env.update(arms[arm])
    env.update({
        "W3D_SEED": str(seed),
        "W3D_WRITE_TIMESERIES": "0",   # keeps parallel arms from colliding
        "W3D_VALIDATION_WINDOW_SEC": str(window_s),
        # Without this the first tripped sentinel exits the simulator before
        # the QUALITY_GATE line is printed, so every arm would be recorded as
        # FAIL for want of a verdict rather than for failing one.
        "W3D_COLLECT_ALL_GATES": "1",
    })
    # A non-zero return code only means a deployed quality gate tripped, which
    # the off-default seeds and off-default variances are not gauged for.  The
    # metrics line is emitted either way and the verdict is recorded with it.
    proc = subprocess.run([str(sim), "--input", record], cwd=sim.parent, env=env,
                          capture_output=True, text=True, check=False)
    try:
        line = next(l for l in proc.stdout.splitlines()
                    if l.startswith("VALIDATION_METRICS"))
    except StopIteration:
        raise SystemExit(f"no metrics from {arm} {record} seed={seed}:\n"
                         f"{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")
    parts = record.split("_")
    row = {"arm": arm, "record": f"{parts[2]}_{parts[3]}", "seed": seed,
           "gate": "PASS" if "QUALITY_GATE: PASS=1" in proc.stdout else "FAIL"}
    for m in METRICS:
        hit = re.search(rf"\b{m}=([-\d.eE+]+)", line)
        row[m] = float(hit.group(1)) if hit else float("nan")
    print(f"  {arm:24} {row['record']:16} seed={seed} "
          f"roll={row['roll_rms_deg']:.4f} pitch={row['pitch_rms_deg']:.4f} "
          f"yaw={row['yaw_rms_deg']:.4f} 3D={row['disp_3d_rms_m']:.4f}", flush=True)
    return row


def short(metric: str) -> str:
    if metric.startswith("disp_"):
        return metric.split("_")[1]
    if metric.startswith("accel_bias"):
        return "accbias"
    if metric.startswith("gyro_bias"):
        return "gyrbias"
    return metric.split("_")[0]


def report(rows: list[dict], arms: dict[str, dict[str, str]]) -> None:
    by: dict[tuple[str, int], dict[str, dict[str, float]]] = defaultdict(dict)
    for r in rows:
        by[(r["record"], int(r["seed"]))][r["arm"]] = {
            m: float(r[m]) for m in REPORTED}
    present = {r["arm"] for r in rows}
    others = [a for a in arms if a != BASE and a in present]
    records = sorted({r["record"] for r in rows})

    def deltas(arm: str, metric: str, record: str | None = None) -> list[float]:
        out = []
        for key, v in by.items():
            if arm not in v or BASE not in v:
                continue
            if record is not None and key[0] != record:
                continue
            b = v[BASE][metric]
            if not (math.isfinite(b) and b > 0.0 and math.isfinite(v[arm][metric])):
                continue
            out.append(100.0 * (v[arm][metric] - b) / b)
        return out

    print("\nPaired per-(record, seed) change against the deployed variances, "
          "percent (negative is better).")
    print("Bracketed range is over records and seeds; * marks one sign on "
          "every pair.\n")
    head = f"{'arm':24} " + " ".join(f"{short(m):>23}" for m in REPORTED)
    print(head)
    print("-" * len(head))
    for arm in others:
        cells = []
        for m in REPORTED:
            d = deltas(arm, m)
            if not d:
                cells.append(" " * 23)
                continue
            star = "*" if all(x > 0 for x in d) or all(x < 0 for x in d) else " "
            cells.append(
                f"{st.fmean(d):+7.2f} [{min(d):+6.1f},{max(d):+6.1f}]{star}".rjust(23))
        print(f"{arm:24} " + " ".join(cells))

    print("\nPer-record pooled means, percent (negative is better).\n")
    head = f"{'arm':24} {'record':16} " + " ".join(f"{short(m):>9}" for m in REPORTED)
    print(head)
    print("-" * len(head))
    for arm in others:
        for record in records:
            cells = []
            for m in REPORTED:
                d = deltas(arm, m, record)
                cells.append(f"{st.fmean(d):+9.2f}" if d else " " * 9)
            print(f"{arm:24} {record:16} " + " ".join(cells))
        print()

    fails = sorted({(r["arm"], r["record"]) for r in rows if r["gate"] == "FAIL"})
    if fails:
        print(f"deterministic quality gate FAIL on {len(fails)} (arm, record) "
              f"pairs, off-default seeds included:")
        for arm, record in fails:
            print(f"  {arm:24} {record}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family", default="ou_iii", choices=sorted(FAMILIES),
                    help="which simulator to sweep (default ou_iii)")
    ap.add_argument("--seeds", default="1,2,3,4,5",
                    help="comma-separated W3D_SEED values (default 1,2,3,4,5)")
    ap.add_argument("--axes", default="",
                    help="comma-separated axes to sweep; default all of the "
                         "family's.  Ignored when --arms is given.")
    ap.add_argument("--arms", nargs="*", default=[],
                    help="explicit arms as 'label:KEY=VALUE,KEY=VALUE'; "
                         "replaces the coordinate sweep")
    ap.add_argument("--base-env", default="",
                    help="KEY=VALUE,... applied to every arm including the "
                         "baseline, to re-centre a refine round")
    ap.add_argument("--families", default="jonswap,pmstokes",
                    help="comma-separated record families (default both)")
    ap.add_argument("--jobs", type=int, default=4, help="concurrent simulator runs")
    ap.add_argument("--window-sec", type=int, default=900,
                    help="trailing scoring window (default 900)")
    ap.add_argument("--output", type=Path,
                    default=None,
                    help="raw CSV destination; defaults to a family-named path "
                         "under reports/results/ou_qmekf_variances/")
    args = ap.parse_args()

    sim = FAMILIES[args.family]["sim"]
    axes = make_axes(args.family)
    if args.output is None:
        args.output = (REPO_ROOT / "reports" / "results" / "ou_qmekf_variances" /
                       f"{args.family}_raw.csv")
    if not sim.is_file():
        raise SystemExit(f"{sim} not built; run make -C {sim.parent} build")

    records: list[str] = []
    for family in (f.strip() for f in args.families.split(",") if f.strip()):
        if family not in RECORDS_BY_FAMILY:
            raise SystemExit(f"unknown family {family!r}; "
                             f"pick from {sorted(RECORDS_BY_FAMILY)}")
        records.extend(RECORDS_BY_FAMILY[family])
    missing = [r for r in records if not (sim.parent / r).is_file()]
    if missing:
        raise SystemExit(f"missing records {missing}; run make fetch-sim-data")

    base_env = {}
    if args.base_env:
        base_env = next(iter(parse_arms([f"base:{args.base_env}"]).values()))

    if args.arms:
        arms = parse_arms(args.arms)
    else:
        selected = ([a.strip() for a in args.axes.split(",") if a.strip()]
                    or list(axes))
        arms = axis_arms(selected, axes)
    if BASE in arms:
        raise SystemExit(f"{BASE!r} is the reserved baseline label")
    arms = {BASE: {}, **arms}

    seeds = [int(s) for s in args.seeds.split(",") if s.strip()]
    jobs = [(a, r, s) for a in arms for r in records for s in seeds]
    print(f"{args.family}: {len(jobs)} runs, {len(arms)} arms x "
          f"{len(records)} records x {len(seeds)} seeds", flush=True)
    if base_env:
        print(f"baseline re-centred on {base_env}", flush=True)
    with cf.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows = list(pool.map(
            lambda j: run_one(j, arms, base_env, args.window_sec, sim), jobs))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["arm", "record", "seed", "gate", *METRICS])
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {args.output}")
    report(rows, arms)
    return 0


if __name__ == "__main__":
    sys.exit(main())
