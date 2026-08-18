#!/usr/bin/env python3
"""Paired multi-seed comparison of the two OU-II pseudo-measurement laws.

doc/kalman_ou_ii/ou2-dual-regularization-mse.tex, Sec. (next-study), asks for
exactly one experiment: run the empirical pair

    r_p0 ~ sigma_aw tau^2,        r_v0 ~ sigma_aw tau

against the joint physical-MSE pair

    r_p0 ~ sigma_a,B^(4/5) tau^(12/5),  r_v0 ~ sigma_a,B^(4/5) tau^(7/5),

on paired multi-seed stationary seas, with the controlled transition reported
separately.  This driver is that experiment.

Both arms are the same binary -- the law is a runtime selection, not a compile
option -- so the comparison isolates the schedule and nothing else.  Each arm is
run through tools/ou_validation.py with identical arguments, so the two see
identical wave-phase, IMU-noise and initialization seeds on identical
scenarios; the per-(scenario, repetition) rows are then paired and differenced.

Typical use:

    python3 tools/ou2_pseudo_law_compare.py --jobs 4
    python3 tools/ou2_pseudo_law_compare.py --arm 'ratio-exact:0.1129,0.5376'

Each --arm adds a PhysicalMSE configuration at an explicit (C_P, C_P/C_V),
paired against the same single Empirical baseline, so alternative coefficient
choices are all measured on identical seeds.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import statistics as st
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATION = REPO_ROOT / "tools" / "ou_validation.py"

METRICS = ("disp_z_pct_hs", "disp_3d_rms_m", "roll_rms_deg", "pitch_rms_deg")


def run_arm(label: str, out_dir: Path, jobs: int, extra: list[str],
            env_extra: dict[str, str]) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(VALIDATION),
        "--mode", "full",
        "--families", "OU_II",
        "--adaptation-modes", "Adaptive",
        "--no-plots",
        "--skip-build",
        "--jobs", str(jobs),
        "--output-dir", str(out_dir),
        *extra,
    ]
    env = os.environ.copy()
    env.update(env_extra)
    detail = " ".join(f"{k}={v}" for k, v in sorted(env_extra.items()))
    print(f"[{label}] {detail} -> {out_dir}", flush=True)
    log = out_dir.with_suffix(".log")
    with log.open("w") as fh:
        proc = subprocess.run(cmd, env=env, stdout=fh,
                              stderr=subprocess.STDOUT, cwd=REPO_ROOT, check=False)
    raw = out_dir / "ou_validation_raw.csv"
    if not raw.exists():
        raise SystemExit(
            f"[{label}] produced no raw rows (exit {proc.returncode}); see {log}")
    # A nonzero exit only means the publication table could not be built; it
    # needs paired FixedNominal/FixedOracle arms this comparison does not run.
    return raw


def load(raw: Path) -> dict[tuple[str, int], dict[str, str]]:
    with raw.open() as fh:
        return {(r["scenario"], int(r["repetition"])): r for r in csv.DictReader(fh)}


def paired(base: dict, rev: dict, keys: list,
           metric: str) -> tuple[float, float, float, float]:
    """(baseline mean, revised mean, mean paired delta, 95% half-width)."""
    b = [float(base[k][metric]) for k in keys]
    r = [float(rev[k][metric]) for k in keys]
    d = [x - y for x, y in zip(r, b)]
    half = (1.96 * st.stdev(d) / math.sqrt(len(d))
            if len(d) > 1 and st.stdev(d) > 0 else 0.0)
    return st.mean(b), st.mean(r), st.mean(d), half


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--arm", action="append", default=[], metavar="NAME:C_P,RATIO",
        help="extra PhysicalMSE arm at an explicit (C_P, C_P/C_V), e.g. "
             "'ratio-exact:0.1129,0.5376'.  Repeatable.  Without any --arm the "
             "deployed compiled-in constants are used.")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--output-dir", type=Path,
                    default=REPO_ROOT / "reports" / "results" / "ou2_pseudo_law")
    ap.add_argument("--reuse", action="store_true",
                    help="reuse existing raw CSVs instead of re-running")
    ap.add_argument("rest", nargs=argparse.REMAINDER,
                    help="extra args forwarded to ou_validation.py")
    args = ap.parse_args(argv)
    extra = [a for a in args.rest if a != "--"]

    revised_specs: list[tuple[str, dict[str, str]]] = []
    for spec in args.arm:
        name, _, pair = spec.partition(":")
        fields = [f.strip() for f in pair.split(",")]
        env = {"OU_II_PSEUDO_LAW": "1", "OU_II_PSEUDO_MSE_COEFF": fields[0]}
        if len(fields) > 1 and fields[1]:
            env["OU_II_PSEUDO_MSE_RATIO"] = fields[1]
        revised_specs.append((name, env))
    if not revised_specs:
        revised_specs = [("physical-mse", {"OU_II_PSEUDO_LAW": "1"})]

    # The baseline is the empirical law, selected explicitly rather than by
    # being the compiled-in default, so this stays a comparison of two named
    # schedules however the default later moves.
    plan = [("empirical", {"OU_II_PSEUDO_LAW": "0"})]
    plan += revised_specs

    arms = {}
    for label, env_extra in plan:
        out = args.output_dir / f"validation_{label}"
        raw = out / "ou_validation_raw.csv"
        if not (args.reuse and raw.exists()):
            raw = run_arm(label, out, args.jobs, extra, env_extra)
        arms[label] = load(raw)

    rows = []
    for label, _env in revised_specs:
        keys = sorted(set(arms["empirical"]) & set(arms[label]))
        if not keys:
            raise SystemExit(f"[{label}] no paired (scenario, repetition) rows in common")
        missing = set(arms["empirical"]) ^ set(arms[label])
        if missing:
            print(f"warning: [{label}] {len(missing)} unpaired rows dropped", flush=True)
        scenarios = sorted({k[0] for k in keys})
        print(f"\n########## arm: {label} ##########")
        print(f"paired over {len(keys)} (scenario, repetition) rows, "
              f"{len(scenarios)} scenarios\n")
        for metric in METRICS:
            print(f"== {metric} ==")
            print(f"{'scenario':<44}{'empirical':>11}{'MSE':>10}{'delta':>10}{'95% hw':>9}")
            for scenario in scenarios:
                sk = [k for k in keys if k[0] == scenario]
                b, r, d, h = paired(arms["empirical"], arms[label], sk, metric)
                print(f"{scenario:<44}{b:11.4f}{r:10.4f}{d:+10.4f}{h:9.4f}")
                rows.append({"arm": label, "metric": metric, "scenario": scenario,
                             "n": len(sk), "empirical": b, "physical_mse": r,
                             "delta": d, "ci95_half": h})
            b, r, d, h = paired(arms["empirical"], arms[label], keys, metric)
            print(f"{'ALL (pooled)':<44}{b:11.4f}{r:10.4f}{d:+10.4f}{h:9.4f}\n")
            rows.append({"arm": label, "metric": metric, "scenario": "ALL",
                         "n": len(keys), "empirical": b, "physical_mse": r,
                         "delta": d, "ci95_half": h})

    out_csv = args.output_dir / "paired_comparison.csv"
    with out_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
