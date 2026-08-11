#!/usr/bin/env python3
"""Amplitude-exponent ablation of the OU-III integral-regularizer schedule.

Section "Amplitude-exponent ablation of the regularizer schedule" of the
OU-III paper asks a single question: does the drift-band regularizer need its
factor of sigma_aw?  The reduced Riccati model, read with the accelerometer
sensor floor, says no --- the deployed filter is deep in the strong
acceleration-observation branch (zeta >> 1), where the pole-preserving
schedule is r_S ~ tau^(5/2) with no leading-order amplitude dependence.  Read
with a sea-state-proportional model-mismatch term instead, it says yes, and
predicts the deployed r_S ~ sigma_aw tau^(5/2).

The filter exposes exactly that one degree of freedom,

    r_S = sqrt(2 r_a) tau^3 / (sqrt(T_S) kappa^3) * (sigma_aw/sigma_ref)^p,

with kappa calibrated so sigma_ref is the sigma_aw of the nominal Hs = 1.5 m
sea.  The common tau^3/sqrt(T_S) factor cancels between family members, so
p = 1 reproduces the deployed schedule exactly, p = 0 is the strong-branch
asymptote, and every member agrees at sigma_ref for every tau.  The comparison
is therefore of the amplitude exponent alone and cannot be confounded by an
overall change of regularizer gain.

This driver runs tools/ou_validation.py once per exponent with the same
scenarios and seeds, then reports the paired differences against p = 1.

Typical use:

    python3 tools/ou_rs_law_ablation.py --exponents 0,0.5,1,1.25 --jobs 4

Add --posterior to include the full transition law
Eq. (riccati-pole-target-rs) as an extra column; it is expected to be
numerically indistinguishable from p = 0 wherever zeta >> 1.
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

REPO_ROOT = Path(__file__).resolve().parent.parent
VALIDATION = REPO_ROOT / "tools" / "ou_validation.py"

# Metrics reported by the ablation table.
METRICS = ("disp_z_pct_hs", "disp_3d_rms_m")

# Law selector understood by tests/kalman_ou_iii/kalman_ou_iii-sim.cpp.
LAW_CUBIC = "0"
LAW_STRONG = "1"
LAW_POSTERIOR = "2"


def run_one(label: str, env_extra: dict[str, str], out_dir: Path,
            jobs: int, extra_args: list[str]) -> Path:
    """Run the validation harness once and return its raw-row CSV."""
    out_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(env_extra)
    cmd = [
        sys.executable, str(VALIDATION),
        "--mode", "full",
        "--families", "OU_III",
        "--adaptation-modes", "Adaptive",
        "--no-plots",
        "--jobs", str(jobs),
        "--output-dir", str(out_dir),
        *extra_args,
    ]
    print(f"[{label}] {' '.join(f'{k}={v}' for k, v in sorted(env_extra.items()))}",
          flush=True)
    log = out_dir.with_suffix(".log")
    with log.open("w") as fh:
        proc = subprocess.run(cmd, env=env, stdout=fh, stderr=subprocess.STDOUT,
                              cwd=REPO_ROOT, check=False)
    raw = out_dir / "ou_validation_raw.csv"
    if not raw.exists():
        raise SystemExit(
            f"[{label}] produced no raw rows (exit {proc.returncode}); see {log}")
    # A nonzero exit only means the publication table could not be built,
    # which needs paired FixedNominal/FixedOracle modes this study does not
    # run.  The per-run rows this study consumes are already written.
    return raw


def load(raw: Path) -> dict[tuple[str, int], dict[str, str]]:
    with raw.open() as fh:
        return {(r["scenario"], int(r["repetition"])): r
                for r in csv.DictReader(fh)}


def paired_delta(runs: dict[str, dict], keys: list, label: str,
                 base: str, metric: str) -> tuple[float, float, float]:
    """Return (mean of label, mean paired delta vs base, 95% half-width)."""
    vals = [float(runs[label][k][metric]) for k in keys]
    diffs = [float(runs[label][k][metric]) - float(runs[base][k][metric])
             for k in keys]
    mean = st.mean(vals)
    dmean = st.mean(diffs)
    half = (1.96 * st.stdev(diffs) / math.sqrt(len(diffs))
            if len(diffs) > 1 and st.stdev(diffs) > 0 else 0.0)
    return mean, dmean, half


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exponents", default="0,0.5,1,1.25",
                    help="comma-separated amplitude exponents p (default 0,0.5,1,1.25)")
    ap.add_argument("--baseline", default="1",
                    help="exponent treated as the deployed baseline (default 1)")
    ap.add_argument("--posterior", action="store_true",
                    help="also run the full posterior transition law")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--output-dir", type=Path,
                    default=REPO_ROOT / "reports" / "results" / "ou_rs_law")
    ap.add_argument("--reuse", action="store_true",
                    help="reuse existing raw CSVs instead of re-running")
    ap.add_argument("rest", nargs=argparse.REMAINDER,
                    help="extra args forwarded to ou_validation.py")
    args = ap.parse_args(argv)

    extra = [a for a in args.rest if a != "--"]
    exponents = [e.strip() for e in args.exponents.split(",") if e.strip()]
    if args.baseline not in exponents:
        raise SystemExit(f"baseline p={args.baseline} must be among --exponents")

    runs: dict[str, dict] = {}
    order: list[str] = []
    for p in exponents:
        label = f"p={p}"
        order.append(label)
        out = args.output_dir / f"p_{p}"
        raw = out / "ou_validation_raw.csv"
        # p == 1 is reproduced exactly by the Riccati path, so the whole family
        # is run through the same code path for a like-for-like comparison.
        if not (args.reuse and raw.exists()):
            raw = run_one(label, {"OU_III_RS_LAW": LAW_STRONG,
                                  "OU_III_RS_SIGMA_EXP": p},
                          out, args.jobs, extra)
        runs[label] = load(raw)

    if args.posterior:
        label = "posterior"
        order.append(label)
        out = args.output_dir / "posterior"
        raw = out / "ou_validation_raw.csv"
        if not (args.reuse and raw.exists()):
            raw = run_one(label, {"OU_III_RS_LAW": LAW_POSTERIOR},
                          out, args.jobs, extra)
        runs[label] = load(raw)

    base = f"p={args.baseline}"
    keys = sorted(set.intersection(*(set(r) for r in runs.values())))
    if not keys:
        raise SystemExit("no (scenario, repetition) units are common to all runs")

    print(f"\nPaired amplitude-exponent ablation, n={len(keys)} units, "
          f"baseline {base}\n")
    width = max(len(x) for x in order) + 2
    for metric in METRICS:
        print(f"{metric}")
        print(f"  {'schedule':{width}} {'mean':>9} {'delta':>9} {'95% CI':>18}")
        for label in order:
            mean, dmean, half = paired_delta(runs, keys, label, base, metric)
            if label == base:
                print(f"  {label:{width}} {mean:9.4f} {'---':>9} {'':>18}")
            else:
                print(f"  {label:{width}} {mean:9.4f} {dmean:+9.4f}"
                      f"   [{dmean - half:+7.4f},{dmean + half:+7.4f}]")
        print()

    print("per-scenario disp_z_pct_hs")
    scenarios = sorted({k[0] for k in keys})
    print(f"  {'scenario':46} " + " ".join(f"{x:>12}" for x in order))
    for scenario in scenarios:
        sub = [k for k in keys if k[0] == scenario]
        cells = " ".join(
            f"{st.mean([float(runs[x][k]['disp_z_pct_hs']) for k in sub]):12.3f}"
            for x in order)
        print(f"  {scenario[:46]:46} {cells}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
