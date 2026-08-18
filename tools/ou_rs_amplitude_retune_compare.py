#!/usr/bin/env python3
"""Paired multi-seed comparison of two OU-III r_S base schedules.

Screening for the (c_sigma, c_R) retune is done on the deterministic
eight-record protocol (tools/ou_rs_amplitude_retune_sweep.py).  The verdict is
not: it is taken on the same paired Monte Carlo harness the paper's primary
endpoint uses, so that the two schedules see identical wave-phase, IMU-noise
and initialization seeds on identical scenarios.

The two arms are two simulator binaries, because the schedule itself differs
in the header rather than in an environment knob:

    --baseline-binary   built from main:  r_S,base = c_R sigma_aw tau^3
    --revised-binary    built here:       r_S,base = c_R tau^3

Each arm is run through tools/ou_validation.py with the same arguments, then
the per-(scenario, repetition) rows are paired and differenced.

Typical use:

    python3 tools/ou_rs_amplitude_retune_compare.py \\
        --baseline-binary tests/kalman_ou_iii/kalman_ou_iii-sim-main \\
        --revised-binary  tests/kalman_ou_iii/kalman_ou_iii-sim-revised \\
        --arm swept:0.32,14.4 --arm theory:1.80,17.11 --jobs 4

Each --arm is a separate revised configuration paired against the same single
baseline run, so the analytical reference points and the sweep optimum are all
measured on identical seeds.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import statistics as st
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATION = REPO_ROOT / "tools" / "ou_validation.py"
LIVE_BINARY = REPO_ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim"

METRICS = ("disp_z_pct_hs", "disp_3d_rms_m", "roll_rms_deg", "pitch_rms_deg")


def run_arm(label: str, binary: Path, out_dir: Path, jobs: int,
            extra: list[str], env_extra: dict[str, str] | None = None) -> Path:
    """Swap `binary` into place and run the validation harness once."""
    out_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(binary, LIVE_BINARY)
    cmd = [
        sys.executable, str(VALIDATION),
        "--mode", "full",
        "--families", "OU_III",
        "--adaptation-modes", "Adaptive",
        "--no-plots",
        "--skip-build",
        "--jobs", str(jobs),
        "--output-dir", str(out_dir),
        *extra,
    ]
    env = os.environ.copy()
    # Coefficient overrides apply to this arm only, so the baseline arm keeps
    # main's own law and coefficients rather than inheriting the revised ones.
    env.update(env_extra or {})
    detail = " ".join(f"{k}={v}" for k, v in sorted((env_extra or {}).items()))
    print(f"[{label}] {binary.name} {detail} -> {out_dir}", flush=True)
    log = out_dir.with_suffix(".log")
    with log.open("w") as fh:
        proc = subprocess.run(cmd, env=env, stdout=fh,
                              stderr=subprocess.STDOUT, cwd=REPO_ROOT, check=False)
    raw = out_dir / "ou_validation_raw.csv"
    if not raw.exists():
        raise SystemExit(f"[{label}] produced no raw rows (exit {proc.returncode}); see {log}")
    # A nonzero exit only means the publication table could not be built; it
    # needs paired FixedNominal/FixedOracle arms this comparison does not run.
    return raw


def load(raw: Path) -> dict[tuple[str, int], dict[str, str]]:
    with raw.open() as fh:
        return {(r["scenario"], int(r["repetition"])): r for r in csv.DictReader(fh)}


def paired(base: dict, rev: dict, keys: list, metric: str) -> tuple[float, float, float, float]:
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
    ap.add_argument("--baseline-binary", type=Path, required=True)
    ap.add_argument("--revised-binary", type=Path, required=True)
    ap.add_argument(
        "--arm", action="append", default=[], metavar="NAME:c_sigma,C_R[,C_J]",
        help="extra revised arm at an explicit coefficient pair, e.g. "
             "'theory:1.80,17.11'.  An optional third field selects the "
             "SpectralMSE law at that C_J, e.g. 'mse:0.9,17.112,0.0538'.  "
             "Repeatable.  Without any --arm the "
             "revised binary is run at its compiled-in defaults.")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--output-dir", type=Path,
                    default=REPO_ROOT / "reports" / "results" / "ou_rs_amplitude_retune")
    ap.add_argument("--reuse", action="store_true",
                    help="reuse existing raw CSVs instead of re-running")
    ap.add_argument("rest", nargs=argparse.REMAINDER,
                    help="extra args forwarded to ou_validation.py")
    args = ap.parse_args(argv)
    extra = [a for a in args.rest if a != "--"]

    # Revised arms: the compiled-in defaults unless explicit pairs are given.
    revised_specs: list[tuple[str, dict[str, str]]] = []
    for spec in args.arm:
        name, _, pair = spec.partition(":")
        fields = [f.strip() for f in pair.split(",")]
        env = {"OU_SIGMA_COEFF": fields[0], "OU_III_R_S_COEFF": fields[1]}
        # Optional third field: C_J, which also selects the SpectralMSE law.
        if len(fields) > 2 and fields[2]:
            env["OU_III_RS_LAW"] = "3"
            env["OU_III_RS_MSE_COEFF"] = fields[2]
        revised_specs.append((name, env))
    if not revised_specs:
        revised_specs = [("revised", {})]

    arms = {}
    plan = [("baseline", args.baseline_binary, {})]
    plan += [(name, args.revised_binary, env) for name, env in revised_specs]
    for label, binary, env_extra in plan:
        out = args.output_dir / f"validation_{label}"
        raw = out / "ou_validation_raw.csv"
        if not (args.reuse and raw.exists()):
            raw = run_arm(label, binary.resolve(), out, args.jobs, extra, env_extra)
        arms[label] = load(raw)

    rows = []
    for label, _env in revised_specs:
        keys = sorted(set(arms["baseline"]) & set(arms[label]))
        if not keys:
            raise SystemExit(f"[{label}] no paired (scenario, repetition) rows in common")
        missing = set(arms["baseline"]) ^ set(arms[label])
        if missing:
            print(f"warning: [{label}] {len(missing)} unpaired rows dropped", flush=True)
        scenarios = sorted({k[0] for k in keys})
        print(f"\n########## arm: {label} ##########")
        print(f"paired over {len(keys)} (scenario, repetition) rows, "
              f"{len(scenarios)} scenarios\n")
        for metric in METRICS:
            print(f"== {metric} ==")
            print(f"{'scenario':<44}{'baseline':>10}{'revised':>10}{'delta':>10}{'95% hw':>9}")
            for scenario in scenarios:
                sk = [k for k in keys if k[0] == scenario]
                b, r, d, h = paired(arms["baseline"], arms[label], sk, metric)
                print(f"{scenario:<44}{b:10.4f}{r:10.4f}{d:+10.4f}{h:9.4f}")
                rows.append({"arm": label, "metric": metric, "scenario": scenario,
                             "n": len(sk), "baseline": b, "revised": r,
                             "delta": d, "ci95_half": h})
            b, r, d, h = paired(arms["baseline"], arms[label], keys, metric)
            print(f"{'ALL (pooled)':<44}{b:10.4f}{r:10.4f}{d:+10.4f}{h:9.4f}\n")
            rows.append({"arm": label, "metric": metric, "scenario": "ALL",
                         "n": len(keys), "baseline": b, "revised": r,
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
