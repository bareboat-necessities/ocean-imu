#!/usr/bin/env python3
"""Run one full-suite cubic effective-rS normalization arm.

This is an experiment-only launcher.  It reuses the established OU-III full
validation harness and keeps the effective law exponent fixed at p_tau=3 and
p_sigma=1.  Only the normalization K in

    rS_eff = K * sigma_aw * tau^3

changes between arms.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import ou_rs_candidate_full_validation_fixed as fixed

exp = fixed.exp

# Coefficients are anchored to the independently measured nominal point
# tau=2.17904091 s, sigma_aw=0.724445343 m/s^2.  Labels encode the nominal
# effective rS target in m*s.
CANDIDATES = {
    "r1161": (0.1548363522, 1.160579),
    "r1300": (0.1734369292, 1.300000),
    "r1400": (0.1867782315, 1.400000),
    "r1500": (0.2001195337, 1.500000),
    "r1600": (0.2134608359, 1.600000),
    "r1700": (0.2268021382, 1.700000),
    "r1860": (0.2481482218, 1.860000),
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", choices=sorted(CANDIDATES), required=True)
    parser.add_argument(
        "--data-dir", type=Path, default=exp.REPO_ROOT / "plots/kalman_ou_ii"
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--shards", type=int, default=2)
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    output_dir = args.output_dir.resolve()
    missing = [name for name in exp.SOURCE_FILES if not (data_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(
            "missing released simulation data: " + ", ".join(missing)
        )
    if args.shards < 2:
        raise ValueError("use at least two shards so ou_validation stops after raw rows")

    coeff, nominal_r = CANDIDATES[args.label]
    print(
        f"CUBIC_NORMALIZATION_ARM label={args.label} K={coeff:.10f} "
        f"nominal_rS_eff={nominal_r:.6f} exponent=3",
        flush=True,
    )

    original = exp.patch_experiment_branch()
    try:
        exp.run_arm(
            args.label,
            coeff,
            3.0,
            data_dir,
            output_dir,
            max(1, args.jobs),
            args.shards,
        )
    finally:
        exp.HEADER.write_text(original, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
