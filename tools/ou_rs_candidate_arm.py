#!/usr/bin/env python3
"""Run one candidate arm of the full OU-III rS validation experiment."""

from __future__ import annotations

import argparse
from pathlib import Path

import ou_rs_candidate_full_validation_fixed as fixed

exp = fixed.exp

CANDIDATES = {
    "cubic_p3": (0.1548363522, 3.0),
    "fitted_p2p9052": (0.1667018769, 2.9052),
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", choices=sorted(CANDIDATES), required=True)
    parser.add_argument("--data-dir", type=Path, default=exp.REPO_ROOT / "plots/kalman_ou_ii")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--shards", type=int, default=2)
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    output_dir = args.output_dir.resolve()
    missing = [name for name in exp.SOURCE_FILES if not (data_dir / name).is_file()]
    if missing:
        raise FileNotFoundError("missing released simulation data: " + ", ".join(missing))
    if args.shards < 2:
        raise ValueError("use at least two shards so ou_validation stops after raw rows")

    coeff, exponent = CANDIDATES[args.label]
    original = exp.patch_experiment_branch()
    try:
        exp.run_arm(
            args.label,
            coeff,
            exponent,
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
