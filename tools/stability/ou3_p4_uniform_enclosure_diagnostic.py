#!/usr/bin/env python3
"""Non-promoting diagnostic for the rejected scalar P3-margin P4 route.

The earlier experiment attempted to pay every nonlinear and binary32 defect
from the source-uniform one-step Riccati margin.  Its H18 admissible radius fell
to about 2.465e-84 on CI.  This producer keeps the architecture rejected using
a stronger invariant argument: even the separately strengthened source-uniform
linear margin is orders of magnitude below binary32 unit roundoff, before any
nonlinear/reset amplification is charged.  Therefore a proof that treats
roundoff as a relative perturbation paid directly from that scalar margin cannot
close a nontrivial shipping tube.

The signed Joseph/reset information ledger remains the intended route.  The old
CI execution is retained only as a human-readable historical label so generated
certificate JSON contains no ephemeral run/workflow identifiers.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_strong_linear_margin as STRONG

UNIT_ROUNDOFF = 5.960464477539063e-08
FAILED_CI_REFERENCE = "historical scalar-route CI 34310904248"
FAILED_H18_RADIUS = 2.465e-84


def build() -> dict:
    strong = STRONG.build()
    failures = STRONG.validate(strong)
    if failures:
        raise RuntimeError(f"strong-margin prerequisite failed: {failures}")
    h = float(strong["H18"]["strong_relative_margin_lower"])
    a = float(strong["A21"]["strong_relative_margin_lower"])
    ratios = {
        "H18_unit_roundoff_over_linear_margin": math.nextafter(UNIT_ROUNDOFF / h, math.inf),
        "A21_unit_roundoff_over_linear_margin": math.nextafter(UNIT_ROUNDOFF / a, math.inf),
    }
    reject = min(ratios.values()) > 1.0e6
    return {
        "qualification": "OU3_P4_FAILED_SCALAR_UNIFORM_ENCLOSURE_DIAGNOSTIC_V2",
        "non_promoting": True,
        "architecture": "one-step P3 moving-Riccati scalar margin pays all nonlinear/reset/rounding defects",
        "strong_linear_margin": {
            "H18": h,
            "A21": a,
        },
        "binary32_unit_roundoff": UNIT_ROUNDOFF,
        "roundoff_to_margin_ratios": ratios,
        "previous_failed_CI_evidence": {
            "historical_reference": FAILED_CI_REFERENCE,
            "H18_allowable_Euclidean_radius_reported": FAILED_H18_RADIUS,
            "role": "historical diagnostic only; not a theorem constant",
        },
        "decision": "ABANDON_THIS_SCALAR_ARCHITECTURE" if reject else "REVIEW",
        "replacement_architecture": "joint signed Joseph/reset information ledger plus explicit finite-precision ISS channel",
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps(d, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
