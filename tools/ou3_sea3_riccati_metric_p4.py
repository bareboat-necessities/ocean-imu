#!/usr/bin/env python3
"""Canonical P4 architecture over the SEA3 moving-Riccati P3 metric.

P4 keeps the validated 0.8-rad Cayley geometry and exact co-rotated a_w
accelerometer coordinate, but retires the 800-endpoint signed-Joseph scan as a
canonical theorem route.  Nonlinear dissipation is to be proved directly in the
same moving shipping-covariance metric produced by P3.

The intended word inequality is

    V(F(x), P_plus) - V(x, P) <= -mu ||x||_P^2 + R_3(x,xi),

where the signed Joseph vector information is retained until the recurrent word
is accumulated and R_3 is a validated higher-order finite-angle remainder.
P4 cannot promote until P3 has a positive Riccati margin and this remainder is
strictly dominated on the declared 0.8-rad sector.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_sea3_riccati_metric_p3 as P3
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_accelerometer_corotated_aw as COROT
import ou3_p4_vector_remainder_sector as REM

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_MOVING_RICCATI_NONLINEAR_P4"


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    p3 = P3.build(path)
    p3f = P3.validate(p3)
    cayley = CAYLEY.build(path)
    cf = CAYLEY.validate(cayley)
    corot = COROT.build(path)
    af = COROT.validate(corot)
    rem = REM.build(path)
    rf = REM.validate(rem)
    prereq_failures = [f"P3: {x}" for x in p3f] + [f"Cayley: {x}" for x in cf]
    prereq_failures += [f"corotated-aw: {x}" for x in af] + [f"remainder: {x}" for x in rf]
    if prereq_failures:
        raise RuntimeError(f"moving-Riccati P4 prerequisites failed: {prereq_failures}")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_P4_architecture": "NONLINEAR_WORD_IN_MOVING_SHIPPING_RICCATI_METRIC",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "SEA3_dynamic_source_used_through_P3": True,
        "old_800_endpoint_signed_Joseph_scan_consumed": False,
        "old_terminal_source_phase_metric_attachment_consumed": False,
        "outer_angle_rad": cayley["outer_angle_rad"],
        "cayley_geometry_validated": True,
        "accelerometer_corotated_aw_coordinate_used": True,
        "accelerometer_aw_nonlinear_eta_eliminated": (
            float(rem["acc_eta_aw_quadratic_coefficient_upper"]) == 0.0
        ),
        "signed_Joseph_directional_forms_retained_to_word_level": True,
        "nonlinear_word_inequality": (
            "V(F(x),P_plus)-V(x,P) <= -mu*||x||_P^2 + R3(x,xi)"
        ),
        "P3_CANONICAL_PASS_consumed": p3["P3_CANONICAL_PASS"],
        "nonlinear_remainder_dominated_on_full_sector": False,
        "P4_CANONICAL_PASS": False,
        "P5_MAY_START": False,
        "P4_CANONICAL_FAIL_REASONS": [
            "canonical moving-Riccati P3 margin is not yet closed",
            "complete H18/A21 recurrent-word nonlinear remainder domination is not yet emitted",
        ],
        "next_obligation": (
            "after P3 emits delta_H,delta_A, accumulate exact vector/S/accelerometer operations in that same metric "
            "and certify the finite-angle higher-order remainder on 0.8 rad; do not return to source-word enumeration"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_P4_architecture") != "NONLINEAR_WORD_IN_MOVING_SHIPPING_RICCATI_METRIC":
        f.append("wrong canonical P4 architecture")
    for key in (
        "source_generated_not_trajectory_fit", "SEA3_dynamic_source_used_through_P3",
        "cayley_geometry_validated", "accelerometer_corotated_aw_coordinate_used",
        "accelerometer_aw_nonlinear_eta_eliminated",
        "signed_Joseph_directional_forms_retained_to_word_level",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "old_800_endpoint_signed_Joseph_scan_consumed",
        "old_terminal_source_phase_metric_attachment_consumed",
        "P3_CANONICAL_PASS_consumed", "nonlinear_remainder_dominated_on_full_sector",
        "P4_CANONICAL_PASS", "P5_MAY_START",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if not d.get("P4_CANONICAL_FAIL_REASONS"):
        f.append("open P4 route does not name remaining obligations")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "architecture": d["canonical_P4_architecture"],
        "P4_CANONICAL_PASS": d["P4_CANONICAL_PASS"],
        "fail_reasons": d["P4_CANONICAL_FAIL_REASONS"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
