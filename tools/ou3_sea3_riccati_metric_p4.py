#!/usr/bin/env python3
"""Canonical P4 architecture over the SEA3 moving-Riccati P3 metric.

P4 is intentionally detached from the retired 800-endpoint P3->P4 metric
attachment and signed-Joseph endpoint scan.  The only retained nonlinear
primitive consumed at this architecture boundary is the validated 0.8-rad
Cayley sector, which is independent of the old source-word metric machinery.

After P3 emits positive H/A moving-Riccati margins, the exact accelerometer and
vector operations will be rebound directly to the *shipping covariance* metric
by covariance congruence, not by the old group-isotropic terminal metric.
The target word inequality is

    V(F(x), P_plus) - V(x, P) <= -mu ||x||_P^2 + R_3(x,xi),

with a validated finite-angle higher-order remainder R_3 over the declared
0.8-rad sector.  Until both the P3 margin and this nonlinear domination are
closed, P4 remains OPEN and P5 may not start.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_sea3_riccati_metric_p3 as P3
import ou3_p4_cayley_sector_certificate as CAYLEY

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_SEA3_MOVING_RICCATI_NONLINEAR_P4"


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    p3 = P3.build(path)
    p3f = P3.validate(p3)
    cayley = CAYLEY.build(path)
    cf = CAYLEY.validate(cayley)
    prereq_failures = [f"P3: {x}" for x in p3f] + [f"Cayley: {x}" for x in cf]
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
        "old_group_isotropic_P3_P4_metric_assumed": False,
        "outer_angle_rad": cayley["outer_angle_rad"],
        "cayley_geometry_validated": True,
        "exact_vector_accelerometer_congruence_rebind_pending": True,
        "moving_covariance_congruence_target": (
            "z_u=T_E z, P_u=T_E P T_E^T; z_u^T P_u^-1 z_u = z^T P^-1 z exactly"
        ),
        "nonlinear_word_inequality": (
            "V(F(x),P_plus)-V(x,P) <= -mu*||x||_P^2 + R3(x,xi)"
        ),
        "P3_CANONICAL_PASS_consumed": p3["P3_CANONICAL_PASS"],
        "nonlinear_remainder_dominated_on_full_sector": False,
        "P4_CANONICAL_PASS": False,
        "P5_MAY_START": False,
        "P4_CANONICAL_FAIL_REASONS": [
            "canonical moving-Riccati P3 margin is not yet closed",
            "exact vector/accelerometer operations must be rebound to the moving covariance metric without the retired group-isotropic attachment",
            "complete H18/A21 recurrent-word nonlinear remainder domination is not yet emitted",
        ],
        "next_obligation": (
            "after P3 emits delta_H,delta_A, bind the exact measurement/reset operations by covariance congruence and certify the finite-angle remainder on 0.8 rad; do not return to endpoint/source-word enumeration"
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
        "cayley_geometry_validated", "exact_vector_accelerometer_congruence_rebind_pending",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "old_800_endpoint_signed_Joseph_scan_consumed",
        "old_terminal_source_phase_metric_attachment_consumed",
        "old_group_isotropic_P3_P4_metric_assumed",
        "P3_CANONICAL_PASS_consumed", "nonlinear_remainder_dominated_on_full_sector",
        "P4_CANONICAL_PASS", "P5_MAY_START",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("outer_angle_rad", 0.0)) < 0.80:
        f.append("declared nonlinear sector fell below 0.8 rad")
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
