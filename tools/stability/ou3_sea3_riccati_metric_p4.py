#!/usr/bin/env python3
"""Canonical, unpromoted P4 over the complete SEA3 moving-Riccati P3 metric.

Prediction/Joseph/reset covariance identities are closed.  Their exact linear
congruence is not a nonlinear storage isometry.  The shipping-H full-shift
residual identity and the full prediction/reset/H-to-A transport identities
still require source-conditioned physical defects and uniform storage bounds.

The next experiment is the non-promoting signed complete-word ratio on legal
same-history SEA3 realizations, with every actual-applied R_S update retained.
Neither structural tests nor an auxiliary H0/congruent H_u identification can
close full nonlinear transport or strict P4 dissipation.  P5 remains blocked.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_sea3_riccati_metric_p3 as P3
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_moving_metric_rebind as REBIND

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 6
QUALIFICATION = "OU3_SEA3_MOVING_RICCATI_NONLINEAR_P4_V6"


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    p3 = P3.build(path)
    p3f = P3.validate(p3)
    cayley = CAYLEY.build(path)
    cf = CAYLEY.validate(cayley)
    rebind = REBIND.build()
    rf = REBIND.validate(rebind)
    prereq_failures = (
        [f"P3: {x}" for x in p3f]
        + [f"Cayley: {x}" for x in cf]
        + [f"rebind: {x}" for x in rf]
    )
    if prereq_failures:
        raise RuntimeError(f"moving-Riccati P4 prerequisites failed: {prereq_failures}")

    p3_pass = bool(p3["P3_CONDITIONAL_SEA3_PASS"])
    h_delta = float(p3["modes"]["H18"]["relative_Riccati_injection_margin_lower"])
    a_delta = float(p3["modes"]["A21"]["relative_Riccati_injection_margin_lower"])
    covariance_closed = bool(rebind["structural_shipping_covariance_identities_closed"])
    transport_closed = bool(rebind["nonlinear_chart_transport_and_storage_closed"])
    remainder_closed = False

    fail_reasons = []
    if not p3_pass:
        fail_reasons.append(
            "canonical moving-Riccati P3 H18/A21 quantitative margin has not met the useful gate"
        )
    if not covariance_closed:
        fail_reasons.append("shipping covariance identities are not closed")
    if not transport_closed:
        fail_reasons.append(
            "full nonlinear shipping transport and uniform storage comparison remain open"
        )
    if not remainder_closed:
        fail_reasons.append(
            "complete H18/A21 recurrent-word nonlinear remainder domination is not yet emitted on the full 0.8-rad geometry sector"
        )

    p4_pass = bool(p3_pass and covariance_closed and transport_closed and remainder_closed)
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
        "exact_vector_accelerometer_congruence_rebind_pending": not covariance_closed,
        "structural_shipping_covariance_identities_closed": covariance_closed,
        "full_nonlinear_measurement_metric_rebind_closed": transport_closed,
        "full_nonlinear_transport_and_storage_closed": transport_closed,
        "moving_metric_rebind_qualification": rebind["qualification"],
        "moving_metric_coordinate_congruence_exact": rebind["moving_metric_coordinate_congruence_exact"],
        "Joseph_nonlinear_injection_metric_closed": rebind["Joseph_nonlinear_injection_metric_closed"],
        "moving_covariance_congruence_target": (
            "z_u=T_E z, P_u=T_E P T_E^T; z_u^T P_u^-1 z_u = z^T P^-1 z exactly"
        ),
        "nonlinear_word_inequality": (
            "V_after(F_W(x)) <= rho_W V_before(x), rho_W < 1, for every admitted complete SEA3 word"
        ),
        "P3_CONDITIONAL_SEA3_PASS_consumed": p3_pass,
        "P3_DEPLOYMENT_PASS_consumed_as_if_closed": False,
        "P3_CANONICAL_PASS_consumed": p3_pass,
        "P3_H18_delta_consumed": h_delta,
        "P3_A21_delta_consumed": a_delta,
        "P3_H_delta_consumed": h_delta,
        "P3_A_delta_consumed": a_delta,
        "nonlinear_coordinate_shipping_binding_closed": False,
        "structural_rebind_does_not_close_nonlinear_coordinate_transport": True,
        "packet_count_remainder_budget_used": False,
        "nonlinear_remainder_dominated_on_full_sector": remainder_closed,
        "P4_FINITE_WINDOW_CLOSED": p4_pass,
        "P4_CANONICAL_PASS": p4_pass,
        "P5_MAY_START": p4_pass,
        "P4_CANONICAL_FAIL_REASONS": fail_reasons,
        "next_obligation": (
            "run the non-promoting legal complete-SEA3 word feasibility experiment before further enclosure; retain full epsilon_aw=(Q_aw-I)*delta_a_w+e_eta transport, actual applied R_S and the separate H-to-A event; follow docs/ou3-proof-research-state.md"
            if p3_pass and covariance_closed
            else "close only the reported prerequisite; do not return to endpoint/source-word enumeration"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_P4_architecture") != "NONLINEAR_WORD_IN_MOVING_SHIPPING_RICCATI_METRIC":
        f.append("wrong canonical P4 architecture")

    for key in (
        "source_generated_not_trajectory_fit",
        "SEA3_dynamic_source_used_through_P3",
        "cayley_geometry_validated",
        "structural_shipping_covariance_identities_closed",
        "moving_metric_coordinate_congruence_exact",
        "Joseph_nonlinear_injection_metric_closed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")

    for key in (
        "nonlinear_coordinate_shipping_binding_closed",
        "full_nonlinear_measurement_metric_rebind_closed",
        "full_nonlinear_transport_and_storage_closed",
        "packet_count_remainder_budget_used",
        "trajectory_replay_used",
        "filter_changed",
        "declared_domain_shrunk",
        "old_800_endpoint_signed_Joseph_scan_consumed",
        "old_terminal_source_phase_metric_attachment_consumed",
        "old_group_isotropic_P3_P4_metric_assumed",
        "exact_vector_accelerometer_congruence_rebind_pending",
        "nonlinear_remainder_dominated_on_full_sector",
        "P4_FINITE_WINDOW_CLOSED",
        "P4_CANONICAL_PASS",
        "P5_MAY_START",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")

    if d.get("P3_CONDITIONAL_SEA3_PASS_consumed") is not True: f.append("P4 did not consume the closed conditional SEA3 P3 verdict")
    if d.get("P3_DEPLOYMENT_PASS_consumed_as_if_closed") is not False: f.append("P4 incorrectly consumed the still-open deployment P3 verdict")
    if d.get("P3_CANONICAL_PASS_consumed") is not True: f.append("deprecated P3 compatibility alias is inconsistent")
    for key in (
        "P3_H18_delta_consumed", "P3_A21_delta_consumed",
        "P3_H_delta_consumed", "P3_A_delta_consumed",
    ):
        x = d.get(key)
        if not isinstance(x, (int, float)) or float(x) < 1.0e-18:
            f.append(f"{key} fell below the useful P3 gate")

    if float(d.get("outer_angle_rad", 0.0)) < 0.80:
        f.append("declared nonlinear sector fell below 0.8 rad")
    reasons = d.get("P4_CANONICAL_FAIL_REASONS", [])
    if len(reasons) != 2 or not any("uniform storage" in x.lower() for x in reasons) or not any(
        "nonlinear remainder" in x.lower() for x in reasons
    ):
        f.append("P4 must report both nonlinear transport/storage and complete-word dissipation blockers")
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
        "P3_CONDITIONAL_SEA3_PASS_consumed": d["P3_CONDITIONAL_SEA3_PASS_consumed"],
        "P3_DEPLOYMENT_PASS_consumed_as_if_closed": d["P3_DEPLOYMENT_PASS_consumed_as_if_closed"],
        "metric_rebind_closed": d["full_nonlinear_measurement_metric_rebind_closed"],
        "P4_CANONICAL_PASS": d["P4_CANONICAL_PASS"],
        "fail_reasons": d["P4_CANONICAL_FAIL_REASONS"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
