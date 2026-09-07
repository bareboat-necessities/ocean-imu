#!/usr/bin/env python3
"""Canonical, unpromoted P4 over the complete SEA3 moving-Riccati P3 metric.

The canonical nonlinear object is now the complete-word endpoint identity, not
an accumulation of packetwise remainder norms.  The full accelerometer-
linearizing shift

    epsilon_aw=(Q_aw-I)delta_a_w+e_eta

is transported through the literal same-history shipping word.  Variation of
constants reduces all prediction/source/hybrid/floor shift transport exactly;
S=0 and magnetometer events have zero interior a_w-shift term.  The remaining
interior shift is the JOINT suffix-weighted accepted-accelerometer operator

    B_W=[M_suffix G K H E_aw]_{acc(W)},

whose suffixes contain every later shipping event, including every due S=0
update with its actual applied per-axis SpectralMSE R_S.

Thus the P4 endpoint has the exact form

    d_W=r_W+E_N epsilon_N-M_W E_0 epsilon_0-B_W epsilon_acc,

and

    Delta V = -Phi_0^T D_W Phi_0
              +2(M_W Phi_0)^T P_N^-1 d_W
              +d_W^T P_N^-1 d_W.

P3 remains frozen at delta=1e-18.  P4 is not closed until the joint nonlinear
endpoint object is enclosed over the SAME complete SEA3 word and the last two
terms are shown to fit inside the full-matrix P3 decrease on one declared
finite-angle candidate.  No correction radius, inverse-metric floor,
packet-count budget, independent R_S schedule, replay, or alternate estimator
is used.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_sea3_riccati_metric_p3 as P3
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_moving_metric_rebind as REBIND
import ou3_p4_complete_word_endpoint_transport as ENDPOINT

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 7
QUALIFICATION = "OU3_SEA3_MOVING_RICCATI_NONLINEAR_P4_V7"


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    p3 = P3.build(path)
    p3f = P3.validate(p3)
    cayley = CAYLEY.build(path)
    cf = CAYLEY.validate(cayley)
    rebind = REBIND.build()
    rf = REBIND.validate(rebind)
    endpoint = ENDPOINT.build(path)
    ef = ENDPOINT.validate(endpoint)
    prereq_failures = (
        [f"P3: {x}" for x in p3f]
        + [f"Cayley: {x}" for x in cf]
        + [f"rebind: {x}" for x in rf]
        + [f"endpoint: {x}" for x in ef]
    )
    if prereq_failures:
        raise RuntimeError(f"moving-Riccati P4 prerequisites failed: {prereq_failures}")

    p3_pass = bool(p3["P3_CONDITIONAL_SEA3_PASS"])
    h_delta = float(p3["modes"]["H18"]["relative_Riccati_injection_margin_lower"])
    a_delta = float(p3["modes"]["A21"]["relative_Riccati_injection_margin_lower"])
    covariance_closed = bool(rebind["structural_shipping_covariance_identities_closed"])
    endpoint_identity_closed = bool(endpoint["master_inequality_object_emitted"])
    joint_endpoint_closed = bool(endpoint["source_uniform_master_endpoint_domination_closed"])

    # Structural coordinate/covariance identities are valid, but the nonlinear
    # Phi storage is not declared isometric to the original physical storage.
    # The endpoint master inequality is the route that must close that gap.
    transport_and_storage_closed = bool(
        endpoint_identity_closed and joint_endpoint_closed
    )
    remainder_closed = joint_endpoint_closed

    fail_reasons = []
    if not p3_pass:
        fail_reasons.append(
            "canonical moving-Riccati P3 H18/A21 quantitative margin has not met the useful gate"
        )
    if not covariance_closed:
        fail_reasons.append("shipping covariance identities are not closed")
    if not endpoint_identity_closed:
        fail_reasons.append("whole-word full-shift endpoint identity is not closed")
    if not joint_endpoint_closed:
        fail_reasons.append(
            "source-uniform joint complete-SEA3 endpoint defect domination is open: enclose r_W, endpoint epsilon terms and the suffix-weighted accelerometer B_W*epsilon history in the full endpoint metric"
        )
    if not transport_and_storage_closed:
        fail_reasons.append(
            "full nonlinear storage comparison is not closed until the joint endpoint inequality is negative on a declared finite-angle candidate"
        )

    p4_pass = bool(
        p3_pass and covariance_closed and endpoint_identity_closed
        and transport_and_storage_closed and remainder_closed
    )
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
        "whole_word_endpoint_transport_qualification": endpoint["qualification"],
        "whole_word_endpoint_transport_consumed": endpoint_identity_closed,
        "exact_full_shift_endpoint_decomposition_closed": endpoint_identity_closed,
        "endpoint_defect_formula": endpoint["endpoint_decomposition_identity"],
        "joint_accelerometer_endpoint_operator": endpoint["accelerometer_joint_operator"],
        "accepted_accelerometer_only_interior_shift_event_class": endpoint[
            "accepted_accelerometer_is_only_interior_epsilon_event_class"
        ],
        "actual_RS_regularization_retained_in_endpoint_operator": endpoint[
            "actual_RS_regularization_enters_every_applicable_suffix"
        ],
        "prediction_source_hybrid_floor_shift_telescoping_closed": endpoint[
            "prediction_source_hybrid_floor_interior_epsilon_terms_cancel_exactly"
        ],
        "S_and_mag_interior_epsilon_terms_zero": bool(
            endpoint["S_zero_interior_epsilon_term_zero_exactly"]
            and endpoint["magnetometer_interior_epsilon_term_zero_exactly"]
        ),
        "master_endpoint_energy_identity": endpoint["master_endpoint_energy_identity"],
        "master_endpoint_D_W_definition": endpoint["D_W_definition"],
        "source_uniform_joint_BW_epsilon_enclosure_closed": endpoint[
            "source_uniform_joint_BW_epsilon_enclosure_closed"
        ],
        "source_uniform_r_word_enclosure_closed": endpoint[
            "source_uniform_r_word_enclosure_closed"
        ],
        "source_uniform_master_endpoint_domination_closed": joint_endpoint_closed,
        "nonlinear_coordinate_shipping_binding_closed": transport_and_storage_closed,
        "full_nonlinear_measurement_metric_rebind_closed": transport_and_storage_closed,
        "full_nonlinear_transport_and_storage_closed": transport_and_storage_closed,
        "structural_rebind_does_not_close_nonlinear_coordinate_transport": True,
        "packet_count_remainder_budget_used": False,
        "packetwise_remainder_norm_sum_used": False,
        "independent_RS_schedule_used": False,
        "correction_radius_claim_used": False,
        "inverse_metric_floor_claim_used": False,
        "nonlinear_remainder_dominated_on_full_sector": remainder_closed,
        "P4_FINITE_WINDOW_CLOSED": p4_pass,
        "P4_CANONICAL_PASS": p4_pass,
        "P5_MAY_START": p4_pass,
        "P4_CANONICAL_FAIL_REASONS": fail_reasons,
        "next_obligation": (
            "enclose the joint suffix-weighted accelerometer operator B_W and r_W over the SAME complete SEA3 H18/A21 execution; retain actual applied per-axis R_S inside every suffix map; prove the two endpoint nonlinear terms fit inside the full-matrix P3 decrease for the widest [30,25,20,15] degree candidate; no packetwise scalarization"
            if p3_pass and covariance_closed and endpoint_identity_closed
            else "close only the named prerequisite; do not return to point-word optimization or source shortcuts"
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
        "whole_word_endpoint_transport_consumed",
        "exact_full_shift_endpoint_decomposition_closed",
        "accepted_accelerometer_only_interior_shift_event_class",
        "actual_RS_regularization_retained_in_endpoint_operator",
        "prediction_source_hybrid_floor_shift_telescoping_closed",
        "S_and_mag_interior_epsilon_terms_zero",
        "structural_rebind_does_not_close_nonlinear_coordinate_transport",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")

    for key in (
        "nonlinear_coordinate_shipping_binding_closed",
        "full_nonlinear_measurement_metric_rebind_closed",
        "full_nonlinear_transport_and_storage_closed",
        "source_uniform_joint_BW_epsilon_enclosure_closed",
        "source_uniform_r_word_enclosure_closed",
        "source_uniform_master_endpoint_domination_closed",
        "packet_count_remainder_budget_used",
        "packetwise_remainder_norm_sum_used",
        "independent_RS_schedule_used",
        "correction_radius_claim_used",
        "inverse_metric_floor_claim_used",
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

    if d.get("P3_CONDITIONAL_SEA3_PASS_consumed") is not True:
        f.append("P4 did not consume the closed conditional SEA3 P3 verdict")
    if d.get("P3_DEPLOYMENT_PASS_consumed_as_if_closed") is not False:
        f.append("P4 incorrectly consumed the still-open deployment P3 verdict")
    if d.get("P3_CANONICAL_PASS_consumed") is not True:
        f.append("deprecated P3 compatibility alias is inconsistent")
    for key in (
        "P3_H18_delta_consumed", "P3_A21_delta_consumed",
        "P3_H_delta_consumed", "P3_A_delta_consumed",
    ):
        x = d.get(key)
        if not isinstance(x, (int, float)) or float(x) != 1.0e-18:
            f.append(f"{key} changed from frozen 1e-18 gate")

    if float(d.get("outer_angle_rad", 0.0)) < 0.80:
        f.append("declared nonlinear sector fell below 0.8 rad")
    reasons = d.get("P4_CANONICAL_FAIL_REASONS", [])
    if len(reasons) != 2 or not any("joint complete-sea3 endpoint" in x.lower() for x in reasons) or not any(
        "storage comparison" in x.lower() for x in reasons
    ):
        f.append("P4 must report the joint endpoint domination/storage blockers")
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
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "architecture": d["canonical_P4_architecture"],
        "P3_CONDITIONAL_SEA3_PASS_consumed": d["P3_CONDITIONAL_SEA3_PASS_consumed"],
        "endpoint_transport_consumed": d["whole_word_endpoint_transport_consumed"],
        "actual_RS_in_endpoint_operator": d["actual_RS_regularization_retained_in_endpoint_operator"],
        "P4_CANONICAL_PASS": d["P4_CANONICAL_PASS"],
        "fail_reasons": d["P4_CANONICAL_FAIL_REASONS"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
