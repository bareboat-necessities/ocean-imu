#!/usr/bin/env python3
"""Canonical unpromoted P4: full-state differential contraction on complete SEA3.

The nonlinear theorem metric is the full-rank pullback

    M(z,zeta)=D Phi(z,zeta)^T P(zeta)^-1 D Phi(z,zeta),

with the exact accelerometer-linearizing Phi coordinate and all active H18/A21
states retained. At zero error D Phi=I, so this is an exact nonlinear
differential continuation of the frozen complete-SEA3 P3 metric.

The complete nonlinear Jacobian is assembled by the literal differential
cocycle. Prediction is differentiated from the SAME source omega/dt/committed
tau using the deployed quaternion and integrated-OU mean. Every Joseph event
derives H,S,K from the SAME source-correlated P/H/R cell; a due S=0 event is
accepted only with ``ACTUAL_APPLIED_SPECTRALMSE_RS`` provenance. In A21 the
shipping residual-bias ball projection is differentiated by its exact outward
Clarke generalized Jacobian and consumes the same-source true-bias cell; it is
not assumed inactive. The H18->A21 homogeneous lift is exactly [I18;0], while
the held b_a error and shipping b_a covariance floor remain explicit separate
hybrid forcing/metric obligations.

The P4 gate is the full matrix inequality

    rho M_0 - D F_W(z)^T M_1 D F_W(z) > 0,  0<rho<1,

uniformly over every admitted complete SEA3 word/state cell. No finite endpoint
energy, packetwise remainder, correction radius, inverse-metric floor, state
elimination, independent K/R_S schedule, replay or longer sampled-word
optimization may promote this gate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_sea3_riccati_metric_p3 as P3
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_complete_sea3_phi_differential_metric as DIFF
import ou3_p4_complete_sea3_differential_prediction as PRED
import ou3_p4_complete_sea3_differential_events as EVENTS
import ou3_p4_complete_sea3_differential_word as DWRD

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 10
QUALIFICATION = "OU3_SEA3_FULL_STATE_DIFFERENTIAL_P4_V10"


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    p3 = P3.build(path)
    cayley = CAYLEY.build(path)
    diff = DIFF.build(path)
    pred = PRED.build(path)
    events = EVENTS.build(path)
    word = DWRD.build(path)
    failures = (
        [f"P3: {x}" for x in P3.validate(p3)]
        + [f"Cayley: {x}" for x in CAYLEY.validate(cayley)]
        + [f"differential metric: {x}" for x in DIFF.validate(diff)]
        + [f"prediction: {x}" for x in PRED.validate(pred)]
        + [f"Joseph events: {x}" for x in EVENTS.validate(events)]
        + [f"differential word: {x}" for x in DWRD.validate(word)]
    )
    if failures:
        raise RuntimeError(f"canonical differential P4 prerequisites failed: {failures}")

    p3_pass = bool(p3["P3_CONDITIONAL_SEA3_PASS"])
    h_delta = float(p3["modes"]["H18"]["relative_Riccati_injection_margin_lower"])
    a_delta = float(p3["modes"]["A21"]["relative_Riccati_injection_margin_lower"])

    pred_closed = bool(pred["source_uniform_prediction_Jacobian_closed"])
    joseph_closed = bool(events["source_uniform_finite_angle_event_Jacobians_closed"])
    projection_differentiated = bool(
        events["A21_bias_projection_generalized_Jacobian_available"]
        and events["A21_bias_projection_same_source_true_bias_required"]
        and word["A21_bias_projection_generalized_Jacobian_available"]
    )
    word_closed = bool(word["source_uniform_complete_word_Jacobian_enclosed"])
    metric_word_closed = bool(diff["source_uniform_complete_word_Jacobian_enclosed"])
    jacobian_closed = bool(
        pred_closed and joseph_closed and projection_differentiated
        and word_closed and metric_word_closed
    )
    contraction_closed = bool(diff["source_uniform_pullback_differential_contraction_closed"])
    p4_pass = bool(p3_pass and jacobian_closed and contraction_closed)

    fail_reasons: list[str] = []
    if not p3_pass:
        fail_reasons.append("frozen conditional complete-SEA3 P3 prerequisite is not closed")
    if not jacobian_closed:
        open_parts = []
        if not pred_closed:
            open_parts.append("same-history prediction")
        if not joseph_closed:
            open_parts.append("same-history P/H/R/true-bias Joseph attachment")
        if not projection_differentiated:
            open_parts.append("A21 projection generalized Jacobian")
        if not word_closed or not metric_word_closed:
            open_parts.append("literal complete-word composition")
        fail_reasons.append(
            "source-uniform outward full-state Jacobian of the literal complete SEA3 word is open: "
            + ", ".join(open_parts)
        )
    if not contraction_closed:
        fail_reasons.append(
            "full interval-LDLT pullback differential contraction is not yet strict on a declared [30,25,20,15] degree cell"
        )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_P4_architecture": "FULL_STATE_COMPLETE_SEA3_DIFFERENTIAL_PULLBACK",
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "source_family_replaced": False,
        "P3_CONDITIONAL_SEA3_PASS_consumed": p3_pass,
        "P3_DEPLOYMENT_PASS_consumed_as_if_closed": False,
        "P3_H18_delta_consumed": h_delta,
        "P3_A21_delta_consumed": a_delta,
        "P3_frozen_not_modified": True,
        "outer_angle_rad": float(cayley["outer_angle_rad"]),
        "candidate_angles_deg": [30.0, 25.0, 20.0, 15.0],

        "differential_metric_qualification": diff["qualification"],
        "differential_metric_type": diff["metric_type"],
        "differential_metric_definition": diff["metric_definition"],
        "prediction_differential_qualification": pred["qualification"],
        "Joseph_differential_qualification": events["qualification"],
        "literal_differential_word_qualification": word["qualification"],

        "finite_Phi_storage_used_as_Lyapunov_function": False,
        "finite_raw_endpoint_storage_used_as_P4_certificate": False,
        "all_active_states_retained": bool(diff["all_active_states_retained"]),
        "H18_full_rank": bool(diff["H18_full_rank"]),
        "A21_full_rank": bool(diff["A21_full_rank"]),
        "Phi_jacobian_determinant_exact": diff["Phi_jacobian_determinant_exact"],
        "metric_reduces_exactly_to_P3_at_zero_error": bool(diff["metric_reduces_to_P3_at_zero_error"]),

        "same_complete_SEA3_word_required": True,
        "same_frontend_tuner_covariance_history_required": True,
        "same_source_omega_h_tau_prediction_required": bool(pred["same_complete_SEA3_omega_h_tau_required"]),
        "prediction_independent_F_forbidden": not bool(pred["independent_F_input_allowed_for_theorem"]),
        "full_prediction_F_Eaw_rows_retained": bool(pred["full_F_Eaw_v_p_S_aw_rows_retained"]),
        "same_P_H_R_cell_required_for_Joseph": bool(events["same_P_H_R_cell_derives_S_and_K"]),
        "independent_K_forbidden": not bool(events["independent_K_input_allowed_for_theorem"]),
        "all_due_S_updates_with_actual_applied_RS_required": True,
        "actual_RS_provenance_token": events["actual_applied_RS_provenance_token"],
        "all_valid_accelerometer_updates_required": True,
        "all_process_Q_floor_reset_events_required": True,

        "A21_bias_projection_generalized_Jacobian_available": bool(
            events["A21_bias_projection_generalized_Jacobian_available"]
        ),
        "A21_bias_projection_same_source_true_bias_required": bool(
            events["A21_bias_projection_same_source_true_bias_required"]
        ),
        "A21_bias_projection_hybrid_differentiated": projection_differentiated,
        "A21_bias_projection_assumed_inactive": False,

        "H_to_A_rectangular_differential_event_required": True,
        "H_to_A_homogeneous_lift": word["H_to_A_homogeneous_lift"],
        "H_to_A_held_ba_error_retained_as_separate_forcing": bool(
            word["H_to_A_held_ba_error_retained_as_separate_forcing"]
        ),
        "H_to_A_covariance_floor_retained_as_separate_metric_event": bool(
            word["H_to_A_covariance_floor_retained_as_separate_metric_event"]
        ),

        "literal_cocycle_same_source_token_required": bool(word["same_complete_SEA3_source_token_required_for_every_event"]),
        "literal_cocycle_same_cell_Joseph_required": bool(word["same_P_H_R_cell_required_for_every_Joseph_event"]),
        "outward_interval_AD_used": True,
        "full_interval_LDLT_required": True,
        "differential_word_inequality": diff["differential_word_inequality"],

        "source_uniform_prediction_Jacobian_closed": pred_closed,
        "source_uniform_Joseph_Jacobians_closed": joseph_closed,
        "source_uniform_literal_cocycle_closed": word_closed,
        "source_uniform_complete_word_Jacobian_enclosed": jacobian_closed,
        "source_uniform_pullback_differential_contraction_closed": contraction_closed,

        "packet_count_remainder_budget_used": False,
        "packetwise_remainder_norm_sum_used": False,
        "state_elimination_used": False,
        "a_w_Schur_final_certificate_used": False,
        "correction_radius_claim_used": False,
        "inverse_metric_floor_claim_used": False,
        "independent_RS_schedule_used": False,
        "point_word_rho_used_to_promote": False,
        "longer_point_window_optimization_used_to_promote": False,

        "P4_FINITE_WINDOW_CLOSED": p4_pass,
        "P4_CANONICAL_PASS": p4_pass,
        "P5_MAY_START": p4_pass,
        "P4_CANONICAL_FAIL_REASONS": fail_reasons,
        "next_obligation": (
            "close the source-uniform same-history prediction and same-cell P/H/R/true-bias Joseph event enclosures, compose them through the literal cocycle with every actual-applied R_S event, then test rho*M0-J^T*M1*J by full interval LDLT from 30 degrees downward; if no cell closes, perform theorem/metric failure analysis rather than retuning or replacing the source"
            if p3_pass else "repair only frozen P3 prerequisite failure"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_P4_architecture") != "FULL_STATE_COMPLETE_SEA3_DIFFERENTIAL_PULLBACK":
        f.append("canonical P4 architecture changed")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")

    for key in (
        "source_generated_not_trajectory_fit", "P3_CONDITIONAL_SEA3_PASS_consumed",
        "P3_frozen_not_modified", "all_active_states_retained", "H18_full_rank", "A21_full_rank",
        "metric_reduces_exactly_to_P3_at_zero_error", "same_complete_SEA3_word_required",
        "same_frontend_tuner_covariance_history_required", "same_source_omega_h_tau_prediction_required",
        "prediction_independent_F_forbidden", "full_prediction_F_Eaw_rows_retained",
        "same_P_H_R_cell_required_for_Joseph", "independent_K_forbidden",
        "all_due_S_updates_with_actual_applied_RS_required", "all_valid_accelerometer_updates_required",
        "all_process_Q_floor_reset_events_required", "A21_bias_projection_generalized_Jacobian_available",
        "A21_bias_projection_same_source_true_bias_required", "A21_bias_projection_hybrid_differentiated",
        "H_to_A_rectangular_differential_event_required", "H_to_A_held_ba_error_retained_as_separate_forcing",
        "H_to_A_covariance_floor_retained_as_separate_metric_event",
        "literal_cocycle_same_source_token_required", "literal_cocycle_same_cell_Joseph_required",
        "outward_interval_AD_used", "full_interval_LDLT_required",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")

    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk", "source_family_replaced",
        "P3_DEPLOYMENT_PASS_consumed_as_if_closed", "finite_Phi_storage_used_as_Lyapunov_function",
        "finite_raw_endpoint_storage_used_as_P4_certificate", "A21_bias_projection_assumed_inactive",
        "source_uniform_prediction_Jacobian_closed", "source_uniform_Joseph_Jacobians_closed",
        "source_uniform_literal_cocycle_closed", "source_uniform_complete_word_Jacobian_enclosed",
        "source_uniform_pullback_differential_contraction_closed", "packet_count_remainder_budget_used",
        "packetwise_remainder_norm_sum_used", "state_elimination_used", "a_w_Schur_final_certificate_used",
        "correction_radius_claim_used", "inverse_metric_floor_claim_used", "independent_RS_schedule_used",
        "point_word_rho_used_to_promote", "longer_point_window_optimization_used_to_promote",
        "P4_FINITE_WINDOW_CLOSED", "P4_CANONICAL_PASS", "P5_MAY_START",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")

    for key in ("P3_H18_delta_consumed", "P3_A21_delta_consumed"):
        if float(d.get(key, 0.0)) != 1.0e-18:
            f.append(f"{key} changed from frozen delta")
    if d.get("actual_RS_provenance_token") != EVENTS.ACTUAL_RS_PROVENANCE:
        f.append("actual R_S provenance token changed")
    if d.get("H_to_A_homogeneous_lift") != "[I18;0]":
        f.append("H->A homogeneous lift changed")
    if float(d.get("Phi_jacobian_determinant_exact", 0.0)) != 1.0:
        f.append("Phi differential coordinate lost exact full rank")
    if float(d.get("outer_angle_rad", 0.0)) < 0.80:
        f.append("outer Cayley sector shrank")
    if d.get("candidate_angles_deg") != [30.0, 25.0, 20.0, 15.0]:
        f.append("declared P4 candidate sequence changed")
    reasons = d.get("P4_CANONICAL_FAIL_REASONS", [])
    if len(reasons) != 2 or not any("jacobian" in x.lower() for x in reasons) or not any(
        "ldlt" in x.lower() or "contraction" in x.lower() for x in reasons
    ):
        f.append("canonical P4 must expose Jacobian and differential-contraction blockers")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "architecture": d["canonical_P4_architecture"],
        "prediction_closed": d["source_uniform_prediction_Jacobian_closed"],
        "Joseph_closed": d["source_uniform_Joseph_Jacobians_closed"],
        "A21_projection_differentiated": d["A21_bias_projection_hybrid_differentiated"],
        "literal_cocycle_closed": d["source_uniform_literal_cocycle_closed"],
        "actual_RS_required": d["all_due_S_updates_with_actual_applied_RS_required"],
        "P4_CANONICAL_PASS": d["P4_CANONICAL_PASS"],
        "fail_reasons": d["P4_CANONICAL_FAIL_REASONS"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
