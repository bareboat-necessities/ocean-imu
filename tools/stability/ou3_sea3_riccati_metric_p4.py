#!/usr/bin/env python3
"""Unpromoted complete-SEA3 P4 gates for motion and full-state targets.

The primary relaxed motion contract is built separately by
ou3_p4_bounded_bias_motion. The full-state construction below is retained as
a stronger unclosed extension; its flags cannot stand in for motion gains.

The stronger P4 object is the source-indexed finite-state quadratic storage

    V(e,zeta) = e^T M(zeta) e,

with two obligations on one complete same-history SEA3 word:

    V_{k+N_W} <= rho V_k + forcing,          0 < rho < 1,
    V_{k+ell} <= kappa_V V_k + forcing,      0 <= ell < N_W,

and every prefix must remain in the certified chart/source domain.

Outward differential/Clarke AD is retained only to enclose the finite physical
map through the generalized mean-value theorem.  It does not redefine P4 as a
state-dependent differential/Riemannian metric.  The finite endpoint algebra
retains the exact full epsilon_aw shift and literal shipping suffix maps, so
every later S=0 update with its actual applied anisotropic SpectralMSE R_S stays
inside the nonlinear defect transport.  Packetwise norm sums are forbidden.

The exact signed-information ledger is also canonical P4 machinery.  Each
accepted Joseph/reset event is organized as

    Delta V = -I_y + E_eta + X_reset + E_reset.

For S=0, eta=0 exactly, so every due S update contributes its favorable
-y^T S^-1 y term with the ACTUAL applied R_S and no nonlinear residual charge.
Accelerometer/vector eta and finite Cayley reset defects remain explicit and
must be dominated jointly over the same complete SEA3 word; they may not be
replaced by an event-count worst-case remainder budget.

Conditional complete-SEA3 P3 is consumed unchanged at delta=1e-18.  P4/P5 stay
open until both the finite endpoint inequality and every finite prefix gain plus
domain-retention condition close on the same source-correlated complete word.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_sea3_riccati_metric_p3 as P3
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_complete_word_endpoint_transport as ENDPOINT
import ou3_p4_complete_sea3_finite_map_mean_value as FINITE
import ou3_p4_complete_sea3_signed_information_ledger as SIGNED
import ou3_p4_complete_sea3_joint_sector_master as JOINT
import ou3_mems_bias_contract as BIAS
import ou3_sea3_response_union as UNION
import ou3_p4_bounded_bias_motion as MOTION

# Exact/outward finite-map differentiation machinery only.
import ou3_p4_complete_sea3_phi_differential_metric as DIFF
import ou3_p4_complete_sea3_differential_prediction as PRED
import ou3_p4_complete_sea3_differential_events as EVENTS
import ou3_p4_complete_sea3_differential_word as DWRD

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 14
QUALIFICATION = "OU3_SEA3_FINITE_STATE_ENDPOINT_PREFIX_P4_V14"
ARCHITECTURE = "FINITE_STATE_COMPLETE_SEA3_QUADRATIC_ENDPOINT_AND_PREFIX"
CANDIDATES = [30.0, 25.0, 20.0, 15.0]


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    p3 = P3.build(path)
    motion = MOTION.build(p3_contract=p3)
    cayley = CAYLEY.build(path)
    endpoint = ENDPOINT.build(path)
    finite = FINITE.build(path)
    signed = SIGNED.build(path)
    joint = JOINT.build(p3_contract=p3, signed_contract=signed)
    diff = DIFF.build(path)
    pred = PRED.build(path)
    events = EVENTS.build(path)
    word = DWRD.build(path)
    failures = (
        [f"P3: {x}" for x in P3.validate(p3)]
        + [f"motion target: {x}" for x in MOTION.validate(motion)]
        + [f"Cayley: {x}" for x in CAYLEY.validate(cayley)]
        + [f"endpoint: {x}" for x in ENDPOINT.validate(endpoint)]
        + [f"finite-map bridge: {x}" for x in FINITE.validate(finite)]
        + [f"signed-information ledger: {x}" for x in SIGNED.validate(signed)]
        + [f"joint-sector master: {x}" for x in JOINT.validate(joint)]
        + [f"differential metric primitive: {x}" for x in DIFF.validate(diff)]
        + [f"prediction AD primitive: {x}" for x in PRED.validate(pred)]
        + [f"event AD primitive: {x}" for x in EVENTS.validate(events)]
        + [f"differential-word primitive: {x}" for x in DWRD.validate(word)]
    )
    if failures:
        raise RuntimeError(f"canonical finite-state P4 prerequisites failed: {failures}")

    bias = p3["mems_bias_preconditions"]
    p3_pass = bool(p3["P3_CONDITIONAL_SEA3_PASS"])
    h_delta = float(p3["modes"]["H18"]["relative_Riccati_injection_margin_lower"])
    a_delta = float(p3["modes"]["A21"]["relative_Riccati_injection_margin_lower"])

    endpoint_master_emitted = bool(endpoint["master_inequality_object_emitted"])
    endpoint_closed = bool(endpoint["source_uniform_master_endpoint_domination_closed"])
    finite_jacobian_closed = bool(finite["source_uniform_complete_word_generalized_Jacobian_enclosed"])
    finite_endpoint_closed = bool(finite["source_uniform_endpoint_finite_map_closed"])
    prefix_gain_closed = bool(finite["source_uniform_all_prefix_gains_closed"])
    prefix_domain_closed = bool(finite["source_uniform_all_prefix_domains_closed"])
    signed_ledger_ready = bool(signed["joint_complete_word_signed_information_composition_available"])
    signed_domination_closed = bool(signed["source_uniform_joint_eta_reset_domination_closed"])
    joint_sector_closed = bool(joint["source_uniform_same_history_joint_sector_closed"])
    joint_ldlt_closed = bool(joint["source_uniform_full_augmented_LDLT_closed"])
    response_coverage = {
        obligation: UNION.branch_mode_coverage({"H18": closed, "A21": closed})
        for obligation, closed in {
            "endpoint": finite_endpoint_closed and joint_ldlt_closed and joint_sector_closed,
            "prefix_gain": prefix_gain_closed,
            "prefix_retention": prefix_domain_closed,
        }.items()
    }

    projection_machinery = bool(
        events["A21_bias_projection_generalized_Jacobian_available"]
        and events["A21_bias_projection_same_source_true_bias_required"]
        and word["A21_bias_projection_generalized_Jacobian_available"]
    )

    # Both endpoint formulations and the signed-information organization must
    # close on the same finite physical map.  None may be replaced by a
    # differential-only, packetwise, or point-word gate.
    p4_pass = bool(
        p3_pass
        and all(UNION.all_branches_modes_closed(c) for c in response_coverage.values())
        and endpoint_master_emitted
        and signed_ledger_ready
        and signed_domination_closed
        and joint_sector_closed
        and joint_ldlt_closed
        and endpoint_closed
        and finite_jacobian_closed
        and finite_endpoint_closed
        and prefix_gain_closed
        and prefix_domain_closed
        and bias["P4_bias_preconditions_closed"]
    )

    fail_reasons: list[str] = []
    if not p3_pass:
        fail_reasons.append("frozen conditional complete-SEA3 P3 prerequisite is not closed")
    if not (
        signed_domination_closed
        and joint_sector_closed
        and joint_ldlt_closed
        and endpoint_closed
        and finite_jacobian_closed
        and finite_endpoint_closed
    ):
        fail_reasons.append(
            "paper finite-state endpoint dissipation is open: materialize the source-correlated same-history nonlinear "
            "graph sectors and close -(L_W+sum lambda_j Pi_j)>0 by full augmented interval LDLT, retaining all A21 "
            "finite-tau_b cross terms and every actual-R_S S event in the complete-word suffixes; "
            "BIAS1 projection compatibility and BIAS2 full corrected-error separation remain unclosed"
        )
    if not (prefix_gain_closed and prefix_domain_closed):
        fail_reasons.append(
            "paper finite-window prefix gain/domain retention is open for every prefix of the same complete SEA3 word"
        )

    return {
        "schema": SCHEMA,
        "primary_proof_target": MOTION.TARGET,
        "bounded_bias_motion_contract": motion,
        "P4_MOTION_PASS": motion["P4_MOTION_PASS"],
        "P5_MOTION_MAY_START": motion["P5_MOTION_MAY_START"],
        "legacy_full_state_gate_scope": "stronger unclosed extension, not the primary motion target",
        "qualification": QUALIFICATION,
        "canonical_P4_architecture": ARCHITECTURE,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "SEA3_response_union": p3["SEA3_response_union"],
        "response_branch_conditional_P3_coverage_consumed": p3["response_branch_conditional_P3_coverage"],
        "response_branch_P4_coverage": response_coverage,
        "mems_bias_preconditions": bias,
        "paper_Lyapunov_function": "V(e,zeta)=e^T M(zeta)e",
        "paper_endpoint_inequality": "V_{k+N_W} <= rho*V_k + gamma_s*D_s,k + gamma_n*D_n,k; 0<rho<1",
        "paper_prefix_gain_inequality": "V_{k+ell} <= kappa_V*V_k + kappa_s*D_s,k + kappa_n*D_n,k; 0<=ell<N_W",
        "finite_state_endpoint_dissipation_required": True,
        "finite_state_prefix_gain_required": True,
        "prefix_chart_and_source_domain_retention_required": True,
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
        "candidate_angles_deg": list(CANDIDATES),

        "same_complete_SEA3_word_required": True,
        "same_frontend_tuner_covariance_history_required": True,
        "same_source_correlated_generalized_Jacobian_required": bool(
            finite["same_source_correlated_generalized_Jacobian_required"]
        ),
        "all_radial_segment_Jacobians_must_be_enclosed": bool(
            finite["all_radial_segment_Jacobians_must_be_enclosed"]
        ),
        "all_due_S_updates_with_actual_applied_RS_required": True,
        "actual_RS_provenance_token": EVENTS.ACTUAL_RS_PROVENANCE,
        "all_valid_accelerometer_updates_required": True,
        "all_process_Q_floor_reset_events_required": True,
        "asynchronous_vector_events_required": True,
        "H_to_A_rectangular_hybrid_event_required": True,
        "H_to_A_homogeneous_lift": word["H_to_A_homogeneous_lift"],
        "H_to_A_held_ba_error_retained_as_separate_forcing": bool(
            word["H_to_A_held_ba_error_retained_as_separate_forcing"]
        ),
        "H_to_A_covariance_floor_retained_as_separate_metric_event": bool(
            word["H_to_A_covariance_floor_retained_as_separate_metric_event"]
        ),
        "full_prediction_F_Eaw_rows_retained": bool(endpoint["full_prediction_F_Eaw_retained"]),
        "full_epsilon_aw_retained": bool(endpoint["full_epsilon_aw_retained"]),

        "signed_information_qualification": signed["qualification"],
        "signed_information_composition_available": signed_ledger_ready,
        "signed_joseph_reset_identity": signed["combined_joseph_reset_identity"],
        "S_zero_nonlinear_eta_exactly_zero": bool(signed["S_zero_nonlinear_eta_exactly_zero"]),
        "actual_RS_S_information_credited_as_favorable_term": bool(
            signed["every_due_S_update_contributes_negative_information_with_actual_RS"]
        ),
        "finite_reset_defect_explicit_in_signed_ledger": bool(signed["finite_reset_defect_remains_explicit"]),
        "source_uniform_joint_eta_reset_domination_closed": signed_domination_closed,

        "joint_sector_master_contract": joint,
        "joint_sector_master_qualification": joint["qualification"],
        "full_state_joint_sector_master_available": bool(
            joint["terminal_full_augmented_interval_LDLT_available"]
        ),
        "A21_finite_taub_full_matrix_P3_required": bool(
            joint["full_A21_finite_taub_P3_margin_required"]
        ),
        "A21_finite_taub_detectability_consumed": bool(
            joint["finite_taub_A21_detectability_consumed"]
        ),
        "joint_sector_all_21_state_cross_terms_retained": bool(
            joint["all_21_state_cross_terms_retained"]
        ),
        "joint_sector_actual_RS_provenance_preserved": bool(
            joint["actual_RS_provenance_preserved"]
        ),
        "source_uniform_same_history_joint_sector_closed": joint_sector_closed,
        "source_uniform_full_augmented_LDLT_closed": joint_ldlt_closed,

        "endpoint_transport_qualification": endpoint["qualification"],
        "endpoint_master_object_emitted": endpoint_master_emitted,
        "endpoint_decomposition_identity": endpoint["endpoint_decomposition_identity"],
        "master_endpoint_energy_identity": endpoint["master_endpoint_energy_identity"],
        "joint_accelerometer_suffix_operator": endpoint["accelerometer_joint_operator"],
        "actual_RS_regularization_enters_suffix_maps": bool(
            endpoint["actual_RS_regularization_enters_every_applicable_suffix"]
        ),
        "source_uniform_joint_BW_epsilon_enclosure_closed": bool(
            endpoint["source_uniform_joint_BW_epsilon_enclosure_closed"]
        ),
        "source_uniform_r_word_enclosure_closed": bool(endpoint["source_uniform_r_word_enclosure_closed"]),
        "source_uniform_master_endpoint_domination_closed": endpoint_closed,

        "finite_map_mean_value_qualification": finite["qualification"],
        "finite_map_mean_value_bridge_validated": True,
        "finite_map_not_differential_metric_replacement": bool(
            finite["finite_map_not_differential_metric_replacement"]
        ),
        "source_uniform_complete_word_generalized_Jacobian_enclosed": finite_jacobian_closed,
        "source_uniform_endpoint_finite_map_closed": finite_endpoint_closed,
        "source_uniform_prefix_gain_closed": prefix_gain_closed,
        "source_uniform_prefix_domain_retention_closed": prefix_domain_closed,

        "differential_AD_used_only_for_finite_map_enclosure": True,
        "differential_pullback_used_as_replacement_P4": False,
        "differential_finite_distance_bridge_closed": bool(
            finite_jacobian_closed and finite_endpoint_closed
        ),
        "prediction_AD_qualification": pred["qualification"],
        "event_AD_qualification": events["qualification"],
        "differential_word_qualification": word["qualification"],
        "same_source_omega_h_tau_prediction_required": bool(pred["same_complete_SEA3_omega_h_tau_required"]),
        "prediction_independent_F_forbidden": not bool(pred["independent_F_input_allowed_for_theorem"]),
        "same_P_H_R_cell_required_for_Joseph": bool(events["same_P_H_R_cell_derives_S_and_K"]),
        "independent_K_forbidden": not bool(events["independent_K_input_allowed_for_theorem"]),
        "A21_bias_projection_generalized_Jacobian_machinery_retained": projection_machinery,
        "outward_interval_AD_machinery_retained": True,

        "packet_count_remainder_budget_used": False,
        "packetwise_remainder_norm_sum_used": False,
        "standalone_eta_Rinv_budget_used": False,
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
            "materialize one source-correlated COMPLETE_SEA3_NORMAL_LIVE_WORD same-history nonlinear graph, including the "
            "A21 bias-projection generalized Jacobian, and close the full augmented joint-sector LDLT at the endpoint and "
            "every prefix; retain finite state, P/H/R, tuner schedule, actual applied R_S and all 21-state cross terms"
            if p3_pass else "repair only the frozen P3 prerequisite failure"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    motion = d.get("bounded_bias_motion_contract", {})
    f.extend(f"motion target: {x}" for x in MOTION.validate(motion))
    if d.get("primary_proof_target") != MOTION.TARGET:
        f.append("primary motion theorem target missing")
    for key in ("P4_MOTION_PASS", "P5_MOTION_MAY_START"):
        if d.get(key) is not False or d.get(key) != motion.get(key):
            f.append(f"{key} missing, detached or falsely promoted")
    f.extend(UNION.validate(d.get("SEA3_response_union", {})))
    if not UNION.all_branches_modes_closed(d.get("response_branch_conditional_P3_coverage_consumed", {})):
        f.append("P4 missing conditional P3 coverage of a response branch/mode")
    coverage = d.get("response_branch_P4_coverage", {})
    expected = UNION.branch_mode_coverage({"H18": False, "A21": False})
    if set(coverage) != {"endpoint", "prefix_gain", "prefix_retention"} or any(
            c != expected for c in coverage.values()):
        f.append("P4 response coverage missing or promoted without branch certificates")
    f.extend(f"MEMS bias: {x}" for x in BIAS.validate(d.get("mems_bias_preconditions", {})))
    if d.get("mems_bias_preconditions") != d.get("joint_sector_master_contract", {}).get("mems_bias_preconditions"):
        f.append("joint sector master detached from canonical bias preconditions")
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_P4_architecture") != ARCHITECTURE:
        f.append("canonical P4 architecture is not the paper finite-state endpoint+prefix theorem")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")

    for key in (
        "finite_state_endpoint_dissipation_required", "finite_state_prefix_gain_required",
        "prefix_chart_and_source_domain_retention_required", "source_generated_not_trajectory_fit",
        "P3_CONDITIONAL_SEA3_PASS_consumed", "P3_frozen_not_modified",
        "same_complete_SEA3_word_required", "same_frontend_tuner_covariance_history_required",
        "same_source_correlated_generalized_Jacobian_required", "all_radial_segment_Jacobians_must_be_enclosed",
        "all_due_S_updates_with_actual_applied_RS_required", "all_valid_accelerometer_updates_required",
        "all_process_Q_floor_reset_events_required", "asynchronous_vector_events_required",
        "H_to_A_rectangular_hybrid_event_required", "H_to_A_held_ba_error_retained_as_separate_forcing",
        "H_to_A_covariance_floor_retained_as_separate_metric_event", "full_prediction_F_Eaw_rows_retained",
        "full_epsilon_aw_retained", "signed_information_composition_available",
        "S_zero_nonlinear_eta_exactly_zero", "actual_RS_S_information_credited_as_favorable_term",
        "finite_reset_defect_explicit_in_signed_ledger", "endpoint_master_object_emitted",
        "full_state_joint_sector_master_available", "A21_finite_taub_full_matrix_P3_required",
        "A21_finite_taub_detectability_consumed",
        "joint_sector_all_21_state_cross_terms_retained", "joint_sector_actual_RS_provenance_preserved",
        "actual_RS_regularization_enters_suffix_maps", "finite_map_mean_value_bridge_validated",
        "finite_map_not_differential_metric_replacement", "differential_AD_used_only_for_finite_map_enclosure",
        "same_source_omega_h_tau_prediction_required", "prediction_independent_F_forbidden",
        "same_P_H_R_cell_required_for_Joseph", "independent_K_forbidden",
        "A21_bias_projection_generalized_Jacobian_machinery_retained", "outward_interval_AD_machinery_retained",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")

    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk", "source_family_replaced",
        "P3_DEPLOYMENT_PASS_consumed_as_if_closed", "source_uniform_joint_eta_reset_domination_closed",
        "source_uniform_same_history_joint_sector_closed", "source_uniform_full_augmented_LDLT_closed",
        "source_uniform_joint_BW_epsilon_enclosure_closed", "source_uniform_r_word_enclosure_closed",
        "source_uniform_master_endpoint_domination_closed", "source_uniform_complete_word_generalized_Jacobian_enclosed",
        "source_uniform_endpoint_finite_map_closed", "source_uniform_prefix_gain_closed",
        "source_uniform_prefix_domain_retention_closed", "differential_pullback_used_as_replacement_P4",
        "differential_finite_distance_bridge_closed",
        "packet_count_remainder_budget_used", "packetwise_remainder_norm_sum_used",
        "standalone_eta_Rinv_budget_used", "state_elimination_used", "a_w_Schur_final_certificate_used",
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
    if float(d.get("outer_angle_rad", 0.0)) < 0.80:
        f.append("outer Cayley sector shrank")
    if d.get("candidate_angles_deg") != CANDIDATES:
        f.append("declared P4 candidate sequence changed")
    if d.get("paper_Lyapunov_function") != "V(e,zeta)=e^T M(zeta)e":
        f.append("paper finite-state Lyapunov function missing")

    reasons = d.get("P4_CANONICAL_FAIL_REASONS", [])
    if len(reasons) != 2:
        f.append("canonical P4 must expose endpoint and prefix blockers")
    else:
        if not any("endpoint" in x.lower() for x in reasons):
            f.append("finite-state endpoint blocker missing")
        if not any("prefix" in x.lower() for x in reasons):
            f.append("finite-state prefix blocker missing")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--joint-sector-output", type=Path)
    ap.add_argument("--motion-output", type=Path)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.motion_output is not None:
        motion = dict(d["bounded_bias_motion_contract"])
        motion_failures = MOTION.validate(motion)
        motion.update(validation_pass=not motion_failures, validation_failures=motion_failures)
        args.motion_output.parent.mkdir(parents=True, exist_ok=True)
        args.motion_output.write_text(json.dumps(motion, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    if args.joint_sector_output is not None:
        joint = dict(d["joint_sector_master_contract"])
        joint_failures = JOINT.validate(joint)
        joint["validation_pass"] = not joint_failures
        joint["validation_failures"] = joint_failures
        args.joint_sector_output.parent.mkdir(parents=True, exist_ok=True)
        args.joint_sector_output.write_text(
            json.dumps(joint, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(json.dumps({
        "primary_proof_target": d["primary_proof_target"],
        "P4_MOTION_PASS": d["P4_MOTION_PASS"],
        "P5_MOTION_MAY_START": d["P5_MOTION_MAY_START"],
        "motion_fail_reasons": d["bounded_bias_motion_contract"]["P4_MOTION_FAIL_REASONS"],
        "architecture": d["canonical_P4_architecture"],
        "signed_information": d["signed_information_composition_available"],
        "signed_joint_domination": d["source_uniform_joint_eta_reset_domination_closed"],
        "finite_map_bridge": d["finite_map_mean_value_bridge_validated"],
        "endpoint_closed": d["source_uniform_endpoint_finite_map_closed"],
        "prefix_gain_closed": d["source_uniform_prefix_gain_closed"],
        "prefix_domain_closed": d["source_uniform_prefix_domain_retention_closed"],
        "actual_RS_required": d["all_due_S_updates_with_actual_applied_RS_required"],
        "P4_CANONICAL_PASS": d["P4_CANONICAL_PASS"],
        "fail_reasons": d["P4_CANONICAL_FAIL_REASONS"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
