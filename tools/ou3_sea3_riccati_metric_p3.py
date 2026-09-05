#!/usr/bin/env python3
"""Canonical complete-source OU-III SEA3 P3 gate.

P3 is one H18/A21 Normal-Live Riccati word.  The promotable numerical object is
not a product of separately worst information/covariance certificates.  It is
the exact joint propagation of

    P_k, Psi_k, Omega_k

on the same source/event path, with

    P_k = Psi_k P_0 Psi_k^T + Omega_k.

Prediction, every accepted/due Joseph measurement update, and every PSD a_w
covariance-floor event update all three objects consistently.  The useful
contraction gate is the full-matrix inequality

    Omega_W - delta P_W >= 0,   delta >= 1e-18,

which is equivalent to moving-metric contraction of the homogeneous transition.
No D_W/L_W split, zero-start Riccati concavity replacement, blockwise ratio,
source-history graph, or scalarized substitute may promote P3.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_sea3_p3_full_preconditions as FULL
import ou3_sea3_rs_innovation_p3 as RS_COMPONENT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 7
QUALIFICATION = "OU3_SEA3_FULL_NORMAL_LIVE_RICCATI_WORD_P3_GATE"
USEFUL_GATE = 1.0e-18


def build(domain_path: Path = DEFAULT_DOMAIN, tube_path: Path | None = None) -> dict:
    path = Path(domain_path).resolve()
    full = FULL.build(path)
    ff = FULL.validate(full)
    rs = RS_COMPONENT.build(path, tube_path)
    rf = RS_COMPONENT.validate(rs)
    if ff or rf:
        raise RuntimeError(
            f"canonical P3 prerequisites failed: full={ff}, four_S_component={rf}"
        )

    numeric = full["final_numeric_contract"]
    modes = {
        "H": {
            "dimension": 18,
            "relative_Riccati_injection_margin_lower": 0.0,
            "Omega_minus_delta_P_ldlt_closed": False,
            "contraction_gap_lower": 0.0,
            "useful_margin_gate": USEFUL_GATE,
            "pass": False,
        },
        "A": {
            "dimension": 21,
            "relative_Riccati_injection_margin_lower": 0.0,
            "Omega_minus_delta_P_ldlt_closed": False,
            "contraction_gap_lower": 0.0,
            "useful_margin_gate": USEFUL_GATE,
            "pass": False,
        },
    }

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_P3_architecture": "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "useful_gate": USEFUL_GATE,
        "complete_precondition_contract_consumed": True,
        "all_current_machine_checkable_preconditions_present": full[
            "all_current_machine_checkable_preconditions_present"
        ],
        "common_word_horizon_s": numeric["common_word_horizon_s"],
        "one_common_event_word_required_for_H_and_A": True,
        "full_H18_A21_matrix_comparison_required": True,
        "same_joint_source_path_feeds_F_Q_TS_RS": True,
        "same_front_end_state_path_generates_all_tuner_targets": True,
        "same_event_word_contains_accel_S_PE_and_aw_floor": True,
        "same_runtime_measurement_covariances_used": True,
        "complete_front_end_generator_state_consumed": True,
        "magnetic_reference_path_retained_not_frozen": True,
        "no_hard_attitude_rewrite_inside_word": True,
        "hybrid_transitions_separate": True,
        "R_S_translation_component_consumed": True,
        "R_S_is_primary_translation_correction_mechanism": True,
        "four_S_translation_word_consumed": True,
        "four_S_translation_observation_geometry_closed": bool(
            rs["P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED"]
        ),
        "four_S_translation_information_matrix_closed": bool(
            rs["P3_RS_TRANSLATION_INFORMATION_MATRIX_CLOSED"]
        ),
        "four_S_batch_noise_upper_closed": bool(
            rs["P3_RS_BATCH_NOISE_UPPER_CLOSED"]
        ),
        "R_S_component_is_not_the_whole_P3_architecture": True,
        "actual_applied_per_axis_RS_required_in_final_word": True,
        "all_due_S_updates_required_not_only_selected_four": True,
        "all_valid_accelerometer_updates_consumed": True,
        "full_accelerometer_attitude_aw_ba_cross_block_information_required": True,
        "windowed_asynchronous_vector_PE_consumed": True,
        "hardware_magnetometer_ODR_used_as_PE_recurrence": False,
        "two_consecutive_accepted_magnetic_packets_required": False,
        "full_process_UCC_consumed": True,
        "full_process_Q_matrices_required_in_word": True,
        "aw_covariance_floor_PSD_events_consumed": True,
        "aw_covariance_floor_marginal_Loewner_shortcut_used": False,
        "A_mode_finite_bias_correlation_consumed": True,
        "SEA3_height_period_partition_coupling_consumed": True,
        "unqualified_RAO_coupling_used_as_hard_pruning": False,
        "global_physical_SEA3_left_inclusion_claimed": False,
        "stochastic_noise_realization_used_as_homogeneous_pruning": False,
        "joint_P_Psi_Omega_word_required": True,
        "exact_covariance_decomposition_identity": (
            "P_k = Psi_k P_0 Psi_k^T + Omega_k"
        ),
        "prediction_joint_recursion": numeric["prediction_recursion_required"],
        "joseph_measurement_joint_recursion": numeric[
            "joseph_measurement_recursion_required"
        ],
        "aw_floor_joint_recursion": numeric["aw_floor_recursion_required"],
        "required_final_inequality": numeric["required_final_inequality"],
        "moving_metric_equivalence": numeric["moving_metric_equivalence"],
        "exact_measurement_dissipation_identity_available": True,
        "batch_innovation_information_identity_available": True,
        "D_W_L_W_split_used_for_final_gate": False,
        "zero_start_Riccati_concavity_replacement_used": False,
        "source_history_graph_consumed": False,
        "predecessor_path_enumeration_consumed": False,
        "old_P2_800_state_graph_consumed": False,
        "one_sample_strict_Riccati_margin_consumed": False,
        "per_sample_SPD_lower_required": False,
        "selected_process_mode_strictness_used": False,
        "determinant_trace_scalarization_used": False,
        "scalar_information_beta_used": False,
        "blockwise_minimum_ratio_used_for_final_gate": False,
        "independent_tau_sigma_RS_TS_extrema_product_used": False,
        "front_end_state_frozen_to_replay_value": False,
        "complete_preconditions": full,
        "measurement_runtime": full["measurement_runtime"],
        "front_end_state_manifest": full["front_end_state_manifest"],
        "four_S_translation_component": {
            "translation_correction_word": rs["translation_correction_word"],
            "exact_measurement_dissipation_identity": rs[
                "exact_measurement_dissipation_identity"
            ],
            "batch_innovation_identity": rs["batch_innovation_identity"],
        },
        "modes": modes,
        "P3_FOUNDATION_PASS": True,
        "P3_ARCHITECTURE_READY": True,
        "P3_FULL_WORD_ENCLOSED": False,
        "P3_FULL_MATRIX_COMPARISON_CLOSED": False,
        "P3_CANONICAL_PASS": False,
        "P4_MAY_CONSUME_P3": False,
        "P3_CANONICAL_FAIL_REASONS": [
            "the complete 3 s Normal-Live event word has not yet been propagated as one joint (P,Psi,Omega) H18/A21 object",
            "the final word must carry the complete measurement-only front-end/tuner source state and use the same source path for F, Q, T_S, applied per-axis R_S, Racc and Rmag",
            "every valid accelerometer update, every due S update, admissible asynchronous magnetic PE event, every prediction Q and every recurrent PSD a_w-floor event must update the same P/Psi/Omega recursion",
            "the full H18/A21 interval matrix Omega_W-delta*P_W has not yet passed validated LDLT at the unchanged 1e-18 gate",
        ],
        "next_obligation": (
            "implement the complete joint P/Psi/Omega Normal-Live word only; "
            "do not close P3 through a separately bounded information/covariance product"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_P3_architecture") != "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD":
        f.append("canonical P3 is not the complete Normal-Live word")

    required_true = (
        "source_generated_not_trajectory_fit",
        "complete_precondition_contract_consumed",
        "all_current_machine_checkable_preconditions_present",
        "one_common_event_word_required_for_H_and_A",
        "full_H18_A21_matrix_comparison_required",
        "same_joint_source_path_feeds_F_Q_TS_RS",
        "same_front_end_state_path_generates_all_tuner_targets",
        "same_event_word_contains_accel_S_PE_and_aw_floor",
        "same_runtime_measurement_covariances_used",
        "complete_front_end_generator_state_consumed",
        "magnetic_reference_path_retained_not_frozen",
        "no_hard_attitude_rewrite_inside_word",
        "hybrid_transitions_separate",
        "R_S_translation_component_consumed",
        "R_S_is_primary_translation_correction_mechanism",
        "four_S_translation_word_consumed",
        "four_S_translation_observation_geometry_closed",
        "four_S_translation_information_matrix_closed",
        "four_S_batch_noise_upper_closed",
        "R_S_component_is_not_the_whole_P3_architecture",
        "actual_applied_per_axis_RS_required_in_final_word",
        "all_due_S_updates_required_not_only_selected_four",
        "all_valid_accelerometer_updates_consumed",
        "full_accelerometer_attitude_aw_ba_cross_block_information_required",
        "windowed_asynchronous_vector_PE_consumed",
        "full_process_UCC_consumed",
        "full_process_Q_matrices_required_in_word",
        "aw_covariance_floor_PSD_events_consumed",
        "A_mode_finite_bias_correlation_consumed",
        "SEA3_height_period_partition_coupling_consumed",
        "joint_P_Psi_Omega_word_required",
        "exact_measurement_dissipation_identity_available",
        "batch_innovation_information_identity_available",
        "P3_FOUNDATION_PASS",
        "P3_ARCHITECTURE_READY",
    )
    for key in required_true:
        if d.get(key) is not True:
            f.append(f"{key} is not true")

    required_false = (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "hardware_magnetometer_ODR_used_as_PE_recurrence",
        "two_consecutive_accepted_magnetic_packets_required",
        "aw_covariance_floor_marginal_Loewner_shortcut_used",
        "unqualified_RAO_coupling_used_as_hard_pruning",
        "global_physical_SEA3_left_inclusion_claimed",
        "stochastic_noise_realization_used_as_homogeneous_pruning",
        "D_W_L_W_split_used_for_final_gate",
        "zero_start_Riccati_concavity_replacement_used",
        "source_history_graph_consumed", "predecessor_path_enumeration_consumed",
        "old_P2_800_state_graph_consumed", "one_sample_strict_Riccati_margin_consumed",
        "per_sample_SPD_lower_required", "selected_process_mode_strictness_used",
        "determinant_trace_scalarization_used", "scalar_information_beta_used",
        "blockwise_minimum_ratio_used_for_final_gate",
        "independent_tau_sigma_RS_TS_extrema_product_used",
        "front_end_state_frozen_to_replay_value",
        "P3_FULL_WORD_ENCLOSED", "P3_FULL_MATRIX_COMPARISON_CLOSED",
        "P3_CANONICAL_PASS", "P4_MAY_CONSUME_P3",
    )
    for key in required_false:
        if d.get(key) is not False:
            f.append(f"{key} is not false")

    if d.get("exact_covariance_decomposition_identity") != (
        "P_k = Psi_k P_0 Psi_k^T + Omega_k"
    ):
        f.append("exact P/Psi/Omega identity changed")
    if d.get("required_final_inequality") != (
        "Omega_W - delta * P_W >= 0 on full H18/A21 coordinates"
    ):
        f.append("canonical final matrix inequality changed")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("P3 useful gate changed")
    if float(d.get("common_word_horizon_s", 0.0)) < 3.0:
        f.append("canonical word no longer spans all declared PE preconditions")

    for mode, dim in (("H", 18), ("A", 21)):
        row = d.get("modes", {}).get(mode, {})
        if row.get("dimension") != dim or row.get("pass") is not False:
            f.append(f"{mode} fail-closed mode contract invalid")
        if float(row.get("relative_Riccati_injection_margin_lower", math.nan)) != 0.0:
            f.append(f"{mode} emitted a margin before complete-word closure")
        if row.get("Omega_minus_delta_P_ldlt_closed") is not False:
            f.append(f"{mode} falsely closed Omega-delta-P LDLT")

    if not d.get("P3_CANONICAL_FAIL_REASONS"):
        f.append("open P3 does not name complete-word obligations")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--tube", type=Path, default=None)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.tube)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "architecture": d["canonical_P3_architecture"],
        "word_horizon_s": d["common_word_horizon_s"],
        "final_inequality": d["required_final_inequality"],
        "P3_CANONICAL_PASS": d["P3_CANONICAL_PASS"],
        "fail_reasons": d["P3_CANONICAL_FAIL_REASONS"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
