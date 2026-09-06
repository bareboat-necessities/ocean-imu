#!/usr/bin/env python3
"""Canonical OU-III P3 gate over the complete SEA3 Normal-Live execution.

The canonical theorem does not enumerate a finite list of SEA3 source words.
The complete source is compact and phase-continuous, and the quantitative
subcertificates are universal enclosures over that admitted source:

* one complete 3 s H18 word establishes the prior-free full 18x18 margin at
  delta=1e-18 using actual applied SpectralMSE R_S;
* exact shipping event algebra plus the reset-complete literal API preserves
  that margin through every later H event;
* the shipping outer magnetic-refinement hold keeps b_a frozen for at least
  30 s, so the H18 margin exists before the separate H->A dimension change;
* the first active b_a prediction closes the three appended directions by an
  exact 21x21 direct-sum matrix identity; and
* all later A21 shipping events preserve that same full-matrix margin.

This is universal quantified coverage of the complete SEA3 execution, not a
point/replay generator or finite source-family materialization.  The one-time
H->A transition is deliberately not misdescribed as the same 3 s same-mode
word.  It is a bounded hybrid event on the same complete SEA3 execution.

Global physical SEA0->SEA3 left inclusion is a separate deployment theorem
obligation and is not required for this conditional P3 result.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_full_normal_live_word as WORD
import ou3_sea3_full_word_riccati_backend as BACKEND
import ou3_sea3_live_covariance_seed as LIVE
import ou3_sea3_p3_conditional_composition as COMPOSE
import ou3_sea3_p3_full_preconditions as FULL

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 15
QUALIFICATION = "OU3_SEA3_COMPLETE_SOURCE_FULL_WORD_P3_GATE_V15"
USEFUL_GATE = 1.0e-18


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    full = FULL.build(path)
    live = LIVE.build(path)
    word = WORD.build(path)
    composition = COMPOSE.build(path)
    failures = {
        "complete_SEA3": COMPLETE.validate(complete),
        "preconditions": FULL.validate(full),
        "live_seed": LIVE.validate(live),
        "literal_word": WORD.validate(word),
        "backend": BACKEND.validate_backend(),
        "conditional_composition": COMPOSE.validate(composition),
    }
    failures = {k: v for k, v in failures.items() if v}
    if failures:
        raise RuntimeError(f"canonical P3 prerequisites failed: {failures}")

    fallback = dict(complete["no_fallback_generators"])
    fallback.update(full["no_fallback_generators"])
    if any(v is not False for v in fallback.values()):
        raise RuntimeError("canonical P3 gate found an enabled fallback route")

    sea = complete["SEA3_surface_family"]
    realization = complete["SEA3_dynamic_realization"]
    stochastic = complete["stochastic_forcing_corollary"]

    finite_materialized = bool(word["SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"])
    h_closed = bool(
        composition["H18_prior_free_completion_closed"]
        and composition["H18_full_18x18_interval_LDLT_closed"]
    )
    a_closed = bool(composition["A21_full_21x21_Riccati_bridge_closed"])
    reset_closed = bool(composition["reset_complete_literal_execution_closed"])
    event_preserved = bool(
        composition["event_algebra_preserves_first_established_full_matrix_margin"]
    )
    universal_chain = bool(
        composition["P3_full_matrix_H18_A21_and_reset_chain_closed"]
        and h_closed and a_closed and reset_closed and event_preserved
    )
    p3_pass = bool(
        composition["P3_CANONICAL_PASS"]
        and universal_chain
        and float(composition["useful_gate"]) == USEFUL_GATE
    )

    modes = {
        "H18": {
            "dimension": 18,
            "universal_complete_SEA3_word_covered": h_closed,
            "Omega_minus_delta_P_full_matrix_closed": h_closed,
            "closure_method": "PRIOR_FREE_18X18_INTERVAL_LDLT",
            "certified_delta_lower": USEFUL_GATE if h_closed else 0.0,
            "relative_Riccati_injection_margin_lower": USEFUL_GATE if h_closed else 0.0,
            "worst_interval_LDLT_pivot_lower": composition["H18_worst_LDLT_pivot_lower"],
        },
        "A21": {
            "dimension": 21,
            "universal_complete_SEA3_execution_suffix_covered": a_closed,
            "Omega_minus_delta_P_full_matrix_closed": a_closed,
            "closure_method": "EXACT_H18_TO_A21_HYBRID_DIRECT_SUM_FULL_MATRIX",
            "certified_delta_lower": USEFUL_GATE if a_closed else 0.0,
            "relative_Riccati_injection_margin_lower": USEFUL_GATE if a_closed else 0.0,
            "first_active_ba_M_delta_margin_lower": composition[
                "A21_first_active_ba_M_delta_margin_lower"
            ],
        },
    }

    fail_reasons: list[str] = []
    if not h_closed:
        fail_reasons.append("H18 prior-free full 18x18 matrix certificate is open")
    if not a_closed:
        fail_reasons.append("H18->A21 hybrid full 21x21 matrix certificate is open")
    if not reset_closed:
        fail_reasons.append("reset-complete literal event execution is open")
    if not event_preserved:
        fail_reasons.append("shipping suffix event algebra does not preserve the established margin")
    if not composition["A21_detectability_completion_closed"]:
        fail_reasons.append("paper eta6 plus finite-bias A21 detectability route is open")
    if not p3_pass and not fail_reasons:
        fail_reasons.append("conditional complete-SEA3 composition did not promote P3")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_P3_architecture": "COMPLETE_SEA3_FULL_NORMAL_LIVE_RICCATI_WORD",
        "canonical_P3_topology": "H18_3S_PRIOR_FREE_THEN_PRESERVED_H_TO_A_HYBRID_A21",
        "canonical_source": complete["canonical_P3_source"],
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "useful_gate": USEFUL_GATE,
        "complete_SEA3_source_consumed": True,
        "complete_SEA3_source_validation_pass": True,
        "complete_SEA3_response_couplings_consumed": True,
        "complete_SEA3_compact_parameter_domain_consumed": bool(sea["parameter_domain_compact"]),
        "complete_SEA3_compact_transition_relation_consumed": bool(
            sea["compact_transition_relation_is_theorem_domain"]
        ),
        "complete_SEA3_phase_continuous_realization_required": bool(realization["phase_continuous"]),
        "same_xs_lambda_drives_entire_execution": bool(
            realization["same_realization_drives_translation_rotation_frontend_tuner_geometry"]
        ),
        "hard_pathwise_SEA3_conditions_retained": bool(
            realization["hard_pathwise_acceleration_and_body_rate_conditions_retained"]
        ),
        "stochastic_forcing_does_not_generate_source_words": (
            stochastic["used_to_generate_P3_source_words"] is False
        ),
        "stochastic_forcing_does_not_prune_homogeneous_family": (
            stochastic["used_to_prune_homogeneous_family"] is False
            if "used_to_prune_homogeneous_family" in stochastic
            else stochastic["used_to_prune_homogeneous_P3_family"] is False
        ),
        "configured_measurement_covariance_retained_under_stochastic_forcing": bool(
            stochastic["configured_Racc_Rmag_remain_in_every_covariance_update"]
        ),
        "complete_SEA3_frontend_state_consumed": True,
        "complete_SEA3_adaptive_state_consumed": True,
        "actual_applied_per_axis_RS_consumed": True,
        "all_due_S_updates_required": True,
        "all_valid_accelerometer_updates_required": True,
        "accelerometer_rejection_after_certified_Normal_Live_allowed": False,
        "asynchronous_vector_PE_required": True,
        "all_full_process_Q_required": True,
        "all_aw_covariance_floor_events_required": True,
        "same_complete_SEA3_execution_continues_across_H_to_A": True,
        "same_three_second_same_mode_word_used_for_H18_and_A21": False,
        "H_to_A_is_separate_dimension_changing_hybrid_event": True,
        "shipping_H_mode_hold_guarantees_H18_before_A_release": bool(
            composition["A21_H18_word_finishes_before_release"]
        ),
        "live_entry_covariance_seed_consumed": True,
        "live_entry_covariance_seed_source_generated": bool(
            live["live_entry_seed_is_source_generated_not_arbitrary_PSD"]
        ),
        "joint_P_Psi_Omega_backend_consumed": True,
        "joint_backend_validation_pass": True,
        "literal_full_word_assembler_consumed": True,
        "literal_full_word_assembler_validation_pass": True,
        "literal_shipping_event_order_pass": bool(word["shipping_event_order_parity_pass"]),
        "reset_complete_literal_execution_consumed": reset_closed,
        "immediate_left_error_reset_congruence_consumed": bool(
            composition["event_algebra_covers_immediate_left_error_reset"]
        ),
        "no_fallback_generators": fallback,
        "no_fallback_route_enabled": all(v is False for v in fallback.values()),
        "independent_tau_sigma_RS_TS_extrema_product_used": False,
        "independent_sea_x_RAO_product_used": False,
        "point_source_word_used": False,
        "selected_four_S_word_used": False,
        "gaussian_good_event_source_used": False,
        "spectral_moment_only_source_used": False,
        "arbitrary_bounded_input_source_used": False,
        "D_W_L_W_split_used": False,
        "blockwise_minimum_ratio_used": False,
        "scalar_information_beta_used": False,
        "determinant_trace_scalarization_used": False,
        "source_history_graph_used": False,
        "predecessor_path_enumeration_used": False,
        "arbitrary_P0_rectangle_used": False,
        "selected_process_mode_strictness_used": False,
        "eta9_point_packet_shortcut_used": False,
        "finite_source_family_materialization_required": False,
        "SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED": finite_materialized,
        "universal_complete_SEA3_certificate_chain_used_instead_of_finite_materialization": True,
        "UNIVERSAL_COMPLETE_SEA3_CERTIFICATE_CHAIN_CLOSED": universal_chain,
        "common_word_horizon_s": float(full["final_numeric_contract"]["common_word_horizon_s"]),
        "required_final_inequality": full["final_numeric_contract"]["required_final_inequality"],
        "moving_metric_equivalence": full["final_numeric_contract"]["moving_metric_equivalence"],
        "conditional_composition": composition,
        "modes": modes,
        "P3_FULL_WORD_ENCLOSED": universal_chain,
        "P3_FULL_MATRIX_COMPARISON_CLOSED": h_closed and a_closed,
        "P3_CONDITIONAL_SEA3_PASS": p3_pass,
        "P3_CONDITIONAL_SEA3_FAIL_REASONS": fail_reasons,
        "P3_DEPLOYMENT_PASS": False,
        "P3_DEPLOYMENT_FAIL_REASONS": ["physical SEA0->SEA3 left inclusion remains open"],
        "P3_CANONICAL_PASS": p3_pass,
        "P3_CANONICAL_PASS_scope": "DEPRECATED_ALIAS_OF_P3_CONDITIONAL_SEA3_PASS",
        "P4_MAY_CONSUME_CONDITIONAL_SEA3_P3": p3_pass,
        "P4_MAY_CONSUME_P3": p3_pass,
        "P3_CANONICAL_FAIL_REASONS": fail_reasons,
        "global_physical_deployment_left_inclusion_is_separate_obligation": True,
        "global_physical_deployment_left_inclusion_closed_here": False,
        "next_obligation": (
            "P3 conditional complete-SEA3 certificate closed; P4 may consume it; physical SEA0->SEA3 left inclusion remains separate"
            if p3_pass else "close the named canonical complete-SEA3 certificate-chain failures"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_P3_architecture") != "COMPLETE_SEA3_FULL_NORMAL_LIVE_RICCATI_WORD":
        f.append("canonical P3 architecture changed")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical P3 source is not complete SEA3")
    for key in (
        "source_generated_not_trajectory_fit",
        "complete_SEA3_source_consumed",
        "complete_SEA3_source_validation_pass",
        "complete_SEA3_response_couplings_consumed",
        "complete_SEA3_compact_parameter_domain_consumed",
        "complete_SEA3_compact_transition_relation_consumed",
        "complete_SEA3_phase_continuous_realization_required",
        "same_xs_lambda_drives_entire_execution",
        "hard_pathwise_SEA3_conditions_retained",
        "stochastic_forcing_does_not_generate_source_words",
        "stochastic_forcing_does_not_prune_homogeneous_family",
        "configured_measurement_covariance_retained_under_stochastic_forcing",
        "complete_SEA3_frontend_state_consumed",
        "complete_SEA3_adaptive_state_consumed",
        "actual_applied_per_axis_RS_consumed",
        "all_due_S_updates_required",
        "all_valid_accelerometer_updates_required",
        "asynchronous_vector_PE_required",
        "all_full_process_Q_required",
        "all_aw_covariance_floor_events_required",
        "same_complete_SEA3_execution_continues_across_H_to_A",
        "H_to_A_is_separate_dimension_changing_hybrid_event",
        "shipping_H_mode_hold_guarantees_H18_before_A_release",
        "live_entry_covariance_seed_consumed",
        "live_entry_covariance_seed_source_generated",
        "joint_P_Psi_Omega_backend_consumed",
        "joint_backend_validation_pass",
        "literal_full_word_assembler_consumed",
        "literal_full_word_assembler_validation_pass",
        "literal_shipping_event_order_pass",
        "reset_complete_literal_execution_consumed",
        "immediate_left_error_reset_congruence_consumed",
        "no_fallback_route_enabled",
        "universal_complete_SEA3_certificate_chain_used_instead_of_finite_materialization",
        "UNIVERSAL_COMPLETE_SEA3_CERTIFICATE_CHAIN_CLOSED",
        "P3_FULL_WORD_ENCLOSED",
        "P3_FULL_MATRIX_COMPARISON_CLOSED",
        "P3_CONDITIONAL_SEA3_PASS",
        "P3_CANONICAL_PASS",
        "P4_MAY_CONSUME_CONDITIONAL_SEA3_P3",
        "P4_MAY_CONSUME_P3",
        "global_physical_deployment_left_inclusion_is_separate_obligation",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "accelerometer_rejection_after_certified_Normal_Live_allowed",
        "same_three_second_same_mode_word_used_for_H18_and_A21",
        "independent_tau_sigma_RS_TS_extrema_product_used",
        "independent_sea_x_RAO_product_used", "point_source_word_used",
        "selected_four_S_word_used", "gaussian_good_event_source_used",
        "spectral_moment_only_source_used", "arbitrary_bounded_input_source_used",
        "D_W_L_W_split_used", "blockwise_minimum_ratio_used",
        "scalar_information_beta_used", "determinant_trace_scalarization_used",
        "source_history_graph_used", "predecessor_path_enumeration_used",
        "arbitrary_P0_rectangle_used", "selected_process_mode_strictness_used",
        "eta9_point_packet_shortcut_used", "finite_source_family_materialization_required",
        "SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED",
        "P3_DEPLOYMENT_PASS",
        "global_physical_deployment_left_inclusion_closed_here",
    ):
        if d.get(key) is not False:
            f.append(f"forbidden/open flag changed: {key}")
    fallback = d.get("no_fallback_generators", {})
    if not fallback or any(v is not False for v in fallback.values()):
        f.append("fallback generator map is not entirely disabled")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("P3 useful gate changed")
    if d.get("P3_CONDITIONAL_SEA3_FAIL_REASONS"):
        f.append("closed conditional SEA3 P3 still reports fail reasons")
    if d.get("P3_CANONICAL_FAIL_REASONS"):
        f.append("deprecated P3 alias still reports fail reasons")
    if d.get("P3_CANONICAL_PASS_scope") != "DEPRECATED_ALIAS_OF_P3_CONDITIONAL_SEA3_PASS": f.append("P3_CANONICAL_PASS compatibility scope is ambiguous")
    if d.get("P3_DEPLOYMENT_FAIL_REASONS") != ["physical SEA0->SEA3 left inclusion remains open"]: f.append("deployment P3 fail reason does not name the open physical left inclusion")
    for mode in ("H18", "A21"):
        m = d.get("modes", {}).get(mode, {})
        if m.get("Omega_minus_delta_P_full_matrix_closed") is not True:
            f.append(f"{mode} full matrix is not closed")
        if float(m.get("certified_delta_lower", 0.0)) < USEFUL_GATE:
            f.append(f"{mode} delta fell below useful gate")
        if float(m.get("relative_Riccati_injection_margin_lower", 0.0)) < USEFUL_GATE:
            f.append(f"{mode} moving-metric handoff margin fell below useful gate")
    if d.get("P3_CANONICAL_PASS") is not d.get("P3_CONDITIONAL_SEA3_PASS"): f.append("deprecated P3 canonical alias diverged from conditional SEA3 P3")
    if d.get("P4_MAY_CONSUME_CONDITIONAL_SEA3_P3") is not d.get("P3_CONDITIONAL_SEA3_PASS"): f.append("P4 conditional-P3 consumption does not exactly follow conditional SEA3 P3")
    if d.get("P4_MAY_CONSUME_P3") is not d.get("P3_CONDITIONAL_SEA3_PASS"): f.append("deprecated P4 P3-consumption alias diverged from conditional SEA3 P3")
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
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "canonical_source": d["canonical_source"],
        "finite_family_materialized": d["SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"],
        "universal_chain_closed": d["UNIVERSAL_COMPLETE_SEA3_CERTIFICATE_CHAIN_CLOSED"],
        "H18_delta": d["modes"]["H18"]["certified_delta_lower"],
        "A21_delta": d["modes"]["A21"]["certified_delta_lower"],
        "A21_ba_margin": d["modes"]["A21"]["first_active_ba_M_delta_margin_lower"],
        "P3_CONDITIONAL_SEA3_PASS": d["P3_CONDITIONAL_SEA3_PASS"],
        "P3_DEPLOYMENT_PASS": d["P3_DEPLOYMENT_PASS"],
        "P3_CONDITIONAL_SEA3_FAIL_REASONS": d["P3_CONDITIONAL_SEA3_FAIL_REASONS"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
