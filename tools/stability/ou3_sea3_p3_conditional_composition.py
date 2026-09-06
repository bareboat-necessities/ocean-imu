#!/usr/bin/env python3
"""Conditional complete-SEA3 P3 composition contract.

The theorem keeps two quantified statements separate:

1. Conditional P3: every admitted complete phase-continuous SEA3 Normal-Live
   word satisfies the shipping H18/A21 Riccati stability inequality at the
   useful gate delta=1e-18.
2. Physical SEA0 -> SEA3 left inclusion: the intended physical/random-sea
   realization family is contained in that admitted deterministic SEA3 set.

The second statement remains a separate deployment obligation and is not used
to manufacture or prune words for the first statement.

The conditional P3 chain is theorem-equivalent to explicit enumeration of a
finite source list and does not require such a list to be materialized:

* complete SEA3 provides the universal compact phase-continuous source;
* H18 uses actual applied SpectralMSE R_S, all due S updates, asynchronous eta6
  PE and the exact prior-free 18x18 interval-LDLT completion;
* the paper's eta6+finite-tau_b A21 detectability/UES route remains independently
  certified;
* shipping holds b_a long enough that H18 is established before the separate
  H->A dimension change, after which the first active GM prediction closes the
  three new directions by an exact 21x21 direct-sum matrix identity;
* the reset-complete literal API applies the exact left-error covariance
  congruence immediately after every accepted S/accelerometer/magnetometer
  Joseph correction, and the event algebra preserves the established margin.

No eta9 packet shortcut, alternate estimator, replay word, independent
parameter box, blockwise contraction ratio, scalar beta, or one-step full-state
Q minimum is a promotable route.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_sea3_a21_detectability_completion as ADET
import ou3_sea3_a21_prior_free_completion as A21PF
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_four_s_translation_information as FOUR_S
import ou3_sea3_full_normal_live_word_reset as RESETWORD
import ou3_sea3_full_word_event_algebra as EVENT
import ou3_sea3_h18_information_composition as H18INFO
import ou3_sea3_h18_prior_free_completion as H18PF
import ou3_sea3_hard_shaping_state as SHAPING
import ou3_sea3_windowed_vector_pe as PE

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 3
QUALIFICATION = "OU3_COMPLETE_SEA3_CONDITIONAL_P3_COMPOSITION_V3"
USEFUL_GATE = 1.0e-18


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    four = FOUR_S.build(path)
    event = EVENT.build()
    h18 = H18INFO.build(path)
    h18pf = H18PF.build(path)
    adet = ADET.build(path)
    a21pf = A21PF.build(path)
    resetword = RESETWORD.build(path)
    shaping = SHAPING.build()
    pe = PE.build(path)
    process = PROCESS.build()
    bad = {
        "complete": COMPLETE.validate(complete),
        "four_S": FOUR_S.validate(four),
        "event_algebra": EVENT.validate(event),
        "H18_information": H18INFO.validate(h18),
        "H18_prior_free": H18PF.validate(h18pf),
        "A21_detectability": ADET.validate(adet),
        "A21_hybrid_full_matrix": A21PF.validate(a21pf),
        "reset_complete_literal_word": RESETWORD.validate(resetword),
        "hard_shaping_state": SHAPING.validate(shaping),
        "PE": PE.validate(pe),
        "process": PROCESS.validate(process),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"conditional complete-SEA3 P3 prerequisites failed: {bad}")

    source = complete["canonical_P3_source"]
    if source != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("conditional P3 source changed")
    for name, value in (
        ("H18 prior-free", h18pf["canonical_source"]),
        ("A21 detectability", adet["canonical_source"]),
        ("A21 hybrid", a21pf["canonical_source"]),
        ("reset-complete literal word", resetword["canonical_source"]),
    ):
        if value != source:
            raise RuntimeError(f"{name} detached from canonical complete SEA3")

    h18_info = float(h18["triangular_information_composition"]["D_H18_lambda_min_lower"])
    h_q = float(process["modes"]["H"]["prediction_Q_lambda_min_lower"])
    a_q = float(process["modes"]["A"]["prediction_Q_lambda_min_lower"])
    ba_gap = float(pe["A_mode_bias_route"]["homogeneous_bias_contraction_gap_lower"])
    a21_gap = float(adet["A21_detectability_asymptotic_word_energy_gap_lower"])
    hard_left_inclusion = bool(shaping["executable_ingredients"]["complete_SEA3_left_inclusion_closed"])

    reset = event["left_error_reset"]
    preservation = event["full_matrix_margin_preservation"]
    event_preserved = bool(
        preservation["covers_prediction"]
        and preservation["covers_every_due_S_update"]
        and preservation["covers_every_Normal_Live_accelerometer_update"]
        and preservation["covers_asynchronous_magnetometer_update"]
        and preservation["covers_immediate_left_error_reset"]
        and preservation["covers_aw_covariance_floor"]
        and preservation["covers_not_due_or_rejected_identity_branches"]
    )
    reset_complete = bool(
        resetword["literal_reset_execution_complete"]
        and resetword["shipping_reset_source_parity_pass"]
        and resetword["S_Joseph_immediately_followed_by_left_reset"]
        and resetword["accelerometer_Joseph_immediately_followed_by_left_reset"]
        and resetword["magnetometer_Joseph_immediately_followed_by_left_reset"]
    )
    h18_closed = bool(h18pf["H18_prior_free_completion_closed"])
    a21_detectability_closed = bool(adet["A21_finite_bias_detectability_closed"])
    a21_full_riccati_bridge = bool(
        a21pf["A21_prior_free_completion_closed"]
        and a21pf["full_21x21_Omega_minus_delta_P_closed"]
    )
    conditional_pass = bool(
        h18_closed
        and a21_detectability_closed
        and a21_full_riccati_bridge
        and reset_complete
        and event_preserved
        and float(h18pf["useful_gate"]) == USEFUL_GATE
        and float(a21pf["useful_gate"]) == USEFUL_GATE
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": source,
        "theorem_conditional_on_admitted_complete_SEA3_word": bool(
            complete["theorem_conditional_on_admitted_complete_SEA3_word"]
        ),
        "conditional_P3_quantifier": "FOR_EVERY_ADMITTED_COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "finite_source_family_materialization_required_for_conditional_theorem": False,
        "global_physical_deployment_left_inclusion_is_separate_obligation": True,
        "global_physical_deployment_left_inclusion_closed_here": bool(
            complete["global_physical_deployment_left_inclusion_closed_here"]
        ),
        "global_physical_left_inclusion_required_before_conditional_P3_math": False,
        "global_physical_left_inclusion_required_before_full_deployment_theorem": True,
        "hard_shaping_physical_left_inclusion_currently_closed": hard_left_inclusion,
        "source_family_replaced": False,
        "gaussian_good_event_used_as_source": False,
        "spectral_moments_used_as_source": False,
        "trajectory_replay_used_as_source": False,
        "independent_bounded_input_box_used_as_source": False,
        "four_S_selected_events_replace_complete_word": False,
        "four_S_information_component_consumed": True,
        "actual_applied_SpectralMSE_R_S_consumed": True,
        "all_due_S_updates_remain_in_complete_word": True,
        "all_valid_accelerometer_updates_required": True,
        "accelerometer_rejection_after_certified_Normal_Live_allowed": False,
        "H18_full_information_matrix_lower_closed": True,
        "H18_information_lambda_min_lower": h18_info,
        "H18_prior_free_completion_closed": h18_closed,
        "H18_full_18x18_interval_LDLT_closed": bool(
            h18pf["full_H18_prior_free_matrix_condition_closed"]
        ),
        "H18_worst_LDLT_pivot_lower": h18pf["worst_full_H18_LDLT_pivot_lower"],
        "H18_delta_squared_completion_penalty": h18pf["delta_squared_over_four_penalty_physical"],
        "eta6_information_lambda_min_lower": float(
            pe["eta6_information"]["alpha_6_information_lower"]
        ),
        "A21_uses_eta9_packet_shortcut": False,
        "A21_finite_bias_correlation_route_consumed": True,
        "A21_bias_homogeneous_contraction_gap_lower": ba_gap,
        "A21_detectability_completion_closed": a21_detectability_closed,
        "A21_paper_UES_hypotheses_closed": bool(adet["A21_paper_UES_hypotheses_closed"]),
        "A21_detectability_asymptotic_word_energy_gap_lower": a21_gap,
        "A21_hybrid_release_full_matrix_consumed": True,
        "A21_H18_word_finishes_before_release": bool(a21pf["H18_word_finishes_before_A_release"]),
        "A21_first_active_ba_M_delta_margin_lower": a21pf["first_active_ba_M_delta_margin_lower"],
        "A21_exact_direct_sum_full_matrix_hybrid_proof_used": bool(
            a21pf["exact_direct_sum_full_matrix_hybrid_proof_used"]
        ),
        "A21_full_21x21_Riccati_bridge_closed": a21_full_riccati_bridge,
        "A21_comparison_observer_is_proof_only_not_alternate_estimator": bool(
            adet["triangular_detectability_observer"]["comparison_observer_only_not_alternate_estimator"]
        ),
        "shipping_process_UCC_consumed": True,
        "H18_one_step_process_Q_lambda_min_lower_diagnostic_only": h_q,
        "A21_one_step_process_Q_lambda_min_lower_diagnostic_only": a_q,
        "one_step_process_Q_used_as_contraction_strictness": False,
        "exact_complete_word_event_algebra_consumed": True,
        "event_algebra_covers_immediate_left_error_reset": bool(
            reset["same_full_matrix_margin_preserved_by_congruence"]
            and reset["determinant_lower"] == 1.0
            and reset["small_angle_needed_for_nonsingularity"] is False
        ),
        "event_algebra_preserves_first_established_full_matrix_margin": event_preserved,
        "reset_complete_literal_word_consumed": True,
        "reset_complete_literal_execution_closed": reset_complete,
        "reset_injection_supplied_by_same_source_word": bool(
            resetword["reset_injection_supplied_by_same_source_word"]
        ),
        "reset_small_angle_bound_required": bool(resetword["reset_small_angle_bound_required"]),
        "useful_gate": USEFUL_GATE,
        "P3_math_subobligations_H18_and_A21_detectability_closed": (
            h18_closed and a21_detectability_closed
        ),
        "P3_full_matrix_H18_A21_and_reset_chain_closed": conditional_pass,
        "P3_remaining_canonical_obligation": "NONE" if conditional_pass else "CERTIFICATE_CHAIN_NOT_CLOSED",
        "P3_CANONICAL_PASS": conditional_pass,
        "P4_MAY_CONSUME_P3": conditional_pass,
        "next_obligation": (
            "conditional complete-SEA3 P3 closed; physical SEA0->SEA3 left inclusion remains a separate deployment obligation"
            if conditional_pass
            else "close H18, hybrid A21, and reset-complete full-matrix chain without weakening delta=1e-18"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "theorem_conditional_on_admitted_complete_SEA3_word",
        "global_physical_deployment_left_inclusion_is_separate_obligation",
        "global_physical_left_inclusion_required_before_full_deployment_theorem",
        "four_S_information_component_consumed",
        "actual_applied_SpectralMSE_R_S_consumed",
        "all_due_S_updates_remain_in_complete_word",
        "all_valid_accelerometer_updates_required",
        "H18_full_information_matrix_lower_closed",
        "H18_prior_free_completion_closed",
        "H18_full_18x18_interval_LDLT_closed",
        "A21_finite_bias_correlation_route_consumed",
        "A21_detectability_completion_closed",
        "A21_paper_UES_hypotheses_closed",
        "A21_hybrid_release_full_matrix_consumed",
        "A21_H18_word_finishes_before_release",
        "A21_exact_direct_sum_full_matrix_hybrid_proof_used",
        "A21_full_21x21_Riccati_bridge_closed",
        "A21_comparison_observer_is_proof_only_not_alternate_estimator",
        "shipping_process_UCC_consumed",
        "exact_complete_word_event_algebra_consumed",
        "event_algebra_covers_immediate_left_error_reset",
        "event_algebra_preserves_first_established_full_matrix_margin",
        "reset_complete_literal_word_consumed",
        "reset_complete_literal_execution_closed",
        "reset_injection_supplied_by_same_source_word",
        "P3_math_subobligations_H18_and_A21_detectability_closed",
        "P3_full_matrix_H18_A21_and_reset_chain_closed",
        "P3_CANONICAL_PASS",
        "P4_MAY_CONSUME_P3",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "finite_source_family_materialization_required_for_conditional_theorem",
        "global_physical_deployment_left_inclusion_closed_here",
        "global_physical_left_inclusion_required_before_conditional_P3_math",
        "hard_shaping_physical_left_inclusion_currently_closed",
        "source_family_replaced",
        "gaussian_good_event_used_as_source",
        "spectral_moments_used_as_source",
        "trajectory_replay_used_as_source",
        "independent_bounded_input_box_used_as_source",
        "four_S_selected_events_replace_complete_word",
        "accelerometer_rejection_after_certified_Normal_Live_allowed",
        "A21_uses_eta9_packet_shortcut",
        "one_step_process_Q_used_as_contraction_strictness",
        "reset_small_angle_bound_required",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    for key in (
        "H18_information_lambda_min_lower",
        "H18_worst_LDLT_pivot_lower",
        "H18_delta_squared_completion_penalty",
        "eta6_information_lambda_min_lower",
        "A21_bias_homogeneous_contraction_gap_lower",
        "A21_detectability_asymptotic_word_energy_gap_lower",
        "A21_first_active_ba_M_delta_margin_lower",
        "H18_one_step_process_Q_lambda_min_lower_diagnostic_only",
        "A21_one_step_process_Q_lambda_min_lower_diagnostic_only",
    ):
        x = d.get(key)
        if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
            f.append(f"invalid positive quantitative field {key}")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("useful gate changed")
    if float(d.get("A21_detectability_asymptotic_word_energy_gap_lower", 0.0)) < USEFUL_GATE:
        f.append("A21 detectability energy gap below useful gate")
    if d.get("P3_remaining_canonical_obligation") != "NONE":
        f.append("conditional P3 still reports a canonical math obligation")
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
        "H18_prior_free_closed": d["H18_prior_free_completion_closed"],
        "H18_worst_LDLT_pivot": d["H18_worst_LDLT_pivot_lower"],
        "A21_detectability_closed": d["A21_detectability_completion_closed"],
        "A21_hybrid_full_matrix_closed": d["A21_full_21x21_Riccati_bridge_closed"],
        "A21_ba_margin": d["A21_first_active_ba_M_delta_margin_lower"],
        "reset_literal_closed": d["reset_complete_literal_execution_closed"],
        "conditional_P3_pass": d["P3_CANONICAL_PASS"],
        "physical_left_inclusion_closed": d["hard_shaping_physical_left_inclusion_currently_closed"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
