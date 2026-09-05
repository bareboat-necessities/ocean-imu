#!/usr/bin/env python3
"""Conditional complete-SEA3 P3 composition contract.

The OU-III theorem has two distinct quantified statements that must not be
collapsed into one gate:

1. Conditional P3: for every *admitted* complete phase-continuous SEA3
   Normal-Live word zeta=(x^s,lambda,z^t,q) satisfying the declared hard
   SEA3/Normal-Live conditions, prove the H18/A21 Riccati stability statement.
2. Physical SEA0 left inclusion: prove that the intended physical/random-sea
   realization family is contained in that admitted deterministic SEA3 set.

The second statement is required before promoting a deployment theorem about
all physical random seas, but it is not a logical prerequisite for proving the
first conditional theorem.  No Gaussian confidence set, spectral moment,
seeded replay, independent bounded-input box or alternate source language may
be invented to make either statement pass.

The complete-source quantitative chain now has two closed mathematical
sub-obligations at the canonical 1e-18 usefulness gate:

* H18: actual-applied-R_S four-S translation information + asynchronous eta6
  PE + shipping process UCC + the exact prior-free full 18x18 matrix completion;
* A21 theorem route: the same H18 complete-word stability plus the implemented
  finite residual accelerometer-bias correlation, bounded one-time H->A mode
  jump, and active A21 process UCC.  This is the paper's eta6+finite-tau_b
  detectability route, not an eta9 point-packet shortcut.

Canonical P3 nevertheless remains fail-closed here because the repository's
implementation-word promotion contract still asks for the stronger explicit
A21 Riccati-word bridge (or a formally accepted theorem-equivalent quantitative
replacement).  Closing H18 and A21 detectability is therefore recorded without
mislabeling the remaining implementation bridge as complete.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_sea3_a21_detectability_completion as A21
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_four_s_translation_information as FOUR_S
import ou3_sea3_full_word_event_algebra as EVENT
import ou3_sea3_h18_information_composition as H18INFO
import ou3_sea3_h18_prior_free_completion as H18PF
import ou3_sea3_hard_shaping_state as SHAPING
import ou3_sea3_windowed_vector_pe as PE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_COMPLETE_SEA3_CONDITIONAL_P3_COMPOSITION"
USEFUL_GATE = 1.0e-18


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    four = FOUR_S.build(path)
    event = EVENT.build()
    h18 = H18INFO.build(path)
    h18pf = H18PF.build(path)
    a21 = A21.build(path)
    shaping = SHAPING.build()
    pe = PE.build(path)
    process = PROCESS.build()
    bad = {
        "complete": COMPLETE.validate(complete),
        "four_S": FOUR_S.validate(four),
        "event_algebra": EVENT.validate(event),
        "H18_information": H18INFO.validate(h18),
        "H18_prior_free": H18PF.validate(h18pf),
        "A21_detectability": A21.validate(a21),
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
        ("A21 detectability", a21["canonical_source"]),
    ):
        if value != source:
            raise RuntimeError(f"{name} detached from the canonical complete SEA3 source")

    h18_info = float(h18["triangular_information_composition"]["D_H18_lambda_min_lower"])
    h_q = float(process["modes"]["H"]["prediction_Q_lambda_min_lower"])
    a_q = float(process["modes"]["A"]["prediction_Q_lambda_min_lower"])
    ba_gap = float(pe["A_mode_bias_route"]["homogeneous_bias_contraction_gap_lower"])
    a21_gap = float(a21["A21_detectability_asymptotic_word_energy_gap_lower"])
    hard_left_inclusion = bool(
        shaping["executable_ingredients"]["complete_SEA3_left_inclusion_closed"]
    )
    reset = event["left_error_reset"]
    preservation = event["full_matrix_margin_preservation"]
    h18_closed = bool(h18pf["H18_prior_free_completion_closed"])
    a21_detectability_closed = bool(a21["A21_finite_bias_detectability_closed"])
    a21_full_riccati_bridge = bool(a21["full_21x21_Omega_minus_delta_P_LDLT_closed_here"])

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": source,
        "theorem_conditional_on_admitted_complete_SEA3_word": bool(
            complete["theorem_conditional_on_admitted_complete_SEA3_word"]
        ),
        "conditional_P3_quantifier": "FOR_EVERY_ADMITTED_COMPLETE_SEA3_NORMAL_LIVE_WORD",
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
        "H18_full_information_matrix_lower_closed": True,
        "H18_information_lambda_min_lower": h18_info,
        "H18_prior_free_completion_closed": h18_closed,
        "H18_full_18x18_interval_LDLT_closed": bool(h18pf["full_H18_prior_free_matrix_condition_closed"]),
        "H18_worst_LDLT_pivot_lower": h18pf["worst_full_H18_LDLT_pivot_lower"],
        "H18_delta_squared_completion_penalty": h18pf["delta_squared_over_four_penalty_physical"],
        "eta6_information_lambda_min_lower": float(
            pe["eta6_information"]["alpha_6_information_lower"]
        ),
        "A21_uses_eta9_packet_shortcut": False,
        "A21_finite_bias_correlation_route_consumed": True,
        "A21_bias_homogeneous_contraction_gap_lower": ba_gap,
        "A21_detectability_completion_closed": a21_detectability_closed,
        "A21_paper_UES_hypotheses_closed": bool(a21["A21_paper_UES_hypotheses_closed"]),
        "A21_detectability_asymptotic_word_energy_gap_lower": a21_gap,
        "A21_full_21x21_Riccati_bridge_closed": a21_full_riccati_bridge,
        "A21_comparison_observer_is_proof_only_not_alternate_estimator": bool(
            a21["triangular_detectability_observer"]["comparison_observer_only_not_alternate_estimator"]
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
        "event_algebra_preserves_first_established_full_matrix_margin": bool(
            preservation["covers_prediction"]
            and preservation["covers_every_due_S_update"]
            and preservation["covers_every_Normal_Live_accelerometer_update"]
            and preservation["covers_asynchronous_magnetometer_update"]
            and preservation["covers_immediate_left_error_reset"]
            and preservation["covers_aw_covariance_floor"]
            and preservation["covers_not_due_or_rejected_identity_branches"]
        ),
        "useful_gate": USEFUL_GATE,
        "P3_math_subobligations_H18_and_A21_detectability_closed": (
            h18_closed and a21_detectability_closed
        ),
        "P3_remaining_canonical_obligation": (
            "FORMAL_A21_RICCATI_WORD_EQUIVALENCE_OR_FULL_21X21_BRIDGE"
            if not a21_full_riccati_bridge else "NONE"
        ),
        "P3_CANONICAL_PASS": False,
        "P4_MAY_CONSUME_P3": False,
        "next_obligation": (
            "prove that the complete-SEA3 eta6+finite-tau_b detectability/UCC certificate is quantitatively "
            "equivalent to the canonical active-mode Riccati-word contraction requirement, or close the full "
            "A21 Omega-delta*P matrix directly; do not weaken delta=1e-18 and do not introduce eta9 packets, "
            "source replay, or an alternate estimator"
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
        "H18_full_information_matrix_lower_closed",
        "H18_prior_free_completion_closed",
        "H18_full_18x18_interval_LDLT_closed",
        "A21_finite_bias_correlation_route_consumed",
        "A21_detectability_completion_closed",
        "A21_paper_UES_hypotheses_closed",
        "A21_comparison_observer_is_proof_only_not_alternate_estimator",
        "shipping_process_UCC_consumed",
        "exact_complete_word_event_algebra_consumed",
        "event_algebra_covers_immediate_left_error_reset",
        "event_algebra_preserves_first_established_full_matrix_margin",
        "P3_math_subobligations_H18_and_A21_detectability_closed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "global_physical_deployment_left_inclusion_closed_here",
        "global_physical_left_inclusion_required_before_conditional_P3_math",
        "hard_shaping_physical_left_inclusion_currently_closed",
        "source_family_replaced",
        "gaussian_good_event_used_as_source",
        "spectral_moments_used_as_source",
        "trajectory_replay_used_as_source",
        "independent_bounded_input_box_used_as_source",
        "four_S_selected_events_replace_complete_word",
        "A21_uses_eta9_packet_shortcut",
        "A21_full_21x21_Riccati_bridge_closed",
        "one_step_process_Q_used_as_contraction_strictness",
        "P3_CANONICAL_PASS",
        "P4_MAY_CONSUME_P3",
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
    if d.get("P3_remaining_canonical_obligation") != "FORMAL_A21_RICCATI_WORD_EQUIVALENCE_OR_FULL_21X21_BRIDGE":
        f.append("canonical P3 remaining obligation is not the A21 bridge")
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
        "H18_information_lower": d["H18_information_lambda_min_lower"],
        "H18_prior_free_closed": d["H18_prior_free_completion_closed"],
        "H18_worst_LDLT_pivot": d["H18_worst_LDLT_pivot_lower"],
        "A21_detectability_closed": d["A21_detectability_completion_closed"],
        "A21_detectability_gap": d["A21_detectability_asymptotic_word_energy_gap_lower"],
        "A21_full_Riccati_bridge_closed": d["A21_full_21x21_Riccati_bridge_closed"],
        "remaining": d["P3_remaining_canonical_obligation"],
        "conditional_P3_pass": d["P3_CANONICAL_PASS"],
        "physical_left_inclusion_closed": d["hard_shaping_physical_left_inclusion_currently_closed"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
