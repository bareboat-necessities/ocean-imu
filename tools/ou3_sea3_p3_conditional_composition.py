#!/usr/bin/env python3
"""Conditional complete-SEA3 P3 composition contract.

The OU-III theorem has two distinct quantified statements that must not be
collapsed into one gate:

1. Conditional P3: for every *admitted* complete phase-continuous SEA3
   Normal-Live word zeta=(x^s,lambda,z^t,q) satisfying the declared hard
   SEA3/Normal-Live conditions, prove the H18/A21 Riccati word contraction.
2. Physical SEA0 left inclusion: prove that the intended physical/random-sea
   realization family is contained in that admitted deterministic SEA3 set.

The second statement is required before promoting a deployment theorem about
all physical random seas, but it is not a logical prerequisite for proving the
first conditional theorem.  In particular, no Gaussian confidence set,
spectral moment, seeded replay, independent bounded-input box or alternate
source language may be invented to make either statement pass.

This module also records the mandatory complete-word quantitative ingredients:
actual-applied-R_S four-S translation information, asynchronous eta6 PE,
shipping process UCC, and exact fixed-dimensional event algebra.  It does not
promote P3 until the remaining prior-free H18 completion and A21 detectability
composition are both closed at the canonical delta=1e-18 gate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_four_s_translation_information as FOUR_S
import ou3_sea3_full_word_event_algebra as EVENT
import ou3_sea3_h18_information_composition as H18INFO
import ou3_sea3_hard_shaping_state as SHAPING
import ou3_sea3_windowed_vector_pe as PE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_COMPLETE_SEA3_CONDITIONAL_P3_COMPOSITION"
USEFUL_GATE = 1.0e-18


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    four = FOUR_S.build(path)
    event = EVENT.build()
    h18 = H18INFO.build(path)
    shaping = SHAPING.build()
    pe = PE.build(path)
    process = PROCESS.build()
    bad = {
        "complete": COMPLETE.validate(complete),
        "four_S": FOUR_S.validate(four),
        "event_algebra": EVENT.validate(event),
        "H18_information": H18INFO.validate(h18),
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

    h18_info = float(h18["triangular_information_composition"]["D_H18_lambda_min_lower"])
    h_q = float(process["modes"]["H"]["prediction_Q_lambda_min_lower"])
    a_q = float(process["modes"]["A"]["prediction_Q_lambda_min_lower"])
    ba_gap = float(pe["A_mode_bias_route"]["homogeneous_bias_contraction_gap_lower"])

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
        "hard_shaping_physical_left_inclusion_currently_closed": bool(
            shaping["complete_SEA3_physical_left_inclusion_closed"]
        ),
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
        "eta6_information_lambda_min_lower": float(
            pe["eta6_information"]["alpha_6_information_lower"]
        ),
        "A21_uses_eta9_packet_shortcut": False,
        "A21_finite_bias_correlation_route_consumed": True,
        "A21_bias_homogeneous_contraction_gap_lower": ba_gap,
        "shipping_process_UCC_consumed": True,
        "H18_one_step_process_Q_lambda_min_lower_diagnostic_only": h_q,
        "A21_one_step_process_Q_lambda_min_lower_diagnostic_only": a_q,
        "one_step_process_Q_used_as_contraction_strictness": False,
        "exact_complete_word_event_algebra_consumed": True,
        "event_algebra_covers_immediate_left_error_reset": bool(
            event["reset"]["generalized_information_margin_congruence_invariant"]
        ),
        "event_algebra_preserves_first_established_full_matrix_margin": bool(
            event["strict_margin_preservation"]["covers_every_operation_class"]
        ),
        "useful_gate": USEFUL_GATE,
        "H18_prior_free_completion_closed": False,
        "A21_detectability_completion_closed": False,
        "P3_CANONICAL_PASS": False,
        "P4_MAY_CONSUME_P3": False,
        "next_obligation": (
            "close the same-word prior-free H18 process completion at delta=1e-18, then "
            "compose the A21 finite-b_a-correlation detectability block; once both full-matrix "
            "conditions pass, bind them into the canonical P3 gate without requiring the separate "
            "global physical SEA0 left-inclusion theorem"
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
        "A21_finite_bias_correlation_route_consumed",
        "shipping_process_UCC_consumed",
        "exact_complete_word_event_algebra_consumed",
        "event_algebra_covers_immediate_left_error_reset",
        "event_algebra_preserves_first_established_full_matrix_margin",
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
        "one_step_process_Q_used_as_contraction_strictness",
        "H18_prior_free_completion_closed",
        "A21_detectability_completion_closed",
        "P3_CANONICAL_PASS",
        "P4_MAY_CONSUME_P3",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    for key in (
        "H18_information_lambda_min_lower",
        "eta6_information_lambda_min_lower",
        "A21_bias_homogeneous_contraction_gap_lower",
        "H18_one_step_process_Q_lambda_min_lower_diagnostic_only",
        "A21_one_step_process_Q_lambda_min_lower_diagnostic_only",
    ):
        x = d.get(key)
        if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
            f.append(f"invalid positive quantitative field {key}")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("useful gate changed")
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
        "H_Q_one_step_diagnostic": d["H18_one_step_process_Q_lambda_min_lower_diagnostic_only"],
        "A_Q_one_step_diagnostic": d["A21_one_step_process_Q_lambda_min_lower_diagnostic_only"],
        "A_bias_contraction_gap": d["A21_bias_homogeneous_contraction_gap_lower"],
        "physical_left_inclusion_closed": d["hard_shaping_physical_left_inclusion_currently_closed"],
        "conditional_P3_pass": d["P3_CANONICAL_PASS"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
