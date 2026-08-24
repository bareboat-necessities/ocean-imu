#!/usr/bin/env python3
"""Bind the OU-III paper's normal-Live theorem hypotheses to a source-word language.

A source-complete fixed-mode tile must cover the declared recurring vector-PE
event and the source-cadence-only four-S observability construction for the
complete [v,p,S,a_w] chain.  The source-valid three-S [v,p,S] detectability plus
stable a_w route may sharpen a Riccati/covariance upper bound; it does not replace
the four-S full-observability qualification and is not a promotion fallback.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import ou3_source_domain_contract as SOURCE
import ou3_translational_uco_ucc as TRANS
import ou3_vector_uco_certificate as VECTOR

SCHEMA = 2


def _finite_positive(x) -> bool:
    try:
        y = float(x)
    except (TypeError, ValueError):
        return False
    return math.isfinite(y) and y > 0.0


def build(pe_recurrence_window_s: float | None = None,
          header: Path = SOURCE.DEFAULT_HEADER.resolve()) -> dict:
    header = header.resolve()
    source = SOURCE.build(header)
    trans = TRANS.build(header)
    vector = VECTOR.build()

    failures: list[str] = []
    if source.get("source_complete_parameter_domain") is not True:
        failures.append("source parameter domain is not complete")
    if TRANS.validate(trans):
        failures.append("translational UCO/UCC contract did not validate")
    if VECTOR.validate(vector):
        failures.append("vector packet UCO contract did not validate")

    runtime = source["configured_runtime_assumption"]
    dt = float(runtime["imu_dt_s"])
    packet_gap = list(vector["operating_envelope"]["packet_gap_s"])
    packet_span_upper = float(packet_gap[1])
    pseudo = trans["S_observation_uco"]
    detect = trans["integrator_detectability"]

    recurrence_supplied = _finite_positive(pe_recurrence_window_s)
    recurrence = float(pe_recurrence_window_s) if recurrence_supplied else None
    if recurrence_supplied and recurrence < packet_span_upper:
        failures.append("PE recurrence window is shorter than one certified consecutive magnetic-packet span")

    ready = bool(not failures and recurrence_supplied)
    word_horizon = word_samples_upper = q_W = selected_spacing_lower = spread_det_factor = None
    if ready:
        four_s_window = float(pseudo["aligned_window_s"])
        word_horizon = max(recurrence, four_s_window)
        word_samples_upper = int(math.ceil(word_horizon / dt)) + 1
        delta_min = float(pseudo["pseudo_gap_min_s"])
        delta_max = float(pseudo["pseudo_gap_max_s"])
        q_W = max(1, int(math.floor(word_horizon / (3.0 * delta_max))))
        selected_spacing_lower = q_W * delta_min
        spread_det_factor = q_W ** 6

    pe = dict(vector["operating_envelope"])
    pe.update({
        "recurrence_window_s": recurrence,
        "recurrence_quantifier": (
            "every normal-Live interval of this duration contains at least one certified two-packet vector-PE event"
        ),
        "accelerometer_required_at_both_vector_times": True,
        "two_consecutive_accepted_magnetic_packets_required": True,
        "arbitrary_rejections_between_required_pe_events_allowed": True,
        "hypothesis_origin": "DEPLOYMENT_THEOREM_ASSUMPTION_NOT_TRAJECTORY_FIT",
    })

    return {
        "schema": SCHEMA,
        "claim": "OU3_CONDITIONAL_SOURCE_COMPLETE_NORMAL_LIVE_WORD_LANGUAGE",
        "qualification": "THEOREM_HYPOTHESIS_CONTRACT_NOT_WORD_ENCLOSURE",
        "source_generated_not_trajectory_fit": True,
        "configured_runtime": runtime,
        "fixed_dimension_modes": {"H": 18, "A": 21},
        "normal_live_scope": {
            "same_mode_only": True,
            "hard_attitude_rewrite_inside_word": False,
            "hybrid_transitions_separate": list(source["hybrid_obligations"]),
            "dimension_change_multiplied_as_square_word": False,
        },
        "source_branch_language": {
            "accelerometer_gate": ["accepted", "rejected"],
            "magnetometer_gate": ["not_due", "accepted", "rejected"],
            "S_zero_pseudo": ["not_due", "due"],
            "aw_covariance_sync": ["not_due", "due_psd_increment"],
            "continuous_parameters": source["validated_parameter_box"]["continuous_parameters"],
            "continuous_parameters_outward_rounded": True,
            "joint_source_reachability_required": True,
            "cartesian_extrema_products_not_a_valid_word": True,
        },
        "translation_recurrence": {
            "full_observability_route": "FOUR_S_SPREAD_COMPLETE_V_P_S_AW_UCO",
            "primary_route": "FOUR_S_SPREAD_COMPLETE_V_P_S_AW_UCO",
            "primary_state_order": ["v", "p", "S", "a_w"],
            "aligned_firing_count": 4,
            "pseudo_gap_min_s": pseudo["pseudo_gap_min_s"],
            "pseudo_gap_max_s": pseudo["pseudo_gap_max_s"],
            "minimum_four_firing_window_s": pseudo["aligned_window_s"],
            "spread_index_q_W": q_W,
            "spread_selected_spacing_lower_s": selected_spacing_lower,
            "determinant_spacing_widening_factor_vs_adjacent": spread_det_factor,
            "three_firing_integrator_detectability_role": "Riccati_covariance_upper_sharpening_only",
            "three_firing_integrator_detectability_is_promotion_fallback": False,
            "three_firing_detectability_window_s": detect["aligned_window_s"],
            "stable_aw_alpha_upper": detect["stable_aw_alpha_upper"],
            "source_complete": trans["translation_source_complete"],
        },
        "vector_persistent_excitation": pe,
        "conditional_word_language": {
            "ready": ready,
            "word_horizon_lower_s": word_horizon,
            "word_samples_upper_at_configured_dt": word_samples_upper,
            "tiling_rule": (
                "tile every fixed-mode normal-Live execution by bounded source-correlated words that each cover one spread four-S complete translation observation and one declared vector-PE recurrence window"
                if ready else None
            ),
            "coverage_rule": (
                "validated backend must enclose every jointly source-reachable branch and continuous-parameter realization; selected accepted packets and independently chosen edge extrema may not replace the source word"
            ),
            "one_sample_decrease_required": False,
            "word_endpoint_decrease_required": True,
        },
        "source_complete_relative_to_theorem_hypotheses": ready,
        "continuous_word_enclosed": False,
        "nonlinear_word_enclosed": False,
        "theorem_promotion": "NOT_ESTABLISHED",
        "failures": failures + ([] if recurrence_supplied else [
            "finite vector-PE recurrence window is an explicit deployment theorem hypothesis and was not supplied"
        ]),
        "next_obligation": (
            "outward-enclose complete jointly source-reachable H/A endpoint maps in group-compatible node metrics; the three-S detectable block may sharpen covariance upper bounds but cannot replace the four-S full-observability word qualification"
        ),
    }


def validate(payload: dict) -> list[str]:
    failures: list[str] = []
    if payload.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if payload.get("source_generated_not_trajectory_fit") is not True:
        failures.append("word-language contract must not use trajectory fitting")
    scope = payload.get("normal_live_scope", {})
    if scope.get("same_mode_only") is not True:
        failures.append("normal-Live word language must be fixed-dimensional")
    if scope.get("dimension_change_multiplied_as_square_word") is not False:
        failures.append("dimension-changing transitions must remain hybrid")
    tr = payload.get("translation_recurrence", {})
    if tr.get("source_complete") is not True:
        failures.append("translation recurrence is not source-complete")
    if tr.get("full_observability_route") != "FOUR_S_SPREAD_COMPLETE_V_P_S_AW_UCO" or tr.get("aligned_firing_count") != 4:
        failures.append("translation full-observability route is not four-S complete-chain UCO")
    if tr.get("three_firing_integrator_detectability_role") != "Riccati_covariance_upper_sharpening_only":
        failures.append("three-S detectability role is not restricted to covariance-upper sharpening")
    if tr.get("three_firing_integrator_detectability_is_promotion_fallback") is not False:
        failures.append("three-S detectability is incorrectly available as promotion fallback")
    alpha = tr.get("stable_aw_alpha_upper")
    if not _finite_positive(alpha) or not float(alpha) < 1.0:
        failures.append("stable a_w tail is not strictly contractive")

    pe = payload.get("vector_persistent_excitation", {})
    gap = pe.get("packet_gap_s")
    if not isinstance(gap, list) or len(gap) != 2 or not all(_finite_positive(x) for x in gap):
        failures.append("vector packet gap is invalid")
    if pe.get("two_consecutive_accepted_magnetic_packets_required") is not True:
        failures.append("consecutive accepted magnetic packets must be explicit")
    if pe.get("accelerometer_required_at_both_vector_times") is not True:
        failures.append("accepted accelerometer packets must be explicit")

    word = payload.get("conditional_word_language", {})
    ready = word.get("ready") is True
    recurrence = pe.get("recurrence_window_s")
    if ready:
        if not _finite_positive(recurrence):
            failures.append("ready word language lacks a finite PE recurrence window")
        elif float(recurrence) < float(gap[1]):
            failures.append("PE recurrence window is shorter than the packet span")
        if payload.get("source_complete_relative_to_theorem_hypotheses") is not True:
            failures.append("ready word language must be conditionally source-complete")
        if not _finite_positive(tr.get("spread_selected_spacing_lower_s")):
            failures.append("spread-selected four-S spacing is not positive")
        if int(tr.get("spread_index_q_W", 0)) < 1:
            failures.append("spread-selected four-S index is invalid")
    elif payload.get("source_complete_relative_to_theorem_hypotheses") is not False:
        failures.append("blocked word language must not claim source completeness")

    if word.get("one_sample_decrease_required") is not False or word.get("word_endpoint_decrease_required") is not True:
        failures.append("word contract reintroduced a one-sample contraction requirement")
    branches = payload.get("source_branch_language", {})
    if branches.get("joint_source_reachability_required") is not True:
        failures.append("word language does not require joint source reachability")
    if payload.get("continuous_word_enclosed") is not False or payload.get("nonlinear_word_enclosed") is not False:
        failures.append("language stage must not masquerade as enclosure")
    if payload.get("theorem_promotion") != "NOT_ESTABLISHED":
        failures.append("language stage must not promote theorem")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--header", type=Path, default=SOURCE.DEFAULT_HEADER)
    ap.add_argument("--pe-recurrence-window-s", type=float, default=None)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.pe_recurrence_window_s, args.header.resolve())
    structural_failures = validate(out)
    out["contract_validation_pass"] = not structural_failures
    out["contract_validation_failures"] = structural_failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    import json
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "ready": out["conditional_word_language"]["ready"],
        "word_horizon_lower_s": out["conditional_word_language"]["word_horizon_lower_s"],
        "translation_recurrence": out["translation_recurrence"],
        "failures": out["failures"],
        "contract_validation_failures": structural_failures,
    }, indent=2, sort_keys=True))
    return 0 if not structural_failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
