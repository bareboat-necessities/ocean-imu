#!/usr/bin/env python3
"""Bind the OU-III paper's normal-Live theorem hypotheses to a source-word language.

The source-domain contract deliberately contains accepted and rejected measurement
branches.  Full-heading stability, however, is conditional on recurring finite-window
vector persistent excitation (PE).  A single accepted accel/mag pair is not enough to
tile an infinite execution.  This producer makes that missing quantifier explicit.

It does not enclose a Riccati or nonlinear word.  Instead it defines the exact
conditional language that a validated backend must cover before it may assert a
source-complete H/A word family:

* configured-runtime timing only;
* fixed-dimensional normal-Live mode H or A;
* source-reachable accepted/rejected branches between required PE packets;
* every PE recurrence window contains two consecutive accepted configured magnetic
  packets, with accepted accelerometer vectors at both packet times and the vector
  geometry/rate bounds of ``ou3_vector_uco_certificate``;
* the source pseudo-measurement scheduler retains its rigorously bounded firing gap;
* hard dimension-changing/reference/reset events are excluded from the same-mode word
  and remain separate hybrid obligations.

The PE recurrence window is a deployment/theorem hypothesis, not an estimator source
constant.  Therefore the default build reports the word language as blocked until a
finite positive recurrence bound is supplied explicitly.  This prevents a later
validated enclosure from silently proving only a favorable subset of source words.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_source_domain_contract as SOURCE
import ou3_translational_uco_ucc as TRANS
import ou3_vector_uco_certificate as VECTOR

SCHEMA = 1


def _finite_positive(x) -> bool:
    try:
        y = float(x)
    except (TypeError, ValueError):
        return False
    return math.isfinite(y) and y > 0.0


def build(
    pe_recurrence_window_s: float | None = None,
    header: Path = SOURCE.DEFAULT_HEADER.resolve(),
) -> dict:
    header = header.resolve()
    source = SOURCE.build(header)
    trans = TRANS.build(header)
    vector = VECTOR.build()

    source_failures = []
    if source.get("source_complete_parameter_domain") is not True:
        source_failures.append("source parameter domain is not complete")
    if TRANS.validate(trans):
        source_failures.append("translational UCO/UCC contract did not validate")
    if VECTOR.validate(vector):
        source_failures.append("vector packet UCO contract did not validate")

    runtime = source["configured_runtime_assumption"]
    dt = float(runtime["imu_dt_s"])
    packet_gap = list(vector["operating_envelope"]["packet_gap_s"])
    packet_span_upper = float(packet_gap[1])
    pseudo = trans["S_observation_uco"]
    detect = trans["integrator_detectability"]

    recurrence_supplied = _finite_positive(pe_recurrence_window_s)
    recurrence = float(pe_recurrence_window_s) if recurrence_supplied else None
    recurrence_failures: list[str] = []
    if recurrence_supplied and recurrence < packet_span_upper:
        recurrence_failures.append(
            "PE recurrence window is shorter than one certified consecutive magnetic-packet span"
        )

    language_ready = bool(
        not source_failures
        and recurrence_supplied
        and not recurrence_failures
    )

    word_horizon = None
    word_samples_upper = None
    if language_ready:
        # A tile must be long enough both to contain the declared vector-PE recurrence
        # window and to realize the source-uniform three-firing translational
        # detectability window.  The backend may choose a longer word.
        word_horizon = max(recurrence, float(detect["aligned_window_s"]))
        word_samples_upper = int(math.ceil(word_horizon / dt)) + 1

    pe = dict(vector["operating_envelope"])
    pe.update({
        "recurrence_window_s": recurrence,
        "recurrence_quantifier": (
            "every normal-Live interval of this duration contains at least one certified "
            "two-packet vector-PE event"
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
        },
        "translation_recurrence": {
            "pseudo_gap_min_s": pseudo["pseudo_gap_min_s"],
            "pseudo_gap_max_s": pseudo["pseudo_gap_max_s"],
            "three_firing_detectability_window_s": detect["aligned_window_s"],
            "stable_aw_alpha_upper": detect["stable_aw_alpha_upper"],
            "source_complete": trans["translation_source_complete"],
        },
        "vector_persistent_excitation": pe,
        "conditional_word_language": {
            "ready": language_ready,
            "word_horizon_lower_s": word_horizon,
            "word_samples_upper_at_configured_dt": word_samples_upper,
            "tiling_rule": (
                "tile every fixed-mode normal-Live execution by bounded words that each cover "
                "the translation recurrence window and one declared vector-PE recurrence window"
                if language_ready else None
            ),
            "coverage_rule": (
                "validated backend must enclose every source-reachable branch/continuous-parameter "
                "realization satisfying the declared theorem hypotheses; selected accepted packets "
                "may not be treated as the only admissible branches"
            ),
        },
        "source_complete_relative_to_theorem_hypotheses": language_ready,
        "continuous_word_enclosed": False,
        "nonlinear_word_enclosed": False,
        "theorem_promotion": "NOT_ESTABLISHED",
        "failures": source_failures + recurrence_failures + (
            [] if recurrence_supplied else [
                "finite vector-PE recurrence window is an explicit deployment theorem hypothesis and was not supplied"
            ]
        ),
        "next_obligation": (
            "provide a declared PE recurrence bound for the deployment envelope, then outward-enclose "
            "all H/A information-metric words in this conditional source-complete language, including "
            "finite prefix gain, exact nonlinear SO(3) endpoint decrease, and prefix safety"
        ),
    }


def validate(payload: dict) -> list[str]:
    failures: list[str] = []
    if payload.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if payload.get("source_generated_not_trajectory_fit") is not True:
        failures.append("word-language contract must not use trajectory fitting")
    if payload.get("normal_live_scope", {}).get("same_mode_only") is not True:
        failures.append("normal-Live word language must be fixed-dimensional")
    if payload.get("normal_live_scope", {}).get("dimension_change_multiplied_as_square_word") is not False:
        failures.append("dimension-changing transitions must remain hybrid")
    tr = payload.get("translation_recurrence", {})
    if tr.get("source_complete") is not True:
        failures.append("translation recurrence is not source-complete")
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

    ready = payload.get("conditional_word_language", {}).get("ready") is True
    recurrence = pe.get("recurrence_window_s")
    if ready:
        if not _finite_positive(recurrence):
            failures.append("ready word language lacks a finite PE recurrence window")
        elif isinstance(gap, list) and len(gap) == 2 and float(recurrence) < float(gap[1]):
            failures.append("PE recurrence window is shorter than the packet span")
        if payload.get("source_complete_relative_to_theorem_hypotheses") is not True:
            failures.append("ready word language must be conditionally source-complete")
    else:
        if payload.get("source_complete_relative_to_theorem_hypotheses") is not False:
            failures.append("blocked word language must not claim source completeness")

    if payload.get("continuous_word_enclosed") is not False:
        failures.append("word-language contract must not claim continuous enclosure")
    if payload.get("nonlinear_word_enclosed") is not False:
        failures.append("word-language contract must not claim nonlinear enclosure")
    if payload.get("theorem_promotion") != "NOT_ESTABLISHED":
        failures.append("word-language contract must not promote the theorem")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--header", type=Path, default=SOURCE.DEFAULT_HEADER)
    ap.add_argument(
        "--pe-recurrence-window-s",
        type=float,
        default=None,
        help="deployment theorem hypothesis: every such Live window contains a certified vector-PE packet pair",
    )
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.pe_recurrence_window_s, args.header.resolve())
    structural_failures = validate(out)
    out["contract_validation_pass"] = not structural_failures
    out["contract_validation_failures"] = structural_failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "qualification": out["qualification"],
        "conditional_word_language_ready": out["conditional_word_language"]["ready"],
        "source_complete_relative_to_theorem_hypotheses": out[
            "source_complete_relative_to_theorem_hypotheses"
        ],
        "word_horizon_lower_s": out["conditional_word_language"]["word_horizon_lower_s"],
        "failures": out["failures"],
        "contract_validation_failures": structural_failures,
    }, indent=2, sort_keys=True))
    return 0 if not structural_failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
