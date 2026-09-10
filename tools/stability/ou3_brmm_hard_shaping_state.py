#!/usr/bin/env python3
"""Fail-closed contract for the deterministic complete-BRMM shaping state.

BRMM is already the compact theorem-domain sea family. The continuum phase
coordinate and its no-reseed propagation are closed. The finite sampled target
B^601_BRMM is compact, and a new validated correlated deterministic outer set
O^601_BRMM now closes the alternative hard finite-window route

    B^601_BRMM subset O^601_BRMM.

That outer set preserves cross-sample primitive recurrence, the all-axis
acceleration-moment IQC, SO(3)/rate constraints and the coupled lambda relation;
it is not an independent bounded-input box. The remaining open source-side
obligation is to materialize this relation through the joint translational /
rotational 601-sample output map consumed by the typed estimator/Riccati
executor. Until that is closed this module remains fail-closed for P3/P4.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_brmm_continuum_phase_state as PHASE
import ou3_brmm_hard_window_behavior as BEHAVIOR

REPO = Path(__file__).resolve().parents[2]
THEOREM = REPO / "doc" / "kalman_ou_iii" / "w3d-marine-reference-models.tex-part"
COMPLETE_SOURCE = REPO / "tools" / "stability" / "ou3_brmm_complete_source.py"
SCHEMA = 5
QUALIFICATION = "OU3_BRMM_HARD_SHAPING_STATE_CONTRACT_V5"

CONTINUUM_PHASE_COORDINATE_SET_CLOSED = PHASE.CONTINUUM_PHASE_COORDINATE_SET_CLOSED
PHASE_CONTINUOUS_PROPAGATION_CLOSED = PHASE.PHASE_CONTINUOUS_PROPAGATION_CLOSED
# Exact continuum spectral-driver membership is not required once a validated
# deterministic correlated outer enclosure is available.
HARD_SPECTRAL_DRIVER_SET_CLOSED = False
CORRELATED_OUTER_ENCLOSURE_CLOSED = True
COMPLETE_BRMM_LEFT_INCLUSION_CLOSED = True
JOINT_SOURCE_OUTPUT_MAP_CLOSED = False

HARD_SHAPING_STATE_OR_EXCITATION_BOUND_CLOSED = all((
    CONTINUUM_PHASE_COORDINATE_SET_CLOSED,
    PHASE_CONTINUOUS_PROPAGATION_CLOSED,
    CORRELATED_OUTER_ENCLOSURE_CLOSED,
    COMPLETE_BRMM_LEFT_INCLUSION_CLOSED,
    JOINT_SOURCE_OUTPUT_MAP_CLOSED,
))


def _normalized_text(text: str) -> str:
    return " ".join(text.split())


def build() -> dict:
    theorem = THEOREM.read_text(encoding="utf-8")
    theorem_flat = _normalized_text(theorem)
    complete = COMPLETE_SOURCE.read_text(encoding="utf-8")
    phase = PHASE.build()
    behavior = BEHAVIOR.build()
    phase_failures = PHASE.validate(phase)
    behavior_failures = BEHAVIOR.validate(behavior)
    if phase_failures or behavior_failures:
        raise RuntimeError(
            f"BRMM shaping prerequisites failed: phase={phase_failures}, behavior={behavior_failures}"
        )

    theorem_has_shaping_system = (
        "x^s_{k+1}&=A_s" in theorem
        and "u^s_k&=C_s" in theorem
        and "oscillator/shaping state or an equivalent hard finite-window" in theorem_flat
    )
    theorem_has_explicit_hard_realization_set = (
        "\\mathcal X^s_{\\rm ref}(\\lambda_{0:N_W})" in theorem
        and "eq:marine-reference-hard-realization-set" in theorem
        and "BRMM itself does not require a finite spectral state" in theorem_flat
        and "machine-readable outward representation" in theorem_flat
    )
    theorem_separates_probabilistic_corollary = (
        "A probabilistic statement for random sea realizations is a later corollary" in theorem_flat
        and "no infinite-time pointwise bound is inferred merely from a Gaussian spectrum" in theorem_flat
    )
    theorem_rejects_finite_or_statistical_surrogates = (
        "neither a Gaussian confidence event" in theorem_flat
        and "spectral moments alone" in theorem_flat
        and "finite seeded harmonic" in theorem_flat
    )
    source_rejects_gaussian_generator = (
        '"used_to_generate_P3_source_words": False' in complete
        and '"used_to_prune_homogeneous_P3_family": False' in complete
    )
    correlated_outer = bool(
        behavior["validated_correlated_outer_enclosure_closed"]
        and behavior["correlated_outer_left_inclusion_closed"]
        and behavior["correlated_outer_retains_cross_sample_and_axis_dependence"]
    )
    if not correlated_outer:
        raise RuntimeError("validated correlated BRMM outer enclosure disappeared")

    executable = {
        "continuum_phase_coordinate_set_closed": CONTINUUM_PHASE_COORDINATE_SET_CLOSED,
        "phase_continuous_propagation_closed": PHASE_CONTINUOUS_PROPAGATION_CLOSED,
        "hard_spectral_driver_set_closed": HARD_SPECTRAL_DRIVER_SET_CLOSED,
        "correlated_outer_enclosure_closed": correlated_outer,
        "complete_BRMM_left_inclusion_closed": COMPLETE_BRMM_LEFT_INCLUSION_CLOSED,
        "joint_source_output_map_closed": JOINT_SOURCE_OUTPUT_MAP_CLOSED,
    }
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "reference_parameter_domain_compact": True,
        "compactness_is_not_an_open_obligation": True,
        "theorem_has_deterministic_shaping_contract": theorem_has_shaping_system,
        "theorem_has_explicit_hard_realization_set": theorem_has_explicit_hard_realization_set,
        "theorem_rejects_statistical_or_seeded_surrogates": theorem_rejects_finite_or_statistical_surrogates,
        "theorem_separates_probabilistic_random_sea_corollary": theorem_separates_probabilistic_corollary,
        "complete_source_rejects_gaussian_word_generator": source_rejects_gaussian_generator,
        "hard_realization_set_symbol": "X^s_ref(lambda_{0:N_W})",
        "continuum_phase_certificate": {
            "qualification": phase["qualification"],
            "phase_state_set": phase["phase_state_set"],
            "continuum_index_set_retained": phase["continuum_index_set_retained"],
            "finite_frequency_grid_used": phase["finite_frequency_grid_used"],
            "finite_direction_grid_used": phase["finite_direction_grid_used"],
            "phase_reset_on_lambda_transition_allowed": phase["phase_reset_on_lambda_transition_allowed"],
            "continuum_phase_coordinate_set_closed": phase["continuum_phase_coordinate_set_closed"],
            "phase_continuous_propagation_closed": phase["phase_continuous_propagation_closed"],
        },
        "sampled_behavior_target": {
            "qualification": behavior["qualification"],
            "symbol": behavior["behavior_set_symbol"],
            "sample_count": behavior["sample_count"],
            "sampled_projection_dimension": behavior["sampled_projection_dimension"],
            "compact": behavior["sampled_behavior_set_compact"],
            "membership_requires_common_BRMM_witness": behavior["membership_requires_common_BRMM_witness"],
            "normal_live_caps_are_membership_sufficient": behavior["normal_live_caps_are_membership_sufficient"],
            "independent_sample_boxes_define_behavior_set": behavior["independent_sample_boxes_define_behavior_set"],
            "validated_membership_or_separation_oracle_closed": behavior["validated_membership_or_separation_oracle_closed"],
            "validated_correlated_outer_enclosure_closed": behavior["validated_correlated_outer_enclosure_closed"],
            "correlated_outer_set_symbol": behavior["correlated_outer_set_symbol"],
            "correlated_outer_left_inclusion_closed": behavior["correlated_outer_left_inclusion_closed"],
        },
        "exact_spectral_membership_oracle_required_for_P4": False,
        "correlated_outer_enclosure_route_used": True,
        "power_spectrum_alone_is_hard_pathwise_bound": False,
        "spectral_moments_alone_may_close_xs": False,
        "gaussian_good_event_may_close_xs": False,
        "replay_may_close_xs": False,
        "seeded_128_frequency_generator_may_close_xs": False,
        "finite_RAO_grid_may_close_xs": False,
        "arbitrary_bounded_input_box_may_close_xs": False,
        "allowed_closure_forms": [
            "validated_compact_oscillator_or_shaping_state_with_hard_driver_set",
            "validated_equivalent_hard_finite_window_dynamic_constraint",
        ],
        "executable_ingredients": executable,
        "hard_shaping_state_or_excitation_bound_closed": HARD_SHAPING_STATE_OR_EXCITATION_BOUND_CLOSED,
        "complete_BRMM_family_materialized_here": False,
        "P3_promoted": False,
        "next_obligation": (
            "left inclusion into the correlated hard finite-window outer relation is closed; materialize that same relation through the joint source-output/typed 601-sample executor without Cartesianizing samples or axes"
        ),
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        failures.append("schema/qualification mismatch")
    for key in (
        "reference_parameter_domain_compact", "compactness_is_not_an_open_obligation",
        "theorem_has_deterministic_shaping_contract", "theorem_has_explicit_hard_realization_set",
        "theorem_rejects_statistical_or_seeded_surrogates",
        "theorem_separates_probabilistic_random_sea_corollary",
        "complete_source_rejects_gaussian_word_generator",
        "correlated_outer_enclosure_route_used",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    if d.get("exact_spectral_membership_oracle_required_for_P4") is not False:
        failures.append("exact spectral membership was incorrectly made mandatory")
    if d.get("hard_realization_set_symbol") != "X^s_ref(lambda_{0:N_W})":
        failures.append("hard realization set symbol drifted")
    phase = d.get("continuum_phase_certificate", {})
    for key in ("continuum_index_set_retained", "continuum_phase_coordinate_set_closed", "phase_continuous_propagation_closed"):
        if phase.get(key) is not True:
            failures.append(f"continuum phase certificate lost {key}")
    for key in ("finite_frequency_grid_used", "finite_direction_grid_used", "phase_reset_on_lambda_transition_allowed"):
        if phase.get(key) is not False:
            failures.append(f"continuum phase certificate reintroduced {key}")
    behavior = d.get("sampled_behavior_target", {})
    for key in ("compact", "membership_requires_common_BRMM_witness", "validated_correlated_outer_enclosure_closed", "correlated_outer_left_inclusion_closed"):
        if behavior.get(key) is not True:
            failures.append(f"sampled behavior target lost {key}")
    for key in ("normal_live_caps_are_membership_sufficient", "independent_sample_boxes_define_behavior_set", "validated_membership_or_separation_oracle_closed"):
        if behavior.get(key) is not False:
            failures.append(f"sampled behavior target falsely closes/reintroduces {key}")
    for key in (
        "power_spectrum_alone_is_hard_pathwise_bound", "spectral_moments_alone_may_close_xs",
        "gaussian_good_event_may_close_xs", "replay_may_close_xs",
        "seeded_128_frequency_generator_may_close_xs", "finite_RAO_grid_may_close_xs",
        "arbitrary_bounded_input_box_may_close_xs", "hard_shaping_state_or_excitation_bound_closed",
        "complete_BRMM_family_materialized_here", "P3_promoted",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    expected = {
        "continuum_phase_coordinate_set_closed": True,
        "phase_continuous_propagation_closed": True,
        "hard_spectral_driver_set_closed": False,
        "correlated_outer_enclosure_closed": True,
        "complete_BRMM_left_inclusion_closed": True,
        "joint_source_output_map_closed": False,
    }
    if d.get("executable_ingredients") != expected:
        failures.append("hard shaping executable ingredient gates drifted")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "BRMM_compact": d["reference_parameter_domain_compact"],
        "correlated_outer": d["sampled_behavior_target"]["validated_correlated_outer_enclosure_closed"],
        "left_inclusion": d["sampled_behavior_target"]["correlated_outer_left_inclusion_closed"],
        "executable_ingredients": d["executable_ingredients"],
        "hard_shaping_closed": d["hard_shaping_state_or_excitation_bound_closed"],
        "next_obligation": d["next_obligation"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
