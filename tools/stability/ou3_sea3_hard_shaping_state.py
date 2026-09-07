#!/usr/bin/env python3
"""Fail-closed contract for the deterministic complete-SEA3 shaping state.

SEA3 is already the compact theorem-domain sea family. This module does not
re-open that compactness question. It assembles the executable SEA0 pieces
needed to propagate the complete compact family over the canonical 3 s word.

The continuum phase coordinate and its no-reseed rotation are closed in
``ou3_sea3_continuum_phase_state``. The hard deterministic driver coordinate is
the single common continuum Hilbert-ball representation in
``ou3_sea3_continuum_driver_gram``. That representation preserves the complete
three-partition spectral family and the joint six-DOF response; it is not a
finite harmonic source, replay, Gaussian event, or independent sample box.

The finite-dimensional sampled physical target is the compact correlated
behavior ``B^601_SEA3`` defined by ``ou3_sea3_hard_window_behavior``. Numerical
closure still requires an outward validated complete-family Gram enclosure and
the same-history nonlinear physical output map consumable by the executor.
Physical SEA0->SEA3 left inclusion remains a separate deployment obligation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_sea3_continuum_phase_state as PHASE
import ou3_sea3_continuum_driver_gram as DRIVER
import ou3_sea3_hard_window_behavior as BEHAVIOR

REPO = Path(__file__).resolve().parents[2]
THEOREM = REPO / "doc" / "kalman_ou_iii" / "w3d-sea3-stability-theorem.tex-part"
COMPLETE_SOURCE = REPO / "tools" / "stability" / "ou3_sea3_complete_source.py"
SCHEMA = 5
QUALIFICATION = "OU3_SEA3_HARD_SHAPING_STATE_CONTRACT_V5"

CONTINUUM_PHASE_COORDINATE_SET_CLOSED = PHASE.CONTINUUM_PHASE_COORDINATE_SET_CLOSED
PHASE_CONTINUOUS_PROPAGATION_CLOSED = PHASE.PHASE_CONTINUOUS_PROPAGATION_CLOSED
HARD_SPECTRAL_DRIVER_SET_CLOSED = True
COMPLETE_SEA3_LEFT_INCLUSION_CLOSED = False
JOINT_SOURCE_OUTPUT_MAP_CLOSED = False

HARD_SHAPING_STATE_OR_EXCITATION_BOUND_CLOSED = all((
    CONTINUUM_PHASE_COORDINATE_SET_CLOSED,
    PHASE_CONTINUOUS_PROPAGATION_CLOSED,
    HARD_SPECTRAL_DRIVER_SET_CLOSED,
))
CONDITIONAL_SEA3_SOURCE_OUTPUT_CLOSED = all((
    HARD_SHAPING_STATE_OR_EXCITATION_BOUND_CLOSED,
    JOINT_SOURCE_OUTPUT_MAP_CLOSED,
))
DEPLOYMENT_LEFT_INCLUSION_CLOSED = COMPLETE_SEA3_LEFT_INCLUSION_CLOSED


def _normalized_text(text: str) -> str:
    return " ".join(text.split())


def build() -> dict:
    theorem = THEOREM.read_text(encoding="utf-8")
    theorem_flat = _normalized_text(theorem)
    complete = COMPLETE_SOURCE.read_text(encoding="utf-8")
    phase = PHASE.build()
    driver = DRIVER.build()
    behavior = BEHAVIOR.build()
    bad = {
        "phase": PHASE.validate(phase),
        "driver": DRIVER.validate(driver),
        "behavior": BEHAVIOR.validate(behavior),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"SEA3 shaping prerequisites failed: {bad}")

    theorem_has_shaping_system = (
        "x^s_{k+1}&=A_s" in theorem
        and "u^s_k&=C_s" in theorem
        and "oscillator/shaping state or an equivalent hard finite-window" in theorem_flat
    )
    theorem_has_explicit_hard_realization_set = (
        "\\mathcal X^s_{\\rm SEA3}(\\lambda_{0:N_W})" in theorem
        and "eq:sea3-hard-realization-set" in theorem
        and "does not reopen SEA3 compactness" in theorem_flat
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

    executable = {
        "continuum_phase_coordinate_set_closed": CONTINUUM_PHASE_COORDINATE_SET_CLOSED,
        "phase_continuous_propagation_closed": PHASE_CONTINUOUS_PROPAGATION_CLOSED,
        "hard_spectral_driver_set_closed": HARD_SPECTRAL_DRIVER_SET_CLOSED,
        "complete_SEA3_left_inclusion_closed": COMPLETE_SEA3_LEFT_INCLUSION_CLOSED,
        "joint_source_output_map_closed": JOINT_SOURCE_OUTPUT_MAP_CLOSED,
    }
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "SEA3_parameter_domain_compact": True,
        "compactness_is_not_an_open_obligation": True,
        "theorem_has_deterministic_shaping_contract": theorem_has_shaping_system,
        "theorem_has_explicit_hard_realization_set": theorem_has_explicit_hard_realization_set,
        "theorem_rejects_statistical_or_seeded_surrogates": theorem_rejects_finite_or_statistical_surrogates,
        "theorem_separates_probabilistic_random_sea_corollary": theorem_separates_probabilistic_corollary,
        "complete_source_rejects_gaussian_word_generator": source_rejects_gaussian_generator,
        "hard_realization_set_symbol": "X^s_SEA3(lambda_{0:N_W})",
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
        "hard_driver_certificate": {
            "qualification": driver["qualification"],
            "driver_space": driver["driver_space"],
            "driver_set": driver["driver_set"],
            "same_driver_field_entire_window": driver["same_driver_field_entire_window"],
            "same_driver_field_translation_and_rotation": driver["same_driver_field_translation_and_rotation"],
            "fixed_history_gram": driver["fixed_history_gram"],
            "fixed_history_membership": driver["fixed_history_membership"],
            "hard_spectral_driver_set_closed": driver["hard_spectral_driver_set_closed"],
            "exact_fixed_history_correlated_oracle_formula_closed": driver["exact_fixed_history_correlated_oracle_formula_closed"],
            "validated_complete_family_gram_enclosure_closed": driver["validated_complete_family_gram_enclosure_closed"],
        },
        "sampled_behavior_target": {
            "qualification": behavior["qualification"],
            "symbol": behavior["behavior_set_symbol"],
            "sample_count": behavior["sample_count"],
            "sampled_projection_dimension": behavior["sampled_projection_dimension"],
            "compact": behavior["sampled_behavior_set_compact"],
            "membership_requires_common_SEA3_witness": behavior["membership_requires_common_SEA3_witness"],
            "normal_live_caps_are_membership_sufficient": behavior["normal_live_caps_are_membership_sufficient"],
            "independent_sample_boxes_define_behavior_set": behavior["independent_sample_boxes_define_behavior_set"],
            "validated_membership_or_separation_oracle_closed": behavior["validated_membership_or_separation_oracle_closed"],
            "validated_correlated_outer_enclosure_closed": behavior["validated_correlated_outer_enclosure_closed"],
        },
        "power_spectrum_alone_is_hard_pathwise_bound": False,
        "spectral_moments_alone_may_close_xs": False,
        "gaussian_good_event_may_close_xs": False,
        "replay_may_close_xs": False,
        "seeded_128_frequency_generator_may_close_xs": False,
        "finite_RAO_grid_may_close_xs": False,
        "arbitrary_bounded_input_box_may_close_xs": False,
        "executable_ingredients": executable,
        "hard_shaping_state_or_excitation_bound_closed": HARD_SHAPING_STATE_OR_EXCITATION_BOUND_CLOSED,
        "conditional_SEA3_source_output_closed": CONDITIONAL_SEA3_SOURCE_OUTPUT_CLOSED,
        "deployment_left_inclusion_closed": DEPLOYMENT_LEFT_INCLUSION_CLOSED,
        "conditional_execution_requires_deployment_left_inclusion": False,
        "complete_SEA3_family_materialized_here": False,
        "P3_promoted": False,
        "next_obligation": "enclose the complete coupled 601-sample Gram family and propagate its same-history translation/rotation core through attitude/specific-force kinematics and the exact frontend/tuner/scheduler; do not add another source surrogate",
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        failures.append("schema/qualification mismatch")
    for key in (
        "SEA3_parameter_domain_compact",
        "compactness_is_not_an_open_obligation",
        "theorem_has_deterministic_shaping_contract",
        "theorem_has_explicit_hard_realization_set",
        "theorem_rejects_statistical_or_seeded_surrogates",
        "theorem_separates_probabilistic_random_sea_corollary",
        "complete_source_rejects_gaussian_word_generator",
        "hard_shaping_state_or_excitation_bound_closed",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    if d.get("hard_realization_set_symbol") != "X^s_SEA3(lambda_{0:N_W})":
        failures.append("hard realization set symbol drifted")
    phase = d.get("continuum_phase_certificate", {})
    for key in ("continuum_index_set_retained", "continuum_phase_coordinate_set_closed", "phase_continuous_propagation_closed"):
        if phase.get(key) is not True:
            failures.append(f"continuum phase certificate lost {key}")
    for key in ("finite_frequency_grid_used", "finite_direction_grid_used", "phase_reset_on_lambda_transition_allowed"):
        if phase.get(key) is not False:
            failures.append(f"continuum phase certificate reintroduced {key}")
    driver = d.get("hard_driver_certificate", {})
    for key in (
        "same_driver_field_entire_window",
        "same_driver_field_translation_and_rotation",
        "hard_spectral_driver_set_closed",
        "exact_fixed_history_correlated_oracle_formula_closed",
    ):
        if driver.get(key) is not True:
            failures.append(f"hard driver certificate lost {key}")
    if driver.get("validated_complete_family_gram_enclosure_closed") is not False:
        failures.append("hard driver falsely claims complete-family Gram enclosure")
    behavior = d.get("sampled_behavior_target", {})
    for key in ("compact", "membership_requires_common_SEA3_witness"):
        if behavior.get(key) is not True:
            failures.append(f"sampled behavior target lost {key}")
    for key in (
        "normal_live_caps_are_membership_sufficient",
        "independent_sample_boxes_define_behavior_set",
        "validated_membership_or_separation_oracle_closed",
        "validated_correlated_outer_enclosure_closed",
    ):
        if behavior.get(key) is not False:
            failures.append(f"sampled behavior target falsely closes/reintroduces {key}")
    for key in (
        "power_spectrum_alone_is_hard_pathwise_bound",
        "spectral_moments_alone_may_close_xs",
        "gaussian_good_event_may_close_xs",
        "replay_may_close_xs",
        "seeded_128_frequency_generator_may_close_xs",
        "finite_RAO_grid_may_close_xs",
        "arbitrary_bounded_input_box_may_close_xs",
        "conditional_SEA3_source_output_closed",
        "deployment_left_inclusion_closed",
        "conditional_execution_requires_deployment_left_inclusion",
        "complete_SEA3_family_materialized_here",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    expected = {
        "continuum_phase_coordinate_set_closed": True,
        "phase_continuous_propagation_closed": True,
        "hard_spectral_driver_set_closed": True,
        "complete_SEA3_left_inclusion_closed": False,
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
        "SEA3_compact": d["SEA3_parameter_domain_compact"],
        "hard_realization_set": d["hard_realization_set_symbol"],
        "executable_ingredients": d["executable_ingredients"],
        "hard_shaping_closed": d["hard_shaping_state_or_excitation_bound_closed"],
        "conditional_source_output_closed": d["conditional_SEA3_source_output_closed"],
        "deployment_left_inclusion_closed": d["deployment_left_inclusion_closed"],
        "next_obligation": d["next_obligation"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
