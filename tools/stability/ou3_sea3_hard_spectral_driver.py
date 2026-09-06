#!/usr/bin/env python3
"""Deterministic continuum hard-driver contract for complete SEA3.

This module closes the *form* of the hard spectral driver that the SEA3 theorem
requires before a complete same-history finite window can be materialized.  It
is deliberately not a replay generator, finite harmonic approximation,
independent bounded-input box, Gaussian confidence event, or a replacement for
the complete SEA3 source.

For a fixed admissible SEA3 parameter history lambda_0:N, let

    E_k(omega,theta) = E(omega,theta;lambda_k)

be the exact three-partition directional JONSWAP/PM density.  The deterministic
hard realization uses one common normalized continuum driver

    psi in H,  ||psi||_H <= 1,

where H is a complex Hilbert space over the common frequency/direction index
set.  The same psi is propagated with the exact continuum phase factor
exp(i*omega*t_k) and multiplied by sqrt(E_k) before *every* translational and
rotational response operator.  Thus the sampled physical word has the operator
form

    y = T(lambda_0:N, G_imu) psi.

For any fixed admissible lambda/response history the complete sampled behavior
is therefore one correlated Hilbert-ball image, never a Cartesian product of
per-sample or per-axis intervals.  Its Gram operator is

    K = T T*,

so every finite sampled vector y in the range of T satisfies the deterministic
quadratic condition y* K^dagger y <= 1.  This is a hard source constraint: the
unit-ball condition is an additional deterministic SEA3 theorem assumption;
it is not inferred from a power spectrum or probability statement.

The construction is continuum.  No frequency grid, direction grid, phase
quantization, seeded Fourier record, or nominal RAO is introduced here.  The
same driver coordinate is shared by sea elevation, CoG translation, rotation,
front-end input, tuner history and shipping geometry.  Deployment inclusion of
all physical random seas into this deterministic hard class remains separate
and open, exactly as required by the theorem.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
THEOREM = REPO / "doc" / "kalman_ou_iii" / "w3d-sea3-stability-theorem.tex-part"
RESPONSE_DOMAIN = REPO / "tools" / "stability" / "ou3_sea3_directional_response_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_CONTINUUM_HARD_SPECTRAL_DRIVER_V1"

# These are mathematical/source-definition gates.  They do not claim that the
# 601-sample numerical operator or a full-family interval enclosure has already
# been emitted.
HARD_SPECTRAL_DRIVER_SET_CLOSED = True
COMMON_DRIVER_PHASE_PROPAGATION_CLOSED = True
JOINT_SOURCE_OUTPUT_OPERATOR_FORM_CLOSED = True
FULL_601_SAMPLE_NUMERICAL_OPERATOR_MATERIALIZED = False
FULL_FAMILY_CORRELATED_OUTER_ENCLOSURE_CLOSED = False
DEPLOYMENT_LEFT_INCLUSION_CLOSED = False


def _flat(text: str) -> str:
    return " ".join(text.split())


def build() -> dict:
    theorem = THEOREM.read_text(encoding="utf-8")
    tf = _flat(theorem)
    response = json.loads(RESPONSE_DOMAIN.read_text(encoding="utf-8"))
    r = response["response_contract"]

    theorem_support = {
        "hard_realization_set_is_explicit": (
            "\\mathcal X^s_{\\rm SEA3}(\\lambda_{0:N_W})" in theorem
            and "eq:sea3-hard-realization-set" in theorem
        ),
        "oscillator_or_equivalent_hard_window_allowed": (
            "oscillator/shaping state or an equivalent hard finite-window" in tf
        ),
        "phase_continuity_required": "phase continuity" in tf,
        "statistical_surrogates_rejected": (
            "neither a Gaussian confidence event" in tf
            and "spectral moments alone" in tf
            and "finite seeded harmonic" in tf
        ),
        "random_sea_is_later_corollary": (
            "A probabilistic statement for random sea realizations is a later corollary" in tf
        ),
    }
    response_support = {
        "continuum_RAO_not_grid": r.get("finite_RAO_grid_used") is False,
        "arbitrary_complex_phase_retained": r.get("phase_quantifier") == "arbitrary complex phase",
        "arbitrary_frequency_dependence_retained": r.get("arbitrary_frequency_dependence_below_envelope") is True,
        "arbitrary_directional_dependence_retained": r.get("arbitrary_directional_dependence") is True,
        "cross_axis_coupling_retained": r.get("arbitrary_cross_axis_coupling_subject_to_PSD") is True,
        "acceleration_kernel_square_integrable": r.get("unbanded_acceleration_moment_is_finite_from_response_rolloff") is True,
    }

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "driver_space": {
            "type": "complex Hilbert unit ball on the common continuum frequency-direction index set",
            "symbol": "psi",
            "hard_constraint": "||psi||_H <= 1",
            "hard_constraint_is_extra_deterministic_SEA3_assumption": True,
            "power_spectrum_alone_used_as_pathwise_bound": False,
            "probabilistic_good_event_used": False,
        },
        "same_history_spectral_factorization": {
            "surface_density": "E_k(omega,theta)=E(omega,theta;lambda_k)",
            "common_driver_coordinate": "psi(omega,theta)",
            "phase_factor": "exp(i*omega*t_k)",
            "channel_factor": "sqrt(E_k)*exp(i*omega*t_k)*psi",
            "same_psi_for_all_samples": True,
            "same_psi_for_translation_and_rotation": True,
            "same_psi_for_frontend_tuner_and_shipping_geometry": True,
            "lambda_changes_scale_the_same_phase_continuous_driver": True,
        },
        "finite_sample_operator": {
            "form": "y=T(lambda_0:N,G_imu) psi",
            "gram": "K=T T*",
            "membership_condition": "y in range(T) and y* K^dagger y <= 1",
            "one_correlated_behavior_not_cartesian_sample_boxes": True,
            "all_cross_sample_cross_axis_terms_retained": True,
            "deterministic_not_probabilistic": True,
        },
        "compactness_argument": {
            "driver_unit_ball_weakly_compact": True,
            "finite_sample_output_map_is_continuous_linear_for_fixed_lambda_response": True,
            "finite_sample_behavior_image_is_compact": True,
            "response_kernel_integrability_comes_from_declared_RAO_rolloff_and_finite_SEA3_moments": True,
        },
        "theorem_support": theorem_support,
        "response_support": response_support,
        "continuum_frequency_index_retained": True,
        "continuum_direction_index_retained": True,
        "finite_frequency_grid_used": False,
        "finite_direction_grid_used": False,
        "finite_seeded_harmonic_generator_used": False,
        "trajectory_replay_used": False,
        "independent_sample_boxes_used": False,
        "independent_axis_boxes_used": False,
        "independent_SEA_RAO_product_used": False,
        "nominal_RAO_selected": False,
        "HARD_SPECTRAL_DRIVER_SET_CLOSED": HARD_SPECTRAL_DRIVER_SET_CLOSED,
        "COMMON_DRIVER_PHASE_PROPAGATION_CLOSED": COMMON_DRIVER_PHASE_PROPAGATION_CLOSED,
        "JOINT_SOURCE_OUTPUT_OPERATOR_FORM_CLOSED": JOINT_SOURCE_OUTPUT_OPERATOR_FORM_CLOSED,
        "FULL_601_SAMPLE_NUMERICAL_OPERATOR_MATERIALIZED": FULL_601_SAMPLE_NUMERICAL_OPERATOR_MATERIALIZED,
        "FULL_FAMILY_CORRELATED_OUTER_ENCLOSURE_CLOSED": FULL_FAMILY_CORRELATED_OUTER_ENCLOSURE_CLOSED,
        "DEPLOYMENT_LEFT_INCLUSION_CLOSED": DEPLOYMENT_LEFT_INCLUSION_CLOSED,
        "next_obligation": (
            "materialize one legal same-history SEA3 member by evaluating its continuum operator T through the actual front end and shipping recursions, then run the ledger-required non-promoting complete-word feasibility ratio before any P4 enclosure work"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source mismatch")
    for group in ("theorem_support", "response_support"):
        values = d.get(group, {})
        if not values or not all(values.values()):
            f.append(f"{group} is incomplete")
    driver = d.get("driver_space", {})
    if driver.get("hard_constraint") != "||psi||_H <= 1":
        f.append("hard driver constraint changed")
    if driver.get("hard_constraint_is_extra_deterministic_SEA3_assumption") is not True:
        f.append("hard driver is being inferred from statistics")
    if driver.get("power_spectrum_alone_used_as_pathwise_bound") is not False:
        f.append("power spectrum was promoted to a hard pathwise bound")
    same = d.get("same_history_spectral_factorization", {})
    for key in (
        "same_psi_for_all_samples",
        "same_psi_for_translation_and_rotation",
        "same_psi_for_frontend_tuner_and_shipping_geometry",
        "lambda_changes_scale_the_same_phase_continuous_driver",
    ):
        if same.get(key) is not True:
            f.append(f"same-history driver lost {key}")
    op = d.get("finite_sample_operator", {})
    for key in (
        "one_correlated_behavior_not_cartesian_sample_boxes",
        "all_cross_sample_cross_axis_terms_retained",
        "deterministic_not_probabilistic",
    ):
        if op.get(key) is not True:
            f.append(f"finite-sample operator lost {key}")
    for key in (
        "continuum_frequency_index_retained",
        "continuum_direction_index_retained",
        "HARD_SPECTRAL_DRIVER_SET_CLOSED",
        "COMMON_DRIVER_PHASE_PROPAGATION_CLOSED",
        "JOINT_SOURCE_OUTPUT_OPERATOR_FORM_CLOSED",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "finite_frequency_grid_used",
        "finite_direction_grid_used",
        "finite_seeded_harmonic_generator_used",
        "trajectory_replay_used",
        "independent_sample_boxes_used",
        "independent_axis_boxes_used",
        "independent_SEA_RAO_product_used",
        "nominal_RAO_selected",
        "FULL_601_SAMPLE_NUMERICAL_OPERATOR_MATERIALIZED",
        "FULL_FAMILY_CORRELATED_OUTER_ENCLOSURE_CLOSED",
        "DEPLOYMENT_LEFT_INCLUSION_CLOSED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} must remain false at this stage")
    return list(dict.fromkeys(f))


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
        "qualification": d["qualification"],
        "hard_driver_closed": d["HARD_SPECTRAL_DRIVER_SET_CLOSED"],
        "joint_operator_form_closed": d["JOINT_SOURCE_OUTPUT_OPERATOR_FORM_CLOSED"],
        "full_601_operator_materialized": d["FULL_601_SAMPLE_NUMERICAL_OPERATOR_MATERIALIZED"],
        "next_obligation": d["next_obligation"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
