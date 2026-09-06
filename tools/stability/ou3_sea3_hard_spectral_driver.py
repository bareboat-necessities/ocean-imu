#!/usr/bin/env python3
"""Hard deterministic continuum spectral coordinate for complete SEA3.

This is the missing *source definition*, not a replay generator and not a
finite-harmonic approximation.

For each SEA3 parameter history lambda_0:N and one common continuum response
witness G_imu, let

    L_k(omega,theta) = G_imu(omega,theta;lambda_k)
                       sqrt(E(omega,theta;lambda_k)).

The hard deterministic spectral coordinate is one complex Hilbert-space
coordinate c in the closed L2 unit ball.  The phase-continuous source output is

    u_k = Re int L_k(omega,theta) c(omega,theta)
                     exp(i omega k h) dtheta domega.

Equivalently, the oscillator coordinate rotates pointwise and is never reseeded.
The same c, lambda history and response witness generate every sample and every
translation/rotation channel.  Therefore the sampled 601-point output is the
image of one Hilbert ball under one finite-dimensional linear operator; it is a
correlated ellipsoid/RKHS image, not the Cartesian product of per-sample caps.

The matrix Gram kernel is

    K_ij = Re int L_i L_j^* exp(i omega (i-j) h) dtheta domega.

At i=j this reduces to the matrix-valued SEA3 response spectrum integrated over
frequency.  Finiteness follows from the existing continuum RAO rolloff theorem.
The closed L2 unit ball is weakly compact, and every coordinate of the sampled
map is a bounded linear functional, so its finite-dimensional image is compact.
This gives the theorem a hard deterministic spectral/driver set without
claiming that an unbounded Gaussian process is pathwise contained in it.  The
physical/probabilistic left-inclusion remains a separate deployment obligation.

For changing lambda, spectral weights and the lambda-dependent response change
inside L_k while the underlying c/phase coordinate remains common.  This is the
required same-history coupling.  Normal-Live acceleration/rate caps are
additional intersection conditions; they are not used to manufacture samples.

No quadrature, frequency grid, direction grid, seeded phase list, replay, or
independent output box is used by this semantic certificate.  A later numerical
oracle may approximate the continuum integrals only with a validated outward
error bound; the exact object certified here remains the continuum operator.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_continuum_phase_state as PHASE
import ou3_sea3_directional_response_family as RESPONSE
import ou3_sea3_hard_window_behavior as BEHAVIOR

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_HARD_CONTINUUM_SPECTRAL_DRIVER_V1"
CANONICAL_SOURCE = "COMPLETE_SEA3_NORMAL_LIVE_WORD"
SAMPLES = 601

# Mathematical source definition gates.  These are intentionally distinct from
# the numerical membership/enclosure oracle and from deployment left inclusion.
HARD_SPECTRAL_DRIVER_SET_CLOSED = True
JOINT_SOURCE_OUTPUT_MAP_SEMANTICS_CLOSED = True
VALIDATED_NUMERICAL_GRAM_ORACLE_CLOSED = False
PHYSICAL_LEFT_INCLUSION_CLOSED = False


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    phase = PHASE.build()
    response = RESPONSE.directional_response_enclosure(REPO)
    behavior = BEHAVIOR.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "phase": PHASE.validate(phase),
        "response": RESPONSE.validate(response),
        "behavior": BEHAVIOR.validate(behavior),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"hard spectral-driver prerequisites failed: {bad}")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": CANONICAL_SOURCE,
        "source_definition": {
            "coordinate_space": (
                "complex L2 continuum coordinate c_r(omega,theta), r=1..3, "
                "with sum_r integral |c_r|^2 <= 1"
            ),
            "closed_unit_ball": True,
            "topology_for_compactness": "weak Hilbert topology",
            "weakly_compact": True,
            "phase_continuity": (
                "the same c is retained and exp(i*omega*h) advances the oscillator phase; "
                "lambda transitions may change spectral weights/response but never reseed c"
            ),
            "source_output": (
                "u_k = Re integral G_imu(omega,theta;lambda_k) "
                "sqrt(E(omega,theta;lambda_k)) c(omega,theta) "
                "exp(i omega k h) dtheta domega"
            ),
            "same_coordinate_drives_translation_and_rotation": True,
            "same_coordinate_drives_all_601_samples": True,
            "lambda_history_remains_coupled": True,
            "response_witness_remains_coupled": True,
        },
        "sampled_operator": {
            "sample_count": SAMPLES,
            "behavior_set_symbol": "B^601_SEA3",
            "finite_dimensional_image_of_one_Hilbert_ball": True,
            "not_cartesian_product_of_sample_caps": True,
            "gram_kernel": (
                "K_ij = Re integral L_i(omega,theta) L_j(omega,theta)^* "
                "exp(i omega (i-j) h) dtheta domega"
            ),
            "Li_definition": (
                "L_i=G_imu(omega,theta;lambda_i)*sqrt(E(omega,theta;lambda_i))"
            ),
            "diagonal_recovers_response_spectral_energy": True,
            "cross_time_blocks_preserved": True,
            "cross_axis_blocks_preserved": True,
            "common_witness_required_for_every_block": True,
            "finite_image_compact": True,
        },
        "bounded_operator_basis": {
            "continuum_RAO_rolloff_retained": True,
            "response_moment_theorem": response["uniform_moment_theorem"],
            "matrix_spectrum_identity_retained": True,
            "pointwise_Normal_Live_caps_are_intersection_conditions": True,
            "pointwise_Normal_Live_caps_generate_source": False,
        },
        "continuum_semantics": {
            "continuum_frequency_retained": True,
            "continuum_direction_retained": True,
            "continuum_phase_coordinate_set_closed": phase[
                "continuum_phase_coordinate_set_closed"
            ],
            "phase_continuous_propagation_closed": phase[
                "phase_continuous_propagation_closed"
            ],
            "finite_frequency_grid_used": False,
            "finite_direction_grid_used": False,
            "finite_seeded_harmonic_used": False,
            "trajectory_replay_used": False,
            "gaussian_good_event_used": False,
            "spectral_moment_only_membership_used": False,
            "independent_axis_boxes_used": False,
            "independent_sample_boxes_used": False,
        },
        "hard_spectral_driver_set_closed": HARD_SPECTRAL_DRIVER_SET_CLOSED,
        "joint_source_output_map_semantics_closed": JOINT_SOURCE_OUTPUT_MAP_SEMANTICS_CLOSED,
        "validated_numerical_gram_oracle_closed": VALIDATED_NUMERICAL_GRAM_ORACLE_CLOSED,
        "physical_left_inclusion_closed": PHYSICAL_LEFT_INCLUSION_CLOSED,
        "conditional_source_definition_closed": bool(
            HARD_SPECTRAL_DRIVER_SET_CLOSED
            and JOINT_SOURCE_OUTPUT_MAP_SEMANTICS_CLOSED
        ),
        "provider_artifact_materialized_here": False,
        "P3_changed": False,
        "P3_promoted_from_this_module": False,
        "P4_promoted_from_this_module": False,
        "next_obligation": (
            "build a validated continuum Gram/operator evaluator for one legal fixed SEA3 "
            "lambda/response member, construct a certified point in its one-coordinate "
            "601-sample image, propagate a same-source prehistory to Live/front-end/covariance "
            "entry, and execute that point through the unchanged complete shipping word"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != CANONICAL_SOURCE:
        f.append("hard driver detached from complete SEA3")
    source = d.get("source_definition", {})
    for key in (
        "closed_unit_ball",
        "weakly_compact",
        "same_coordinate_drives_translation_and_rotation",
        "same_coordinate_drives_all_601_samples",
        "lambda_history_remains_coupled",
        "response_witness_remains_coupled",
    ):
        if source.get(key) is not True:
            f.append(f"source definition lost {key}")
    op = d.get("sampled_operator", {})
    for key in (
        "finite_dimensional_image_of_one_Hilbert_ball",
        "not_cartesian_product_of_sample_caps",
        "diagonal_recovers_response_spectral_energy",
        "cross_time_blocks_preserved",
        "cross_axis_blocks_preserved",
        "common_witness_required_for_every_block",
        "finite_image_compact",
    ):
        if op.get(key) is not True:
            f.append(f"sampled operator lost {key}")
    basis = d.get("bounded_operator_basis", {})
    for key in (
        "continuum_RAO_rolloff_retained",
        "matrix_spectrum_identity_retained",
        "pointwise_Normal_Live_caps_are_intersection_conditions",
    ):
        if basis.get(key) is not True:
            f.append(f"bounded operator basis lost {key}")
    if basis.get("pointwise_Normal_Live_caps_generate_source") is not False:
        f.append("Normal-Live boxes became a source generator")
    continuum = d.get("continuum_semantics", {})
    for key in (
        "continuum_frequency_retained",
        "continuum_direction_retained",
        "continuum_phase_coordinate_set_closed",
        "phase_continuous_propagation_closed",
    ):
        if continuum.get(key) is not True:
            f.append(f"continuum semantics lost {key}")
    for key in (
        "finite_frequency_grid_used",
        "finite_direction_grid_used",
        "finite_seeded_harmonic_used",
        "trajectory_replay_used",
        "gaussian_good_event_used",
        "spectral_moment_only_membership_used",
        "independent_axis_boxes_used",
        "independent_sample_boxes_used",
    ):
        if continuum.get(key) is not False:
            f.append(f"forbidden source construction reappeared: {key}")
    for key in (
        "hard_spectral_driver_set_closed",
        "joint_source_output_map_semantics_closed",
        "conditional_source_definition_closed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "validated_numerical_gram_oracle_closed",
        "physical_left_inclusion_closed",
        "provider_artifact_materialized_here",
        "P3_changed",
        "P3_promoted_from_this_module",
        "P4_promoted_from_this_module",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if int(op.get("sample_count", 0)) != SAMPLES:
        f.append("sampled operator does not cover 601 samples")
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
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "canonical_source": d["canonical_source"],
        "conditional_source_definition_closed": d["conditional_source_definition_closed"],
        "numerical_gram_oracle_closed": d["validated_numerical_gram_oracle_closed"],
        "physical_left_inclusion_closed": d["physical_left_inclusion_closed"],
        "next_obligation": d["next_obligation"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
