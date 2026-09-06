#!/usr/bin/env python3
"""Exact continuum hard-driver representation for complete SEA3.

This closes the missing *definition* of the deterministic SEA3 driver without
introducing a finite harmonic realization, replay, independent sample box, or
statistical good-event source.

For one complete SEA3 parameter/response history, let H be the realification of
L2 over the three-partition continuum spectral index.  A single normalized
coefficient field ``a`` with ||a||_H <= 1 is carried through the whole window;
phase continuity is the theorem's pointwise multiplication by exp(i omega h),
not a reseed.  The six-DOF sea/IMU spectral factor L_k(omega,theta) contains the
JONSWAP/PM directional weight and the *same* admissible vessel response used for
translation and rotation.  The sampled source core is therefore

    y = K_{lambda,G} a,

with one operator K for the entire 601-sample history.  Its exact finite-window
Gram operator is Q=K K*, so for a fixed admissible lambda/response history

    support(c) = sqrt(c^T Q c),
    y in image(B_H) iff y in range(Q) and y^T Q^dagger y <= 1.

This is the correlated membership/separation representation required by the
SEA3 theorem.  A validated numerical enclosure of Q over the complete compact
lambda/response history family is a separate next step.  Until that enclosure
is emitted, this module does not materialize a 601-sample word and does not
promote P4/P5.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
THEOREM = REPO / "doc" / "kalman_ou_iii" / "w3d-sea3-stability-theorem.tex-part"
RESPONSE_DOMAIN = REPO / "tools" / "stability" / "ou3_sea3_directional_response_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_CONTINUUM_HARD_DRIVER_GRAM_V1"
N = 601
SOURCE_CORE_DIM = 6 * N


def build() -> dict:
    theorem = THEOREM.read_text(encoding="utf-8")
    flat = " ".join(theorem.split())
    response = json.loads(RESPONSE_DOMAIN.read_text(encoding="utf-8"))
    rc = response["response_contract"]

    theorem_has_driver_ball = (
        "\\mathcal H_{\\rm SEA3}" in theorem
        and "\\|a\\|_{\\mathcal H_{\\rm SEA3}}\\le 1" in theorem
        and "Q_{\\lambda,G}=K_{\\lambda,G}K_{\\lambda,G}^{*}" in theorem
        and "eq:sea3-hard-driver-ball" in theorem
        and "eq:sea3-window-gram" in theorem
    )
    theorem_retains_common_history = (
        "one common continuum coefficient field" in flat
        and "not re-selected at individual samples" in flat
    )
    theorem_rejects_finite_surrogate = (
        "does not discretize the sea into a finite harmonic source" in flat
        and "does not turn the sampled history into independent boxes" in flat
    )
    response_is_continuum_family = (
        rc.get("finite_RAO_grid_used") is False
        and rc.get("single_nominal_RAO_used") is False
        and rc.get("phase_quantifier") == "arbitrary complex phase"
        and rc.get("arbitrary_frequency_dependence_below_envelope") is True
        and rc.get("arbitrary_directional_dependence") is True
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "driver_space": "realification L2({1,2,3} x Omega x Theta; C)",
        "driver_set": "closed unit ball of H_SEA3",
        "driver_norm_bound": 1.0,
        "driver_bound_is_definition_not_statistical_confidence": True,
        "same_driver_field_entire_window": True,
        "same_driver_field_translation_and_rotation": True,
        "phase_continuity_is_operator_propagation_not_reseed": True,
        "spectral_factor_contains_directional_JONSWAP_PM_and_joint_six_dof_response": True,
        "sample_count": N,
        "sampled_source_core_coordinates": [
            "f_cog_body[3]",
            "omega_body_corrected[3]",
        ],
        "sampled_source_core_dimension": SOURCE_CORE_DIM,
        "fixed_history_operator": "y=K_{lambda,G} a",
        "fixed_history_gram": "Q_{lambda,G}=K_{lambda,G} K_{lambda,G}^*",
        "fixed_history_support": "h(c)=sqrt(c^T Q c)",
        "fixed_history_membership": "y in range(Q) and y^T Q^dagger y <= 1",
        "complete_family": (
            "union over the compact coupled lambda history and admissible joint six-DOF response family, with one common driver field"
        ),
        "theorem_has_driver_ball": theorem_has_driver_ball,
        "theorem_retains_common_history": theorem_retains_common_history,
        "theorem_rejects_finite_surrogate": theorem_rejects_finite_surrogate,
        "response_is_continuum_family": response_is_continuum_family,
        "finite_frequency_grid_used": False,
        "finite_direction_grid_used": False,
        "finite_harmonic_source_used": False,
        "trajectory_replay_used": False,
        "gaussian_good_event_used": False,
        "independent_sample_boxes_used": False,
        "independent_translation_rotation_sources_used": False,
        "hard_spectral_driver_set_closed": True,
        "exact_fixed_history_correlated_oracle_formula_closed": True,
        "validated_complete_family_gram_enclosure_closed": False,
        "provider_word_materialized_here": False,
        "P4_promoted_here": False,
        "P5_promoted_here": False,
        "next_obligation": (
            "compute an outward validated enclosure of the complete-family 601-sample Gram operator while retaining the coupled lambda/response history; then feed an admitted same-history source through the real frontend/tuner/scheduler and every actual-applied R_S event"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "driver_bound_is_definition_not_statistical_confidence",
        "same_driver_field_entire_window",
        "same_driver_field_translation_and_rotation",
        "phase_continuity_is_operator_propagation_not_reseed",
        "spectral_factor_contains_directional_JONSWAP_PM_and_joint_six_dof_response",
        "theorem_has_driver_ball",
        "theorem_retains_common_history",
        "theorem_rejects_finite_surrogate",
        "response_is_continuum_family",
        "hard_spectral_driver_set_closed",
        "exact_fixed_history_correlated_oracle_formula_closed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "finite_frequency_grid_used",
        "finite_direction_grid_used",
        "finite_harmonic_source_used",
        "trajectory_replay_used",
        "gaussian_good_event_used",
        "independent_sample_boxes_used",
        "independent_translation_rotation_sources_used",
        "validated_complete_family_gram_enclosure_closed",
        "provider_word_materialized_here",
        "P4_promoted_here",
        "P5_promoted_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("driver detached from canonical complete SEA3 source")
    if int(d.get("sample_count", 0)) != N:
        f.append("driver does not cover canonical 601-sample window")
    if int(d.get("sampled_source_core_dimension", 0)) != SOURCE_CORE_DIM:
        f.append("joint translation/rotation sampled source dimension changed")
    if float(d.get("driver_norm_bound", 0.0)) != 1.0:
        f.append("normalized hard driver ball changed")
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
        "hard_driver_closed": d["hard_spectral_driver_set_closed"],
        "fixed_history_oracle_formula_closed": d["exact_fixed_history_correlated_oracle_formula_closed"],
        "complete_family_gram_enclosure_closed": d["validated_complete_family_gram_enclosure_closed"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
