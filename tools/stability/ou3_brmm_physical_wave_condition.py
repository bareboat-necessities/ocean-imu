#!/usr/bin/env python3
"""Primary physical COMPLETE-BRMM wave condition.

The theorem is about *physical wave histories*, not about a preferred signal
construction.  A wave displacement p_wave is defined relative to a local
slowly varying equilibrium/reference.  Current, propulsion, leeway, secular
translation and arbitrary global-position origin offsets belong to the
reference coordinate, not to p_wave.

Every admitted complete wave history must satisfy the hard pathwise condition

    p_wave is an oscillatory bounded displacement about that local equilibrium,
    v_wave = d p_wave / dt,
    a_wave = d v_wave / dt,
    and there is a finite family constant D_S such that

        || integral_u^t p_wave(s) ds || <= D_S

    for every relevant u,t on every successive admissible continuation.

The bounded primitive is part of the physical BRMM condition because it rules
out a permanent displacement DC component without banning finite flat segments
or quiet water.  It is not an arbitrary S cap and it is not chosen from the P4
working radius.

Harmonic/spectral amplitude measures and bounded shaping-state realizations are
only sufficient *certificate methods* for proving that a concrete wave family
satisfies this physical condition and for deriving D_S.  They do not define
what wave motion is.  Other deterministic representations are admissible when
they prove the same hard physical condition and preserve the same-history
p/v/a/primitive ancestry.
"""
from __future__ import annotations

import json
from fractions import Fraction as F
from pathlib import Path
from typing import Sequence

import ou3_brmm_physical_wave_source as CERT

QUALIFICATION = "OU3_BRMM_PHYSICAL_WAVE_CONDITION_V1"
CANONICAL_SOURCE = "COMPLETE_BRMM_NORMAL_LIVE_WORD"


def _r(x) -> F:
    if isinstance(x, (bool, float)):
        raise TypeError("physical hard bounds use exact integer/Fraction/rational strings")
    return F(x)


def physical_condition() -> dict:
    """Machine-readable theorem-level physics, independent of certificate type."""
    return {
        "qualification": QUALIFICATION,
        "canonical_source": CANONICAL_SOURCE,
        "physics_is_primary": True,
        "wave_coordinate": "p_wave=x_CoG-x_local_equilibrium",
        "absolute_global_position_is_not_wave_coordinate": True,
        "slow_translation_outside_wave_coordinate": [
            "current", "propulsion", "leeway", "secular_drift",
            "arbitrary_global_origin", "slow_reference_translation",
        ],
        "required_jet_relation": "v_wave=dp_wave/dt; a_wave=dv_wave/dt",
        "required_pathwise_property": (
            "bounded oscillatory wave displacement about local equilibrium with "
            "uniformly bounded centered primitive"
        ),
        "bounded_primitive": "for all u,t: ||integral_u^t p_wave(s) ds|| <= D_S",
        "D_S_must_be_derived_from_physical_family": True,
        "D_S_may_be_chosen_from_P4_working_radius": False,
        "zero_mean_alone_is_sufficient": False,
        "power_spectrum_alone_is_sufficient": False,
        "finite_replay_alone_is_sufficient": False,
        "quiet_zero_wave_is_admissible": True,
        "finite_constant_position_segment_is_blanket_forbidden": False,
        "indefinite_nonzero_position_DC_is_admissible": False,
        "same_history_ancestry_required": True,
        "wordwise_rezero_of_S_allowed": False,
        "position_reanchor_allowed": False,
        "certificate_methods_are_sufficient_not_definitional": True,
        "currently_implemented_certificate_methods": [
            "hard_positive_frequency_amplitude_measure",
            "bounded_shaping_state_with_exact_potential_identity",
        ],
    }


def qualify_with_certificate(cert: dict) -> dict:
    """Prove the primary physical condition using one supported certificate method."""
    failures = CERT.verify_certificate(cert)
    if failures:
        raise ValueError("wave-family certificate does not prove physical condition: " + repr(failures))
    ds = _r(cert["D_S_upper_m_s"])
    if ds < 0:
        raise ValueError("derived D_S must be nonnegative")
    return {
        "physical_condition": physical_condition(),
        "certificate_kind": cert["kind"],
        "derived_D_S_m_s": str(ds),
        "bounded_primitive_proved": True,
        "certificate_is_definition_of_wave_motion": False,
        "same_history_p_v_a_primitive_required": True,
        "P4_working_radius_used_to_choose_D_S": False,
    }


def constant_history_admitted(displacement: Sequence, cert: dict) -> bool:
    """Regression for PR #515: derive rejection from the bounded-wave condition."""
    q = qualify_with_certificate(cert)
    d = tuple(_r(x) for x in displacement)
    if len(d) != 3:
        raise ValueError("displacement must have three components")
    if all(x == 0 for x in d):
        return True
    # A permanent d != 0 has ||int_0^t p ds|| growing linearly and violates
    # every finite derived D_S.  CERT computes the exact contradiction horizon.
    return bool(CERT.constant_history_admission(d, cert)["admitted_by_physical_wave_condition"])


def build() -> dict:
    example = CERT.spectral_certificate((CERT.SpectralBand(F(1, 10), F(2), F(3)),))
    qualification = qualify_with_certificate(example)
    return {
        **physical_condition(),
        "certificate_example": qualification,
        "constant_nonzero_position_zero_velocity_history_admitted":
            constant_history_admitted((F(1, 8), 0, 0), example),
        "historical_old_finite_window_classification": "B",
        "physical_source_specification_omission_classification": "E",
        "shipping_filter_instability_claimed": False,
        "physical_D_S_numeric_qualification_closed_for_complete_family": False,
        "P3_delta": 1e-18,
        "P4_PASS": False,
        "P5_MAY_START": False,
    }


def validate(d: dict) -> list[str]:
    expected = build()
    return [k + " differs from physical BRMM condition" for k, v in expected.items()
            if d.get(k) != v]


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "physics_is_primary": d["physics_is_primary"],
        "dc_history_admitted": d["constant_nonzero_position_zero_velocity_history_admitted"],
        "P4_PASS": d["P4_PASS"],
        "failures": failures,
    }, sort_keys=True))
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
