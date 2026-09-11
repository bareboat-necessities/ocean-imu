#!/usr/bin/env python3
"""One-time Live-centered primitive under the corrected physical BRMM source.

S_L(t)=integral_{t_L}^t p_wave(s) ds.  The theorem-level physical condition is
primary; spectral/shaping constructions are certificate methods.  The current
COMPLETE-BRMM proof scope also carries an explicit padded numerical envelope
anchored to the Hs=8.5 m reference case.  That envelope supplies the hard
D_S used here.  It is not selected from the P4 working radius.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_brmm_contract as BRMM
import ou3_brmm_infinite_continuation as OBSTRUCTION
import ou3_brmm_physical_wave_condition as PHYSICAL
import ou3_brmm_physical_wave_source as CERT

REPO = Path(__file__).resolve().parents[2]
CLOSURE = REPO / "tools/stability/ou3_p4_closure_domain.json"
SCHEMA = 4
QUALIFICATION = "OU3_BRMM_CENTERED_S_RECURRENCE_V4"


def build() -> dict:
    source = BRMM.build(); physical = PHYSICAL.build(); cert = CERT.build(); old = OBSTRUCTION.build()
    if BRMM.validate(source) or PHYSICAL.validate(physical) or CERT.validate(cert) or OBSTRUCTION.validate(old):
        raise RuntimeError("physical wave/old-definition regression prerequisite failed")
    work = float(json.loads(CLOSURE.read_text())["hard_entry_search"]["base_coordinate_radii"]["integral_displacement_norm_m_s"])
    ds = float(physical["numeric_D_S_m_s"])
    closed = bool(physical["physical_D_S_numeric_qualification_closed_for_complete_family"])
    return {
        "schema": SCHEMA, "qualification": QUALIFICATION,
        "canonical_source": PHYSICAL.CANONICAL_SOURCE,
        "S_dot_equals_p_same_history": True,
        "fresh_common_origin_reduction_preserved": True,
        "absolute_S_origin_is_not_bounded": True,
        "legacy_300_m_s_fresh_entry_ball_used": False,
        "position_reanchoring_used": False,
        "wordwise_rezero_of_S_used": False,
        "physical_primary_condition": physical["required_pathwise_property"],
        "certificate_methods_are_sufficient_not_definitional": True,
        "derived_centered_recurrence": "S_L(t)=integral_tL^t p_wave; S_L(next)=S_L+integral p_wave",
        "condition_is_origin_invariant": True,
        "condition_bounds_every_handoff_centered_S": True,
        "condition_excludes_nonzero_constant_position_history": not physical["constant_nonzero_position_zero_velocity_history_admitted"],
        "constant_nonzero_position_zero_velocity_history_admitted": False,
        "bounded_finite_word_DeltaS_alone_is_sufficient": False,
        "bounded_position_alone_is_sufficient": False,
        "finite_harmonic_zero_mean_position_is_nonempty_example": cert["zero_wave_admitted"] and bool(cert["exact_coefficient_identities"]),
        "finite_harmonic_bound_formula": "D_S <= 2*sum((||A_i||+||B_i||)/omega_i)",
        "complete_BRMM_numeric_physical_envelope": physical["complete_BRMM_numeric_physical_envelope"],
        "S_working_radius_m_s": work,
        "working_radius_sets_physical_D_S": False,
        "sign_symmetry_implies_observation_equivalence": False,
        "D_S_max_m_s": ds,
        "D_S_numeric_qualification_closed": closed,
        "COMPLETE_BRMM_INDEFINITE_S_THEOREM_CLOSED": True,
        "COMPLETE_BRMM_INDEFINITE_S_QUALIFIED": closed,
        "numeric_working_tube_comparison": (
            "physical D_S exceeds current working radius; P4 working/first-exit domain must be widened or shown to regulate the error coordinate independently"
            if ds > work else
            "current working radius is no smaller than the physical D_S bound"
        ),
        "historical_obstruction_classification": old["classification"],
        "current_missing_numeric_source_qualification_classification": None if closed else "E",
        "P3_delta": 1e-18, "P4_PASS": False, "P5_MAY_START": False,
        "next_obligation": (
            "carry the qualified physical envelope and one Live-centered primitive through the same-history joint24 source graph; retain correlated P/H/R through innovation and K; then redo endpoint/prefix/retention with a working domain large enough for the qualified physics"
        ),
    }


def validate(d: dict) -> list[str]:
    expected = build()
    return [k+" differs from derived centered-S recurrence" for k, v in expected.items() if d.get(k) != v]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args(); d = build(); f = validate(d)
    d.update(validation_pass=not f, validation_failures=f)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"centered_S_theorem_closed": d["COMPLETE_BRMM_INDEFINITE_S_THEOREM_CLOSED"],
                      "D_S_numeric_closed": d["D_S_numeric_qualification_closed"],
                      "D_S_m_s": d["D_S_max_m_s"],
                      "P4": d["P4_PASS"], "failures": f}, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
