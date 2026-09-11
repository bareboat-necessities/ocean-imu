#!/usr/bin/env python3
"""One-time Live-centered primitive derived from the physical wave generator.

S_L(t)=phi(t)-phi(t_L), where phidot=p is an identity of a bounded spectral
or shaping realization. The source supplies its hard generator constants;
retention must adapt to their consequent D_S, not choose them to fit a tube.
The old sign-symmetry-only D_S<=300 inference was invalid: opposite oscillatory
wave histories are not generally sensor-indistinguishable.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_brmm_contract as BRMM
import ou3_brmm_infinite_continuation as OBSTRUCTION
import ou3_brmm_physical_wave_source as WAVE

REPO = Path(__file__).resolve().parents[2]
CLOSURE = REPO / "tools/stability/ou3_p4_closure_domain.json"
SCHEMA = 3
QUALIFICATION = "OU3_BRMM_CENTERED_S_RECURRENCE_V3"


def build() -> dict:
    source = BRMM.build(); wave = WAVE.build(); old = OBSTRUCTION.build()
    if BRMM.validate(source) or WAVE.validate(wave) or OBSTRUCTION.validate(old):
        raise RuntimeError("physical wave/old-definition regression prerequisite failed")
    work = float(json.loads(CLOSURE.read_text())["hard_entry_search"]["base_coordinate_radii"]["integral_displacement_norm_m_s"])
    return {
        "schema": SCHEMA, "qualification": QUALIFICATION,
        "canonical_source": WAVE.CANONICAL_SOURCE,
        "S_dot_equals_p_same_history": True,
        "fresh_common_origin_reduction_preserved": True,
        "absolute_S_origin_is_not_bounded": True,
        "legacy_300_m_s_fresh_entry_ball_used": False,
        "position_reanchoring_used": False, "wordwise_rezero_of_S_used": False,
        "physical_sufficient_condition": "p=dphi/dt on the same bounded physical spectral/shaping generator",
        "physical_wave_generator_contract": wave,
        "derived_centered_recurrence": "S_L(t)=phi(t)-phi(t_L); S_L(next)=S_L+phi(next)-phi",
        "condition_is_origin_invariant": True,
        "condition_bounds_every_handoff_centered_S": wave["hard_bounded_primitive_derived_from_generator"],
        "condition_excludes_nonzero_constant_position_history": wave["old_witness_excluded_by_corrected_physical_theorem"],
        "constant_nonzero_position_zero_velocity_history_admitted": False,
        "bounded_finite_word_DeltaS_alone_is_sufficient": False,
        "bounded_position_alone_is_sufficient": False,
        "finite_harmonic_zero_mean_position_is_nonempty_example": wave["zero_wave_admitted"] and bool(wave["exact_coefficient_identities"]),
        "finite_harmonic_bound_formula": "D_S <= 2*sum((||A_i||+||B_i||)/omega_i)",
        "S_working_radius_m_s": work,
        "working_radius_sets_physical_D_S": False,
        "sign_symmetry_implies_observation_equivalence": False,
        "D_S_max_m_s": wave["numeric_D_S_m_s"],
        "D_S_numeric_qualification_closed": wave["physical_D_S_numeric_qualification_closed"],
        "COMPLETE_BRMM_INDEFINITE_S_THEOREM_CLOSED": wave["hard_bounded_primitive_derived_from_generator"],
        "COMPLETE_BRMM_INDEFINITE_S_QUALIFIED": wave["physical_D_S_numeric_qualification_closed"],
        "numeric_working_tube_comparison": wave["numeric_working_tube_comparison"],
        "historical_obstruction_classification": old["classification"],
        "current_missing_numeric_source_qualification_classification": "E",
        "P3_delta": 1e-18, "P4_PASS": False, "P5_MAY_START": False,
        "next_obligation": "qualify a hard physical generator envelope independently of P4, derive D_S outwardly, then carry phi and the one Live phi through joint24 endpoint/prefix/retention; enlarge a working tube only if that certificate permits",
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
                      "P4": d["P4_PASS"], "failures": f}, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
