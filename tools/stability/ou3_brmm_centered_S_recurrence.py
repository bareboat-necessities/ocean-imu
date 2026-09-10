#!/usr/bin/env python3
"""Origin-invariant indefinite S recurrence qualification for COMPLETE-BRMM.

This is the positive continuation of the exact quiet-source obstruction in
``ou3_brmm_infinite_continuation.py``.  It does not restore the legacy 300 m*s
fresh-entry factor and it does not re-zero S at word boundaries.

For one physical history let S_dot=p and fix the real Live handoff time t_L.
PR #513 permits the proof-coordinate change S_L(t)=S(t)-S(t_L), so S_L(t_L)=0.
A finite all-time P4 tube cannot follow from bounded p and bounded finite-word
Delta-S alone.  The source must additionally control the *indefinite* displacement
primitive.  A simple physical sufficient condition is

    || integral_u^t p(s) ds || <= D_S       for all t,u >= t_L,

for one uniform finite D_S.  This is invariant to the absolute S origin and gives
||S_L(t)||<=D_S for every possible handoff.  It is stronger than the exact
necessary observation-equivalence-class condition, but unlike an absolute S ball
it is a genuine recurrence property of the physical history.

A numerical D_S is not guessed here.  The existing 300 m*s number is read only as
the current *working/retention* radius.  Sign symmetry of the norm-defined BRMM
source gives a necessary feasibility ceiling D_S<=R_work for the strengthened
symmetric source family: an indistinguishable p/-p pair can separate centered S
by 2 D_S while two errors inside the same radius-R_work tube can separate by at
most 2 R_work.  The real admissible D_S may be smaller after Joseph/reset, prefix
and storage margins are included.  It may be frozen only by that calculation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_brmm_contract as BRMM
import ou3_brmm_infinite_continuation as OBSTRUCTION

REPO = Path(__file__).resolve().parents[2]
CLOSURE = REPO / "tools" / "stability" / "ou3_p4_closure_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_BRMM_CENTERED_S_RECURRENCE_V2"


def build() -> dict:
    source = BRMM.build()
    obstruction = OBSTRUCTION.build()
    closure = json.loads(CLOSURE.read_text())
    if BRMM.validate(source):
        raise RuntimeError("BRMM source declaration failed")
    if not obstruction["bounded_all18_indefinite_target_refuted_under_finite_window_definition"]:
        raise RuntimeError("indefinite-S obstruction prerequisite lost")
    if obstruction["P3_delta"] != 1e-18:
        raise RuntimeError("canonical P3 changed")
    work = float(closure["hard_entry_search"]["base_coordinate_radii"]["integral_displacement_norm_m_s"])
    if not (work > 0.0):
        raise RuntimeError("positive S working radius required")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "S_dot_equals_p_same_history": True,
        "fresh_common_origin_reduction_preserved": True,
        "absolute_S_origin_is_not_bounded": True,
        "legacy_300_m_s_fresh_entry_ball_used": False,
        "position_reanchoring_used": False,
        "wordwise_rezero_of_S_used": False,
        "necessary_observation_class_condition": (
            "every admitted estimator-observation equivalence class must have "
            "uniformly bounded handoff-centered S diameter"
        ),
        "physical_sufficient_condition": (
            "exists uniform finite D_S such that norm(integral_u^t p(s) ds) "
            "<= D_S for every admitted history and all t,u>=t_L"
        ),
        "condition_is_origin_invariant": True,
        "condition_bounds_every_handoff_centered_S": True,
        "condition_excludes_nonzero_constant_position_history": True,
        "strengthened_source_is_sign_symmetric": True,
        "bounded_finite_word_DeltaS_alone_is_sufficient": False,
        "bounded_position_alone_is_sufficient": False,
        "finite_harmonic_zero_mean_position_is_nonempty_example": True,
        "finite_harmonic_bound_formula": (
            "for p=sum(A_i cos(w_i t)+B_i sin(w_i t)), "
            "D_S <= 2*sum((||A_i||+||B_i||)/w_i)"
        ),
        "S_working_radius_m_s": work,
        "D_S_retention_search_lower_m_s": 0.0,
        "D_S_retention_search_upper_m_s": work,
        "search_upper_is_working_tube_necessary_ceiling_not_source_assumption": True,
        "D_S_max_m_s": None,
        "D_S_numeric_qualification_closed": False,
        "COMPLETE_BRMM_INDEFINITE_S_QUALIFIED": False,
        "source_uniform_transition_cover_may_resume_after_D_S_qualification": True,
        "P3_delta": 1e-18,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "search D_S in [0,S_working_radius] with the same-history joint24 endpoint, "
            "literal-prefix and first-exit construction; freeze only the largest rigorously "
            "retained D_S, then resume universal source-cover promotion"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for k in (
        "S_dot_equals_p_same_history",
        "fresh_common_origin_reduction_preserved",
        "absolute_S_origin_is_not_bounded",
        "condition_is_origin_invariant",
        "condition_bounds_every_handoff_centered_S",
        "condition_excludes_nonzero_constant_position_history",
        "strengthened_source_is_sign_symmetric",
        "finite_harmonic_zero_mean_position_is_nonempty_example",
        "search_upper_is_working_tube_necessary_ceiling_not_source_assumption",
        "source_uniform_transition_cover_may_resume_after_D_S_qualification",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "legacy_300_m_s_fresh_entry_ball_used",
        "position_reanchoring_used",
        "wordwise_rezero_of_S_used",
        "bounded_finite_word_DeltaS_alone_is_sufficient",
        "bounded_position_alone_is_sufficient",
        "D_S_numeric_qualification_closed",
        "COMPLETE_BRMM_INDEFINITE_S_QUALIFIED",
        "P4_PASS",
        "P5_MAY_START",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    work = float(d.get("S_working_radius_m_s", -1.0))
    if work <= 0.0:
        f.append("working radius invalid")
    if d.get("D_S_retention_search_lower_m_s") != 0.0:
        f.append("D_S search lower bound changed")
    if d.get("D_S_retention_search_upper_m_s") != work:
        f.append("D_S search ceiling detached from working radius")
    if d.get("D_S_max_m_s") is not None:
        f.append("D_S_max must remain unfrozen until retention computes it")
    if d.get("P3_delta") != 1e-18:
        f.append("P3 delta changed")
    return f


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build()
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "centered_S_recurrence_materialized": True,
        "D_S_search_upper_m_s": d["D_S_retention_search_upper_m_s"],
        "D_S_numeric_closed": d["D_S_numeric_qualification_closed"],
        "P4": d["P4_PASS"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
