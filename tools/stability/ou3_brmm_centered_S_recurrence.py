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

The executable object below deliberately leaves D_S_max unfrozen.  A numerical
value may be promoted only after the source-cover/retention calculation computes
the largest value compatible with the certified P4 tube.  Thus this module
materializes the missing theorem coordinate and its exact implications without
silently shrinking the source family for proof convenience.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_brmm_contract as BRMM
import ou3_brmm_infinite_continuation as OBSTRUCTION

SCHEMA = 1
QUALIFICATION = "OU3_BRMM_CENTERED_S_RECURRENCE_V1"


def build() -> dict:
    source = BRMM.build()
    obstruction = OBSTRUCTION.build()
    if BRMM.validate(source):
        raise RuntimeError("BRMM source declaration failed")
    if not obstruction["bounded_all18_indefinite_target_refuted_under_finite_window_definition"]:
        raise RuntimeError("indefinite-S obstruction prerequisite lost")
    if obstruction["P3_delta"] != 1e-18:
        raise RuntimeError("canonical P3 changed")

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
        "bounded_finite_word_DeltaS_alone_is_sufficient": False,
        "bounded_position_alone_is_sufficient": False,
        "finite_harmonic_zero_mean_position_is_nonempty_example": True,
        "finite_harmonic_bound_formula": (
            "for p=sum(A_i cos(w_i t)+B_i sin(w_i t)), "
            "D_S <= 2*sum((||A_i||+||B_i||)/w_i)"
        ),
        "D_S_max_m_s": None,
        "D_S_numeric_qualification_closed": False,
        "COMPLETE_BRMM_INDEFINITE_S_QUALIFIED": False,
        "source_uniform_transition_cover_may_resume_after_D_S_qualification": True,
        "P3_delta": 1e-18,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "compute the largest admissible D_S from the same-history joint24 "
            "endpoint/prefix first-exit calculation; only then freeze COMPLETE-BRMM(D_S) "
            "and resume universal source-cover promotion"
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
        "finite_harmonic_zero_mean_position_is_nonempty_example",
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
        "D_S_numeric_closed": d["D_S_numeric_qualification_closed"],
        "P4": d["P4_PASS"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
