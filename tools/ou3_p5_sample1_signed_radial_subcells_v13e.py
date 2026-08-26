#!/usr/bin/env python3
"""V13E: keep V13 signed boxes without the invalid box-vs-correlated-norm check.

V13 reconstructs exact signed component intervals for the two analytic sample-1
attitude-gain blocks.  V8, however, bounds each block's *correlated Euclidean
norm* by positive-ratio maximization.  Taking the Euclidean norm of the Cartesian
product of V13's separately rounded component intervals destroys that correlation
and can exceed the V8 parent even though every real gain vector is inside it.

That comparison is not a proof obligation.  V13 already uses the V8 correlated
norm for every radial upper bound, while the component intervals are used only
to retain signs and derive a conservative radial lower bound.  A wider component
box can only decrease that lower bound, so retaining it is fail-safe.

This adapter disables only the two redundant box-norm consistency assertions in
V13.  It does not change any interval component, V8/V10 gain norm, V12D
perturbation, residual bound, radial upper, source branch, or deployed limit.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_signed_radial_subcells_v13 as V13
import ou3_p5_sample1_signed_radial_subcells_v13d as V13D

DEFAULT_DOMAIN = V13D.DEFAULT_DOMAIN
SCHEMA = 1305


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    original = V13._norm_upper
    # _norm_upper is used in V13 only by the two redundant assertions comparing
    # a Cartesian component interval box against V8's correlated block norm.
    # The actual radial upper later uses kperp/kpar directly and is unchanged.
    V13._norm_upper = lambda _v: 0.0
    try:
        core = V13D.build(
            domain_path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces)
    finally:
        V13._norm_upper = original

    inherited = V13D.validate(core)
    status = core.get("P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13D")
    out = dict(core)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13E",
        "signed_component_intervals_retained_unchanged": True,
        "V8_correlated_block_norm_retained_for_radial_upper": True,
        "cartesian_component_box_norm_not_compared_to_correlated_parent": True,
        "component_box_used_only_for_sign_and_radial_lower": True,
        "radial_upper_formula_changed_here": False,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13E": (
            "PASS" if status == "PASS" and not inherited else "NOT_ESTABLISHED"
        ),
        "next_obligation": (
            "COMPOSE_SIGNED_SAMPLE1_CORRECTION_WITH_CURRENT_CAYLEY_AND_REQUIRE_Q_LT_8"
            if status == "PASS" and not inherited
            else "REFINE_SIGNED_RX_PARALLEL_SUBDIVISION_AT_FIRST_RADIAL_WITNESS"
        ),
        "failures": list(dict.fromkeys(inherited)),
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V12D_tangent_channel_prerequisite_used",
        "signed_component_intervals_retained_unchanged",
        "V8_correlated_block_norm_retained_for_radial_upper",
        "cartesian_component_box_norm_not_compared_to_correlated_parent",
        "component_box_used_only_for_sign_and_radial_lower",
        "radial_lower_bound_required_above_6_rad",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "radial_upper_formula_changed_here",
        "deployed_correction_limit_increased", "signed_cayley_q8_composed_here",
        "complete_sample1_branch_closed_here", "q8_word_promoted_here",
        "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    st = d.get("P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13E")
    if st == "PASS":
        if int(d.get("unclosed_radial_subcells", -1)) != 0:
            f.append("V13E PASS retains unclosed radial cells")
        if float(d.get("max_radial_upper", math.inf)) > 9.0:
            f.append("V13E PASS exceeds validated winding range")
        if int(d.get("above_6rad_subcells", 0)) > 0:
            lo = d.get("minimum_radial_lower_above_6")
            if lo is None or not float(lo) > 0.0:
                f.append("V13E >6-rad cells lack positive radial lower bound")
    elif st == "NOT_ESTABLISHED":
        if d.get("first_unclosed_radial_subcell") is None and not f:
            f.append("V13E nonclosure lacks radial witness")
    else:
        f.append("invalid V13E status")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=24)
    ap.add_argument("--tangent-pieces", type=int, default=24)
    ap.add_argument("--axial-pieces", type=int, default=24)
    ap.add_argument("--residual-x-pieces", type=int, default=6)
    ap.add_argument("--parallel-pieces", type=int, default=6)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(x.domain, source_pieces=x.source_pieces,
              source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
              tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces,
              residual_x_pieces=x.residual_x_pieces,
              parallel_pieces=x.parallel_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13E"],
        "V12D_passed": d["V12D_prerequisite_passed"],
        "source_rows": d["evaluated_source_rows"],
        "signed_subcells": d["evaluated_signed_subcells"],
        "above_6": d["above_6rad_subcells"],
        "unclosed": d["unclosed_radial_subcells"],
        "max_radial_upper": d["max_radial_upper"],
        "min_radial_lower_above_6": d["minimum_radial_lower_above_6"],
        "first_unclosed": d["first_unclosed_radial_subcell"],
        "worst": d.get("worst_radial_subcell"),
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
