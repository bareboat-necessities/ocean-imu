#!/usr/bin/env python3
"""V14B fix for the OU-III P5 sample-1 signed Cayley q<8 composition.

The authoritative V14 run reached the signed product calculation only to fail
because its helper delegated a two-component yz vector to a legacy norm routine
that intentionally accepts exactly three components.  This wrapper preserves
V14's source family, signed/radial correction cells, quaternion algebra, q<8
target, and all no-promotion guards.  It replaces only that helper by the
outward-rounded Euclidean norm valid for any finite interval vector.

No proof bound is reduced beyond evaluating the intended two-dimensional norm;
no filter, source domain, deployed six-radian correction limit, or theorem gate
is changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_signed_cayley_q8_v14 as V14

DEFAULT_DOMAIN = V14.DEFAULT_DOMAIN
SCHEMA = 1401
FULL = V14.FULL


def interval_euclidean_norm_upper(v) -> float:
    """Rigorous Euclidean upper for any nonempty finite interval vector."""
    if not v:
        raise ValueError("interval vector must be nonempty")
    s = 0.0
    for x in v:
        a = x.abs_upper()
        if not math.isfinite(a):
            raise ValueError("interval vector must be finite")
        s = FULL.up(s + FULL.up(a * a))
    return FULL.up(math.sqrt(max(0.0, s)))


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    original = V14._norm2_upper
    V14._norm2_upper = interval_euclidean_norm_upper
    try:
        core = V14.build(
            domain_path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces)
    finally:
        V14._norm2_upper = original

    failures = list(V14.validate(core))
    out = dict(core)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14B",
        "V14_signed_cayley_algebra_retained": True,
        "yz_interval_norm_uses_outward_euclidean_two_component_bound": True,
        "legacy_three_vector_norm_misuse_retired": True,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14B": (
            "PASS" if core.get("P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14") == "PASS" and not failures
            else "NOT_ESTABLISHED"
        ),
        "failures": failures,
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V14_signed_cayley_algebra_retained",
        "yz_interval_norm_uses_outward_euclidean_two_component_bound",
        "legacy_three_vector_norm_misuse_retired",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "q8_word_promoted_here", "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    st = d.get("P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14B")
    if st == "PASS":
        if d.get("signed_cayley_q8_composed_here") is not True:
            f.append("V14B PASS lacks signed q8 composition")
        if int(d.get("unclosed_q8_cells", -1)) != 0:
            f.append("V14B PASS retains unclosed q8 cells")
        if not float(d.get("max_post_sample1_cayley_norm_upper", math.inf)) < V14.Q_TARGET:
            f.append("V14B PASS does not satisfy strict q<8")
    elif st == "NOT_ESTABLISHED":
        if d.get("complete_sample1_branch_closed_here") is not False:
            f.append("V14B nonclosure claims sample1 closure")
        if d.get("first_unclosed_q8_cell") is None and not f:
            f.append("V14B nonclosure lacks witness")
    else:
        f.append("invalid V14B status")
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
    d = build(
        x.domain, source_pieces=x.source_pieces,
        source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
        tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces,
        residual_x_pieces=x.residual_x_pieces,
        parallel_pieces=x.parallel_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14B"],
        "cells": d["evaluated_signed_cayley_cells"],
        "radial_not_ready": d["radial_not_ready_cells"],
        "antipode_cells": d["product_scalar_antipode_cells"],
        "unclosed": d["unclosed_q8_cells"],
        "min_abs_W": d["minimum_abs_product_scalar_lower"],
        "max_q_plus": d["max_post_sample1_cayley_norm_upper"],
        "first_unclosed": d["first_unclosed_q8_cell"],
        "worst": d["worst_q8_cell"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
