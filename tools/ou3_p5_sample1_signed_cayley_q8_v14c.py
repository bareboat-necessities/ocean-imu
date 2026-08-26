#!/usr/bin/env python3
"""V14C branch-local radial clipping for OU-III P5 sample-1 q<8.

V14B fixed the accidental use of a three-vector norm routine on the yz block.
Its full run then showed an impossible product-scalar interval: for q1<1 and a
unit correction quaternion, W=2w-v^T c cannot span tens in magnitude.  The
source was another Cartesian dependency artifact in V14's shipping-quaternion
enclosure.  A signed component box was reused unchanged for both the tiny
series branch and the much larger axis-angle branch, even though each branch is
also conditioned on a radial interval.

For every real vector d with ||d||<=r one has |d_i|<=r.  V14C therefore
intersects each signed component interval with [-r,r] separately for each
radial branch before evaluating exactly the same shipping series/axis-angle
formula.  This is a set intersection implied by the already certified radial
cell; it does not narrow the source family or alter any filter/proof limit.

V14's algebra, V14B's outward Euclidean norm fix, V12D/V13E prerequisites,
strict q<8 target, and no-promotion guards are otherwise unchanged.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, hull
import ou3_p5_sample1_signed_cayley_q8_v14 as V14
import ou3_p5_sample1_signed_cayley_q8_v14b as V14B

DEFAULT_DOMAIN = V14.DEFAULT_DOMAIN
SCHEMA = 1402
FULL = V14.FULL


def _clip_component_box(dbox, radial_upper: float):
    r = float(radial_upper)
    if not (math.isfinite(r) and r >= 0.0):
        raise ValueError("finite nonnegative radial upper required")
    out = []
    for x in dbox:
        lo = max(float(x.lo), -r)
        hi = min(float(x.hi), r)
        if lo > hi:
            return None
        out.append(Interval.outward_bounds(lo, hi))
    return out


def branch_local_normalized_shipping_quaternion(dbox, *, radial_lower: float,
                                                 radial_upper: float):
    hi = float(radial_upper)
    lo = max(0.0, float(radial_lower))
    if not (math.isfinite(lo) and math.isfinite(hi) and 0.0 <= lo <= hi <= 9.0):
        raise RuntimeError("invalid source correction radial cell")
    box_norm = V14B.interval_euclidean_norm_upper(dbox)
    if hi > box_norm + 64.0 * math.ulp(max(1.0, box_norm)):
        raise RuntimeError("radial upper exceeds signed component box")

    parts = []
    if lo < V14.SERIES:
        ser_hi = min(hi, V14.SERIES)
        if ser_hi > 0.0:
            ser_box = _clip_component_box(dbox, ser_hi)
            if ser_box is not None:
                p = V14._normalized_series_part(ser_box, ser_hi)
                if p is not None:
                    parts.append((p[0], p[1], "SERIES_NORMALIZED_RADIAL_CLIP"))
    if hi >= V14.SERIES:
        axis_lo = max(lo, V14.SERIES)
        axis_box = _clip_component_box(dbox, hi)
        if axis_box is not None:
            p = V14._axis_part(axis_box, axis_lo, hi)
            if p is not None:
                parts.append((p[0], p[1], "AXIS_ANGLE_UNIT_RADIAL_CLIP"))
    if not parts:
        if hi == 0.0:
            return FULL.I(1.0), [FULL.I(0.0) for _ in range(3)], ["ZERO"]
        raise RuntimeError("radial/component intersection is empty")

    w = hull(*(p[0] for p in parts))
    v = [hull(*(p[1][i] for p in parts)) for i in range(3)]
    return w, v, [p[2] for p in parts]


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    original_norm = V14._norm2_upper
    original_quat = V14._normalized_shipping_quaternion
    V14._norm2_upper = V14B.interval_euclidean_norm_upper
    V14._normalized_shipping_quaternion = branch_local_normalized_shipping_quaternion
    try:
        core = V14.build(
            domain_path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces)
    finally:
        V14._normalized_shipping_quaternion = original_quat
        V14._norm2_upper = original_norm

    failures = list(V14.validate(core))
    out = dict(core)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14C",
        "V14_signed_cayley_algebra_retained": True,
        "V14B_outward_yz_norm_retained": True,
        "branch_local_radial_component_intersection_used": True,
        "component_radial_intersection_is_source_implied": True,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14C": (
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
        "source_generated_not_trajectory_fit", "V14_signed_cayley_algebra_retained",
        "V14B_outward_yz_norm_retained",
        "branch_local_radial_component_intersection_used",
        "component_radial_intersection_is_source_implied",
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
    st = d.get("P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14C")
    if st == "PASS":
        if d.get("signed_cayley_q8_composed_here") is not True:
            f.append("V14C PASS lacks signed q8 composition")
        if int(d.get("unclosed_q8_cells", -1)) != 0:
            f.append("V14C PASS retains unclosed q8 cells")
        if not float(d.get("max_post_sample1_cayley_norm_upper", math.inf)) < V14.Q_TARGET:
            f.append("V14C PASS does not satisfy strict q<8")
    elif st == "NOT_ESTABLISHED":
        if d.get("complete_sample1_branch_closed_here") is not False:
            f.append("V14C nonclosure claims sample1 closure")
        if d.get("first_unclosed_q8_cell") is None and not f:
            f.append("V14C nonclosure lacks witness")
    else:
        f.append("invalid V14C status")
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
        "status": d["P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14C"],
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
