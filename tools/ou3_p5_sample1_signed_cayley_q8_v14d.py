#!/usr/bin/env python3
"""V14D: source-faithful radial sinc evaluation for sample-1 Cayley q<8.

V14C correctly intersected each signed correction component with the already
certified radial branch, but its axis-angle helper still evaluated
``sin(r/2)/r`` by dividing two independent intervals.  On a perfectly benign
radial cell such as [0.01,0.75] this can pair the largest sine with the smallest
radius and manufacture quaternion-vector components far larger than one.

The deployed source coefficient is

    k(r) = sin(r/2)/r = 0.5 sinc(r/2).

For r<=6 the half-angle lies in [0,3]<pi and the existing V1 proof backend has a
validated monotone sinc enclosure.  V14D uses that same validated sinc interval
on the *actual radial branch*.  For r>6, where monotonicity across the winding is
not available, it uses the V2 centered nonmonotone trig enclosure and divides
only by the strictly positive half-angle interval [r_lo/2,r_hi/2].

The <1e-2 source polynomial branch remains separately normalized exactly as in
V14.  The signed component boxes, V12D/V13E bounds, q<8 target, source family,
shipping six-radian proof limit, and all no-promotion guards are unchanged.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, hull
import ou3_p5_sample1_signed_cayley_q8_v14 as V14
import ou3_p5_sample1_signed_cayley_q8_v14b as V14B
import ou3_p5_sample1_signed_cayley_q8_v14c as V14C

DEFAULT_DOMAIN = V14.DEFAULT_DOMAIN
SCHEMA = 1404
FULL = V14.FULL
VT = V14.CAYLEY1.VT


def _axis_monotone_part(dbox, lo: float, hi: float):
    lo = max(V14.SERIES, float(lo)); hi = min(6.0, float(hi))
    if hi < lo:
        return None
    clipped = V14C._clip_component_box(dbox, hi)
    if clipped is None:
        return None
    half = Interval(FULL.down(0.5 * lo), FULL.up(0.5 * hi))
    sinc = VT.sinc_interval(half)
    k = FULL.I(0.5) * sinc
    clo = VT.cos_point(half.hi)
    chi = VT.cos_point(half.lo)
    w = Interval(clo.lo, chi.hi)
    return w, [k * x for x in clipped], "AXIS_ANGLE_MONOTONE_SINC_RADIAL"


def _axis_winding_part(dbox, lo: float, hi: float):
    lo = max(6.0, float(lo)); hi = float(hi)
    if hi < lo:
        return None
    if not (6.0 <= lo <= hi <= 9.0):
        raise RuntimeError("invalid winding radial branch")
    clipped = V14C._clip_component_box(dbox, hi)
    if clipped is None:
        return None
    half = Interval(FULL.down(0.5 * lo), FULL.up(0.5 * hi))
    sin_half, cos_half = V14.CAYLEY2._trig_interval(half.lo, half.hi)
    # Here half.lo >= 3, so interval division cannot manufacture the tiny-radius
    # dependency failure that V14C exhibited.  It is the same V2 winding proof.
    k = FULL.I(0.5) * sin_half / half
    return cos_half, [k * x for x in clipped], "AXIS_ANGLE_V2_WINDING_RADIAL"


def radial_sinc_normalized_shipping_quaternion(dbox, *, radial_lower: float,
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
            ser_box = V14C._clip_component_box(dbox, ser_hi)
            if ser_box is not None:
                p = V14._normalized_series_part(ser_box, ser_hi)
                if p is not None:
                    parts.append((p[0], p[1], "SERIES_NORMALIZED_RADIAL_CLIP"))

    if hi >= V14.SERIES and lo <= 6.0:
        p = _axis_monotone_part(dbox, max(lo, V14.SERIES), min(hi, 6.0))
        if p is not None:
            parts.append(p)
    if hi > 6.0:
        p = _axis_winding_part(dbox, max(lo, 6.0), hi)
        if p is not None:
            parts.append(p)

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
    V14._normalized_shipping_quaternion = radial_sinc_normalized_shipping_quaternion
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
        "qualification": "OU3_P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14D",
        "V14_signed_cayley_algebra_retained": True,
        "V14B_outward_yz_norm_retained": True,
        "V14C_branch_local_radial_component_intersection_retained": True,
        "validated_monotone_sinc_used_through_6_rad": True,
        "V2_nonmonotone_winding_trig_used_above_6_rad": True,
        "independent_small_radius_sine_division_retired": True,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14D": (
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
        "V14C_branch_local_radial_component_intersection_retained",
        "validated_monotone_sinc_used_through_6_rad",
        "V2_nonmonotone_winding_trig_used_above_6_rad",
        "independent_small_radius_sine_division_retired",
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
    st = d.get("P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14D")
    if st == "PASS":
        if d.get("signed_cayley_q8_composed_here") is not True:
            f.append("V14D PASS lacks signed q8 composition")
        if int(d.get("unclosed_q8_cells", -1)) != 0:
            f.append("V14D PASS retains unclosed q8 cells")
        if not float(d.get("max_post_sample1_cayley_norm_upper", math.inf)) < V14.Q_TARGET:
            f.append("V14D PASS does not satisfy strict q<8")
    elif st == "NOT_ESTABLISHED":
        if d.get("complete_sample1_branch_closed_here") is not False:
            f.append("V14D nonclosure claims sample1 closure")
        if d.get("first_unclosed_q8_cell") is None and not f:
            f.append("V14D nonclosure lacks witness")
    else:
        f.append("invalid V14D status")
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
        residual_x_pieces=x.residual_x_pieces, parallel_pieces=x.parallel_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14D"],
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
