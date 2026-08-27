#!/usr/bin/env python3
"""V16: retain the signed correction-axis cone in the sample-1 q<8 proof.

V15 adds an independent SO(3) geodesic closure route, but the remaining
high-principal-angle cells are still evaluated by V14D's signed quaternion
product.  In that route the deployed axis-angle quaternion vector is formed as

    v_d = sin(r/2) d / r.

V14D evaluates ``sin(r/2)/r`` over the certified radial interval and then
multiplies it by each signed Cartesian component interval.  That is rigorous,
but it forgets a strong correlation when one component of ``d`` dominates the
other two.  For example, the completed #416 worst witness has a positive
x-component above 3.22 rad while the y/z components remain below about 0.39 rad
jointly.  Every real correction axis in that box is therefore tightly aligned
with +x even though independent ``k(r)*d_x`` interval multiplication permits a
much smaller quaternion x component.

For a Cartesian correction box d=(d_x,d_y,d_z), the unit axis u=d/||d|| obeys,
for each sign-definite component i,

    |u_i| >= a_i,min / sqrt(a_i,min^2 + sum_{j!=i} a_j,max^2),

and for every component

    |u_i| <= a_i,max / sqrt(a_i,max^2 + sum_{j!=i} a_j,min^2).

The already certified radial branch also gives

    a_i,min/r_max <= |u_i| <= a_i,max/r_min.

V16 intersects these source-complete cone bounds with V14D's existing
quaternion component enclosure.  The radial interval is split at 2*pi before
multiplying by sin(r/2), so the quaternion-vector sign flip across a complete
2*pi winding is not silently treated as monotone.  The polynomial source branch
below 1e-2 rad is left exactly unchanged.

V16 is only a tightening of the same deployed quaternion image.  It does not
change the estimator, source domain, six-radian shipping correction limit,
q<8 target, source family, or any theorem-promotion flag.  V15's geodesic bound
and V14D's signed-product bound remain independent valid parents.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, hull
import ou3_p5_sample1_signed_cayley_q8_v14 as V14
import ou3_p5_sample1_signed_cayley_q8_v14d as V14D
import ou3_p5_sample1_signed_cayley_q8_v15 as V15

DEFAULT_DOMAIN = V15.DEFAULT_DOMAIN
SCHEMA = 1600
FULL = V14.FULL
SERIES = V14.SERIES
TWO_PI = 2.0 * math.pi
Q_TARGET = V14.Q_TARGET


def _minimum_abs(x: Interval) -> float:
    if x.lo <= 0.0 <= x.hi:
        return 0.0
    return min(abs(x.lo), abs(x.hi))


def _axis_component_interval(dbox, index: int, radial_lower: float,
                             radial_upper: float) -> Interval:
    """Outward enclosure of d_i/||d|| using the signed Cartesian cone."""
    if len(dbox) != 3 or index not in (0, 1, 2):
        raise ValueError("three-component correction box required")
    rlo = float(radial_lower)
    rhi = float(radial_upper)
    if not (math.isfinite(rlo) and math.isfinite(rhi) and 0.0 < rlo <= rhi):
        raise ValueError("positive radial branch required for axis cone")

    amin = [_minimum_abs(x) for x in dbox]
    amax = [x.abs_upper() for x in dbox]
    ai0 = amin[index]
    ai1 = amax[index]
    others = [j for j in range(3) if j != index]

    if ai0 > 0.0:
        den_hi2 = FULL.up(FULL.up(ai0 * ai0) + FULL.up(
            FULL.up(amax[others[0]] * amax[others[0]]) +
            FULL.up(amax[others[1]] * amax[others[1]])))
        cone_lo = FULL.down(ai0 / math.sqrt(max(ai0 * ai0, den_hi2)))
    else:
        cone_lo = 0.0

    if ai1 > 0.0:
        den_lo2 = FULL.down(FULL.down(ai1 * ai1) + FULL.down(
            FULL.down(amin[others[0]] * amin[others[0]]) +
            FULL.down(amin[others[1]] * amin[others[1]])))
        den_lo = math.sqrt(max(0.0, den_lo2))
        cone_hi = 1.0 if den_lo == 0.0 else FULL.up(ai1 / den_lo)
    else:
        cone_hi = 0.0

    radial_lo = FULL.down(ai0 / rhi) if ai0 > 0.0 else 0.0
    radial_hi = FULL.up(ai1 / rlo) if ai1 > 0.0 else 0.0
    lo_abs = max(0.0, cone_lo, radial_lo)
    hi_abs = min(1.0, cone_hi, radial_hi)
    if hi_abs < lo_abs:
        # This can only arise from mutually inconsistent outward parent boxes.
        raise RuntimeError("axis-cone/radial intersection is empty")

    x = dbox[index]
    if x.lo >= 0.0:
        return Interval(FULL.down(lo_abs), FULL.up(hi_abs))
    if x.hi <= 0.0:
        return Interval(FULL.down(-hi_abs), FULL.up(-lo_abs))
    return Interval(FULL.down(-hi_abs), FULL.up(hi_abs))


def _sin_half_interval(radial_lower: float, radial_upper: float) -> Interval:
    half_lo = FULL.down(0.5 * float(radial_lower))
    half_hi = FULL.up(0.5 * float(radial_upper))
    sin_half, _cos_half = V14.CAYLEY2._trig_interval(half_lo, half_hi)
    return sin_half


def _axis_vector_cone(dbox, radial_lower: float, radial_upper: float):
    """Quaternion-vector cone for an axis-angle radial branch >= SERIES."""
    lo = float(radial_lower)
    hi = float(radial_upper)
    if not (SERIES <= lo <= hi <= 9.0):
        raise ValueError("axis cone requires an axis-angle radial branch")

    # Preserve the sign of sin(r/2) on either side of one complete 2*pi turn.
    cuts = [lo]
    if lo < TWO_PI < hi:
        cuts.append(TWO_PI)
    cuts.append(hi)

    parts = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        if b < a:
            continue
        clipped = V14D.V14C._clip_component_box(dbox, b)
        if clipped is None:
            continue
        s = _sin_half_interval(a, b)
        u = [_axis_component_interval(clipped, i, a, b) for i in range(3)]
        parts.append([s * x for x in u])
    if not parts:
        raise RuntimeError("axis cone has no nonempty radial branch")
    return [hull(*(p[i] for p in parts)) for i in range(3)]


def _intersect(a: Interval, b: Interval) -> Interval:
    lo = max(a.lo, b.lo)
    hi = min(a.hi, b.hi)
    if hi < lo:
        raise RuntimeError("independent quaternion enclosures do not intersect")
    return Interval(lo, hi)


def axis_cone_normalized_shipping_quaternion(dbox, *, radial_lower: float,
                                              radial_upper: float,
                                              parent=None):
    """Intersect V14D's quaternion vector with the source Cartesian axis cone."""
    base = V14D.radial_sinc_normalized_shipping_quaternion if parent is None else parent
    w, v, branches = base(
        dbox, radial_lower=radial_lower, radial_upper=radial_upper)
    lo = max(0.0, float(radial_lower))
    hi = float(radial_upper)
    if lo < SERIES or hi < SERIES:
        return w, v, branches, False

    cone = _axis_vector_cone(dbox, max(lo, SERIES), hi)
    refined = [_intersect(v[i], cone[i]) for i in range(3)]
    narrowed = any(
        refined[i].lo > v[i].lo or refined[i].hi < v[i].hi
        for i in range(3)
    )
    return w, refined, branches, narrowed


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    context = {"calls": 0, "refined": 0}
    original = V14D.radial_sinc_normalized_shipping_quaternion

    def refined_quat(dbox, *, radial_lower: float, radial_upper: float):
        context["calls"] += 1
        w, v, branches, narrowed = axis_cone_normalized_shipping_quaternion(
            dbox, radial_lower=radial_lower, radial_upper=radial_upper,
            parent=original)
        context["refined"] += int(narrowed)
        return w, v, branches

    V14D.radial_sinc_normalized_shipping_quaternion = refined_quat
    try:
        core = V15.build(
            Path(domain_path).resolve(),
            source_pieces=source_pieces,
            source_cell_index=source_cell_index,
            p_pieces=p_pieces,
            tangent_pieces=tangent_pieces,
            axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces,
        )
    finally:
        V14D.radial_sinc_normalized_shipping_quaternion = original

    inherited = V15.validate(core)
    parent_status = core.get("P5_SAMPLE1_GEODESIC_SIGNED_CAYLEY_Q8_V15")
    out = dict(core)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_AXIS_CONE_GEODESIC_SIGNED_CAYLEY_Q8_V16",
        "V15_geodesic_signed_product_parent_retained": True,
        "signed_cartesian_correction_axis_cone_used": True,
        "axis_cone_radial_interval_intersection_used": True,
        "axis_cone_split_at_two_pi_winding": True,
        "series_branch_axis_cone_refined": False,
        "quaternion_axis_cone_calls": int(context["calls"]),
        "quaternion_axis_cone_refined_cells": int(context["refined"]),
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_AXIS_CONE_GEODESIC_SIGNED_CAYLEY_Q8_V16": (
            "PASS" if parent_status == "PASS" and not inherited else "NOT_ESTABLISHED"
        ),
        "next_obligation": (
            "LIFT_CLOSED_SAMPLE1_PREFIX_TO_ALL_SOURCE_PHASE_CELLS_AND_CONTINUE_SAMPLE2_PREFIX"
            if parent_status == "PASS" and not inherited else
            "REFINE_REMAINING_Q8_CELLS_WITH_SOURCE_CORRELATED_CURRENT_YZ_SUPPORT"
        ),
        "failures": list(dict.fromkeys(inherited)),
    })
    return out


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_AXIS_CONE_GEODESIC_SIGNED_CAYLEY_Q8_V16":
        failures.append("qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "V15_geodesic_signed_product_parent_retained",
        "signed_cartesian_correction_axis_cone_used",
        "axis_cone_radial_interval_intersection_used",
        "axis_cone_split_at_two_pi_winding",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "series_branch_axis_cone_refined",
        "deployed_correction_limit_increased", "q8_word_promoted_here",
        "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    calls = int(d.get("quaternion_axis_cone_calls", -1))
    refined = int(d.get("quaternion_axis_cone_refined_cells", -1))
    if not (calls >= refined >= 0):
        failures.append("invalid axis-cone accounting")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        failures.append("deployed correction limit changed")
    st = d.get("P5_SAMPLE1_AXIS_CONE_GEODESIC_SIGNED_CAYLEY_Q8_V16")
    if st == "PASS":
        if int(d.get("unclosed_q8_cells", -1)) != 0:
            failures.append("V16 PASS retains unclosed q8 cells")
        if not float(d.get("max_post_sample1_cayley_norm_upper", math.inf)) < Q_TARGET:
            failures.append("V16 PASS does not satisfy strict q<8")
    elif st == "NOT_ESTABLISHED":
        if d.get("complete_sample1_branch_closed_here") is True:
            failures.append("V16 nonclosure claims sample1 closure")
        if d.get("first_unclosed_q8_cell") is None and not failures:
            failures.append("V16 numerical nonclosure lacks q8 witness")
    else:
        failures.append("invalid V16 status")
    return list(dict.fromkeys(failures))


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
    args = ap.parse_args()
    out = build(
        args.domain,
        source_pieces=args.source_pieces,
        source_cell_index=args.source_cell_index,
        p_pieces=args.p_pieces,
        tangent_pieces=args.tangent_pieces,
        axial_pieces=args.axial_pieces,
        residual_x_pieces=args.residual_x_pieces,
        parallel_pieces=args.parallel_pieces,
    )
    vf = validate(out)
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_SAMPLE1_AXIS_CONE_GEODESIC_SIGNED_CAYLEY_Q8_V16"],
        "cells": out.get("evaluated_signed_cayley_cells"),
        "unclosed": out.get("unclosed_q8_cells"),
        "axis_cone_calls": out["quaternion_axis_cone_calls"],
        "axis_cone_refined": out["quaternion_axis_cone_refined_cells"],
        "geodesic_newly_closed": out.get("geodesic_bound_newly_closed_cells"),
        "first_unclosed": out.get("first_unclosed_q8_cell"),
        "worst": out.get("worst_q8_cell"),
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
