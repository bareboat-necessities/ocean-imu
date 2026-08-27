#!/usr/bin/env python3
"""V17: tighten the sample-1 current Cayley radius from the exact sample-0 product scalar.

V14 already evaluates the exact source-correlated sample-0 quaternion product
scalar W0 before propagating the sample-1 prefix.  It then discards that scalar
for the radius calculation and falls back to the SO(3) triangle bound
``_post_correction_q_upper(qpre,d0_hi)``.  This is safe but loses the same
current/correction correlation that W0 was introduced to preserve.

For the normalized shipping correction quaternion and the unnormalized current
Cayley quaternion [2,c], quaternion multiplication preserves norm, hence

    W0^2 + ||V0||^2 = 4 + ||c||^2.

Therefore any strict lower bound |W0| >= wmin > 0 gives the independent exact
radius enclosure

    q0 <= 2 sqrt((4 + qpre^2)/wmin^2 - 1).

V17 takes the minimum of that source-correlated bound and the existing geodesic
triangle bound, propagates it through the *same* first prediction and possible
sample-1 S injection, and intersects V14's already-valid signed c_x interval
with the tighter norm ball.  V16's correction-axis cone, V15's geodesic closure,
and V14D's signed-product route remain unchanged.

No estimator, source domain, shipping correction limit, q<8 target, source
branch, or theorem-promotion flag is changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_signed_cayley_q8_v14 as V14
import ou3_p5_sample1_signed_cayley_q8_v16 as V16

DEFAULT_DOMAIN = V16.DEFAULT_DOMAIN
SCHEMA = 1700
FULL = V14.FULL
Q_TARGET = V14.Q_TARGET


def _q_upper_from_product_scalar(q_before: float, W: Interval) -> float:
    """Exact Cayley-radius upper from a separated quaternion product scalar."""
    q = float(q_before)
    if not (math.isfinite(q) and q >= 0.0):
        raise ValueError("invalid current Cayley radius")
    if W.lo <= 0.0 <= W.hi:
        return math.inf
    wmin = min(abs(W.lo), abs(W.hi))
    if not wmin > 0.0:
        return math.inf
    num = FULL.up(4.0 + FULL.up(q * q))
    den = FULL.down(wmin * wmin)
    if not den > 0.0:
        return math.inf
    ratio = FULL.up(num / den)
    r = max(0.0, FULL.up(ratio - 1.0))
    return FULL.up(2.0 * math.sqrt(r))


def _clip_component_to_radius(x: Interval, q_upper: float) -> Interval:
    q = float(q_upper)
    if not (math.isfinite(q) and q >= 0.0):
        raise ValueError("invalid Cayley radius for component clipping")
    lo = max(x.lo, FULL.down(-q))
    hi = min(x.hi, FULL.up(q))
    if hi < lo:
        raise RuntimeError("independent Cayley component/radius enclosures are disjoint")
    return Interval(lo, hi)


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    context = {
        "charts": 0,
        "sample0_radius_improved": 0,
        "sample1_radius_improved": 0,
        "max_sample1_q_before": 0.0,
        "max_sample1_q_after": 0.0,
        "first_improved_chart": None,
    }
    original_chart = V14._sample1_current_chart

    def tightened_chart(*, first: dict, base: dict, vr: dict,
                        dom: dict, src: dict, sample1_s_angle: float):
        chart = dict(original_chart(
            first=first, base=base, vr=vr, dom=dom, src=src,
            sample1_s_angle=sample1_s_angle))
        context["charts"] += 1

        qpre = float(first["post_prediction_full_cayley_norm_upper"])
        d0_hi = float(chart["sample0_correction_radial_upper"])
        W0 = chart["sample0_product_scalar"]
        q0_triangle = V14.PREFIX2._post_correction_q_upper(qpre, d0_hi)
        q0_product = _q_upper_from_product_scalar(qpre, W0)
        q0_best = min(q0_triangle, q0_product)
        if not math.isfinite(q0_best):
            raise RuntimeError("sample-0 product-scalar tightening left Cayley chart")
        if q0_product < q0_triangle:
            context["sample0_radius_improved"] += 1

        qpred_best = V14.RG._q_after_first_prediction(
            q0_best, dom, float(src["dt_s"]))
        ds = max(0.0, float(sample1_s_angle))
        q1_best = V14.PREFIX2._post_correction_q_upper(qpred_best, ds)
        q1_parent = float(chart["q1"])
        q1_best = min(q1_parent, q1_best)
        if not math.isfinite(q1_best):
            raise RuntimeError("sample-1 product-scalar tightening left Cayley chart")

        cx_parent = chart["cx"]
        cx = _clip_component_to_radius(cx_parent, q1_best)
        cx_min = V14._minimum_abs(cx)
        yz2 = max(0.0, FULL.up(q1_best * q1_best) - FULL.down(cx_min * cx_min))
        cyz = FULL.up(math.sqrt(yz2))

        context["max_sample1_q_before"] = max(context["max_sample1_q_before"], q1_parent)
        context["max_sample1_q_after"] = max(context["max_sample1_q_after"], q1_best)
        if q1_best < q1_parent:
            context["sample1_radius_improved"] += 1
            if context["first_improved_chart"] is None:
                context["first_improved_chart"] = {
                    "qpre": qpre,
                    "sample0_product_scalar": W0.as_list(),
                    "q0_triangle_upper": q0_triangle,
                    "q0_product_upper": q0_product,
                    "q0_selected_upper": q0_best,
                    "q1_parent_upper": q1_parent,
                    "q1_selected_upper": q1_best,
                    "cx_parent": cx_parent.as_list(),
                    "cx_selected": cx.as_list(),
                }

        chart.update({
            "cx": cx,
            "q1": q1_best,
            "cyz_norm_upper": cyz,
            "sample0_q_triangle_upper": q0_triangle,
            "sample0_q_product_scalar_upper": q0_product,
            "sample0_q_selected_upper": q0_best,
            "sample1_q_parent_upper": q1_parent,
            "sample1_q_product_tightened_upper": q1_best,
            "sample0_product_scalar_radius_tightening_used": True,
        })
        return chart

    V14._sample1_current_chart = tightened_chart
    try:
        core = V16.build(
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
        V14._sample1_current_chart = original_chart

    inherited = V16.validate(core)
    parent_status = core.get("P5_SAMPLE1_AXIS_CONE_GEODESIC_SIGNED_CAYLEY_Q8_V16")
    out = dict(core)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_SAMPLE0_PRODUCT_RADIUS_SIGNED_CAYLEY_Q8_V17",
        "V16_axis_cone_geodesic_parent_retained": True,
        "sample0_exact_product_scalar_radius_bound_used": True,
        "sample0_triangle_radius_parent_retained_by_minimum": True,
        "sample1_prediction_and_S_maps_unchanged": True,
        "sample1_signed_cx_intersected_with_tighter_radius": True,
        "sample1_current_chart_calls": int(context["charts"]),
        "sample0_product_radius_improved_charts": int(context["sample0_radius_improved"]),
        "sample1_product_radius_improved_charts": int(context["sample1_radius_improved"]),
        "max_sample1_q_before_product_tightening": float(context["max_sample1_q_before"]),
        "max_sample1_q_after_product_tightening": float(context["max_sample1_q_after"]),
        "first_product_radius_improved_chart": context["first_improved_chart"],
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_SAMPLE0_PRODUCT_RADIUS_SIGNED_CAYLEY_Q8_V17": (
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
    if d.get("qualification") != "OU3_P5_SAMPLE1_SAMPLE0_PRODUCT_RADIUS_SIGNED_CAYLEY_Q8_V17":
        failures.append("qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "V16_axis_cone_geodesic_parent_retained",
        "sample0_exact_product_scalar_radius_bound_used",
        "sample0_triangle_radius_parent_retained_by_minimum",
        "sample1_prediction_and_S_maps_unchanged",
        "sample1_signed_cx_intersected_with_tighter_radius",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed",
        "deployed_correction_limit_increased", "q8_word_promoted_here",
        "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        failures.append("deployed correction limit changed")
    charts = int(d.get("sample1_current_chart_calls", -1))
    s0 = int(d.get("sample0_product_radius_improved_charts", -1))
    s1 = int(d.get("sample1_product_radius_improved_charts", -1))
    if not (charts >= s0 >= 0 and charts >= s1 >= 0):
        failures.append("invalid product-radius accounting")
    if float(d.get("max_sample1_q_after_product_tightening", math.inf)) > \
            float(d.get("max_sample1_q_before_product_tightening", -math.inf)):
        failures.append("product-radius tightening increased sample1 q")
    st = d.get("P5_SAMPLE1_SAMPLE0_PRODUCT_RADIUS_SIGNED_CAYLEY_Q8_V17")
    if st == "PASS":
        if int(d.get("unclosed_q8_cells", -1)) != 0:
            failures.append("V17 PASS retains unclosed q8 cells")
        if not float(d.get("max_post_sample1_cayley_norm_upper", math.inf)) < Q_TARGET:
            failures.append("V17 PASS does not satisfy strict q<8")
    elif st == "NOT_ESTABLISHED":
        if d.get("complete_sample1_branch_closed_here") is True:
            failures.append("V17 nonclosure claims sample1 closure")
        if d.get("first_unclosed_q8_cell") is None and not failures:
            failures.append("V17 numerical nonclosure lacks q8 witness")
    else:
        failures.append("invalid V17 status")
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
        "status": out["P5_SAMPLE1_SAMPLE0_PRODUCT_RADIUS_SIGNED_CAYLEY_Q8_V17"],
        "cells": out.get("evaluated_signed_cayley_cells"),
        "unclosed": out.get("unclosed_q8_cells"),
        "sample1_chart_calls": out["sample1_current_chart_calls"],
        "sample0_radius_improved_charts": out["sample0_product_radius_improved_charts"],
        "sample1_radius_improved_charts": out["sample1_product_radius_improved_charts"],
        "max_q1_before": out["max_sample1_q_before_product_tightening"],
        "max_q1_after": out["max_sample1_q_after_product_tightening"],
        "first_radius_improvement": out["first_product_radius_improved_chart"],
        "first_unclosed": out.get("first_unclosed_q8_cell"),
        "worst": out.get("worst_q8_cell"),
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
