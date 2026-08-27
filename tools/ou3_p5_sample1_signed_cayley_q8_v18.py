#!/usr/bin/env python3
"""V18: retain source-correlated current y/z support in the sample-1 q<8 proof.

V17 tightens the complete sample-1 current Cayley radius with the exact
sample-0 quaternion product scalar, but V14's signed-product route still keeps
only the signed current x component.  The other current components are replaced
by one scalar ||c_yz|| bound and the product

    v_d,yz^T c_yz

is enclosed only by Cauchy.  The V17 fine run therefore still leaves a very
small product-scalar denominator on high-angle cells even though the exact
sample-0 product already contains component information.

V18 keeps all three components of that sample-0 Cayley product.  After the
sample-0 correction it applies the same simultaneous R_x(d)^T proof-gauge
change used by the structured V7/V10 gain proof, propagates each component
through the same bounded first prediction and possible sample-1 S attitude
injection, and intersects every component with V17's tighter current-radius
ball.  The resulting signed y/z box and the existing Euclidean y/z ball are two
independent enclosures of the same current support.

For every V16 signed correction quaternion, V18 therefore forms both

    D_box = v_y c_y + v_z c_z
    D_ball in [-||v_yz|| ||c_yz||, +||v_yz|| ||c_yz||]

and intersects them before constructing

    W = 2 w_d - v_x c_x - D_yz.

The resulting W interval is also intersected with the original V14/V16 W
interval.  Thus V18 can only tighten the existing signed-product certificate;
V15's independent SO(3) geodesic route remains available and V16's correction
axis cone remains unchanged.

No estimator, source domain, source branch, shipping six-radian correction
limit, q<8 target, or theorem-promotion flag is changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_signed_cayley_q8_v14 as V14
import ou3_p5_sample1_signed_cayley_q8_v16 as V16
import ou3_p5_sample1_signed_cayley_q8_v17 as V17

DEFAULT_DOMAIN = V17.DEFAULT_DOMAIN
SCHEMA = 1800
FULL = V14.FULL
Q_TARGET = V14.Q_TARGET


def _intersect(a: Interval, b: Interval) -> Interval:
    lo = max(a.lo, b.lo)
    hi = min(a.hi, b.hi)
    if hi < lo:
        raise RuntimeError("independent interval enclosures do not intersect")
    return Interval(lo, hi)


def _rotate_yz_rx_transpose(cy: Interval, cz: Interval,
                            angle: Interval) -> tuple[Interval, Interval]:
    """Interval image of the y/z pair under R_x(angle)^T."""
    _sin, _cos = V14.CAYLEY2._trig_interval(angle.lo, angle.hi)
    # R_x(a)^T has yz block [[cos(a), sin(a)],[-sin(a), cos(a)]].
    return _cos * cy + _sin * cz, -(_sin * cy) + _cos * cz


def _yz_norm_upper(cy: Interval, cz: Interval) -> float:
    y = cy.abs_upper(); z = cz.abs_upper()
    return FULL.up(math.sqrt(FULL.up(FULL.up(y * y) + FULL.up(z * z))))


def _support_product_scalar(parent_W: Interval, wd: Interval, vd,
                            chart: dict) -> tuple[Interval, Interval, Interval]:
    """Intersect componentwise and Cauchy y/z support, then tighten W."""
    if len(vd) != 3:
        raise ValueError("three-component correction quaternion required")
    cx = chart["cx"]
    cy = chart["cy"]
    cz = chart["cz"]
    cyz = float(chart["cyz_norm_upper"])

    yz_box = vd[1] * cy + vd[2] * cz
    vdyz = _yz_norm_upper(vd[1], vd[2])
    yz_abs = FULL.up(vdyz * cyz)
    yz_ball = Interval.outward_bounds(-yz_abs, yz_abs)
    yz_joint = _intersect(yz_box, yz_ball)

    dot = vd[0] * cx + yz_joint
    support_W = FULL.I(2.0) * wd - dot
    joint_W = _intersect(parent_W, support_W)
    return joint_W, yz_box, yz_joint


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    context = {
        "chart": None,
        "quat": None,
        "charts": 0,
        "sample0_radius_improved": 0,
        "sample1_radius_improved": 0,
        "max_sample1_q_before": 0.0,
        "max_sample1_q_after": 0.0,
        "support_calls": 0,
        "support_refined": 0,
        "support_newly_closed": 0,
        "first_support_newly_closed": None,
        "first_support_refinement": None,
    }

    original_chart = V14._sample1_current_chart
    original_qplus = V14._qplus_from_product_scalar
    original_axis_cone = V16.axis_cone_normalized_shipping_quaternion

    def component_chart(*, first: dict, base: dict, vr: dict,
                        dom: dict, src: dict, sample1_s_angle: float):
        parent = dict(original_chart(
            first=first, base=base, vr=vr, dom=dom, src=src,
            sample1_s_angle=sample1_s_angle))
        context["charts"] += 1

        qpre = float(first["post_prediction_full_cayley_norm_upper"])
        ctan = float(first["post_prediction_cayley_tangent_norm_upper"])
        d0 = Interval.outward_bounds(*map(float, base["first_attitude_correction_rad"]))
        dd = float(vr["first_offaxis_attitude_correction_upper_rad"])
        e = Interval.outward_bounds(-dd, dd)
        dbox0 = [d0 + e, e, e]
        d0_hi = FULL.up(max(0.0, d0.hi) + dd)
        d0_lo = max(0.0, FULL.down(max(0.0, d0.lo) - dd))

        cpre = [
            Interval.outward_bounds(-ctan, ctan),
            Interval.outward_bounds(-ctan, ctan),
            Interval.outward_bounds(-qpre, qpre),
        ]
        w0, v0, branches0 = V14._normalized_shipping_quaternion(
            dbox0, radial_lower=d0_lo, radial_upper=d0_hi)
        dot0 = V14.CAYLEY1._dot(v0, cpre)
        W0 = FULL.I(2.0) * w0 - dot0
        if W0.lo <= 0.0 <= W0.hi:
            raise RuntimeError("sample-0 source-correlated product scalar crosses zero")
        cross0 = V14.CAYLEY1._cross(v0, cpre)
        V0 = [w0 * cpre[i] + FULL.I(2.0) * v0[i] + cross0[i]
              for i in range(3)]
        c0 = [FULL.I(2.0) * V0[i] / W0 for i in range(3)]

        # Apply the exact proof-coordinate change used by the structured gain
        # derivation.  It is a basis rotation, not an estimator operation.
        cy0, cz0 = _rotate_yz_rx_transpose(c0[1], c0[2], d0)
        comps = [c0[0], cy0, cz0]

        q0_triangle = V14.PREFIX2._post_correction_q_upper(qpre, d0_hi)
        q0_product = V17._q_upper_from_product_scalar(qpre, W0)
        q0_best = min(q0_triangle, q0_product)
        if not math.isfinite(q0_best):
            raise RuntimeError("sample-0 product-radius tightening left Cayley chart")
        context["sample0_radius_improved"] += int(q0_product < q0_triangle)

        transport = float(first["first_prediction_transport_angle_upper_rad"])
        comps = [V14._widen_cx_by_unknown_rotation(x, q0_best, transport)
                 for x in comps]
        qpred = V14.RG._q_after_first_prediction(q0_best, dom, float(src["dt_s"]))
        if not math.isfinite(qpred):
            raise RuntimeError("sample-1 prediction q upper left Cayley chart")

        ds = max(0.0, float(sample1_s_angle))
        comps = [V14._widen_cx_by_unknown_rotation(x, qpred, ds) for x in comps]
        q1_new = V14.PREFIX2._post_correction_q_upper(qpred, ds)
        q1_parent = float(parent["q1"])
        q1 = min(q1_parent, q1_new)
        if not math.isfinite(q1):
            raise RuntimeError("sample-1 current component support left Cayley chart")

        comps = [V17._clip_component_to_radius(x, q1) for x in comps]
        comps[0] = _intersect(comps[0], parent["cx"])

        cx_min = V14._minimum_abs(comps[0])
        yz2 = max(0.0, FULL.up(q1 * q1) - FULL.down(cx_min * cx_min))
        sphere_cyz = FULL.up(math.sqrt(yz2))
        box_cyz = _yz_norm_upper(comps[1], comps[2])
        cyz = min(float(parent["cyz_norm_upper"]), sphere_cyz, box_cyz)

        context["max_sample1_q_before"] = max(context["max_sample1_q_before"], q1_parent)
        context["max_sample1_q_after"] = max(context["max_sample1_q_after"], q1)
        context["sample1_radius_improved"] += int(q1 < q1_parent)

        parent.update({
            "cx": comps[0],
            "cy": comps[1],
            "cz": comps[2],
            "q1": q1,
            "cyz_norm_upper": cyz,
            "sample0_q_triangle_upper": q0_triangle,
            "sample0_q_product_scalar_upper": q0_product,
            "sample0_q_selected_upper": q0_best,
            "sample1_q_parent_upper": q1_parent,
            "sample1_q_product_tightened_upper": q1,
            "sample0_product_scalar": W0,
            "sample0_quaternion_branches": branches0,
            "source_correlated_current_yz_support_retained": True,
        })
        context["chart"] = parent
        return parent

    def tracked_axis_cone(dbox, *, radial_lower: float, radial_upper: float,
                          parent=None):
        ans = original_axis_cone(
            dbox, radial_lower=radial_lower, radial_upper=radial_upper,
            parent=parent)
        w, v, branches, narrowed = ans
        context["quat"] = (w, v, float(radial_lower), float(radial_upper))
        return w, v, branches, narrowed

    def support_qplus(q1: float, parent_W: Interval):
        context["support_calls"] += 1
        w_parent, q_parent = original_qplus(q1, parent_W)
        chart = context.get("chart")
        quat = context.get("quat")
        if chart is None or quat is None:
            return w_parent, q_parent
        wd, vd, rlo, rhi = quat
        joint_W, yz_box, yz_joint = _support_product_scalar(
            parent_W, wd, vd, chart)
        w_support, q_support = original_qplus(q1, joint_W)

        refined = joint_W.lo > parent_W.lo or joint_W.hi < parent_W.hi
        if refined:
            context["support_refined"] += 1
            if context["first_support_refinement"] is None:
                context["first_support_refinement"] = {
                    "q_current_upper": float(q1),
                    "correction_radial_lower_rad": rlo,
                    "correction_radial_upper_rad": rhi,
                    "parent_W": parent_W.as_list(),
                    "support_W": joint_W.as_list(),
                    "current_cy": chart["cy"].as_list(),
                    "current_cz": chart["cz"].as_list(),
                    "yz_component_box_dot": yz_box.as_list(),
                    "yz_joint_dot": yz_joint.as_list(),
                }

        parent_closed = math.isfinite(q_parent) and q_parent < Q_TARGET and w_parent > 0.0
        support_closed = math.isfinite(q_support) and q_support < Q_TARGET and w_support > 0.0
        if support_closed and not parent_closed:
            context["support_newly_closed"] += 1
            if context["first_support_newly_closed"] is None:
                context["first_support_newly_closed"] = {
                    "q_current_upper": float(q1),
                    "correction_radial_lower_rad": rlo,
                    "correction_radial_upper_rad": rhi,
                    "parent_q_upper": q_parent,
                    "support_q_upper": q_support,
                    "parent_abs_W_lower": w_parent,
                    "support_abs_W_lower": w_support,
                    "parent_W": parent_W.as_list(),
                    "support_W": joint_W.as_list(),
                }
        return max(w_parent, w_support), min(q_parent, q_support)

    V14._sample1_current_chart = component_chart
    V14._qplus_from_product_scalar = support_qplus
    V16.axis_cone_normalized_shipping_quaternion = tracked_axis_cone
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
        V16.axis_cone_normalized_shipping_quaternion = original_axis_cone
        V14._qplus_from_product_scalar = original_qplus
        V14._sample1_current_chart = original_chart

    inherited = V16.validate(core)
    parent_status = core.get("P5_SAMPLE1_AXIS_CONE_GEODESIC_SIGNED_CAYLEY_Q8_V16")
    out = dict(core)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_CURRENT_YZ_SUPPORT_SIGNED_CAYLEY_Q8_V18",
        "V16_axis_cone_geodesic_parent_retained": True,
        "V17_sample0_product_radius_construction_retained": True,
        "sample0_full_cayley_product_components_retained": True,
        "structured_Rx_proof_gauge_applied_to_current_yz": True,
        "prediction_and_sample1_S_component_widening_unchanged": True,
        "componentwise_yz_dot_intersected_with_cauchy_parent": True,
        "support_W_intersected_with_signed_product_parent": True,
        "sample1_current_chart_calls": int(context["charts"]),
        "sample0_product_radius_improved_charts": int(context["sample0_radius_improved"]),
        "sample1_product_radius_improved_charts": int(context["sample1_radius_improved"]),
        "max_sample1_q_before_product_tightening": float(context["max_sample1_q_before"]),
        "max_sample1_q_after_product_tightening": float(context["max_sample1_q_after"]),
        "current_yz_support_qplus_calls": int(context["support_calls"]),
        "current_yz_support_refined_cells": int(context["support_refined"]),
        "current_yz_support_newly_closed_cells": int(context["support_newly_closed"]),
        "first_current_yz_support_refinement": context["first_support_refinement"],
        "first_current_yz_support_newly_closed": context["first_support_newly_closed"],
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_CURRENT_YZ_SUPPORT_SIGNED_CAYLEY_Q8_V18": (
            "PASS" if parent_status == "PASS" and not inherited else "NOT_ESTABLISHED"
        ),
        "next_obligation": (
            "LIFT_CLOSED_SAMPLE1_PREFIX_TO_ALL_SOURCE_PHASE_CELLS_AND_CONTINUE_SAMPLE2_PREFIX"
            if parent_status == "PASS" and not inherited else
            "REFINE_REMAINING_Q8_CELLS_WITH_JOINT_CURRENT_CORRECTION_YZ_DIRECTION_SUBDIVISION"
        ),
        "failures": list(dict.fromkeys(inherited)),
    })
    return out


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_CURRENT_YZ_SUPPORT_SIGNED_CAYLEY_Q8_V18":
        failures.append("qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "V16_axis_cone_geodesic_parent_retained",
        "V17_sample0_product_radius_construction_retained",
        "sample0_full_cayley_product_components_retained",
        "structured_Rx_proof_gauge_applied_to_current_yz",
        "prediction_and_sample1_S_component_widening_unchanged",
        "componentwise_yz_dot_intersected_with_cauchy_parent",
        "support_W_intersected_with_signed_product_parent",
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
    calls = int(d.get("current_yz_support_qplus_calls", -1))
    refined = int(d.get("current_yz_support_refined_cells", -1))
    newly = int(d.get("current_yz_support_newly_closed_cells", -1))
    if not (calls >= refined >= newly >= 0):
        failures.append("invalid current-yz support accounting")
    charts = int(d.get("sample1_current_chart_calls", -1))
    s0 = int(d.get("sample0_product_radius_improved_charts", -1))
    s1 = int(d.get("sample1_product_radius_improved_charts", -1))
    if not (charts >= s0 >= 0 and charts >= s1 >= 0):
        failures.append("invalid product-radius accounting")
    if float(d.get("max_sample1_q_after_product_tightening", math.inf)) > \
            float(d.get("max_sample1_q_before_product_tightening", -math.inf)):
        failures.append("product-radius tightening increased sample1 q")
    st = d.get("P5_SAMPLE1_CURRENT_YZ_SUPPORT_SIGNED_CAYLEY_Q8_V18")
    if st == "PASS":
        if int(d.get("unclosed_q8_cells", -1)) != 0:
            failures.append("V18 PASS retains unclosed q8 cells")
        if not float(d.get("max_post_sample1_cayley_norm_upper", math.inf)) < Q_TARGET:
            failures.append("V18 PASS does not satisfy strict q<8")
    elif st == "NOT_ESTABLISHED":
        if d.get("complete_sample1_branch_closed_here") is True:
            failures.append("V18 nonclosure claims sample1 closure")
        if d.get("first_unclosed_q8_cell") is None and not failures:
            failures.append("V18 numerical nonclosure lacks q8 witness")
    else:
        failures.append("invalid V18 status")
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
        "status": out["P5_SAMPLE1_CURRENT_YZ_SUPPORT_SIGNED_CAYLEY_Q8_V18"],
        "cells": out.get("evaluated_signed_cayley_cells"),
        "unclosed": out.get("unclosed_q8_cells"),
        "sample1_chart_calls": out["sample1_current_chart_calls"],
        "sample0_radius_improved_charts": out["sample0_product_radius_improved_charts"],
        "sample1_radius_improved_charts": out["sample1_product_radius_improved_charts"],
        "max_q1_before": out["max_sample1_q_before_product_tightening"],
        "max_q1_after": out["max_sample1_q_after_product_tightening"],
        "yz_support_calls": out["current_yz_support_qplus_calls"],
        "yz_support_refined": out["current_yz_support_refined_cells"],
        "yz_support_newly_closed": out["current_yz_support_newly_closed_cells"],
        "first_yz_support_refinement": out["first_current_yz_support_refinement"],
        "first_yz_support_newly_closed": out["first_current_yz_support_newly_closed"],
        "first_unclosed": out.get("first_unclosed_q8_cell"),
        "worst": out.get("worst_q8_cell"),
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
