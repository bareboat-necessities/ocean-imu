#!/usr/bin/env python3
"""V19: adaptive joint current/correction y-z subdivision for the sample-1 q<8 proof.

V18B retains componentwise current y/z support and intersects its signed dot
product with the Cauchy parent, but its first remaining q<8 witness is unchanged.
The remaining loss is inside the Cartesian product of the current and correction
y/z interval boxes: extrema from different directions can still be paired before
the Euclidean constraints are applied.

V19 partitions both y/z boxes into two outward subintervals per component. On
each joint subbox it intersects the componentwise dot-product enclosure with the
local current-radius and unit-quaternion-vector norm bounds, discards only
subboxes whose simultaneous bounds prove them empty, then hulls the surviving
dot intervals. The hull is intersected with the V18B parent before W is formed.
Thus the construction can only tighten an already-valid product-scalar parent.

To keep the diagnostic run bounded, subdivision is applied only to V18B cells
that are not already q<8 and whose parent q upper is at most 64. Higher-q cells
retain the V18B parent unchanged and are counted explicitly. No estimator,
source domain, source branch, six-radian shipping correction limit, q<8 target,
or theorem-promotion flag is changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_signed_cayley_q8_v14 as V14
import ou3_p5_sample1_signed_cayley_q8_v18 as V18
import ou3_p5_sample1_signed_cayley_q8_v18b as V18B

DEFAULT_DOMAIN = V18B.DEFAULT_DOMAIN
SCHEMA = 1900
FULL = V14.FULL
Q_TARGET = V14.Q_TARGET
SUBDIVISION_Q_CEILING = 64.0
BASE_QPLUS = V14._qplus_from_product_scalar


def _parts2(x: Interval) -> list[Interval]:
    """Two outward pieces whose union covers x exactly up to overlap."""
    if x.lo == x.hi:
        return [x]
    mid = 0.5 * (x.lo + x.hi)
    if not math.isfinite(mid):
        raise ValueError("finite interval required for directional subdivision")
    return [Interval.outward_bounds(x.lo, mid), Interval.outward_bounds(mid, x.hi)]


def _hull(xs: list[Interval]) -> Interval:
    if not xs:
        raise RuntimeError("joint y-z subdivision eliminated every subbox")
    return Interval(min(x.lo for x in xs), max(x.hi for x in xs))


def _joint_yz_dot_subdivision(vd, chart: dict, parent: Interval):
    """Bound v_y c_y + v_z c_z by joint 2x2 component subdivision."""
    cy = chart["cy"]
    cz = chart["cz"]
    current_norm_parent = float(chart["cyz_norm_upper"])
    if not (math.isfinite(current_norm_parent) and current_norm_parent >= 0.0):
        raise RuntimeError("invalid current y-z norm parent")

    ys_c = _parts2(cy)
    zs_c = _parts2(cz)
    ys_v = _parts2(vd[1])
    zs_v = _parts2(vd[2])
    survivors = []
    total_pairs = empty_pairs = 0

    for cyp in ys_c:
        for czp in zs_c:
            cn = min(current_norm_parent, V18._yz_norm_upper(cyp, czp))
            for vyp in ys_v:
                for vzp in zs_v:
                    total_pairs += 1
                    vn = min(1.0, V18._yz_norm_upper(vyp, vzp))
                    box = vyp * cyp + vzp * czp
                    rad = FULL.up(vn * cn)
                    ball = Interval.outward_bounds(-rad, rad)
                    lo = max(box.lo, ball.lo)
                    hi = min(box.hi, ball.hi)
                    if hi < lo:
                        empty_pairs += 1
                        continue
                    survivors.append(Interval(lo, hi))

    subdiv = _hull(survivors)
    joint = V18._intersect(parent, subdiv)
    return joint, total_pairs, empty_pairs


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6,
          subdivision_q_ceiling: float = SUBDIVISION_Q_CEILING) -> dict:
    if not (math.isfinite(subdivision_q_ceiling) and subdivision_q_ceiling >= Q_TARGET):
        raise ValueError("subdivision q ceiling must be finite and >= q target")

    context = {
        "calls": 0,
        "attempted": 0,
        "skipped_closed": 0,
        "skipped_high_q": 0,
        "pair_evaluations": 0,
        "empty_pairs": 0,
        "refined": 0,
        "newly_closed": 0,
        "first_refined": None,
        "first_newly_closed": None,
    }
    original_support = V18._support_product_scalar

    def subdivided_support(parent_W: Interval, wd: Interval, vd, chart: dict):
        context["calls"] += 1
        joint_W, yz_box, yz_parent = original_support(parent_W, wd, vd, chart)
        w_parent, q_parent = BASE_QPLUS(float(chart["q1"]), joint_W)
        parent_closed = math.isfinite(q_parent) and q_parent < Q_TARGET and w_parent > 0.0
        if parent_closed:
            context["skipped_closed"] += 1
            return joint_W, yz_box, yz_parent
        if not math.isfinite(q_parent) or q_parent > subdivision_q_ceiling:
            context["skipped_high_q"] += 1
            return joint_W, yz_box, yz_parent

        context["attempted"] += 1
        yz_sub, pairs, empty = _joint_yz_dot_subdivision(vd, chart, yz_parent)
        context["pair_evaluations"] += pairs
        context["empty_pairs"] += empty
        dot = vd[0] * chart["cx"] + yz_sub
        support_W = FULL.I(2.0) * wd - dot
        W = V18._intersect(joint_W, support_W)
        w_new, q_new = BASE_QPLUS(float(chart["q1"]), W)

        refined = W.lo > joint_W.lo or W.hi < joint_W.hi
        if refined:
            context["refined"] += 1
            if context["first_refined"] is None:
                context["first_refined"] = {
                    "parent_q_upper": q_parent,
                    "subdivided_q_upper": q_new,
                    "parent_W": joint_W.as_list(),
                    "subdivided_W": W.as_list(),
                    "parent_yz_dot": yz_parent.as_list(),
                    "subdivided_yz_dot": yz_sub.as_list(),
                    "current_cy": chart["cy"].as_list(),
                    "current_cz": chart["cz"].as_list(),
                    "correction_vy": vd[1].as_list(),
                    "correction_vz": vd[2].as_list(),
                }

        new_closed = math.isfinite(q_new) and q_new < Q_TARGET and w_new > 0.0
        if new_closed:
            context["newly_closed"] += 1
            if context["first_newly_closed"] is None:
                context["first_newly_closed"] = {
                    "parent_q_upper": q_parent,
                    "subdivided_q_upper": q_new,
                    "parent_abs_W_lower": w_parent,
                    "subdivided_abs_W_lower": w_new,
                    "parent_W": joint_W.as_list(),
                    "subdivided_W": W.as_list(),
                }
        return W, yz_box, yz_sub

    V18._support_product_scalar = subdivided_support
    try:
        core = V18B.build(
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
        V18._support_product_scalar = original_support

    inherited = V18B.validate(core)
    unclosed = int(core.get("unclosed_q8_cells", -1))
    focused_closed = unclosed == 0 and not inherited
    out = dict(core)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_JOINT_YZ_DIRECTION_SUBDIVISION_Q8_V19",
        "V18B_signed_full_angle_parent_retained": True,
        "joint_current_correction_yz_subdivision_used": True,
        "joint_yz_subdivision_pieces_per_component": 2,
        "correction_quaternion_vector_unit_ball_retained": True,
        "joint_yz_subdivision_q_ceiling": float(subdivision_q_ceiling),
        "joint_yz_subdivision_calls": int(context["calls"]),
        "joint_yz_subdivision_attempted_cells": int(context["attempted"]),
        "joint_yz_subdivision_skipped_parent_closed_cells": int(context["skipped_closed"]),
        "joint_yz_subdivision_skipped_high_q_cells": int(context["skipped_high_q"]),
        "joint_yz_subdivision_pair_evaluations": int(context["pair_evaluations"]),
        "joint_yz_subdivision_empty_pair_cells": int(context["empty_pairs"]),
        "joint_yz_subdivision_refined_cells": int(context["refined"]),
        "joint_yz_subdivision_newly_closed_cells": int(context["newly_closed"]),
        "first_joint_yz_subdivision_refinement": context["first_refined"],
        "first_joint_yz_subdivision_newly_closed": context["first_newly_closed"],
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_JOINT_YZ_DIRECTION_SUBDIVISION_Q8_V19": "PASS" if focused_closed else "NOT_ESTABLISHED",
        "next_obligation": (
            "LIFT_CLOSED_SAMPLE1_PREFIX_TO_ALL_SOURCE_PHASE_CELLS_AND_CONTINUE_SAMPLE2_PREFIX"
            if focused_closed else
            "REFINE_REMAINING_Q8_CELLS_WITH_SOURCE_CORRELATED_BASE_ROW_DIRECTION_SUBDIVISION"
        ),
        "failures": list(dict.fromkeys(inherited)),
    })
    return out


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_JOINT_YZ_DIRECTION_SUBDIVISION_Q8_V19":
        failures.append("qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "V18B_signed_full_angle_parent_retained",
        "joint_current_correction_yz_subdivision_used",
        "correction_quaternion_vector_unit_ball_retained",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "q8_word_promoted_here", "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        failures.append("deployed correction limit changed")
    calls = int(d.get("joint_yz_subdivision_calls", -1))
    attempted = int(d.get("joint_yz_subdivision_attempted_cells", -1))
    skipped_closed = int(d.get("joint_yz_subdivision_skipped_parent_closed_cells", -1))
    skipped_high = int(d.get("joint_yz_subdivision_skipped_high_q_cells", -1))
    refined = int(d.get("joint_yz_subdivision_refined_cells", -1))
    newly = int(d.get("joint_yz_subdivision_newly_closed_cells", -1))
    if calls != attempted + skipped_closed + skipped_high:
        failures.append("joint y-z subdivision call accounting mismatch")
    if not (attempted >= refined >= newly >= 0):
        failures.append("invalid joint y-z refinement accounting")
    if int(d.get("joint_yz_subdivision_pair_evaluations", -1)) < attempted:
        failures.append("invalid joint y-z pair accounting")
    st = d.get("P5_SAMPLE1_JOINT_YZ_DIRECTION_SUBDIVISION_Q8_V19")
    if st == "PASS":
        if int(d.get("unclosed_q8_cells", -1)) != 0:
            failures.append("V19 PASS retains unclosed q8 cells")
    elif st == "NOT_ESTABLISHED":
        if d.get("first_unclosed_q8_cell") is None and not failures:
            failures.append("V19 numerical nonclosure lacks q8 witness")
    else:
        failures.append("invalid V19 status")
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
    ap.add_argument("--subdivision-q-ceiling", type=float, default=SUBDIVISION_Q_CEILING)
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
        subdivision_q_ceiling=args.subdivision_q_ceiling,
    )
    vf = validate(out)
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_SAMPLE1_JOINT_YZ_DIRECTION_SUBDIVISION_Q8_V19"],
        "cells": out.get("evaluated_signed_cayley_cells"),
        "unclosed": out.get("unclosed_q8_cells"),
        "attempted": out["joint_yz_subdivision_attempted_cells"],
        "refined": out["joint_yz_subdivision_refined_cells"],
        "newly_closed": out["joint_yz_subdivision_newly_closed_cells"],
        "skipped_high_q": out["joint_yz_subdivision_skipped_high_q_cells"],
        "pair_evaluations": out["joint_yz_subdivision_pair_evaluations"],
        "empty_pairs": out["joint_yz_subdivision_empty_pair_cells"],
        "first_refined": out["first_joint_yz_subdivision_refinement"],
        "first_newly_closed": out["first_joint_yz_subdivision_newly_closed"],
        "first_unclosed": out.get("first_unclosed_q8_cell"),
        "worst": out.get("worst_q8_cell"),
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
