#!/usr/bin/env python3
"""V42: adaptive directional-current subdivision over the V41 full q<8 cover.

V41 installs the exact V40 first-PSD Joseph transport globally and substantially
reduces the complete source-cell-0 sample-1 q<8 obstruction, but its remaining
signed-product cells still evaluate the current Cayley component box as one
interval object.  The exact product scalar

    W = 2 w_d - v_d^T c

retains directional dependence on the current Cayley vector c.  V42 refines
only cells for which V18's component-supported product route is still open.  It
bisects the current component with the largest certified contribution

    sup |v_{d,i}| * width(c_i)

to the dot-product uncertainty, recomputes the exact V18 component/Cauchy
intersection on each child, and certifies the union by taking the minimum child
|W| lower bound and the maximum child q upper bound.  Empty children outside
the already certified current q-ball are discarded.  Recursion stops as soon
as every child closes q<8 or the declared finite depth is reached.

The subdivision is a proof partition only.  It does not alter the estimator,
source domain, source branches, six-radian shipping correction limit, q<8
target, source-word language, theorem promotion, or N_H.  V40 remains installed
through V41 exactly as before.  PASS means only that the complete source-cell-0
sample-1 signed-Cayley cover closes under this dependency-preserving refinement.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_v40_full_source_cell0_q8_lift_v41 as V41

DEFAULT_DOMAIN = V41.DEFAULT_DOMAIN
SCHEMA = 4200
Q_TARGET = V41.Q_TARGET
V18B = V41.V18B
V18 = V18B.V18
V14 = V18.V14
FULL = V14.FULL


def _components(chart: dict) -> list[Interval]:
    return [chart["cx"], chart["cy"], chart["cz"]]


def _chart_from_components(parent: dict, comps: list[Interval]) -> dict:
    out = dict(parent)
    out["cx"], out["cy"], out["cz"] = comps
    out["cyz_norm_upper"] = min(
        float(parent["cyz_norm_upper"]),
        V18._yz_norm_upper(comps[1], comps[2]),
    )
    return out


def _split_interval(x: Interval) -> tuple[Interval, Interval] | None:
    if not (math.isfinite(x.lo) and math.isfinite(x.hi) and x.lo < x.hi):
        return None
    mid = x.lo + 0.5 * (x.hi - x.lo)
    if not (x.lo < mid < x.hi):
        return None
    return Interval(x.lo, mid), Interval(mid, x.hi)


def _split_dimension(vd, chart: dict) -> int | None:
    comps = _components(chart)
    scores = []
    for i, c in enumerate(comps):
        width = max(0.0, float(c.hi - c.lo))
        scores.append(float(vd[i].abs_upper()) * width)
    best = max(range(3), key=lambda i: scores[i])
    if scores[best] > 0.0 and _split_interval(comps[best]) is not None:
        return best
    widths = [max(0.0, float(c.hi - c.lo)) for c in comps]
    best = max(range(3), key=lambda i: widths[i])
    return best if widths[best] > 0.0 and _split_interval(comps[best]) is not None else None


def _radially_nonempty(q_parent: float, comps: list[Interval]) -> bool:
    return V14.CAYLEY2._norm_lower(comps) <= FULL.up(float(q_parent))


def _leaf_product(*, q_parent: float, wd: Interval, vd,
                  chart: dict, support_fn, qplus_fn) -> dict | None:
    comps = _components(chart)
    if not _radially_nonempty(q_parent, comps):
        return None
    q_child = min(float(q_parent), float(V14.CAYLEY1._norm_upper(comps)))
    parent_W = FULL.I(2.0) * wd - V14.CAYLEY1._dot(vd, comps)
    joint_W, _yz_box, _yz_joint = support_fn(parent_W, wd, vd, chart)
    w, q = qplus_fn(q_child, joint_W)
    return {
        "abs_W_lower": float(w),
        "q_upper": float(q),
        "q_current_upper": float(q_child),
        "closed": bool(math.isfinite(q) and q < Q_TARGET and w > 0.0),
    }


def _adaptive_product_cover(*, q_parent: float, wd: Interval, vd,
                            chart: dict, max_depth: int, support_fn=None,
                            qplus_fn=None, depth: int = 0) -> dict:
    if max_depth < 0:
        raise ValueError("max_depth must be nonnegative")
    support_fn = V18._support_product_scalar if support_fn is None else support_fn
    qplus_fn = V14._qplus_from_product_scalar if qplus_fn is None else qplus_fn

    leaf = _leaf_product(
        q_parent=q_parent, wd=wd, vd=vd, chart=chart,
        support_fn=support_fn, qplus_fn=qplus_fn)
    if leaf is None:
        return {
            "empty": True, "closed": True, "abs_W_lower": math.inf,
            "q_upper": 0.0, "leaf_evaluations": 0, "split_nodes": 0,
            "max_depth_used": depth, "first_split_dimension": None,
        }
    if leaf["closed"] or depth >= max_depth:
        return {
            "empty": False, "closed": leaf["closed"],
            "abs_W_lower": leaf["abs_W_lower"], "q_upper": leaf["q_upper"],
            "leaf_evaluations": 1, "split_nodes": 0,
            "max_depth_used": depth, "first_split_dimension": None,
        }

    dim = _split_dimension(vd, chart)
    if dim is None:
        return {
            "empty": False, "closed": False,
            "abs_W_lower": leaf["abs_W_lower"], "q_upper": leaf["q_upper"],
            "leaf_evaluations": 1, "split_nodes": 0,
            "max_depth_used": depth, "first_split_dimension": None,
        }
    parts = _split_interval(_components(chart)[dim])
    if parts is None:
        return {
            "empty": False, "closed": False,
            "abs_W_lower": leaf["abs_W_lower"], "q_upper": leaf["q_upper"],
            "leaf_evaluations": 1, "split_nodes": 0,
            "max_depth_used": depth, "first_split_dimension": None,
        }

    children = []
    for part in parts:
        comps = list(_components(chart))
        comps[dim] = part
        child_chart = _chart_from_components(chart, comps)
        child = _adaptive_product_cover(
            q_parent=q_parent, wd=wd, vd=vd, chart=child_chart,
            max_depth=max_depth, support_fn=support_fn, qplus_fn=qplus_fn,
            depth=depth + 1)
        if not child["empty"]:
            children.append(child)

    if not children:
        return {
            "empty": True, "closed": True, "abs_W_lower": math.inf,
            "q_upper": 0.0, "leaf_evaluations": 0, "split_nodes": 1,
            "max_depth_used": depth + 1, "first_split_dimension": dim,
        }
    return {
        "empty": False,
        "closed": all(bool(c["closed"]) for c in children),
        "abs_W_lower": min(float(c["abs_W_lower"]) for c in children),
        "q_upper": max(float(c["q_upper"]) for c in children),
        "leaf_evaluations": sum(int(c["leaf_evaluations"]) for c in children),
        "split_nodes": 1 + sum(int(c["split_nodes"]) for c in children),
        "max_depth_used": max(int(c["max_depth_used"]) for c in children),
        "first_split_dimension": dim,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6,
          adaptive_depth: int = 1) -> dict:
    if adaptive_depth < 1:
        raise ValueError("adaptive_depth must be at least one")
    path = Path(domain_path).resolve()

    original_support = V18._support_product_scalar
    original_qplus = V14._qplus_from_product_scalar
    pending: list[dict] = []
    context = {
        "parent_qplus_calls": 0,
        "eligible_open_support_calls": 0,
        "adaptive_refined_support_calls": 0,
        "adaptive_newly_closed_support_calls": 0,
        "adaptive_leaf_evaluations": 0,
        "adaptive_split_nodes": 0,
        "adaptive_max_depth_used": 0,
        "first_refinement": None,
        "first_newly_closed": None,
        "last_parent_result": None,
    }

    def tracked_support(parent_W, wd, vd, chart):
        ans = original_support(parent_W, wd, vd, chart)
        pending.append({
            "parent_W": parent_W,
            "wd": wd,
            "vd": vd,
            "chart": dict(chart),
            "parent_result": context.get("last_parent_result"),
        })
        context["last_parent_result"] = None
        return ans

    def adaptive_qplus(q1: float, W: Interval):
        if not pending:
            ans = original_qplus(q1, W)
            context["parent_qplus_calls"] += 1
            context["last_parent_result"] = ans
            return ans

        item = pending.pop()
        base_w, base_q = original_qplus(q1, W)
        if math.isfinite(base_q) and base_q < Q_TARGET and base_w > 0.0:
            return base_w, base_q

        context["eligible_open_support_calls"] += 1
        refined = _adaptive_product_cover(
            q_parent=q1, wd=item["wd"], vd=item["vd"], chart=item["chart"],
            max_depth=adaptive_depth, support_fn=original_support,
            qplus_fn=original_qplus)
        context["adaptive_leaf_evaluations"] += int(refined["leaf_evaluations"])
        context["adaptive_split_nodes"] += int(refined["split_nodes"])
        context["adaptive_max_depth_used"] = max(
            int(context["adaptive_max_depth_used"]), int(refined["max_depth_used"]))

        if refined["empty"]:
            return base_w, base_q
        ref_w = float(refined["abs_W_lower"])
        ref_q = float(refined["q_upper"])
        out_w = max(float(base_w), ref_w)
        out_q = min(float(base_q), ref_q)
        improved = out_q < float(base_q) or out_w > float(base_w)
        if improved:
            context["adaptive_refined_support_calls"] += 1
            row = {
                "base_product_abs_W_lower": float(base_w),
                "base_product_q_upper": float(base_q),
                "adaptive_product_abs_W_lower": out_w,
                "adaptive_product_q_upper": out_q,
                "current_q_parent_upper": float(q1),
                "first_split_dimension": refined["first_split_dimension"],
                "leaf_evaluations": int(refined["leaf_evaluations"]),
                "split_nodes": int(refined["split_nodes"]),
                "max_depth_used": int(refined["max_depth_used"]),
            }
            if context["first_refinement"] is None:
                context["first_refinement"] = row
            parent_result = item.get("parent_result")
            parent_closed = bool(
                parent_result is not None
                and math.isfinite(float(parent_result[1]))
                and float(parent_result[1]) < Q_TARGET
                and float(parent_result[0]) > 0.0)
            newly_closed = (
                not parent_closed
                and not (math.isfinite(base_q) and base_q < Q_TARGET and base_w > 0.0)
                and math.isfinite(out_q) and out_q < Q_TARGET and out_w > 0.0)
            if newly_closed:
                context["adaptive_newly_closed_support_calls"] += 1
                if context["first_newly_closed"] is None:
                    context["first_newly_closed"] = row
        return out_w, out_q

    V18._support_product_scalar = tracked_support
    V14._qplus_from_product_scalar = adaptive_qplus
    try:
        core = V41.build(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces)
    finally:
        V14._qplus_from_product_scalar = original_qplus
        V18._support_product_scalar = original_support

    failures = [f"V41: {x}" for x in V41.validate(core)]
    if pending:
        failures.append("adaptive current subdivision left unmatched support contexts")
    if int(context["eligible_open_support_calls"]) <= 0:
        failures.append("adaptive current subdivision was never exercised")
    restored = (
        V18._support_product_scalar is original_support
        and V14._qplus_from_product_scalar is original_qplus)
    if not restored:
        failures.append("temporary adaptive current helpers were not restored")

    cells = int(core.get("evaluated_signed_cayley_cells", -1))
    unclosed = int(core.get("unclosed_q8_cells", -1))
    full_closed = (
        cells > 0 and unclosed == 0
        and core.get("P5_SAMPLE1_V40_FULL_SOURCE_CELL0_Q8_LIFT_V41") == "PASS"
        and not failures)

    out = dict(core)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_V41_ADAPTIVE_DIRECTIONAL_CURRENT_V42",
        "V41_full_source_cell0_parent_retained": True,
        "V40_exact_Joseph_first_PSD_retained_through_V41": True,
        "adaptive_directional_current_subdivision_used": True,
        "adaptive_split_strategy": "MAX_ABS_CORRECTION_QUATERNION_COMPONENT_TIMES_CURRENT_COMPONENT_WIDTH",
        "adaptive_depth": int(adaptive_depth),
        "adaptive_parent_qplus_calls": int(context["parent_qplus_calls"]),
        "adaptive_eligible_open_support_calls": int(context["eligible_open_support_calls"]),
        "adaptive_refined_support_calls": int(context["adaptive_refined_support_calls"]),
        "adaptive_newly_closed_support_calls": int(context["adaptive_newly_closed_support_calls"]),
        "adaptive_leaf_evaluations": int(context["adaptive_leaf_evaluations"]),
        "adaptive_split_nodes": int(context["adaptive_split_nodes"]),
        "adaptive_max_depth_used": int(context["adaptive_max_depth_used"]),
        "first_adaptive_refinement": context["first_refinement"],
        "first_adaptive_newly_closed": context["first_newly_closed"],
        "adaptive_children_outside_parent_q_ball_discarded_only_by_norm_lower": True,
        "adaptive_union_uses_min_abs_W_and_max_q": True,
        "temporary_adaptive_helpers_restored_after_build": restored,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_V41_ADAPTIVE_DIRECTIONAL_CURRENT_V42": (
            "PASS" if full_closed else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_SOURCE_CELL0_SAMPLE1_Q8_CLOSURE_TO_ALL_DUE_SOURCE_PHASE_CHILDREN"
            if full_closed else
            ("INCREASE_ADAPTIVE_DIRECTIONAL_CURRENT_DEPTH_ON_REMAINING_V42_CELLS"
             if int(context["adaptive_refined_support_calls"]) > 0 else
             "COUPLE_CORRECTION_INTERVAL_TO_CURRENT_CHILD_BEFORE_MORE_SUBDIVISION")),
        "failures": list(dict.fromkeys(failures)),
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_V41_ADAPTIVE_DIRECTIONAL_CURRENT_V42":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V41_full_source_cell0_parent_retained",
        "V40_exact_Joseph_first_PSD_retained_through_V41",
        "adaptive_directional_current_subdivision_used",
        "adaptive_children_outside_parent_q_ball_discarded_only_by_norm_lower",
        "adaptive_union_uses_min_abs_W_and_max_q",
        "temporary_adaptive_helpers_restored_after_build",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed",
        "deployed_correction_limit_increased", "q8_word_promoted_here",
        "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if int(d.get("adaptive_depth", 0)) < 1:
        f.append("adaptive depth is not positive")
    if int(d.get("adaptive_eligible_open_support_calls", 0)) <= 0:
        f.append("adaptive current subdivision has no eligible calls")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    st = d.get("P5_SAMPLE1_V41_ADAPTIVE_DIRECTIONAL_CURRENT_V42")
    if st == "PASS":
        if int(d.get("unclosed_q8_cells", -1)) != 0:
            f.append("V42 PASS retains unclosed q8 cells")
    elif st == "NOT_ESTABLISHED":
        if d.get("first_unclosed_q8_cell") is None and not f:
            f.append("V42 numerical nonclosure lacks first q8 witness")
    else:
        f.append("invalid V42 status")
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
    ap.add_argument("--adaptive-depth", type=int, default=1)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(
        x.domain, source_pieces=x.source_pieces,
        source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
        tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces,
        residual_x_pieces=x.residual_x_pieces,
        parallel_pieces=x.parallel_pieces, adaptive_depth=x.adaptive_depth)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_V41_ADAPTIVE_DIRECTIONAL_CURRENT_V42"],
        "cells": d.get("evaluated_signed_cayley_cells"),
        "closed": d.get("closed_q8_cells"),
        "unclosed": d.get("unclosed_q8_cells"),
        "max_q": d.get("max_post_sample1_cayley_norm_upper"),
        "eligible": d.get("adaptive_eligible_open_support_calls"),
        "refined": d.get("adaptive_refined_support_calls"),
        "newly_closed": d.get("adaptive_newly_closed_support_calls"),
        "leaf_evaluations": d.get("adaptive_leaf_evaluations"),
        "first_refinement": d.get("first_adaptive_refinement"),
        "first_newly_closed": d.get("first_adaptive_newly_closed"),
        "first_unclosed": d.get("first_unclosed_q8_cell"),
        "worst": d.get("worst_q8_cell"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
