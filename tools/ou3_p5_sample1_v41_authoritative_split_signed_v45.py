#!/usr/bin/env python3
"""V45: evaluate V44's split-signed correction on the authoritative V41 chart.

V44 proved that the split signed V28/V31 source geometry is nontrivial on the
new V41 first-survivor row (p,t,a)=(0,0,23): two of four parent-open correction
subcells were source-incompatible.  Its copied current-chart reconstruction,
however, produced q=0.659377... instead of V41's q=0.641521..., so V44
correctly failed its provenance contract.

This stage does not copy the V18 current chart.  It executes the actual
V41/V18B stack with V40 installed, wraps V14's parent-chart call only to track
the current source row, and wraps V18's final support call to capture the full
V18 component chart on the first final-product subcell of (0,0,23).  A private
sentinel then aborts the expensive full cover.  V18/V18B/V41 temporary hooks
unwind through their own finally blocks, and V45 restores its wrappers too.

The captured chart is required to reproduce V41's archived first-survivor
q/cx/cyz witness.  Only then is V44's independently source-derived correction
intersection reevaluated with V16/V15/V18.  This is still a focused diagnostic:
no estimator, source domain, six-radian limit, q<8 target, promotion flag, or
N_H changes.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_v40_full_source_cell0_q8_lift_v41 as V41
import ou3_p5_sample1_first_psd_exact_joseph_components_v40 as V40
import ou3_p5_sample1_signed_cayley_q8_v14 as V14
import ou3_p5_sample1_signed_cayley_q8_v18 as V18
import ou3_p5_sample1_signed_cayley_q8_v18b as V18B
import ou3_p5_sample1_v40_split_signed_first_survivor_v44 as V44

DEFAULT_DOMAIN = V41.DEFAULT_DOMAIN
SCHEMA = 4500
Q_TARGET = 8.0
WITNESS = (0, 0, 23)
V41_Q_CURRENT = 0.6415212986499801
V41_CX = (-0.420125769880562, 0.47923885652470244)
V41_CYZ = 0.6415212986499802
V41_Q_POST = 8.344528951460543


class _CapturedAuthoritativeChart(RuntimeError):
    pass


def _I(pair):
    return Interval.outward_bounds(*map(float, pair))


def _matches(x: float, ref: float, *, atol: float = 2.0e-12) -> bool:
    return math.isfinite(float(x)) and abs(float(x) - float(ref)) <= atol


def _capture_authoritative_chart(path: Path, *, source_pieces: int,
                                 source_cell_index: int, p_pieces: int,
                                 tangent_pieces: int, axial_pieces: int,
                                 residual_x_pieces: int,
                                 parallel_pieces: int) -> dict:
    context = {"ids": None, "captured": None, "parent_calls": 0,
               "support_calls": 0}
    V12D = V41.V12D
    root_psd = V12D._first_psd_perturbation_tangent
    root_parent_chart = V14._sample1_current_chart
    root_support = V18._support_product_scalar

    def tracked_parent_chart(*, first, base, vr, dom, src, sample1_s_angle):
        context["ids"] = (
            int(base["p_cell"]), int(base["tangent_residual_cell"]),
            int(base["axial_residual_cell"]))
        context["parent_calls"] += 1
        return root_parent_chart(
            first=first, base=base, vr=vr, dom=dom, src=src,
            sample1_s_angle=sample1_s_angle)

    def capture_support(parent_W, wd, vd, chart):
        context["support_calls"] += 1
        ans = root_support(parent_W, wd, vd, chart)
        if context["ids"] == WITNESS and context["captured"] is None:
            context["captured"] = {
                "q1": float(chart["q1"]),
                "cx": chart["cx"].as_list(),
                "cy": chart["cy"].as_list(),
                "cz": chart["cz"].as_list(),
                "cyz_norm_upper": float(chart["cyz_norm_upper"]),
                "sample0_q_triangle_upper": float(chart.get("sample0_q_triangle_upper", math.nan)),
                "sample0_q_product_scalar_upper": float(chart.get("sample0_q_product_scalar_upper", math.nan)),
                "sample0_q_selected_upper": float(chart.get("sample0_q_selected_upper", math.nan)),
                "sample1_q_parent_upper": float(chart.get("sample1_q_parent_upper", math.nan)),
                "sample1_q_product_tightened_upper": float(chart.get("sample1_q_product_tightened_upper", math.nan)),
            }
            raise _CapturedAuthoritativeChart("captured V41 first-survivor chart")
        return ans

    V12D._first_psd_perturbation_tangent = \
        V40._first_psd_perturbation_exact_joseph_components
    V14._sample1_current_chart = tracked_parent_chart
    V18._support_product_scalar = capture_support
    try:
        try:
            V18B.build(
                path, source_pieces=source_pieces,
                source_cell_index=source_cell_index, p_pieces=p_pieces,
                tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
                residual_x_pieces=residual_x_pieces,
                parallel_pieces=parallel_pieces)
        except _CapturedAuthoritativeChart:
            pass
    finally:
        V18._support_product_scalar = root_support
        V14._sample1_current_chart = root_parent_chart
        V12D._first_psd_perturbation_tangent = root_psd

    if context["captured"] is None:
        raise RuntimeError("authoritative V41 chart capture was never exercised")
    context["hooks_restored"] = (
        V18._support_product_scalar is root_support
        and V14._sample1_current_chart is root_parent_chart
        and V12D._first_psd_perturbation_tangent is root_psd)
    return context


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    path = Path(domain_path).resolve()
    failures = []
    cap = _capture_authoritative_chart(
        path, source_pieces=source_pieces,
        source_cell_index=source_cell_index, p_pieces=p_pieces,
        tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
        residual_x_pieces=residual_x_pieces,
        parallel_pieces=parallel_pieces)
    chart_raw = cap["captured"]
    chart = {
        "q1": float(chart_raw["q1"]),
        "cx": _I(chart_raw["cx"]),
        "cy": _I(chart_raw["cy"]),
        "cz": _I(chart_raw["cz"]),
        "cyz_norm_upper": float(chart_raw["cyz_norm_upper"]),
    }
    chart_matches = (
        _matches(chart["q1"], V41_Q_CURRENT)
        and _matches(chart["cx"].lo, V41_CX[0])
        and _matches(chart["cx"].hi, V41_CX[1])
        and _matches(chart["cyz_norm_upper"], V41_CYZ))
    if not chart_matches:
        failures.append("captured chart does not reproduce archived V41 first survivor")
    if cap.get("hooks_restored") is not True:
        failures.append("authoritative chart-capture hooks were not restored")

    # V44's source correction geometry is independent of its copied current
    # chart.  Reuse that geometry but discard all V44 q values.
    v44 = V44.build(
        path, source_pieces=source_pieces,
        source_cell_index=source_cell_index, p_pieces=p_pieces,
        tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
        residual_x_pieces=residual_x_pieces,
        parallel_pieces=parallel_pieces)
    parent = v44.get("first_parent_open_subcell") or {}
    refined = v44.get("first_parent_open_after_V44") or {}
    if not parent or not refined:
        failures.append("V44 first split-signed subcell data missing")
        parent_eval = refined_eval = None
    else:
        dbox = [_I(x) for x in parent["correction_component_box_rad"]]
        parent_eval = V44._eval_q(
            q=chart["q1"], chart=chart, dbox=dbox,
            radial_lo=float(parent["correction_radial_lower_rad"]),
            radial_hi=float(parent["correction_radial_upper_rad"]))
        parent_reproduces = _matches(float(parent_eval["best_q"]), V41_Q_POST, atol=3.0e-11)
        if not parent_reproduces:
            failures.append("authoritative chart plus V44 parent does not reproduce V41 q")

        if refined.get("source_incompatible") is True:
            refined_eval = {
                "closed": True, "incompatible": True, "best_q": 0.0,
                "geodesic_q": 0.0, "product_q": 0.0,
                "product_w": math.inf,
            }
        else:
            joint = [_I(x) for x in refined["source_joint_correction_box_rad"]]
            refined_eval = V44._eval_q(
                q=chart["q1"], chart=chart, dbox=joint,
                radial_lo=float(refined["source_radial_lower_rad"]),
                radial_hi=float(refined["source_radial_upper_rad"]))

    parent_q = math.inf if parent_eval is None else float(parent_eval["best_q"])
    refined_q = math.inf if refined_eval is None else float(refined_eval["best_q"])
    refined_closed = bool(refined_eval and refined_eval.get("closed"))
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_V41_AUTHORITATIVE_SPLIT_SIGNED_V45",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "actual_V41_V18B_stack_executed_for_chart_capture": True,
        "V40_exact_Joseph_parent_installed_during_capture": True,
        "capture_aborted_after_first_target_final_subcell": True,
        "capture_hooks_restored": bool(cap.get("hooks_restored")),
        "capture_parent_chart_calls": int(cap.get("parent_calls", 0)),
        "capture_support_calls": int(cap.get("support_calls", 0)),
        "V41_first_survivor_row": list(WITNESS),
        "authoritative_current_chart": chart_raw,
        "archived_V41_current_q_reference": V41_Q_CURRENT,
        "archived_V41_current_cx_reference": list(V41_CX),
        "archived_V41_current_cyz_reference": V41_CYZ,
        "authoritative_chart_matches_archived_V41_witness": chart_matches,
        "V44_split_signed_source_geometry_reused_without_V44_q_values": True,
        "authoritative_parent_best_q_upper": parent_q,
        "archived_V41_post_sample1_q_reference": V41_Q_POST,
        "authoritative_parent_reproduces_V41_post_q": (
            math.isfinite(parent_q) and _matches(parent_q, V41_Q_POST, atol=3.0e-11)),
        "authoritative_refined_best_q_upper": refined_q,
        "authoritative_refined_geodesic_q_upper": (
            None if refined_eval is None else float(refined_eval["geodesic_q"])),
        "authoritative_refined_product_q_upper": (
            None if refined_eval is None else float(refined_eval["product_q"])),
        "authoritative_refined_product_abs_W_lower": (
            None if refined_eval is None else float(refined_eval["product_w"])),
        "first_V41_survivor_closed_by_authoritative_V45": refined_closed,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_V41_AUTHORITATIVE_SPLIT_SIGNED_V45": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_AUTHORITATIVE_V45_SPLIT_SIGNED_REFINEMENT_OVER_FULL_V41_COVER"
            if refined_closed and not failures else
            "PARTITION_AUTHORITATIVE_V41_CURRENT_COMPONENT_BOX_AND_INTERSECT_V43_AND_V44_CORRECTION_IMAGES"
        ),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "actual_V41_V18B_stack_executed_for_chart_capture",
        "V40_exact_Joseph_parent_installed_during_capture",
        "capture_aborted_after_first_target_final_subcell",
        "capture_hooks_restored",
        "authoritative_chart_matches_archived_V41_witness",
        "V44_split_signed_source_geometry_reused_without_V44_q_values",
        "authoritative_parent_reproduces_V41_post_q",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "q8_composed_here", "q8_word_promoted_here", "whole_word_promoted_here",
        "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if int(d.get("capture_parent_chart_calls", 0)) <= 0 or int(d.get("capture_support_calls", 0)) <= 0:
        f.append("authoritative chart capture was not exercised")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    if d.get("P5_SAMPLE1_V41_AUTHORITATIVE_SPLIT_SIGNED_V45") not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid V45 status")
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
    vf = validate(d); d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_V41_AUTHORITATIVE_SPLIT_SIGNED_V45"],
        "chart": d["authoritative_current_chart"],
        "parent_q": d["authoritative_parent_best_q_upper"],
        "refined_q": d["authoritative_refined_best_q_upper"],
        "refined_geodesic": d["authoritative_refined_geodesic_q_upper"],
        "refined_product": d["authoritative_refined_product_q_upper"],
        "closed": d["first_V41_survivor_closed_by_authoritative_V45"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
