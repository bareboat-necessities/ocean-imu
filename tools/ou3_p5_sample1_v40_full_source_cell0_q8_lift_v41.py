#!/usr/bin/env python3
"""V41: lift V40 through the complete source-cell-0 signed-Cayley q<8 cover.

The focused V39 construction closes its 64 current subboxes, but the historical
V18B source-cell-0 cover contains 461376 signed-Cayley cells and leaves a large
family open.  Applying the focused 4^3 current subdivision independently to
every open V18B cell would multiply work before using the strongest global
refinement already available.

V40 rigorously tightens the first omitted-PSD Joseph covariance transport while
retaining V38's exact first-PSD mean-correction geometry.  V41 therefore installs
that V40 helper at the common V12D first-PSD proof hook and reruns the complete
V18B partition once.  Unlike focused V39, V41 intentionally lets V18 recompute
its source-correlated current component chart from the refined correction bound:
that chart is part of the newly refined proof construction, not the frozen V18B
reference witness used by V21B's focused provenance equality check.  V18/V18B
still intersect every new support enclosure with their existing parents.

This file is an inventory/lift step only.  It does not change the estimator,
source domain, six-radian shipping correction limit, q<8 target, source-word
language, whole-word promotion, or finite N_H.  PASS means only that the entire
source-cell-0 sample-1 signed-Cayley partition closes under the refined parent.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p5_sample1_first_psd_exact_joseph_components_v40 as V40
import ou3_p5_sample1_signed_cayley_q8_v18b as V18B

DEFAULT_DOMAIN = V18B.DEFAULT_DOMAIN
SCHEMA = 4100
Q_TARGET = V18B.Q_TARGET
V12D = V40.V12D


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    path = Path(domain_path).resolve()
    context = {"calls": 0, "first": None}
    original = V12D._first_psd_perturbation_tangent

    def tracked_v40(**kwargs):
        context["calls"] += 1
        ans = V40._first_psd_perturbation_exact_joseph_components(**kwargs)
        if context["first"] is None:
            context["first"] = {
                "first_offaxis_attitude_correction_upper_rad": float(
                    ans["first_offaxis_attitude_correction_upper_rad"]),
                "first_posterior_covariance_perturbation_upper": float(
                    ans["first_posterior_covariance_perturbation_upper"]),
                "sample1_reduced_covariance_PSD_perturbation_upper": float(
                    ans["sample1_reduced_covariance_PSD_perturbation_upper"]),
                "first_PSD_Joseph_tangent_column_norm_upper": float(
                    ans["first_PSD_Joseph_tangent_column_norm_upper"]),
            }
        return ans

    V12D._first_psd_perturbation_tangent = tracked_v40
    try:
        core = V18B.build(
            path,
            source_pieces=source_pieces,
            source_cell_index=source_cell_index,
            p_pieces=p_pieces,
            tangent_pieces=tangent_pieces,
            axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces,
        )
    finally:
        V12D._first_psd_perturbation_tangent = original

    failures = [f"V18B: {x}" for x in V18B.validate(core)]
    if context["calls"] <= 0:
        failures.append("V40 exact-Joseph first-PSD helper was not exercised")
    if V12D._first_psd_perturbation_tangent is not original:
        failures.append("V12D first-PSD helper was not restored")

    cells = int(core.get("evaluated_signed_cayley_cells", -1))
    unclosed = int(core.get("unclosed_q8_cells", -1))
    full_closed = (
        cells > 0 and unclosed == 0
        and core.get("P5_SAMPLE1_SIGNED_FULL_ANGLE_CURRENT_YZ_Q8_V18B") == "PASS"
        and not failures
    )

    out = dict(core)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_V40_FULL_SOURCE_CELL0_Q8_LIFT_V41",
        "V18B_full_signed_angle_current_yz_parent_retained": True,
        "V40_exact_Joseph_first_PSD_installed_globally": True,
        "V40_exact_Joseph_first_PSD_helper_calls": int(context["calls"]),
        "V40_first_refined_witness": context["first"],
        "refined_V18_current_chart_recomputed_from_V40_bound": True,
        "V18_support_intersections_with_parent_retained": True,
        "temporary_V12D_helper_restored_after_build": (
            V12D._first_psd_perturbation_tangent is original),
        "source_cell_index_certified_here": int(source_cell_index),
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_V40_FULL_SOURCE_CELL0_Q8_LIFT_V41": (
            "PASS" if full_closed else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_SOURCE_CELL0_SAMPLE1_Q8_CLOSURE_TO_ALL_DUE_SOURCE_PHASE_CHILDREN"
            if full_closed else
            "REFINE_REMAINING_V41_FULL_COVER_CELLS_WITH_ADAPTIVE_DIRECTIONAL_CURRENT_SUBDIVISION"
        ),
        "failures": list(dict.fromkeys(failures)),
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_V40_FULL_SOURCE_CELL0_Q8_LIFT_V41":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V18B_full_signed_angle_current_yz_parent_retained",
        "V40_exact_Joseph_first_PSD_installed_globally",
        "refined_V18_current_chart_recomputed_from_V40_bound",
        "V18_support_intersections_with_parent_retained",
        "temporary_V12D_helper_restored_after_build",
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
    if int(d.get("V40_exact_Joseph_first_PSD_helper_calls", 0)) <= 0:
        f.append("no V40 exact-Joseph helper calls")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    st = d.get("P5_SAMPLE1_V40_FULL_SOURCE_CELL0_Q8_LIFT_V41")
    if st == "PASS":
        if int(d.get("unclosed_q8_cells", -1)) != 0:
            f.append("V41 PASS retains unclosed q8 cells")
    elif st == "NOT_ESTABLISHED":
        if d.get("first_unclosed_q8_cell") is None and not f:
            f.append("V41 numerical nonclosure lacks first q8 witness")
    else:
        f.append("invalid V41 status")
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
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_V40_FULL_SOURCE_CELL0_Q8_LIFT_V41"],
        "cells": d.get("evaluated_signed_cayley_cells"),
        "closed": d.get("closed_q8_cells"),
        "unclosed": d.get("unclosed_q8_cells"),
        "max_q": d.get("max_post_sample1_cayley_norm_upper"),
        "helper_calls": d.get("V40_exact_Joseph_first_PSD_helper_calls"),
        "first_refined_witness": d.get("V40_first_refined_witness"),
        "first_unclosed": d.get("first_unclosed_q8_cell"),
        "worst": d.get("worst_q8_cell"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
