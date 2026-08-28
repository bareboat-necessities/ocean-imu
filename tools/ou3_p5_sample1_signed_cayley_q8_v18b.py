#!/usr/bin/env python3
"""V18B: audited accounting wrapper for V18's signed full-angle proof gauge.

V18's first implementation used the deployed-quaternion nonnegative half-angle
radial helper on a full signed proof-gauge angle.  The shared
``ou3_p5_signed_full_angle_trig`` primitive now supplies the validated signed
full-angle enclosure directly to V18.  V18B retains the historical accounting
contract: it counts every signed-angle call and every deliberate broad fallback
while leaving all V18 current-y/z algebra and earlier parents unchanged.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_signed_cayley_q8_v14 as V14
import ou3_p5_sample1_signed_cayley_q8_v18 as V18
import ou3_p5_signed_full_angle_trig as SIGNED_TRIG

DEFAULT_DOMAIN = V18.DEFAULT_DOMAIN
SCHEMA = 1801
FULL = V14.FULL
Q_TARGET = V14.Q_TARGET
AUDITED_POINT_ABS_MAX = SIGNED_TRIG.AUDITED_POINT_ABS_MAX
_signed_full_angle_trig_interval = SIGNED_TRIG.signed_full_angle_trig_interval


def _rotate_yz_rx_transpose(cy: Interval, cz: Interval,
                            angle: Interval) -> tuple[Interval, Interval]:
    sinx, cosx, _broad = _signed_full_angle_trig_interval(angle)
    # R_x(a)^T has yz block [[cos(a), sin(a)],[-sin(a), cos(a)]].
    return cosx * cy + sinx * cz, -(sinx * cy) + cosx * cz


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    original_rotate = V18._rotate_yz_rx_transpose
    context = {"calls": 0, "broad": 0}

    def tracked_rotate(cy: Interval, cz: Interval, angle: Interval):
        context["calls"] += 1
        sinx, cosx, broad = _signed_full_angle_trig_interval(angle)
        context["broad"] += int(broad)
        return cosx * cy + sinx * cz, -(sinx * cy) + cosx * cz

    V18._rotate_yz_rx_transpose = tracked_rotate
    try:
        core = V18.build(
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
        V18._rotate_yz_rx_transpose = original_rotate

    inherited = V18.validate(core)
    parent_status = core.get("P5_SAMPLE1_CURRENT_YZ_SUPPORT_SIGNED_CAYLEY_Q8_V18")
    out = dict(core)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_SIGNED_FULL_ANGLE_CURRENT_YZ_Q8_V18B",
        "V18_current_yz_support_parent_retained": True,
        "signed_full_angle_proof_gauge_trig_used": True,
        "signed_full_angle_trig_uses_validated_point_backend": True,
        "signed_full_angle_trig_uses_no_range_reduction": True,
        "signed_full_angle_trig_uses_no_libm_proof_call": True,
        "signed_full_angle_trig_shared_with_V18": True,
        "signed_full_angle_trig_calls": int(context["calls"]),
        "signed_full_angle_trig_broad_fallback_calls": int(context["broad"]),
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_SIGNED_FULL_ANGLE_CURRENT_YZ_Q8_V18B": (
            "PASS" if parent_status == "PASS" and not inherited else "NOT_ESTABLISHED"
        ),
        "next_obligation": (
            "LIFT_CLOSED_SAMPLE1_PREFIX_TO_ALL_SOURCE_PHASE_CELLS_AND_CONTINUE_SAMPLE2_PREFIX"
            if parent_status == "PASS" and not inherited else
            core.get("next_obligation", "REFINE_REMAINING_Q8_CELLS_WITH_JOINT_CURRENT_CORRECTION_YZ_DIRECTION_SUBDIVISION")
        ),
        "failures": list(dict.fromkeys(inherited)),
    })
    return out


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_SIGNED_FULL_ANGLE_CURRENT_YZ_Q8_V18B":
        failures.append("qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "V18_current_yz_support_parent_retained",
        "signed_full_angle_proof_gauge_trig_used",
        "signed_full_angle_trig_uses_validated_point_backend",
        "signed_full_angle_trig_uses_no_range_reduction",
        "signed_full_angle_trig_uses_no_libm_proof_call",
        "signed_full_angle_trig_shared_with_V18",
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
    calls = int(d.get("signed_full_angle_trig_calls", -1))
    broad = int(d.get("signed_full_angle_trig_broad_fallback_calls", -1))
    if not (calls >= broad >= 0):
        failures.append("invalid signed full-angle trig accounting")
    st = d.get("P5_SAMPLE1_SIGNED_FULL_ANGLE_CURRENT_YZ_Q8_V18B")
    if st == "PASS":
        if int(d.get("unclosed_q8_cells", -1)) != 0:
            failures.append("V18B PASS retains unclosed q8 cells")
        if not float(d.get("max_post_sample1_cayley_norm_upper", math.inf)) < Q_TARGET:
            failures.append("V18B PASS does not satisfy strict q<8")
    elif st == "NOT_ESTABLISHED":
        if d.get("complete_sample1_branch_closed_here") is True:
            failures.append("V18B nonclosure claims sample1 closure")
        if d.get("first_unclosed_q8_cell") is None and not failures:
            failures.append("V18B numerical nonclosure lacks q8 witness")
    else:
        failures.append("invalid V18B status")
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
        "status": out["P5_SAMPLE1_SIGNED_FULL_ANGLE_CURRENT_YZ_Q8_V18B"],
        "cells": out.get("evaluated_signed_cayley_cells"),
        "unclosed": out.get("unclosed_q8_cells"),
        "signed_trig_calls": out["signed_full_angle_trig_calls"],
        "signed_trig_broad": out["signed_full_angle_trig_broad_fallback_calls"],
        "yz_support_calls": out.get("current_yz_support_qplus_calls"),
        "yz_support_refined": out.get("current_yz_support_refined_cells"),
        "yz_support_newly_closed": out.get("current_yz_support_newly_closed_cells"),
        "first_yz_support_refinement": out.get("first_current_yz_support_refinement"),
        "first_yz_support_newly_closed": out.get("first_current_yz_support_newly_closed"),
        "first_unclosed": out.get("first_unclosed_q8_cell"),
        "worst": out.get("worst_q8_cell"),
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
