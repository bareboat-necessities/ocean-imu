#!/usr/bin/env python3
"""V37: lift the V36 PSD off-diagonal cone through the V34 64-subbox cover.

V36 proves that after absorbing the diagonal PSD remainder into V12D's existing
t/Y intervals, the omitted symmetric off-diagonal attitude covariance has
operator norm <=eps rather than the historical 2eps.  V37 changes only that
first-PSD helper, then reuses V34 unchanged:

  V32 exact theta-x Delta-C,
  V34 first-row Delta-S / sparse theta-yz gain rows,
  V31/V23 current-Cayley 4x4x4 cover,
  V16/V15/V18 q<8 tests.

The V36 helper is installed temporarily and restored after the build.  No
estimator parameter, source domain, correction limit, q target, or promotion
criterion changes.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_first_psd_offdiagonal_cone_v36 as V36
import ou3_p5_sample1_directional_innovation_row_lift_v34 as V34

DEFAULT_DOMAIN = V34.DEFAULT_DOMAIN
SCHEMA = 3700
Q_TARGET = V34.Q_TARGET
V12D = V34.V31.V23.V22.V21B.V21.V12D


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6,
          current_component_pieces: int = 4) -> dict:
    context = {"psd_cone_calls": 0}
    original = V12D._first_psd_perturbation_tangent

    def tracked_psd_cone(**kwargs):
        context["psd_cone_calls"] += 1
        return V36._first_psd_perturbation_psd_cone(**kwargs)

    V12D._first_psd_perturbation_tangent = tracked_psd_cone
    try:
        parent = V34.build(
            Path(domain_path).resolve(),
            source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces,
            current_component_pieces=current_component_pieces)
    finally:
        V12D._first_psd_perturbation_tangent = original

    failures = [f"V34: {x}" for x in V34.validate(parent)]
    if parent.get("P5_SAMPLE1_DIRECTIONAL_INNOVATION_ROW_LIFT_V34") != "PASS":
        failures.append("V34 parent did not pass under V36 PSD cone")
    if context["psd_cone_calls"] <= 0:
        failures.append("V36 PSD-cone helper was not exercised")

    open_count = int(parent.get("open_current_subboxes", -1))
    out = dict(parent)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_V36_PSD_CONE_CURRENT_SUBBOX_LIFT_V37",
        "V34_directional_innovation_construction_retained": True,
        "V36_PSD_offdiagonal_cone_installed": True,
        "V36_PSD_cone_helper_calls": int(context["psd_cone_calls"]),
        "V36_changes_only_first_PSD_offdiagonal_operator_bound": True,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_V36_PSD_CONE_CURRENT_SUBBOX_LIFT_V37": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_V37_ROW_REFINEMENT_INTO_FULL_SOURCE_CELL0_Q8_COVER"
            if open_count == 0 and not failures else
            "REFINE_FIRST_REMAINING_V37_SUBBOX_WITH_EXACT_FIRST_PSD_RESET_COMPONENT_MATRIX"
        ),
        "failures": list(dict.fromkeys(failures)),
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_V36_PSD_CONE_CURRENT_SUBBOX_LIFT_V37":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V34_directional_innovation_construction_retained",
        "V36_PSD_offdiagonal_cone_installed",
        "V36_changes_only_first_PSD_offdiagonal_operator_bound",
        "all_candidate_current_subboxes_accounted",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed",
        "deployed_correction_limit_increased", "q8_composed_here",
        "q8_word_promoted_here", "whole_word_promoted_here",
        "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if int(d.get("V36_PSD_cone_helper_calls", 0)) <= 0:
        f.append("no V36 PSD-cone helper calls")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    if d.get("P5_SAMPLE1_V36_PSD_CONE_CURRENT_SUBBOX_LIFT_V37") not in (
            "PASS", "NOT_ESTABLISHED"):
        f.append("invalid V37 status")
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
    ap.add_argument("--current-component-pieces", type=int, default=4)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(
        x.domain, source_pieces=x.source_pieces,
        source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
        tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces,
        residual_x_pieces=x.residual_x_pieces,
        parallel_pieces=x.parallel_pieces,
        current_component_pieces=x.current_component_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_V36_PSD_CONE_CURRENT_SUBBOX_LIFT_V37"],
        "psd_cone_calls": d.get("V36_PSD_cone_helper_calls"),
        "directional_innovation": d.get("directional_innovation_detail"),
        "candidate": d.get("candidate_current_subboxes"),
        "closed": d.get("closed_current_subboxes"),
        "open": d.get("open_current_subboxes"),
        "minimum_best_q": d.get("minimum_best_q_upper"),
        "maximum_best_q": d.get("maximum_best_q_upper"),
        "first_open": d.get("first_open_current_subbox"),
        "worst_open": d.get("worst_open_current_subbox"),
        "witness_closed": d.get("focused_first_witness_signed_subcell_closed_by_V30_lift"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
