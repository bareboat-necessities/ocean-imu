#!/usr/bin/env python3
"""V39: lift V38 while freezing the authoritative V18B current chart.

V38 gives a large, rigorous tightening of the first PSD-induced correction
perturbation.  A naive global installation of the V36/V38 helper also changes
the V12D row passed into V21's current-chart reconstruction.  That chart is a
*different proof role*: V21B deliberately requires it to reproduce the already
measured V18B first-witness current q exactly.  V37 therefore closed all 64
artificial current subboxes numerically but correctly failed its provenance
contract because the refined V12D row leaked into that authoritative chart.

V39 separates the two roles.  Before installing V38 it builds the baseline
V12D witness row.  During the refined V34/V31 lift:

  * V21._current_component_chart is wrapped so its ``vr`` argument is replaced
    by that baseline V12D witness, preserving the authoritative V18B current
    Cayley geometry and V21B equality check;
  * every other V12D consumer receives V38's exact canonical first-PSD
    correction geometry, so residual/covariance/gain perturbation bounds are
    tightened on the correction side;
  * the existing V32/V34/V31/V23/V16/V15/V18 composition is otherwise
    unchanged.

Both temporary patches are restored after the build.  No estimator setting,
source domain, six-radian shipping correction limit, q<8 target, or theorem
promotion state changes.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_first_psd_exact_correction_geometry_v38 as V38
import ou3_p5_sample1_directional_innovation_row_lift_v34 as V34

DEFAULT_DOMAIN = V34.DEFAULT_DOMAIN
SCHEMA = 3900
Q_TARGET = V34.Q_TARGET
V31 = V34.V31
V30 = V31.V30
V23 = V31.V23
V22 = V23.V22
V21B = V22.V21B
V21 = V21B.V21
V12D = V21.V12D
FULL = V12D.FULL


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6,
          current_component_pieces: int = 4) -> dict:
    path = Path(domain_path).resolve()

    # Baseline row is the current-state parent.  It is intentionally computed
    # before installing V38 and is used only by V21's current-chart helper.
    baseline_v12 = V12D.build(
        path, source_pieces=source_pieces,
        source_cell_index=source_cell_index, p_pieces=p_pieces,
        tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    failures = [f"baseline V12D: {x}" for x in V12D.validate(baseline_v12)]
    baseline_vr = V30._witness_row(baseline_v12)
    baseline_dd = float(baseline_vr["first_offaxis_attitude_correction_upper_rad"])

    context = {"exact_psd_calls": 0, "frozen_chart_calls": 0}
    original_psd = V12D._first_psd_perturbation_tangent
    original_chart = V21._current_component_chart

    def tracked_exact_psd(**kwargs):
        context["exact_psd_calls"] += 1
        return V38._first_psd_perturbation_exact_correction(**kwargs)

    def frozen_current_chart(*, first, base, vr, dom, src, sample1_s_angle):
        context["frozen_chart_calls"] += 1
        return original_chart(
            first=first, base=base, vr=baseline_vr, dom=dom, src=src,
            sample1_s_angle=sample1_s_angle)

    V12D._first_psd_perturbation_tangent = tracked_exact_psd
    V21._current_component_chart = frozen_current_chart
    try:
        parent = V34.build(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces,
            current_component_pieces=current_component_pieces)
    finally:
        V21._current_component_chart = original_chart
        V12D._first_psd_perturbation_tangent = original_psd

    failures += [f"V34: {x}" for x in V34.validate(parent)]
    if parent.get("P5_SAMPLE1_DIRECTIONAL_INNOVATION_ROW_LIFT_V34") != "PASS":
        failures.append("V34 parent did not pass under isolated V38 refinement")
    if context["exact_psd_calls"] <= 0:
        failures.append("V38 exact first-PSD helper was not exercised")
    if context["frozen_chart_calls"] <= 0:
        failures.append("authoritative V21 current-chart freeze was not exercised")

    q_current = float(parent.get("sample1_current_cayley_norm_upper", math.inf))
    q_ref = float(V21B.V18B_FIRST_WITNESS_CURRENT_Q)
    q_matches = V21B._matches_reference(q_current)
    if not q_matches:
        failures.append("V39 did not preserve authoritative V18B first-witness current q")

    open_count = int(parent.get("open_current_subboxes", -1))
    out = dict(parent)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_V38_ISOLATED_EXACT_PSD_CURRENT_SUBBOX_LIFT_V39",
        "V34_directional_innovation_construction_retained": True,
        "V38_exact_first_PSD_correction_installed": True,
        "V38_exact_first_PSD_helper_calls": int(context["exact_psd_calls"]),
        "V38_exact_canonical_tangent_geometry_retained": True,
        "V36_full_Joseph_gain_operator_parent_retained": True,
        "authoritative_V18B_current_chart_frozen_to_baseline_V12D_witness": True,
        "authoritative_V18B_current_chart_freeze_calls": int(context["frozen_chart_calls"]),
        "baseline_V12D_first_offaxis_correction_upper_rad": baseline_dd,
        "authoritative_V18B_first_witness_current_q_reference": q_ref,
        "current_q_matches_authoritative_V18B_reference": q_matches,
        "refined_PSD_used_only_outside_authoritative_current_chart": True,
        "temporary_helpers_restored_after_build": (
            V12D._first_psd_perturbation_tangent is original_psd
            and V21._current_component_chart is original_chart),
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_V38_ISOLATED_EXACT_PSD_CURRENT_SUBBOX_LIFT_V39": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_V39_REFINEMENT_INTO_FULL_SOURCE_CELL0_Q8_COVER"
            if open_count == 0 and not failures else
            "REFINE_FIRST_REMAINING_V39_SUBBOX_WITH_EXACT_FIRST_PSD_JOSEPH_COMPONENT_MATRIX"
        ),
        "failures": list(dict.fromkeys(failures)),
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_V38_ISOLATED_EXACT_PSD_CURRENT_SUBBOX_LIFT_V39":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V34_directional_innovation_construction_retained",
        "V38_exact_first_PSD_correction_installed",
        "V38_exact_canonical_tangent_geometry_retained",
        "V36_full_Joseph_gain_operator_parent_retained",
        "authoritative_V18B_current_chart_frozen_to_baseline_V12D_witness",
        "current_q_matches_authoritative_V18B_reference",
        "refined_PSD_used_only_outside_authoritative_current_chart",
        "temporary_helpers_restored_after_build",
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
    if int(d.get("V38_exact_first_PSD_helper_calls", 0)) <= 0:
        f.append("no V38 exact first-PSD helper calls")
    if int(d.get("authoritative_V18B_current_chart_freeze_calls", 0)) <= 0:
        f.append("no authoritative current-chart freeze calls")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    q = d.get("sample1_current_cayley_norm_upper")
    if not isinstance(q, (int, float)) or not V21B._matches_reference(float(q)):
        f.append("authoritative first-witness current q mismatch")
    if d.get("P5_SAMPLE1_V38_ISOLATED_EXACT_PSD_CURRENT_SUBBOX_LIFT_V39") not in (
            "PASS", "NOT_ESTABLISHED"):
        f.append("invalid V39 status")
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
        "status": d["P5_SAMPLE1_V38_ISOLATED_EXACT_PSD_CURRENT_SUBBOX_LIFT_V39"],
        "exact_psd_calls": d.get("V38_exact_first_PSD_helper_calls"),
        "frozen_chart_calls": d.get("authoritative_V18B_current_chart_freeze_calls"),
        "q_current": d.get("sample1_current_cayley_norm_upper"),
        "q_reference": d.get("authoritative_V18B_first_witness_current_q_reference"),
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
