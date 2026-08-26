#!/usr/bin/env python3
"""Bind V13 signed radial cells to the closed V12D sample-1 family.

V12D closes all 24^3 source cells after replacing the spurious first-PSD axial
noise-floor amplification by the exact gravity-tangent innovation algebra.  The
existing V13 signed-cell construction is otherwise the right next step: it
reconstructs the signed one-plus-two attitude correction components and retains
an additive norm ball for all omitted PSD/S effects.

This module adapts V12D to V13's established input contract without duplicating
V13's signed geometry.  The only compatibility aliases are

* V12D's attitude-only gain perturbation -> V13's perturbation field;
* V12D's certified correction upper -> V13's radial parent upper;
* V12D PASS -> the historical V12 prerequisite status key consumed by V13.

No numerical quantity is reduced by the adapter.  No source branch, filter
parameter, six-radian deployed proof limit, or theorem gate is changed.  The
stage still does not compose the current Cayley state; it only proves that every
sample-1 correction belongs to a signed/radial cell accepted by the existing
winding-aware deployed-quaternion primitive.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p5_sample1_signed_radial_subcells_v13 as V13
import ou3_p5_sample1_structured_full_gain_v12d as V12D

DEFAULT_DOMAIN = V12D.DEFAULT_DOMAIN
SCHEMA = 1301


class _V12DAdapter:
    DEFAULT_DOMAIN = V12D.DEFAULT_DOMAIN
    FULL = V12D.FULL
    V11 = V12D.V11

    @staticmethod
    def build(domain_path=DEFAULT_DOMAIN, **kwargs):
        d = V12D.build(domain_path, **kwargs)
        out = dict(d)
        out["P5_SAMPLE1_PSD_S_ACTUAL_INNOVATION_V12"] = (
            "PASS" if d.get("P5_SAMPLE1_FIRST_PSD_TANGENT_REFINEMENT_V12D") == "PASS"
            else "NOT_ESTABLISHED"
        )
        rows = []
        for row0 in d.get("rows", []):
            row = dict(row0)
            row["sample1_gain_operator_perturbation_upper"] = float(
                row0["sample1_attitude_gain_operator_perturbation_upper"])
            row["V12_correction_norm_upper_rad"] = float(
                row0["V12C_correction_norm_upper_rad"])
            rows.append(row)
        out["rows"] = rows
        return out

    @staticmethod
    def validate(d):
        return V12D.validate(d)


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 12, parallel_pieces: int = 12) -> dict:
    original = V13.V12
    V13.V12 = _V12DAdapter
    try:
        core = V13.build(
            domain_path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces)
    finally:
        V13.V12 = original

    inherited = V13.validate(core)
    status = core.get("P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13")
    out = dict(core)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13D",
        "V12D_tangent_channel_prerequisite_used": True,
        "V12D_prerequisite_passed": core.get("V12_prerequisite_passed") is True,
        "V12D_attitude_gain_perturbation_used": True,
        "V12D_radial_parent_upper_used": True,
        "P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13D": (
            "PASS" if status == "PASS" and not inherited else "NOT_ESTABLISHED"
        ),
        "next_obligation": (
            "COMPOSE_SIGNED_SAMPLE1_CORRECTION_WITH_CURRENT_CAYLEY_AND_REQUIRE_Q_LT_8"
            if status == "PASS" and not inherited
            else "REFINE_SIGNED_RX_PARALLEL_SUBDIVISION"
        ),
        "failures": list(dict.fromkeys(inherited)),
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit", "V12D_tangent_channel_prerequisite_used",
        "V12D_attitude_gain_perturbation_used", "V12D_radial_parent_upper_used",
        "radial_lower_bound_required_above_6_rad",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "signed_cayley_q8_composed_here", "complete_sample1_branch_closed_here",
        "q8_word_promoted_here", "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if d.get("V12D_prerequisite_passed") is not True:
        f.append("V12D prerequisite did not pass")
    st = d.get("P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13D")
    if st == "PASS":
        if int(d.get("unclosed_radial_subcells", -1)) != 0:
            f.append("V13D PASS retains unclosed radial cells")
        if float(d.get("max_radial_upper", 1e300)) > 9.0:
            f.append("V13D PASS exceeds validated winding range")
        if int(d.get("above_6rad_subcells", 0)) > 0 and not float(
                d.get("minimum_radial_lower_above_6", 0.0)) > 0.0:
            f.append("V13D >6-rad cells lack positive radial lower bound")
    elif st == "NOT_ESTABLISHED":
        if d.get("first_unclosed_radial_subcell") is None and not f:
            f.append("V13D nonclosure lacks radial witness")
    else:
        f.append("invalid V13D status")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=24)
    ap.add_argument("--tangent-pieces", type=int, default=24)
    ap.add_argument("--axial-pieces", type=int, default=24)
    ap.add_argument("--residual-x-pieces", type=int, default=12)
    ap.add_argument("--parallel-pieces", type=int, default=12)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(x.domain, source_pieces=x.source_pieces,
              source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
              tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces,
              residual_x_pieces=x.residual_x_pieces,
              parallel_pieces=x.parallel_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13D"],
        "V12D_passed": d["V12D_prerequisite_passed"],
        "source_rows": d["evaluated_source_rows"],
        "signed_subcells": d["evaluated_signed_subcells"],
        "above_6": d["above_6rad_subcells"],
        "unclosed": d["unclosed_radial_subcells"],
        "max_radial_upper": d["max_radial_upper"],
        "min_radial_lower_above_6": d["minimum_radial_lower_above_6"],
        "first_unclosed": d["first_unclosed_radial_subcell"],
        "worst": d.get("worst_radial_subcell"),
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
