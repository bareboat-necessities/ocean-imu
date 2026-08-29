#!/usr/bin/env python3
"""V52: lift V51's exact monotone gain enclosure over the full source-cell-0 cover.

V51 proves the exact monotone corner enclosure of the eight first-accelerometer
block rationals at the authoritative V41/V45 first survivor, and shows the
refined V10 correction magnitude closes that witness on the SO(3) triangle.  A
single cell is not a cover, so V52 installs the same refinement into the shared
`ou3_p5_sample1_structured_full_gain_v8._first_block_quantities` hook and
re-evaluates the complete V10 source-cell-0 cover.

The hook is the only place the block parameters enter V8's row loop; installing
V51's exact path there refines every cell by construction and is removed again
in a `finally`, so no other producer sees a changed backend.

V52 reports three things.

1. **Containment.**  Every refined cell correction is required to be at or below
   its parent, cell by cell over the whole cover, not only in aggregate.  A
   single cell that widened fails the producer closed.

2. **Cover narrowing.**  The parent and refined maxima of the sample-1 residual,
   the two gain blocks, and the V10 directional correction, together with the
   best and worst per-cell narrowing ratio.

3. **Witness closure carried from V51.**  The authoritative survivor's refined
   correction composed with V41's archived sample-0 chart.

What V52 does *not* do: it does not re-run V41's signed/chart q<8 composition
with the refined block, so it does not by itself establish `q<8` over the cover.
The cover-wide correction maximum is a magnitude bound over every cell, not the
per-cell composed q that V41 evaluates, and the two are not interchangeable.
Re-running that composition is the next obligation.  Nothing here composes q<8,
promotes sample 1 or P5, or sets `N_H_words`, and no filter setting, source
domain, six-radian correction limit, `q<8` target, or source language changes.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_structured_full_gain_v8 as V8
import ou3_p5_sample1_structured_full_gain_v10 as V10
import ou3_p5_sample1_exact_monotone_source_gain_v51 as V51

DEFAULT_DOMAIN = V51.DEFAULT_DOMAIN
SCHEMA = 5200
Q_TARGET = V51.Q_TARGET
WITNESS = V51.WITNESS
CORRECTION_KEY = "combined_directional_correction_norm_upper_rad"

#: Archived cover maxima of the unrefined V10 source-cell-0 evaluation.  V52
#: reproduces them from source and refuses to continue otherwise.
PARENT_COVER = {
    "evaluated_joint_cells": 12816,
    "max_full_residual_norm_upper_mps2": 46.01460061569009,
    "max_Ktheta_perpendicular_block_upper": 1.3986770467171177,
    "max_Ktheta_parallel_block_upper": 0.19289137244367335,
    "max_combined_directional_correction_norm_upper_rad": 7.016940736774492,
}

_SUMMARY_KEYS = (
    "max_full_residual_norm_upper_mps2",
    "max_Ktheta_perpendicular_block_upper",
    "max_Ktheta_parallel_block_upper",
    "max_combined_directional_correction_norm_upper_rad",
)


def _exact_block(*, t, p, r, g: float) -> dict:
    """V51's exact monotone corner enclosure, in V8's hook signature."""
    return V51._first_block(t=t, p=p, r=r, g=g, exact=True)


def _cell_key(row: dict) -> tuple:
    return (int(row["p_cell"]), int(row["tangent_residual_cell"]),
            int(row["axial_residual_cell"]))


def _corrections(cover: dict) -> dict:
    out = {}
    for row in cover.get("rows", []):
        out[_cell_key(row)] = float(row[CORRECTION_KEY])
    return out


def _cover_summary(cover: dict) -> dict:
    d = {k: float(cover[k]) for k in _SUMMARY_KEYS}
    d["evaluated_joint_cells"] = int(cover["evaluated_joint_cells"])
    d["unclosed_joint_cells"] = int(cover["unclosed_joint_cells"])
    d["status"] = cover.get("P5_SAMPLE1_COMBINED_PERPENDICULAR_RESIDUAL_V10")
    return d


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24) -> dict:
    """Evaluate the V10 source-cell-0 cover with and without the refinement."""
    path = Path(domain_path).resolve()
    failures: list[str] = []
    parent_summary = refined_summary = comparison = None
    witness = None
    hook_restored = False

    root = V8._first_block_quantities
    try:
        parent = V10.build(path, source_pieces=source_pieces,
                           source_cell_index=source_cell_index,
                           p_pieces=p_pieces, tangent_pieces=tangent_pieces,
                           axial_pieces=axial_pieces)
        failures += [f"parent V10: {x}" for x in V10.validate(parent)]

        V8._first_block_quantities = _exact_block
        try:
            refined = V10.build(path, source_pieces=source_pieces,
                                source_cell_index=source_cell_index,
                                p_pieces=p_pieces, tangent_pieces=tangent_pieces,
                                axial_pieces=axial_pieces)
        finally:
            V8._first_block_quantities = root
        hook_restored = V8._first_block_quantities is root
        if not hook_restored:
            failures.append("V52 temporary V8 block hook was not restored")
        failures += [f"refined V10: {x}" for x in V10.validate(refined)]

        parent_summary = _cover_summary(parent)
        refined_summary = _cover_summary(refined)
        for key, want in PARENT_COVER.items():
            got = parent_summary.get(key)
            if got is None or float(got) != float(want):
                failures.append(f"parent cover {key} is {got}, archived {want}")

        pc = _corrections(parent)
        rc = _corrections(refined)
        if set(pc) != set(rc):
            failures.append("refined cover does not cover the same cells")
        shared = sorted(set(pc) & set(rc))
        widened = [k for k in shared if rc[k] > pc[k]]
        narrowed = [k for k in shared if rc[k] < pc[k]]
        if widened:
            failures.append(
                f"{len(widened)} cells widened, first {list(widened[0])}")
        ratios = [pc[k] / rc[k] for k in shared if rc[k] > 0.0]
        for key in _SUMMARY_KEYS:
            if refined_summary[key] > parent_summary[key]:
                failures.append(f"refined cover {key} exceeded its parent")
        if int(refined_summary["unclosed_joint_cells"]) != 0:
            failures.append("refined cover left an unclosed joint cell")

        comparison = {
            "cells_evaluated": len(shared),
            "cells_narrowed": len(narrowed),
            "cells_widened": len(widened),
            "every_cell_inside_parent": not widened,
            "narrowing_ratio_min": min(ratios) if ratios else math.inf,
            "narrowing_ratio_max": max(ratios) if ratios else math.inf,
            "parent_cover": parent_summary,
            "refined_cover": refined_summary,
            "cover_correction_reduction_rad": V51.FULL.up(
                parent_summary["max_combined_directional_correction_norm_upper_rad"]
                - refined_summary["max_combined_directional_correction_norm_upper_rad"]),
        }

        w = V51.build(path, source_pieces=source_pieces,
                      source_cell_index=source_cell_index, p_pieces=p_pieces,
                      tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
        failures += [f"V51: {x}" for x in V51.validate(w)]
        witness = {
            "V51_status": w.get("P5_SAMPLE1_EXACT_MONOTONE_SOURCE_GAIN_V51"),
            "authoritative_witness_closed": bool(
                w.get("authoritative_witness_closed")),
            "comparison": w.get("witness_comparison"),
        }
        wc = (refined.get("rows") or [])
        wrow = next((r for r in wc if _cell_key(r) == WITNESS), None)
        if wrow is None:
            failures.append("refined cover does not contain the V41 witness")
        elif float(wrow[CORRECTION_KEY]) != float(
                (w.get("exact_monotone_witness_chain") or {}).get(
                    CORRECTION_KEY, math.nan)):
            failures.append(
                "refined cover and V51 disagree on the witness correction")
    except Exception as exc:
        V8._first_block_quantities = root
        failures.append(f"V52 exact monotone cover lift: {exc}")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_EXACT_MONOTONE_COVER_LIFT_V52",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "exact_monotone_corner_enclosure_used": True,
        "parent_enclosure_retained_as_intersection": True,
        "temporary_V8_block_hook_restored": hook_restored,
        "complete_source_cell0_cover_evaluated": bool(comparison),
        "V41_signed_chart_q8_composition_rerun_here": False,
        "V41_first_survivor_row": list(WITNESS),
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_established_here": False,
        "cover_comparison": comparison,
        "authoritative_witness": witness,
        "P5_SAMPLE1_EXACT_MONOTONE_COVER_LIFT_V52": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "RERUN_V41_SIGNED_CHART_Q8_COMPOSITION_WITH_EXACT_MONOTONE_BLOCK"
            if not failures else
            "REPAIR_V52_EXACT_MONOTONE_COVER_LIFT"),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    """Validate the fail-closed V52 proof-artifact contract."""
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_EXACT_MONOTONE_COVER_LIFT_V52":
        f.append("qualification mismatch")
    for k in ("source_generated_not_trajectory_fit",
              "exact_monotone_corner_enclosure_used",
              "parent_enclosure_retained_as_intersection",
              "temporary_V8_block_hook_restored",
              "complete_source_cell0_cover_evaluated"):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed",
              "deployed_correction_limit_increased",
              "V41_signed_chart_q8_composition_rerun_here",
              "q8_composed_here", "q8_word_promoted_here",
              "whole_word_promoted_here", "N_H_words_set_here",
              "P5_established_here"):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if tuple(d.get("V41_first_survivor_row", ())) != tuple(WITNESS):
        f.append("V41 witness changed")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")

    cmp_ = d.get("cover_comparison") or {}
    if cmp_.get("every_cell_inside_parent") is not True:
        f.append("a refined cell escaped its parent")
    if int(cmp_.get("cells_widened", 1)) != 0:
        f.append("a refined cell widened")
    if int(cmp_.get("cells_evaluated", 0)) != PARENT_COVER["evaluated_joint_cells"]:
        f.append("cover cell count changed")
    parent_cover = cmp_.get("parent_cover") or {}
    refined_cover = cmp_.get("refined_cover") or {}
    for key in _SUMMARY_KEYS:
        pv = float(parent_cover.get(key, -math.inf))
        rv = float(refined_cover.get(key, math.inf))
        if not (math.isfinite(pv) and math.isfinite(rv) and 0.0 <= rv <= pv):
            f.append(f"invalid cover comparison for {key}")
    if int(refined_cover.get("unclosed_joint_cells", 1)) != 0:
        f.append("refined cover left an unclosed joint cell")

    w = d.get("authoritative_witness") or {}
    if w.get("V51_status") != "PASS":
        f.append("V51 witness stage did not pass")
    if w.get("authoritative_witness_closed") is not True:
        f.append("authoritative witness is not closed")
    if d.get("P5_SAMPLE1_EXACT_MONOTONE_COVER_LIFT_V52") not in (
            "PASS", "NOT_ESTABLISHED"):
        f.append("invalid V52 status")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=24)
    ap.add_argument("--tangent-pieces", type=int, default=24)
    ap.add_argument("--axial-pieces", type=int, default=24)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(x.domain, source_pieces=x.source_pieces,
              source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
              tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_EXACT_MONOTONE_COVER_LIFT_V52"],
        "cover": d.get("cover_comparison"),
        "witness": d.get("authoritative_witness"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
