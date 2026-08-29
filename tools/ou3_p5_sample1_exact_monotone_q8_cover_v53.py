#!/usr/bin/env python3
"""V53: re-run V41's signed-chart q<8 cover with the exact monotone block.

V51 proves the exact monotone corner enclosure of the eight first-accelerometer
block rationals and closes the authoritative V41/V45 first survivor on the SO(3)
triangle.  V52 lifts the same enclosure over the complete V10 source-cell-0
cover and shows every one of the 12816 cells narrows and none widens.  Neither
re-runs V41's signed-chart composition, so neither says what the refinement does
to the `q<8` cover itself.  That is V53's only job.

V53 evaluates `ou3_p5_sample1_v40_full_source_cell0_q8_lift_v41.build` twice in
one process: once with the shipping backend, and once with V51's exact monotone
enclosure installed in the shared
`ou3_p5_sample1_structured_full_gain_v8._first_block_quantities` hook.  Both runs
are validated by V41's own contract, the hook is removed in a `finally`, and the
parent run is required to reproduce V45's archived `V41_Q_POST` first-survivor q
so the comparison is anchored to the same cover the rest of the chain cites.

Running the parent here rather than quoting an archived count is deliberate.
The refinement claim is a difference between two covers, and this producer
regenerates both sides rather than trusting a number recorded elsewhere.

The refinement is required to be an improvement in every direction that is
comparable: the same cell count, no increase in open cells, and no increase in
the worst composed q.  Any regression fails the producer closed.

V53 does not compose `q<8`, promote sample 1 or P5, or set `N_H_words`, and it
changes no filter setting, source domain, six-radian correction limit, `q<8`
target, or source language.  It reports where the cover stands after the
refinement and names the first cell that still does not close.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_structured_full_gain_v8 as V8
import ou3_p5_sample1_exact_monotone_source_gain_v51 as V51
import ou3_p5_sample1_exact_monotone_cover_lift_v52 as V52
import ou3_p5_sample1_v40_full_source_cell0_q8_lift_v41 as V41
import ou3_p5_sample1_v41_authoritative_split_signed_v45 as V45

DEFAULT_DOMAIN = V51.DEFAULT_DOMAIN
SCHEMA = 5300
Q_TARGET = V51.Q_TARGET
WITNESS = V51.WITNESS
V41_Q_POST = V45.V41_Q_POST

_CELL_KEYS = ("p_cell", "tangent_residual_cell", "axial_residual_cell")


def _cell(row: dict | None) -> dict | None:
    """Return the comparable fields of a V41 q8 cell record."""
    if not row:
        return None
    out = {k: int(row[k]) for k in _CELL_KEYS if k in row}
    for k in ("post_sample1_cayley_norm_upper",
              "sample1_current_cayley_norm_upper",
              "correction_radial_lower_rad", "correction_radial_upper_rad"):
        if k in row:
            out[k] = float(row[k])
    return out


def _summary(cover: dict) -> dict:
    cells = int(cover["evaluated_signed_cayley_cells"])
    unclosed = int(cover["unclosed_q8_cells"])
    return {
        "status": cover.get("P5_SAMPLE1_V40_FULL_SOURCE_CELL0_Q8_LIFT_V41"),
        "evaluated_signed_cayley_cells": cells,
        "unclosed_q8_cells": unclosed,
        "closed_q8_cells": cells - unclosed,
        "geodesic_bound_newly_closed_cells": int(
            cover.get("geodesic_bound_newly_closed_cells", 0)),
        "current_yz_support_newly_closed_cells": int(
            cover.get("current_yz_support_newly_closed_cells", 0)),
        "max_sample1_q_after_product_tightening": float(
            cover.get("max_sample1_q_after_product_tightening", math.nan)),
        "first_unclosed_q8_cell": _cell(cover.get("first_unclosed_q8_cell")),
        "worst_q8_cell": _cell(cover.get("worst_q8_cell")),
    }


def _worst_q(summary: dict) -> float:
    worst = summary.get("worst_q8_cell") or {}
    return float(worst.get("post_sample1_cayley_norm_upper", math.inf))


def _first_q(summary: dict) -> float:
    first = summary.get("first_unclosed_q8_cell") or {}
    return float(first.get("post_sample1_cayley_norm_upper", math.inf))


def _first_key(summary: dict) -> tuple | None:
    first = summary.get("first_unclosed_q8_cell") or {}
    if not all(k in first for k in _CELL_KEYS):
        return None
    return tuple(int(first[k]) for k in _CELL_KEYS)


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    """Run V41's q<8 cover with and without the exact monotone block."""
    path = Path(domain_path).resolve()
    failures: list[str] = []
    parent_summary = refined_summary = comparison = None
    hook_restored = False
    kwargs = dict(source_pieces=source_pieces,
                  source_cell_index=source_cell_index, p_pieces=p_pieces,
                  tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
                  residual_x_pieces=residual_x_pieces,
                  parallel_pieces=parallel_pieces)

    root = V8._first_block_quantities
    try:
        parent = V41.build(path, **kwargs)
        failures += [f"parent V41: {x}" for x in V41.validate(parent)]

        V8._first_block_quantities = V52._exact_block
        try:
            refined = V41.build(path, **kwargs)
        finally:
            V8._first_block_quantities = root
        hook_restored = V8._first_block_quantities is root
        if not hook_restored:
            failures.append("V53 temporary V8 block hook was not restored")
        failures += [f"refined V41: {x}" for x in V41.validate(refined)]

        parent_summary = _summary(parent)
        refined_summary = _summary(refined)

        if _first_key(parent_summary) != tuple(WITNESS):
            failures.append(
                f"parent first survivor is {_first_key(parent_summary)}, "
                f"archived {tuple(WITNESS)}")
        if _first_q(parent_summary) != V41_Q_POST:
            failures.append(
                f"parent first survivor q is {_first_q(parent_summary)}, "
                f"archived {V41_Q_POST}")

        pc = parent_summary["evaluated_signed_cayley_cells"]
        rc = refined_summary["evaluated_signed_cayley_cells"]
        if pc != rc:
            failures.append(f"cover size changed from {pc} to {rc}")
        if refined_summary["unclosed_q8_cells"] > parent_summary["unclosed_q8_cells"]:
            failures.append("refined cover left more open cells than its parent")
        if _worst_q(refined_summary) > _worst_q(parent_summary):
            failures.append("refined cover worsened the worst composed q")

        archived_closed = _first_key(refined_summary) != tuple(WITNESS)
        newly_closed = (parent_summary["unclosed_q8_cells"]
                        - refined_summary["unclosed_q8_cells"])
        comparison = {
            "parent": parent_summary,
            "refined": refined_summary,
            "additional_cells_closed": newly_closed,
            "open_cell_reduction_fraction": (
                0.0 if parent_summary["unclosed_q8_cells"] == 0
                else newly_closed / parent_summary["unclosed_q8_cells"]),
            "archived_first_survivor_closed_by_refinement": archived_closed,
            "parent_worst_q_upper": _worst_q(parent_summary),
            "refined_worst_q_upper": _worst_q(refined_summary),
            "q_target": Q_TARGET,
            "cover_fully_closed": refined_summary["unclosed_q8_cells"] == 0,
        }
        if not archived_closed:
            failures.append(
                "refined cover did not close the archived first survivor")
    except Exception as exc:
        V8._first_block_quantities = root
        failures.append(f"V53 exact monotone q8 cover: {exc}")

    fully_closed = bool(comparison and comparison["cover_fully_closed"])
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_EXACT_MONOTONE_Q8_COVER_V53",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "exact_monotone_corner_enclosure_used": True,
        "both_covers_regenerated_here": bool(comparison),
        "archived_parent_cover_quoted_instead_of_regenerated": False,
        "temporary_V8_block_hook_restored": hook_restored,
        "V41_first_survivor_row": list(WITNESS),
        "archived_V41_post_sample1_q_reference": V41_Q_POST,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_established_here": False,
        "cover_comparison": comparison,
        "source_cell0_q8_cover_fully_closed": fully_closed,
        "P5_SAMPLE1_EXACT_MONOTONE_Q8_COVER_V53": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "COMPOSE_SIGNED_CAYLEY_Q8_OVER_REMAINING_SOURCE_CELLS"
            if fully_closed and not failures else
            "REFINE_NOMINAL_GEOMETRY_AT_THE_REMAINING_FIRST_UNCLOSED_Q8_CELL"),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    """Validate the fail-closed V53 proof-artifact contract."""
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_EXACT_MONOTONE_Q8_COVER_V53":
        f.append("qualification mismatch")
    for k in ("source_generated_not_trajectory_fit",
              "exact_monotone_corner_enclosure_used",
              "both_covers_regenerated_here",
              "temporary_V8_block_hook_restored"):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed",
              "archived_parent_cover_quoted_instead_of_regenerated",
              "deployed_correction_limit_increased", "q8_composed_here",
              "q8_word_promoted_here", "whole_word_promoted_here",
              "N_H_words_set_here", "P5_established_here"):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if tuple(d.get("V41_first_survivor_row", ())) != tuple(WITNESS):
        f.append("V41 witness changed")
    if float(d.get("archived_V41_post_sample1_q_reference", 0.0)) != V41_Q_POST:
        f.append("archived V41 q reference changed")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")

    cmp_ = d.get("cover_comparison") or {}
    parent = cmp_.get("parent") or {}
    refined = cmp_.get("refined") or {}
    if not parent or not refined:
        f.append("cover summaries missing")
    else:
        if parent["evaluated_signed_cayley_cells"] != refined[
                "evaluated_signed_cayley_cells"]:
            f.append("cover size changed between the two runs")
        if int(refined["unclosed_q8_cells"]) > int(parent["unclosed_q8_cells"]):
            f.append("refined cover left more open cells")
        if int(cmp_.get("additional_cells_closed", -1)) != (
                int(parent["unclosed_q8_cells"]) - int(refined["unclosed_q8_cells"])):
            f.append("inconsistent additional closed-cell count")
        if float(cmp_.get("refined_worst_q_upper", math.inf)) > float(
                cmp_.get("parent_worst_q_upper", -math.inf)):
            f.append("refined cover worsened the worst composed q")
    if cmp_.get("archived_first_survivor_closed_by_refinement") is not True:
        f.append("archived first survivor is not closed")
    if bool(d.get("source_cell0_q8_cover_fully_closed")) != bool(
            cmp_.get("cover_fully_closed")):
        f.append("inconsistent cover closure verdict")
    if d.get("P5_SAMPLE1_EXACT_MONOTONE_Q8_COVER_V53") not in (
            "PASS", "NOT_ESTABLISHED"):
        f.append("invalid V53 status")
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
    d = build(x.domain, source_pieces=x.source_pieces,
              source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
              tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces,
              residual_x_pieces=x.residual_x_pieces,
              parallel_pieces=x.parallel_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_EXACT_MONOTONE_Q8_COVER_V53"],
        "fully_closed": d.get("source_cell0_q8_cover_fully_closed"),
        "comparison": d.get("cover_comparison"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
