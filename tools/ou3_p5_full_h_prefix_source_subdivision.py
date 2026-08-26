#!/usr/bin/env python3
"""Source-complete subdivision wrapper for the OU-III P5 full-H prefix map.

The V3 full-matrix producer carries the shipping 18x18 covariance, recomputes
H/R/S/K/r/d_eff in the same prefix cell, applies the Joseph update and immediate
reset, and composes the deployed quaternion before returning to Cayley
coordinates. Its remaining numerical obstruction is therefore a cell-width
problem, not missing shipping algebra.

This stage keeps that V3 map unchanged and replaces its one broad tuner-source
box by a finite Cartesian cover in (tau, sigma_aw, R_S). The pseudo cadence is
recomputed from each tau child with the same source schedule, so every child
still carries a joint source cell. No child is hulled with another child.

The CLI supports one-cell evaluation for CI matrix jobs and a separate aggregate
mode. The aggregate is source complete only when every expected child artifact
is present and validates. Numerical nonclosure remains a valid fail-closed
result; this producer never promotes P5 or sets N_H_words by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_first_accel_subdivision as SUB
import ou3_p5_full_h_prefix_cells as V1
import ou3_p5_full_h_prefix_cells_v3 as V3
import ou3_source_reachable_matrix_p3 as P3CELL

DEFAULT_DOMAIN = V1.DEFAULT_DOMAIN
CHILD_SCHEMA = 4
AGGREGATE_SCHEMA = 5


def _source_children(source_pieces: int) -> list[dict]:
    if source_pieces < 1:
        raise ValueError("source_pieces must be positive")
    sched = P3CELL.source_schedule()
    src0 = V1._source_cell()
    taus = SUB._geom_split(src0["tau_s"], source_pieces)
    sigmas = SUB._geom_split(src0["sigma_aw_mps2"], source_pieces)
    rss = SUB._geom_split(src0["R_S_filter_std"], source_pieces)
    out: list[dict] = []
    for tau in taus:
        plo, phi = P3CELL.cadence_bounds(tau, sched)
        period = Interval.outward_bounds(float(plo), float(phi))
        for sigma in sigmas:
            for rs in rss:
                src = dict(src0)
                src["tau_s"] = tau
                src["sigma_aw_mps2"] = sigma
                src["R_S_filter_std"] = rs
                src["pseudo_period_s"] = period
                out.append(src)
    return out


def _serialize_source(src: dict) -> dict:
    return {
        "dt_s": float(src["dt_s"]),
        "tau_s": src["tau_s"].as_list(),
        "sigma_aw_mps2": src["sigma_aw_mps2"].as_list(),
        "R_S_filter_std": src["R_S_filter_std"].as_list(),
        "pseudo_period_s": src["pseudo_period_s"].as_list(),
    }


def build_cell(
    domain_path: Path = DEFAULT_DOMAIN,
    *,
    source_pieces: int = 2,
    source_cell_index: int = 0,
) -> dict:
    children = _source_children(source_pieces)
    if not 0 <= source_cell_index < len(children):
        raise IndexError(
            f"source_cell_index {source_cell_index} outside [0,{len(children)})"
        )
    src = children[source_cell_index]
    original_source_cell = V1._source_cell
    V1._source_cell = lambda: dict(src)
    try:
        core = V3.build(Path(domain_path).resolve())
        core_failures = V3.validate(core)
    finally:
        V1._source_cell = original_source_cell

    closed = bool(
        not core_failures
        and core.get("P5_FULL_H_PREFIX_MATRIX_CERTIFICATE") == "PASS"
        and core.get("complete_q_le_8_prefix_family_closed") is True
    )
    out = dict(core)
    out.update({
        "schema": CHILD_SCHEMA,
        "qualification": "OU3_P5_FULL_H_PREFIX_SOURCE_SUBDIVISION_CHILD",
        "source_partition_method": "GEOMETRIC_CARTESIAN_TAU_SIGMA_AW_RS_WITH_TAU_DEPENDENT_CADENCE",
        "source_partition_is_finite_cover": True,
        "source_partition_cells_hulled_together": False,
        "source_partition_pieces_per_axis": int(source_pieces),
        "source_partition_total_cells": len(children),
        "source_partition_cell_index": int(source_cell_index),
        "selected_source_cell": _serialize_source(src),
        "v3_full_matrix_schema": V3.SCHEMA,
        "v3_validation_failures": list(core_failures),
        "same_full_18x18_Joseph_reset_quaternion_map_as_v3": True,
        "filter_changed": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_FULL_H_PREFIX_SOURCE_CELL_CERTIFICATE": "PASS" if closed else "NOT_ESTABLISHED",
        "next_obligation": (
            "AGGREGATE_ALL_SOURCE_PARTITION_CELLS_AND_CHECK_COMPLETE_Q8_CLOSURE"
            if closed else
            "REFINE_THIS_SOURCE_CELL_ALONG_VECTOR_DIRECTION_OR_BRANCH_PHASE_WITH_THE_SAME_FULL_MATRIX_CALCULUS"
        ),
    })
    return out


def validate_cell(d: dict) -> list[str]:
    failures = list(d.get("v3_validation_failures", []))
    if d.get("schema") != CHILD_SCHEMA:
        failures.append("child schema mismatch")
    if d.get("qualification") != "OU3_P5_FULL_H_PREFIX_SOURCE_SUBDIVISION_CHILD":
        failures.append("child qualification mismatch")
    for key in (
        "source_partition_is_finite_cover",
        "same_full_18x18_Joseph_reset_quaternion_map_as_v3",
        "full_18x18_covariance_propagated",
        "H_R_S_K_r_d_eff_recomputed_in_same_prefix_cell",
        "shipping_Joseph_update_used",
        "immediate_left_error_reset_congruence_used",
        "deployed_quaternion_composed_before_result_cayley",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "source_partition_cells_hulled_together",
        "filter_changed",
        "whole_word_promoted_here",
        "N_H_words_set_here",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    pieces = int(d.get("source_partition_pieces_per_axis", 0))
    total = int(d.get("source_partition_total_cells", 0))
    idx = int(d.get("source_partition_cell_index", -1))
    if pieces < 1 or total != pieces ** 3:
        failures.append("source partition cardinality mismatch")
    if not 0 <= idx < total:
        failures.append("source partition index invalid")
    src = d.get("selected_source_cell", {})
    for key in ("tau_s", "sigma_aw_mps2", "R_S_filter_std", "pseudo_period_s"):
        x = src.get(key)
        if not (isinstance(x, list) and len(x) == 2 and float(x[0]) <= float(x[1])):
            failures.append(f"invalid selected source interval {key}")
    closed = d.get("P5_FULL_H_PREFIX_SOURCE_CELL_CERTIFICATE") == "PASS"
    if closed:
        if d.get("complete_q_le_8_prefix_family_closed") is not True:
            failures.append("child PASS without complete q<=8 prefix closure")
        if d.get("first_failure") is not None:
            failures.append("child PASS retains first failure")
    else:
        if d.get("complete_q_le_8_prefix_family_closed") is True and not failures:
            failures.append("closed V3 child not promoted to child PASS")
        if d.get("complete_q_le_8_prefix_family_closed") is False and d.get("first_failure") is None:
            failures.append("nonclosed child lacks first failure witness")
    return list(dict.fromkeys(failures))


def aggregate(children: list[dict], *, source_pieces: int) -> dict:
    expected = source_pieces ** 3
    by_index: dict[int, dict] = {}
    failures: list[str] = []
    for child in children:
        idx = int(child.get("source_partition_cell_index", -1))
        if idx in by_index:
            failures.append(f"duplicate source cell {idx}")
            continue
        by_index[idx] = child
        for f in validate_cell(child):
            failures.append(f"cell {idx}: {f}")
    missing = sorted(set(range(expected)) - set(by_index))
    extra = sorted(set(by_index) - set(range(expected)))
    if missing:
        failures.append(f"missing source cells: {missing}")
    if extra:
        failures.append(f"unexpected source cells: {extra}")

    ordered = [by_index[i] for i in range(expected) if i in by_index]
    unclosed = [
        c for c in ordered
        if c.get("P5_FULL_H_PREFIX_SOURCE_CELL_CERTIFICATE") != "PASS"
    ]
    all_closed = bool(not failures and len(ordered) == expected and not unclosed)
    finite_q = [
        float(c.get("max_reached_cayley_norm_upper", math.inf))
        for c in ordered
        if math.isfinite(float(c.get("max_reached_cayley_norm_upper", math.inf)))
    ]
    max_q = max(finite_q, default=math.inf)
    first_unclosed = None
    if unclosed:
        c = unclosed[0]
        first_unclosed = {
            "source_partition_cell_index": c["source_partition_cell_index"],
            "selected_source_cell": c["selected_source_cell"],
            "max_reached_cayley_norm_upper": c.get("max_reached_cayley_norm_upper"),
            "first_failure": c.get("first_failure"),
            "next_obligation": c.get("next_obligation"),
        }

    return {
        "schema": AGGREGATE_SCHEMA,
        "qualification": "OU3_P5_FULL_H_PREFIX_SOURCE_SUBDIVISION_AGGREGATE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "source_partition_method": "GEOMETRIC_CARTESIAN_TAU_SIGMA_AW_RS_WITH_TAU_DEPENDENT_CADENCE",
        "source_partition_pieces_per_axis": int(source_pieces),
        "expected_source_partition_cells": expected,
        "evaluated_source_partition_cells": len(ordered),
        "source_partition_cell_indices": [int(c["source_partition_cell_index"]) for c in ordered],
        "complete_source_partition_covered": len(ordered) == expected and not missing and not extra,
        "source_partition_cells_hulled_together": False,
        "same_full_18x18_Joseph_reset_quaternion_map_in_every_child": True,
        "all_partition_cells_numerically_closed": all_closed,
        "complete_q_le_8_prefix_family_closed_over_source_partition": all_closed,
        "max_reached_cayley_norm_upper_over_partition": max_q,
        "first_unclosed_source_cell": first_unclosed,
        "child_statuses": [
            {
                "source_partition_cell_index": int(c["source_partition_cell_index"]),
                "status": c["P5_FULL_H_PREFIX_SOURCE_CELL_CERTIFICATE"],
                "max_reached_cayley_norm_upper": c.get("max_reached_cayley_norm_upper"),
                "first_failure": c.get("first_failure"),
            }
            for c in ordered
        ],
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_FULL_H_PREFIX_SOURCE_SUBDIVISION_CERTIFICATE": "PASS" if all_closed else "NOT_ESTABLISHED",
        "next_obligation": (
            "COMPUTE_CERTIFIED_INNER_FUNNEL_OVERLAP_AND_SET_N_H_WORDS"
            if all_closed else
            "REFINE_FIRST_UNCLOSED_SOURCE_CELL_ALONG_VECTOR_DIRECTION_OR_BRANCH_PHASE_WITH_THE_SAME_FULL_MATRIX_CALCULUS"
        ),
        "failures": failures,
    }


def validate_aggregate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != AGGREGATE_SCHEMA:
        failures.append("aggregate schema mismatch")
    if d.get("qualification") != "OU3_P5_FULL_H_PREFIX_SOURCE_SUBDIVISION_AGGREGATE":
        failures.append("aggregate qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "complete_source_partition_covered",
        "same_full_18x18_Joseph_reset_quaternion_map_in_every_child",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "source_replay_used",
        "filter_changed",
        "source_partition_cells_hulled_together",
        "whole_word_promoted_here",
        "N_H_words_set_here",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    expected = int(d.get("expected_source_partition_cells", 0))
    evaluated = int(d.get("evaluated_source_partition_cells", -1))
    if expected < 1 or evaluated != expected:
        failures.append("aggregate does not cover every expected source cell")
    status = d.get("P5_FULL_H_PREFIX_SOURCE_SUBDIVISION_CERTIFICATE")
    if status == "PASS":
        if d.get("all_partition_cells_numerically_closed") is not True:
            failures.append("aggregate PASS without all child closures")
        if d.get("complete_q_le_8_prefix_family_closed_over_source_partition") is not True:
            failures.append("aggregate PASS without q<=8 closure")
        if d.get("first_unclosed_source_cell") is not None:
            failures.append("aggregate PASS retains unclosed witness")
    elif status == "NOT_ESTABLISHED":
        if d.get("all_partition_cells_numerically_closed") is True and not failures:
            failures.append("all children closed but aggregate not promoted")
        if not failures and d.get("first_unclosed_source_cell") is None:
            failures.append("numerical nonclosure lacks source-cell witness")
    else:
        failures.append("invalid aggregate status")
    return list(dict.fromkeys(failures))


def _load_children(directory: Path) -> list[dict]:
    paths = sorted(directory.glob("p5-source-cell-*.json"))
    return [json.loads(p.read_text(encoding="utf-8")) for p in paths]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--source-cell-index", type=int)
    group.add_argument("--aggregate-dir", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    if args.source_cell_index is not None:
        out = build_cell(
            args.domain.resolve(),
            source_pieces=args.source_pieces,
            source_cell_index=args.source_cell_index,
        )
        vf = validate_cell(out)
        out["validation_pass"] = not vf
        out["validation_failures"] = vf
        summary = {
            "mode": "source-cell",
            "cell": out["source_partition_cell_index"],
            "status": out["P5_FULL_H_PREFIX_SOURCE_CELL_CERTIFICATE"],
            "q8_closed": out["complete_q_le_8_prefix_family_closed"],
            "max_q": out["max_reached_cayley_norm_upper"],
            "first_failure": out["first_failure"],
            "next": out["next_obligation"],
            "validation_failures": vf,
        }
    else:
        out = aggregate(_load_children(args.aggregate_dir), source_pieces=args.source_pieces)
        vf = validate_aggregate(out)
        out["validation_pass"] = not vf
        out["validation_failures"] = vf
        summary = {
            "mode": "aggregate",
            "status": out["P5_FULL_H_PREFIX_SOURCE_SUBDIVISION_CERTIFICATE"],
            "evaluated": out["evaluated_source_partition_cells"],
            "expected": out["expected_source_partition_cells"],
            "q8_closed": out["complete_q_le_8_prefix_family_closed_over_source_partition"],
            "max_q": out["max_reached_cayley_norm_upper_over_partition"],
            "first_unclosed": out["first_unclosed_source_cell"],
            "next": out["next_obligation"],
            "validation_failures": vf,
        }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
