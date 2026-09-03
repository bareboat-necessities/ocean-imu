#!/usr/bin/env python3
"""Finite-cost matched-history source frontier for OU-III P3.

The 49-step gap-forgotten quotient is a valid source superset, but CI showed that
it reaches the single global adverse label and therefore recreates the rejected
cross-tau P3 comparison.  This producer restores the finite 635-sample budget
without returning to the expensive exact elapsed-sample dynamic program.

For each exact P2 V1 source transition s -> t, keep only the *minimum* certified
gap g in {13,...,26} for which t is a labelled successor of s.  This is exact
for reachability of the four path-max statistics used by the covariance upper:
transition legality depends only on (s,g,t), the source statistics depend only
on s, and replacing a larger supporting gap by the minimum supporting gap keeps
the same source sequence and path maxima while consuming no more samples.

At one current source, retain pairs

    (adverse path-max label, minimum elapsed samples).

A pair (qa,ca) safely dominates (qb,cb) when qa is componentwise at least as
adverse as qb and ca <= cb.  Both pairs have the same future transition set; qa
can only increase the covariance upper and enlarge the source set admitted by
the conservative lower, while ca leaves at least as much finite-word budget for
every continuation.  Dominated pairs can therefore be removed without losing a
worst matched lower/upper class.

When one minimum-gap edge crosses the 635-sample target, the terminal endpoint
source is folded into the label before global adverse Pareto reduction.  This
covers the canonical 0..25-sample terminal phase.  The emitted labels are then
mapped to their same-history covariance upper and same-tau covariance-lower
dominator source nodes.

This module computes no covariance margin and cannot promote P3/P4/P5.
"""
from __future__ import annotations

import argparse
import heapq
import json
import math
from pathlib import Path

import ou3_p2_correlation_path_memory as CORR
import ou3_p3_correlated_translation_covariance_upper as CUPPER
import ou3_p3_matched_history_label_frontier as UNION
import ou3_p3_p2_v1_history_frontier as HIST
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
MAX_FRONTIER_PER_SOURCE = 50000
MAX_PROCESSED_ENTRIES = 5000000

Label = tuple[int, int, int, int]


def min_gap_successors(rt: dict) -> list[tuple[tuple[int, int], ...]]:
    """Return exact source transitions with the cheapest certified gap per edge."""
    gaps = tuple(map(int, rt["gaps"]))
    labelled = rt["labelled_successors"]
    out: list[tuple[tuple[int, int], ...]] = []
    for s in range(len(rt["nodes"])):
        best: dict[int, int] = {}
        for gi, gap in enumerate(gaps):
            for t0 in labelled[s][gi]:
                t = int(t0)
                old = best.get(t)
                if old is None or gap < old:
                    best[t] = gap
        if not best:
            raise RuntimeError(f"P2 source {s} has no finite-clock successor")
        row = tuple(sorted((t, g) for t, g in best.items()))
        for t, g in row:
            gi = gaps.index(g)
            if t not in labelled[s][gi]:
                raise RuntimeError("minimum-gap edge is not an exact labelled P2 transition")
        out.append(row)
    return out


def cost_dominates(label_a: Label, cost_a: int, label_b: Label, cost_b: int) -> bool:
    """Whether A is no less adverse and no more expensive than B."""
    return int(cost_a) <= int(cost_b) and HIST.dominates(label_a, label_b)


def insert_cost_frontier(front: dict[Label, int], label: Label, cost: int) -> bool:
    """Insert one (label,min-cost) state, removing only safely dominated states."""
    q = tuple(map(int, label))
    c = int(cost)
    old_same = front.get(q)
    if old_same is not None and old_same <= c:
        return False
    for y, yc in front.items():
        if y != q and cost_dominates(y, yc, q, c):
            return False
    dead = [
        y for y, yc in front.items()
        if y != q and cost_dominates(q, c, y, yc)
    ]
    for y in dead:
        del front[y]
    front[q] = c
    if len(front) > MAX_FRONTIER_PER_SOURCE:
        raise RuntimeError(
            f"matched-history cost frontier exceeded {MAX_FRONTIER_PER_SOURCE} labels at one source"
        )
    return True


def _terminal_insert(front: set[Label], label: Label) -> bool:
    """Pareto-insert a terminal adverse label after finite-word reachability is paid."""
    return HIST.pareto_insert(front, tuple(map(int, label)))


def finite_cost_frontier(rt: dict, ranks: list[Label], target_samples: int) -> dict:
    """Propagate the exact minimum-cost source-label quotient to the word target."""
    N = int(target_samples)
    if N <= 0:
        raise ValueError("positive target sample count required")
    edges = min_gap_successors(rt)
    starts = HIST._start_nodes(rt)
    if not starts:
        raise RuntimeError("P2 V1 has no legal staged source start")

    state: list[dict[Label, int]] = [dict() for _ in rt["nodes"]]
    heap: list[tuple[int, int, Label]] = []
    for s in starts:
        q = HIST.EMPTY_LABEL
        state[int(s)][q] = 0
        heapq.heappush(heap, (0, int(s), q))

    terminal: set[Label] = set()
    processed = 0
    generated = 0
    dominated = 0
    terminal_candidates = 0
    max_source_frontier = 1
    max_heap = len(heap)

    while heap:
        cost, s, label = heapq.heappop(heap)
        if state[s].get(label) != cost:
            continue
        processed += 1
        if processed > MAX_PROCESSED_ENTRIES:
            raise RuntimeError(
                f"matched-history cost frontier exceeded {MAX_PROCESSED_ENTRIES} processed states"
            )
        q_after_s = HIST.update_label(label, ranks[s])
        for t, gap in edges[s]:
            generated += 1
            c2 = cost + int(gap)
            if c2 >= N:
                terminal_candidates += 1
                q_terminal = HIST.update_label(q_after_s, ranks[int(t)])
                if not _terminal_insert(terminal, q_terminal):
                    dominated += 1
                continue
            dst = state[int(t)]
            if insert_cost_frontier(dst, q_after_s, c2):
                heapq.heappush(heap, (c2, int(t), q_after_s))
                max_source_frontier = max(max_source_frontier, len(dst))
                max_heap = max(max_heap, len(heap))
            else:
                dominated += 1

    if not terminal:
        raise RuntimeError("finite-cost source quotient produced no target-crossing history")
    active_entries = sum(len(x) for x in state)
    return {
        "min_gap_successors": edges,
        "start_sources": tuple(map(int, starts)),
        "terminal_labels": terminal,
        "processed_entries": processed,
        "generated_candidates": generated,
        "dominated_candidates": dominated,
        "terminal_candidates": terminal_candidates,
        "active_nonterminal_frontier_entries": active_entries,
        "maximum_frontier_entries_at_one_source": max_source_frontier,
        "maximum_heap_entries": max_heap,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    """Build the non-promoting finite-cost matched-history classes."""
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("matched-history cost frontier must not be trajectory fitted")

    # CORR.runtime is the frozen non-JSON P2 V1 interface.  Do not call
    # CORR.build here: that would re-materialize the ~963k-edge JSON table only
    # to validate a structure runtime() has already checked.  The dedicated P2
    # workflow remains the promotion authority for the interface certificate.
    rt = CORR.runtime(path)
    if CORR.INTERFACE_VERSION != "OU3_P2_CORRELATED_STAGE_TRANSFER_V1":
        raise RuntimeError("P2 correlation interface version changed")
    if len(rt["nodes"]) != 800 or tuple(map(int, rt["gaps"])) != tuple(range(13, 27)):
        raise RuntimeError("frozen P2 V1 source partition or gap alphabet changed")
    if any(not rt["labelled_successors"][s][gi]
           for s in range(len(rt["nodes"])) for gi in range(len(rt["gaps"]))):
        raise RuntimeError("frozen P2 V1 interface has a dead finite-clock kernel")

    sched = BASE.source_schedule()
    h = float(rt["clock"]["dt_binary32_s"])
    target = HIST._global_word_target(domain, sched, h)
    N = int(target["target_samples"])
    stats = HIST._stat_tables(rt, sched)
    ranks = stats["node_ranks"]
    fr = finite_cost_frontier(rt, ranks, N)

    Tpe = BASE.pos(domain["normal_live"]["vector_pe_recurrence_window_s"], "PE recurrence")
    rows = []
    for label in sorted(fr["terminal_labels"]):
        allowed, dominators = UNION._tau_dominators(label, rt, ranks)
        summary = HIST.label_summary(label, {"stats": stats, "target": target})
        upper, timing = CUPPER.translation_upper_from_summary(
            summary, Tpe, sched, require_history_cover=True
        )
        rows.append({
            "adverse_label": list(label),
            "allowed_physical_source_nodes": len(allowed),
            "allowed_tau_cells": [int(rt["nodes"][d]["tau_index"]) for d in dominators],
            "same_tau_covariance_lower_dominator_nodes": dominators,
            "Sigma_translation_diagonal_upper": list(map(float, upper)),
            "timing": timing,
        })

    min_edges = fr["min_gap_successors"]
    edge_count = sum(len(row) for row in min_edges)
    min_gap_histogram: dict[str, int] = {}
    for row in min_edges:
        for _, gap in row:
            key = str(int(gap))
            min_gap_histogram[key] = min_gap_histogram.get(key, 0) + 1

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_MATCHED_HISTORY_FINITE_COST_ADVERSE_FRONTIER",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "canonical_gate_changed": False,
        "P2_correlation_interface_consumed": True,
        "P2_correlation_interface_version": CORR.INTERFACE_VERSION,
        "exact_gap_labelled_successors_consumed": True,
        "minimum_supporting_gap_per_source_edge_used": True,
        "minimum_gap_replacement_preserves_source_sequence_and_can_only_reduce_cost": True,
        "finite_word_sample_budget_retained": True,
        "cost_adverse_dominance_continuation_safe": True,
        "terminal_endpoint_rank_included_before_global_reduction": True,
        "global_cartesian_source_extrema_used": False,
        "target_samples": N,
        "global_covariance_word_upper_s": float(target["global_covariance_word_upper_s"]),
        "min_gap_source_edges": edge_count,
        "min_gap_histogram": min_gap_histogram,
        "processed_cost_frontier_entries": fr["processed_entries"],
        "generated_cost_frontier_candidates": fr["generated_candidates"],
        "dominated_cost_frontier_candidates": fr["dominated_candidates"],
        "target_crossing_candidates": fr["terminal_candidates"],
        "active_nonterminal_frontier_entries": fr["active_nonterminal_frontier_entries"],
        "maximum_frontier_entries_at_one_source": fr["maximum_frontier_entries_at_one_source"],
        "maximum_heap_entries": fr["maximum_heap_entries"],
        "global_pareto_adverse_labels": len(fr["terminal_labels"]),
        "matched_history_classes": rows,
        "same_tau_covariance_lower_dominators_verified": True,
        "matched_margin_computed": False,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "P5_PROMOTED": False,
        "next_obligation": (
            "if the finite-cost labels remain source-correlated, group them by admitted same-tau dominator set and "
            "compute one rigorous whole-word covariance lower per group in that group's own covariance-upper metric; "
            "if they still collapse to the global label, the path-max upper representation rather than elapsed-state enumeration is the next structural limiter"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    """Fail closed on any change to the finite-cost quotient contract."""
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_MATCHED_HISTORY_FINITE_COST_ADVERSE_FRONTIER":
        f.append("wrong qualification")
    for key in (
        "source_only", "P2_correlation_interface_consumed",
        "exact_gap_labelled_successors_consumed",
        "minimum_supporting_gap_per_source_edge_used",
        "minimum_gap_replacement_preserves_source_sequence_and_can_only_reduce_cost",
        "finite_word_sample_budget_retained", "cost_adverse_dominance_continuation_safe",
        "terminal_endpoint_rank_included_before_global_reduction",
        "same_tau_covariance_lower_dominators_verified",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "canonical_gate_changed", "global_cartesian_source_extrema_used",
        "matched_margin_computed", "P3_PROMOTED", "P4_PROMOTED", "P5_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("P2_correlation_interface_version") != CORR.INTERFACE_VERSION:
        f.append("lost frozen P2 V1 binding")
    if int(d.get("target_samples", 0)) != 635:
        f.append("global covariance-word target changed")
    if not int(d.get("min_gap_source_edges", 0)) > 0:
        f.append("minimum-gap source graph is empty")
    if not int(d.get("global_pareto_adverse_labels", 0)) > 0:
        f.append("finite-cost quotient emitted no global label")
    rows = d.get("matched_history_classes", [])
    if len(rows) != int(d.get("global_pareto_adverse_labels", -1)):
        f.append("matched-history row count does not equal terminal label count")
    for row in rows:
        if not row.get("same_tau_covariance_lower_dominator_nodes"):
            f.append("matched history class has no covariance-lower source dominator")
        u = row.get("Sigma_translation_diagonal_upper", [])
        if len(u) != 4 or any(not (isinstance(x, (int, float)) and math.isfinite(float(x)) and float(x) > 0.0) for x in u):
            f.append("matched history class has invalid covariance upper")
        if row.get("timing", {}).get("summarized_history_covers_covariance_word") is not True:
            f.append("matched history label does not cover its covariance word")
    return list(dict.fromkeys(f))


def main() -> int:
    """Write the auditable finite-cost source quotient artifact."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "target_samples": d["target_samples"],
        "min_gap_source_edges": d["min_gap_source_edges"],
        "processed_entries": d["processed_cost_frontier_entries"],
        "max_frontier_per_source": d["maximum_frontier_entries_at_one_source"],
        "global_labels": d["global_pareto_adverse_labels"],
        "classes": [
            {
                "label": row["adverse_label"],
                "allowed_nodes": row["allowed_physical_source_nodes"],
                "tau_dominators": row["same_tau_covariance_lower_dominator_nodes"],
                "upper": row["Sigma_translation_diagonal_upper"],
            }
            for row in d["matched_history_classes"]
        ],
        "validation_failures": vf,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
