#!/usr/bin/env python3
"""Shortest legal witness for collapse of the four-max P3 history summary.

The current same-history covariance upper remembers only four path maxima:
maximum pseudo-update cadence, sigma^2, q_c, and S-measurement variance.  If one
legal P2 V1 history can attain the global maximum of all four statistics before
the 635-sample word target, then that global label is an *actual* history label.
It componentwise dominates every other four-max label, so no Pareto/cost
enumeration built from these same four maxima can retain useful source
correlation at the terminal word.

This producer tests exactly that condition without enumerating the full history
frontier.  Each P2 source transition s->t is assigned its minimum certified
supporting gap.  Replacing a larger supporting gap by the minimum keeps the same
source sequence and can only make a witness cheaper.  A 4-bit mask records which
of the four global maxima have been attained by completed source segments.
Dijkstra over only 800*16 states gives the minimum sample cost of attaining the
global four-max label.  Once attained, the label is absorbing under max-update;
P2 has no dead finite-clock states, so any witness at cost <= target can be
extended until the full word is covered without losing the label.

This is a falsification/selection producer.  It computes no covariance margin
and cannot promote P3/P4/P5.
"""
from __future__ import annotations

import argparse
import heapq
import json
import math
from pathlib import Path

import ou3_p2_correlation_path_memory as CORR
import ou3_p3_matched_history_cost_frontier as COST
import ou3_p3_p2_v1_history_frontier as HIST
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
FULL_MASK = 0b1111


def global_rank_tuple(ranks: list[tuple[int, int, int, int]]) -> tuple[int, int, int, int]:
    """Return the coordinatewise global adverse rank."""
    if not ranks:
        raise ValueError("nonempty rank table required")
    return tuple(max(int(r[j]) for r in ranks) for j in range(4))


def source_max_mask(rank: tuple[int, int, int, int], global_rank: tuple[int, int, int, int]) -> int:
    """Bit mask of global adverse coordinates attained by one physical source."""
    return sum((1 << j) for j in range(4) if int(rank[j]) == int(global_rank[j]))


def shortest_global_label_witness(rt: dict, ranks: list[tuple[int, int, int, int]],
                                  target_samples: int) -> dict:
    """Find the minimum-cost legal source history attaining all four global maxima.

    State `(s,mask)` means the execution is at a P2 stage boundary with staged
    source `s`; `mask` contains maxima attained by *completed following
    segments*.  Traversing exact edge `s->t` costs its minimum certified gap and
    folds source `s` into the mask, exactly matching HIST.frontier_runtime's
    label update order.  If FULL_MASK is reached at cost <= target, the actual
    terminal history family contains the global four-max label because the path
    can be legally extended until it crosses the target.
    """
    N = int(target_samples)
    if N <= 0:
        raise ValueError("positive target sample count required")
    edges = COST.min_gap_successors(rt)
    starts = HIST._start_nodes(rt)
    if not starts:
        raise RuntimeError("P2 V1 has no legal staged start source")
    g = global_rank_tuple(ranks)
    node_masks = [source_max_mask(r, g) for r in ranks]

    inf = 10**18
    dist = [[inf] * 16 for _ in rt["nodes"]]
    prev: dict[tuple[int, int], tuple[int, int, int]] = {}
    heap: list[tuple[int, int, int]] = []
    for s0 in starts:
        s = int(s0)
        dist[s][0] = 0
        heapq.heappush(heap, (0, s, 0))

    goal: tuple[int, int] | None = None
    while heap:
        cost, s, mask = heapq.heappop(heap)
        if cost != dist[s][mask]:
            continue
        if mask == FULL_MASK:
            goal = (s, mask)
            break
        if cost > N:
            break
        mask2 = mask | node_masks[s]
        for t, gap in edges[s]:
            c2 = cost + int(gap)
            if c2 > N:
                continue
            if c2 < dist[int(t)][mask2]:
                dist[int(t)][mask2] = c2
                prev[(int(t), mask2)] = (s, mask, int(gap))
                heapq.heappush(heap, (c2, int(t), mask2))

    if goal is None:
        best_cost = min(dist[s][FULL_MASK] for s in range(len(rt["nodes"])))
        reachable = best_cost <= N
        if reachable:
            end = min(range(len(rt["nodes"])), key=lambda s: dist[s][FULL_MASK])
            goal = (end, FULL_MASK)
    else:
        best_cost = dist[goal[0]][FULL_MASK]
        reachable = best_cost <= N

    path = []
    if reachable and goal is not None:
        s, mask = goal
        chain = []
        while (s, mask) in prev:
            ps, pm, gap = prev[(s, mask)]
            chain.append((ps, s, gap, pm, mask))
            s, mask = ps, pm
        chain.reverse()
        cumulative = 0
        running_mask = 0
        for ps, t, gap, pm, qm in chain:
            if pm != running_mask:
                raise RuntimeError("witness reconstruction lost mask continuity")
            cumulative += int(gap)
            new_mask = running_mask | node_masks[int(ps)]
            if new_mask != qm:
                raise RuntimeError("witness reconstruction lost source-max update")
            path.append({
                "source": int(ps),
                "successor": int(t),
                "gap_samples": int(gap),
                "cumulative_samples": cumulative,
                "source_global_max_mask": int(node_masks[int(ps)]),
                "path_global_max_mask": int(new_mask),
                "source_rank": list(map(int, ranks[int(ps)])),
            })
            running_mask = new_mask
        if cumulative != int(best_cost) or running_mask != FULL_MASK:
            raise RuntimeError("global-label witness reconstruction is inconsistent")

    max_sources = {
        str(j): [i for i, r in enumerate(ranks) if int(r[j]) == int(g[j])]
        for j in range(4)
    }
    return {
        "global_rank": g,
        "node_masks": node_masks,
        "minimum_cost_samples": int(best_cost) if reachable else None,
        "reachable_within_target": bool(reachable),
        "witness_path": path,
        "global_max_source_nodes_by_coordinate": max_sources,
        "min_gap_source_edges": sum(len(row) for row in edges),
        "start_source_count": len(starts),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    """Build the non-promoting four-max-summary collapse witness."""
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("global-label witness must not be trajectory fitted")
    rt = CORR.runtime(path)
    if len(rt["nodes"]) != 800 or tuple(map(int, rt["gaps"])) != tuple(range(13, 27)):
        raise RuntimeError("frozen P2 V1 source partition or gap alphabet changed")
    if any(not rt["labelled_successors"][s][gi]
           for s in range(len(rt["nodes"])) for gi in range(len(rt["gaps"]))):
        raise RuntimeError("P2 V1 finite-clock kernel has a dead labelled transition set")

    sched = BASE.source_schedule()
    h = float(rt["clock"]["dt_binary32_s"])
    target = HIST._global_word_target(domain, sched, h)
    stats = HIST._stat_tables(rt, sched)
    ranks = stats["node_ranks"]
    witness = shortest_global_label_witness(rt, ranks, int(target["target_samples"]))
    reachable = bool(witness["reachable_within_target"])

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_FOUR_MAX_HISTORY_SUMMARY_GLOBAL_LABEL_WITNESS",
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
        "label_update_order_matches_HIST_frontier_runtime": True,
        "target_samples": int(target["target_samples"]),
        "global_covariance_word_upper_s": float(target["global_covariance_word_upper_s"]),
        "global_adverse_rank_label": list(map(int, witness["global_rank"])),
        "global_label_reachable_within_word": reachable,
        "global_label_minimum_cost_samples": witness["minimum_cost_samples"],
        "global_label_witness_path": witness["witness_path"],
        "global_max_source_nodes_by_coordinate": witness["global_max_source_nodes_by_coordinate"],
        "min_gap_source_edges": witness["min_gap_source_edges"],
        "start_source_count": witness["start_source_count"],
        "four_max_summary_retains_terminal_source_correlation": not reachable,
        "four_max_summary_structurally_collapses_to_global_if_reachable": True,
        "matched_margin_computed": False,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "P5_PROMOTED": False,
        "classification": (
            "FOUR_MAX_HISTORY_SUMMARY_COLLAPSES_TO_ACTUAL_GLOBAL_LABEL_WITHIN_WORD"
            if reachable else
            "GLOBAL_FOUR_MAX_LABEL_NOT_REACHABLE_WITHIN_WORD_KEEP_FINITE_COST_FRONTIER"
        ),
        "next_obligation": (
            "replace the four independent path maxima in the covariance upper by a time-ordered/segmented source-history representation; "
            "do not run another covariance lower against the global four-max envelope"
            if reachable else
            "finish the finite-cost matched-history frontier and compute grouped matched lower/upper margins"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    """Fail closed if the witness no longer matches the frozen source contract."""
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_FOUR_MAX_HISTORY_SUMMARY_GLOBAL_LABEL_WITNESS":
        f.append("wrong qualification")
    for key in (
        "source_only", "P2_correlation_interface_consumed",
        "exact_gap_labelled_successors_consumed",
        "minimum_supporting_gap_per_source_edge_used",
        "minimum_gap_replacement_preserves_source_sequence_and_can_only_reduce_cost",
        "label_update_order_matches_HIST_frontier_runtime",
        "four_max_summary_structurally_collapses_to_global_if_reachable",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "canonical_gate_changed", "matched_margin_computed",
        "P3_PROMOTED", "P4_PROMOTED", "P5_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("P2_correlation_interface_version") != CORR.INTERFACE_VERSION:
        f.append("lost frozen P2 V1 binding")
    if int(d.get("target_samples", 0)) != 635:
        f.append("global covariance-word target changed")
    label = d.get("global_adverse_rank_label", [])
    if len(label) != 4 or any(not isinstance(x, int) or x < 0 for x in label):
        f.append("invalid global adverse rank label")
    reachable = d.get("global_label_reachable_within_word") is True
    cost = d.get("global_label_minimum_cost_samples")
    if reachable:
        if not isinstance(cost, int) or not 0 <= cost <= int(d["target_samples"]):
            f.append("reachable global label has invalid witness cost")
        path = d.get("global_label_witness_path", [])
        if not path or int(path[-1].get("path_global_max_mask", -1)) != FULL_MASK:
            f.append("reachable global label has no complete reconstructed witness")
        if d.get("four_max_summary_retains_terminal_source_correlation") is not False:
            f.append("reachable global label incorrectly claims four-max correlation survives")
    else:
        if cost is not None:
            f.append("unreachable global label unexpectedly has a finite witness cost")
        if d.get("four_max_summary_retains_terminal_source_correlation") is not True:
            f.append("unreachable global label lost finite-cost correlation claim")
    return list(dict.fromkeys(f))


def main() -> int:
    """Write the auditable four-max-summary witness artifact."""
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
        "global_label": d["global_adverse_rank_label"],
        "reachable": d["global_label_reachable_within_word"],
        "minimum_cost_samples": d["global_label_minimum_cost_samples"],
        "witness_path": d["global_label_witness_path"],
        "classification": d["classification"],
        "validation_failures": vf,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
