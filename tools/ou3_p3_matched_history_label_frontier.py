#!/usr/bin/env python3
"""Finite-word adverse-label frontier for matched-history OU-III P3 work.

PR #475 established that P3 must compare a whole-word translation covariance
lower with a covariance upper built from the *same source-history class*.  The
first whole-word probe deliberately collapsed all legal upper classes into one
global envelope; that sufficient test was rigorous but failed at 1.39e-27 and,
as expected, reintroduced a cross-tau mismatch.

This producer constructs the next source quotient without propagating covariance.
The 635-sample covariance word can contain at most ceil(635/13)=49 complete P2
finite-clock segments.  Every shorter legal finite-clock history can be extended
to 49 segments because P2 V1 has no dead finite-clock source states.  Extending
can only increase the four adverse maxima used by the retained covariance upper:

    (pseudo cadence, sigma^2, q_c, S-measurement variance).

We therefore propagate exactly 49 steps on the *union-successor* P2 V1 graph,
ignoring the particular 13..26 gap label.  Ignoring the gap can only add source
paths.  At each current source we retain the componentwise Pareto-maximal adverse
labels.  Finally we Pareto-reduce across endpoints as well.  Proving a matched
lower/upper margin for every emitted global label is sufficient for every legal
shorter full-word history, because a more-adverse label has a no-smaller upper
and admits a no-smaller source set for the covariance lower.

For each emitted label this module also identifies the physical source nodes
whose four statistics lie below the label.  Within each represented tau cell the
physical (sigma-index 0, R_S-index 0) node is verified to be admitted and to
Loewner-dominate the covariance lower of every other admitted node in that tau
cell.  Those same-tau nodes are the only source kernels the next rigorous
matched-margin recurrence needs.

This module does not compute a P3 margin and cannot promote P3/P4/P5.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

import ou3_p2_correlation_path_memory as CORR
import ou3_p3_correlated_translation_covariance_upper as CUPPER
import ou3_p3_p2_v1_history_frontier as HIST
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
TAU_CELLS = 10
SIGMA_CELLS = 8
RS_CELLS = 10
TAU_STRIDE = SIGMA_CELLS * RS_CELLS


def maximum_complete_segments(target_samples: int, minimum_gap_samples: int) -> int:
    """Maximum complete finite-clock segments needed to cover a target word."""
    n = int(target_samples)
    g = int(minimum_gap_samples)
    if n <= 0 or g <= 0:
        raise ValueError("positive target and gap are required")
    return int(math.ceil(n / g))


def _reduce(labels) -> set[tuple[int, int, int, int]]:
    """Return the componentwise Pareto-maximal adverse labels."""
    out: set[tuple[int, int, int, int]] = set()
    # Adverse-first ordering usually makes the retained frontier small early.
    for label in sorted(set(tuple(map(int, x)) for x in labels), reverse=True):
        HIST.pareto_insert(out, label)
    return out


def _step(front: dict[int, set[tuple[int, int, int, int]]], rt: dict,
          ranks: list[tuple[int, int, int, int]]) -> dict[int, set[tuple[int, int, int, int]]]:
    """One conservative segment step on the gap-forgotten union-successor graph."""
    candidates: dict[int, set[tuple[int, int, int, int]]] = defaultdict(set)
    for s, labels in front.items():
        nr = ranks[int(s)]
        updated = {HIST.update_label(label, nr) for label in labels}
        for t0 in rt["union_successors"][int(s)]:
            candidates[int(t0)].update(updated)
    return {t: _reduce(labels) for t, labels in candidates.items() if labels}


def _allowed_nodes(label: tuple[int, int, int, int],
                   ranks: list[tuple[int, int, int, int]]) -> list[int]:
    q = tuple(map(int, label))
    return [
        i for i, r in enumerate(ranks)
        if all(int(r[j]) <= q[j] for j in range(4))
    ]


def _tau_dominators(label: tuple[int, int, int, int], rt: dict,
                    ranks: list[tuple[int, int, int, int]]) -> tuple[list[int], list[int]]:
    """Return admitted physical nodes and the same-tau min-(sigma,R_S) dominators."""
    allowed = _allowed_nodes(label, ranks)
    if not allowed:
        raise RuntimeError("adverse label admits no physical P2 source node")
    nodes = rt["nodes"]
    allowed_set = set(allowed)
    tau_present = sorted({int(nodes[i]["tau_index"]) for i in allowed})
    dominators = []
    for ti in tau_present:
        d = ti * TAU_STRIDE
        if d not in allowed_set:
            raise RuntimeError(
                f"label admits tau cell {ti} but not its physical sigma0/R_S0 dominator"
            )
        dn = nodes[d]
        dsig = tuple(map(float, dn["sigma_filter_committed_mps2"]))
        drs = tuple(map(float, dn["R_S_filter_std"]))
        for i in allowed:
            n = nodes[i]
            if int(n["tau_index"]) != ti:
                continue
            nsig = tuple(map(float, n["sigma_filter_committed_mps2"]))
            nrs = tuple(map(float, n["R_S_filter_std"]))
            if dsig[0] > nsig[0] or dsig[1] > nsig[1]:
                raise RuntimeError("same-tau sigma0 node is not covariance-lower dominating")
            if drs[0] > nrs[0] or drs[1] > nrs[1]:
                raise RuntimeError("same-tau R_S0 node is not strongest admitted S measurement")
        dominators.append(d)
    return allowed, dominators


def _summary(label: tuple[int, int, int, int], stats: dict, target: dict) -> dict:
    """Build the retained covariance-upper summary directly from one adverse label."""
    fr = {"stats": stats, "target": target}
    return HIST.label_summary(label, fr)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("matched-history frontier must not be trajectory fitted")

    rt = CORR.runtime(path)
    corr = CORR.build(path)
    cf = CORR.validate(corr)
    if cf:
        raise RuntimeError(f"P2 correlation interface failed: {cf}")
    if corr.get("interface_version") != CORR.INTERFACE_VERSION:
        raise RuntimeError("P2 correlation interface version changed")

    sched = BASE.source_schedule()
    h = float(rt["clock"]["dt_binary32_s"])
    target = HIST._global_word_target(domain, sched, h)
    stats = HIST._stat_tables(rt, sched)
    ranks = stats["node_ranks"]
    gaps = tuple(map(int, rt["gaps"]))
    if gaps != tuple(range(13, 27)):
        raise RuntimeError("P2 finite-clock gap alphabet changed")
    max_segments = maximum_complete_segments(int(target["target_samples"]), min(gaps))

    starts = HIST._start_nodes(rt)
    front: dict[int, set[tuple[int, int, int, int]]] = {
        int(s): {HIST.EMPTY_LABEL} for s in starts
    }
    progress = []
    for k in range(1, max_segments + 1):
        front = _step(front, rt, ranks)
        if not front:
            raise RuntimeError(f"union-successor adverse frontier died at segment {k}")
        total = sum(len(v) for v in front.values())
        progress.append({
            "segments": k,
            "endpoint_source_states": len(front),
            "endpoint_pareto_labels": total,
            "maximum_labels_at_one_endpoint": max(len(v) for v in front.values()),
        })

    global_front = _reduce(label for labels in front.values() for label in labels)
    if not global_front:
        raise RuntimeError("49-segment global adverse frontier is empty")

    Tpe = BASE.pos(domain["normal_live"]["vector_pe_recurrence_window_s"], "PE recurrence")
    rows = []
    for label in sorted(global_front):
        allowed, dominators = _tau_dominators(label, rt, ranks)
        summary = _summary(label, stats, target)
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

    min_history_samples = max_segments * min(gaps)
    coverage = BASE.down(min_history_samples * math.nextafter(h, -math.inf))
    if coverage <= float(target["global_covariance_word_upper_s"]):
        raise RuntimeError("49 minimum-gap segments do not cover global covariance word")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_MATCHED_HISTORY_FINITE_WORD_ADVERSE_LABEL_FRONTIER",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "canonical_gate_changed": False,
        "P2_correlation_interface_consumed": True,
        "P2_correlation_interface_version": CORR.INTERFACE_VERSION,
        "exact_gap_labels_forgotten_only_by_source_superset": True,
        "union_successor_graph_used": True,
        "finite_clock_no_dead_states_required": True,
        "shorter_histories_covered_by_adverse_49_segment_extension": True,
        "endpoint_identity_forgotten_only_after_global_adverse_dominance": True,
        "global_cartesian_source_extrema_used": False,
        "target_samples": int(target["target_samples"]),
        "minimum_gap_samples": min(gaps),
        "maximum_complete_segments": max_segments,
        "minimum_49_segment_history_samples": min_history_samples,
        "minimum_49_segment_history_duration_s": coverage,
        "global_covariance_word_upper_s": float(target["global_covariance_word_upper_s"]),
        "frontier_progress": progress,
        "endpoint_source_states_after_49_segments": len(front),
        "endpoint_pareto_labels_after_49_segments": sum(len(v) for v in front.values()),
        "global_pareto_adverse_labels": len(global_front),
        "matched_history_classes": rows,
        "same_tau_covariance_lower_dominators_verified": True,
        "matched_margin_computed": False,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "P5_PROMOTED": False,
        "next_obligation": (
            "for each emitted adverse label, iterate the rigorous whole-word covariance lower only over its "
            "admitted same-tau dominator nodes and compare in that label's own covariance-upper metric"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_MATCHED_HISTORY_FINITE_WORD_ADVERSE_LABEL_FRONTIER":
        f.append("wrong qualification")
    for key in (
        "source_only", "P2_correlation_interface_consumed",
        "exact_gap_labels_forgotten_only_by_source_superset",
        "union_successor_graph_used", "finite_clock_no_dead_states_required",
        "shorter_histories_covered_by_adverse_49_segment_extension",
        "endpoint_identity_forgotten_only_after_global_adverse_dominance",
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
        f.append("lost frozen P2 V1 interface binding")
    if int(d.get("target_samples", 0)) != 635:
        f.append("global covariance-word target changed")
    if int(d.get("minimum_gap_samples", 0)) != 13:
        f.append("minimum finite-clock gap changed")
    if int(d.get("maximum_complete_segments", 0)) != 49:
        f.append("49-segment finite-word extension bound changed")
    if not int(d.get("global_pareto_adverse_labels", 0)) > 0:
        f.append("no global adverse labels emitted")
    rows = d.get("matched_history_classes", [])
    if len(rows) != int(d.get("global_pareto_adverse_labels", -1)):
        f.append("matched-history row count does not equal global label frontier")
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
        "maximum_complete_segments": d["maximum_complete_segments"],
        "endpoint_labels": d["endpoint_pareto_labels_after_49_segments"],
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
