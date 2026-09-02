#!/usr/bin/env python3
"""Endpoint-conditioned source envelopes for the OU-III P3 word.

P3 needs a covariance upper that is valid for every tuner trajectory reaching
an endpoint source node, not a static-node replay and not the full Cartesian
source box.  The shipping tuner changes only at the staged/committed clock
boundaries certified by :mod:`ou3_p2_clock_phase_tuner_graph`.  This producer
therefore computes, for every physical P2 node, the complete reverse-reachable
set over one conservative P3 word horizon.

The result is source evidence only.  It does not claim a covariance bound or
promote P3.  Its purpose is to determine the exact endpoint-conditioned
(tau,sigma,R_S) envelope that the next P3 covariance/process comparison must
cover.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p2_clock_phase_tuner_graph as CLOCK
import ou3_p4_source_node_cells as NODES
import ou3_source_reachable_matrix_p3 as P3

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def _p3_word_horizon_upper_s() -> float:
    """Return the source-uniform horizon used by P3's covariance upper."""
    sched = P3.source_schedule()
    h = float(sched["dt_s"])
    tau_lo, tau_hi = map(float, sched["tau_applied_invariant_s"])
    cadence = P3.cadence_bounds(Interval.outward_bounds(tau_lo, tau_hi), sched)
    gap = P3.up(float(cadence[1]) + h)
    Tpe = 1.0
    spacing = P3.up(max(Tpe, 2.0 * gap))
    Tobs = P3.up(2.0 * spacing + gap)
    return P3.up(Tobs + Tpe)


def _reverse_union(labelled, gaps):
    n = len(labelled)
    rev = [set() for _ in range(n)]
    for s in range(n):
        for gi, _gap in enumerate(gaps):
            for t in labelled[s][gi]:
                rev[t].add(s)
    # The certified floating-clock stagnation branch holds the committed tuple.
    for s in range(n):
        rev[s].add(s)
    return rev


def _ancestor_set(endpoint: int, reverse, stages: int):
    reached = {int(endpoint)}
    frontier = {int(endpoint)}
    for _ in range(int(stages)):
        nxt = set()
        for t in frontier:
            nxt.update(reverse[t])
        nxt.difference_update(reached)
        if not nxt:
            break
        reached.update(nxt)
        frontier = nxt
    return reached


def _envelope(indices, nodes):
    rows = [nodes[i] for i in sorted(indices)]

    def span(key):
        return [
            min(float(r[key][0]) for r in rows),
            max(float(r[key][1]) for r in rows),
        ]

    return {
        "tau_s": span("tau_s"),
        "sigma_tuner_raw_mps2": span("sigma_tuner_raw_mps2"),
        "sigma_filter_committed_mps2": span("sigma_filter_committed_mps2"),
        "R_S_filter_std": span("R_S_filter_std"),
        "pseudo_update_period_s": span("pseudo_update_period_s"),
        "tau_indices": sorted({int(r["tau_index"]) for r in rows}),
        "sigma_raw_indices": sorted({int(r["sigma_raw_index"]) for r in rows}),
        "R_S_indices": sorted({int(r["R_S_index"]) for r in rows}),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P3 ancestor envelope must not be trajectory fitted")

    materialized = NODES.build()
    nf = NODES.validate(materialized)
    if nf:
        raise RuntimeError(f"P2 source-node materialization failed: {nf}")
    nodes = materialized["nodes"]

    states, gaps, labelled, _union, clock = CLOCK._build_labelled_successors(domain_path)
    if len(states) != len(nodes) or len(nodes) != 800:
        raise RuntimeError("clock-phase graph and P2 materialization disagree")
    reverse = _reverse_union(labelled, gaps)

    word_upper = _p3_word_horizon_upper_s()
    h = float(clock["dt_binary32_s"])
    min_gap = int(clock["finite_stage_spacing_valid_samples_lower"])
    # Reverse reachability must admit the largest possible number of stage
    # boundaries inside the word, hence divide by the minimum stage spacing and
    # round upward. One extra boundary covers endpoint phase convention.
    max_stage_transitions = int(math.ceil(word_upper / (min_gap * h))) + 1

    endpoint_rows = []
    full = len(nodes)
    for endpoint in range(full):
        anc = _ancestor_set(endpoint, reverse, max_stage_transitions)
        endpoint_rows.append({
            "endpoint_node": endpoint,
            "ancestor_count": len(anc),
            "ancestor_fraction": len(anc) / full,
            "all_800_nodes_reachable_backward": len(anc) == full,
            "envelope": _envelope(anc, nodes),
        })

    counts = [r["ancestor_count"] for r in endpoint_rows]
    widest = max(endpoint_rows, key=lambda r: r["ancestor_count"])
    narrowest = min(endpoint_rows, key=lambda r: r["ancestor_count"])
    failures = []
    if gaps != list(range(13, 27)):
        failures.append("clock-phase gap alphabet changed")
    if not (word_upper > 0.0 and max_stage_transitions > 0):
        failures.append("invalid P3 word horizon/stage count")
    if any(c <= 0 or c > full for c in counts):
        failures.append("invalid endpoint ancestor count")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_CLOCK_PHASE_ENDPOINT_ANCESTOR_ENVELOPES",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "physical_source_nodes": full,
        "clock_phase_gap_alphabet_samples": gaps,
        "P3_word_horizon_upper_s": word_upper,
        "minimum_stage_gap_samples": min_gap,
        "reverse_stage_transitions_admitted": max_stage_transitions,
        "endpoint_ancestor_count_min": min(counts),
        "endpoint_ancestor_count_max": max(counts),
        "endpoint_ancestor_count_mean": sum(counts) / len(counts),
        "narrowest_endpoint": narrowest,
        "widest_endpoint": widest,
        "diagnostic_endpoints": {
            "0": endpoint_rows[0],
            "729": endpoint_rows[729],
        },
        "endpoints": endpoint_rows,
        "P3_COVARIANCE_UPPER_ESTABLISHED_HERE": False,
        "P3_PROMOTED": False,
        "next_obligation": (
            "evaluate the P3 covariance upper and whole-word process floor on each endpoint-conditioned ancestor envelope; "
            "retain the staged/committed path language rather than replacing these envelopes by the full Cartesian source box"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_CLOCK_PHASE_ENDPOINT_ANCESTOR_ENVELOPES":
        f.append("wrong qualification")
    if d.get("source_only") is not True:
        f.append("source_only is not true")
    for key in (
        "trajectory_replay_used", "filter_changed",
        "P3_COVARIANCE_UPPER_ESTABLISHED_HERE", "P3_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if int(d.get("physical_source_nodes", 0)) != 800:
        f.append("physical source partition is not 800 nodes")
    if d.get("clock_phase_gap_alphabet_samples") != list(range(13, 27)):
        f.append("clock phase alphabet is not 13..26")
    rows = d.get("endpoints", [])
    if len(rows) != 800:
        f.append("missing endpoint ancestor rows")
    for row in rows:
        c = int(row.get("ancestor_count", 0))
        if not 1 <= c <= 800:
            f.append("endpoint has invalid ancestor count")
            break
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
        "validation_pass": not vf,
        "word_horizon_upper_s": d["P3_word_horizon_upper_s"],
        "reverse_stage_transitions": d["reverse_stage_transitions_admitted"],
        "ancestor_count_min": d["endpoint_ancestor_count_min"],
        "ancestor_count_max": d["endpoint_ancestor_count_max"],
        "ancestor_count_mean": d["endpoint_ancestor_count_mean"],
        "node_0": d["diagnostic_endpoints"]["0"],
        "node_729": d["diagnostic_endpoints"]["729"],
        "failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
