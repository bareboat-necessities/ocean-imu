#!/usr/bin/env python3
"""Lossless filter-side quotient for the terminal OU-III P4 source cases.

P2 keeps raw tuner sigma as part of its 800-state physical partition because the
EMA candidate memory needs it.  The nonlinear filter word, however, sees only
the committed filter tuple

    (tau, sigma_aw_after_filter_floor, R_S, pseudo cadence).

Two P2 nodes are therefore interchangeable *inside a fixed terminal cluster*
only when these four already-outward interval cells are bit-for-bit identical.
This producer builds that exact quotient and then maps the reachable P2
stage-boundary pairs (c,s) onto ordered pairs of filter-side equivalence
classes.  No hull, midpoint, or free Cartesian switch is introduced.

The terminal-cluster reduction separately proves that samples 192..201 contain
at most one stage boundary, so an ordered quotient pair plus boundary position
is sufficient for the nonlinear source schedule there.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p2_clock_phase_tuner_graph as P2
import ou3_p4_source_node_cells as NODES
import ou3_p4_terminal_cluster_p2_reduction as RED

DEFAULT_DOMAIN = RED.DEFAULT_DOMAIN
SCHEMA = 1


def _node_key(n: dict):
    return (
        tuple(n["tau_s"]),
        tuple(n["sigma_filter_committed_mps2"]),
        tuple(n["R_S_filter_std"]),
        tuple(n["pseudo_update_period_s"]),
    )


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    nodes = NODES.build()
    nf = NODES.validate(nodes)
    red = RED.build(path)
    rf = RED.validate(red)
    states, gaps, labelled, union, _clock = P2._build_labelled_successors(path)
    failures = [f"nodes: {x}" for x in nf] + [f"reduction: {x}" for x in rf]

    if len(states) != len(nodes["nodes"]):
        failures.append("P2 state count differs from source-node materialization")

    class_of = []
    classes = {}
    for n in nodes["nodes"]:
        key = _node_key(n)
        if key not in classes:
            classes[key] = len(classes)
        class_of.append(classes[key])

    reachable_pairs = set()
    quotient_pairs = set()
    for c, outs in enumerate(union):
        for s in outs:
            reachable_pairs.add((c, s))
            quotient_pairs.add((class_of[c], class_of[s]))

    if len(reachable_pairs) != int(red["reachable_stage_pair_states"]):
        failures.append("reachable pair count differs from certified terminal reduction")
    if not quotient_pairs:
        failures.append("filter-side quotient pair set is empty")

    shapes = int(red["including_no_boundary_schedule_shapes_upper"])
    raw_cases = len(reachable_pairs) * shapes
    quotient_cases = len(quotient_pairs) * shapes
    if quotient_cases > raw_cases:
        failures.append("exact quotient increased terminal source cases")

    members = [[] for _ in range(len(classes))]
    for i, q in enumerate(class_of):
        members[q].append(i)

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_TERMINAL_EXACT_FILTER_SIDE_SOURCE_EQUIVALENCE",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "exact_bitwise_interval_equivalence_only": True,
        "source_cells_hulled": False,
        "cartesian_pair_completion_used": False,
        "raw_P2_physical_nodes": len(nodes["nodes"]),
        "filter_side_equivalence_classes": len(classes),
        "class_of_node": class_of,
        "class_member_nodes": members,
        "reachable_stage_pairs": len(reachable_pairs),
        "reachable_filter_class_pairs": len(quotient_pairs),
        "terminal_schedule_shapes_upper": shapes,
        "raw_terminal_pair_phase_cases_upper": raw_cases,
        "quotiented_terminal_pair_phase_cases_upper": quotient_cases,
        "quotient_reduction_factor": (raw_cases / quotient_cases) if quotient_cases else 0.0,
        "terminal_single_boundary_reduction_consumed": red.get("cluster_span_strictly_below_minimum_stage_gap") is True,
        "reachable_pair_correlation_retained": True,
        "ready_for_terminal_nonlinear_quotient_enumeration": not failures,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_only", "exact_bitwise_interval_equivalence_only",
        "terminal_single_boundary_reduction_consumed",
        "reachable_pair_correlation_retained",
        "ready_for_terminal_nonlinear_quotient_enumeration",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "trajectory_replay_used", "filter_changed", "source_cells_hulled",
        "cartesian_pair_completion_used", "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if int(d.get("raw_P2_physical_nodes", 0)) != 800:
        f.append("raw P2 physical node count changed")
    if not 0 < int(d.get("filter_side_equivalence_classes", 0)) <= 800:
        f.append("invalid filter-side class count")
    if not 0 < int(d.get("reachable_filter_class_pairs", 0)) <= int(d.get("reachable_stage_pairs", 0)):
        f.append("invalid reachable quotient-pair count")
    if int(d.get("quotiented_terminal_pair_phase_cases_upper", 0)) > int(d.get("raw_terminal_pair_phase_cases_upper", 0)):
        f.append("terminal quotient did not preserve/nonincrease case count")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain)
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "raw_nodes": d["raw_P2_physical_nodes"],
        "filter_classes": d["filter_side_equivalence_classes"],
        "reachable_pairs": d["reachable_stage_pairs"],
        "quotient_pairs": d["reachable_filter_class_pairs"],
        "shapes": d["terminal_schedule_shapes_upper"],
        "raw_cases": d["raw_terminal_pair_phase_cases_upper"],
        "quotient_cases": d["quotiented_terminal_pair_phase_cases_upper"],
        "reduction_factor": d["quotient_reduction_factor"],
        "validation_failures": f,
    }, indent=2, sort_keys=True))
    return 0 if not f else 2


if __name__ == "__main__":
    raise SystemExit(main())
