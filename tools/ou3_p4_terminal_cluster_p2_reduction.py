#!/usr/bin/env python3
"""Exact P2 reduction for the terminal nonlinear P4 observation cluster.

The complete H/A word contains a late cluster carrying the mandatory vector-PE
packet and the four spread S=0 observations.  The source-faithful schedule used
by the complete-word backend places its first mandatory vector packet at sample
192 and its final vector/S packet at sample 201, a span of only nine valid-sample
intervals.

P2 proves that two finite-clock tuner staging boundaries are separated by at
least 13 valid samples.  A stage-boundary node is (c,s): the staging sample uses
committed tuple c, and the next valid sample through the next boundary uses the
staged tuple s.  Consequently a 9-interval terminal cluster can contain at most
one stage boundary and therefore at most one applied-tuple change.  No
second-order (c,s)->(s,t) transition is needed inside this cluster.

Every possible terminal-cluster source schedule is exactly represented by:

  * one reachable P2 stage pair (c,s), and
  * either no boundary in the cluster, or one boundary sample b.

If b is present, the staging sample uses c and the next valid sample onward uses
s.  Depending on the incoming clock phase, the cluster may lie wholly on one
side of that boundary.  Equivalently, the applied source sequence is one
constant tuple or one c-to-s step.  This reduction preserves staged/committed
correlation and never permits a Cartesian free switch.

This is a source-language reduction only.  It does not itself establish P4
nonlinear dissipation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p2_clock_phase_tuner_graph as P2
import ou3_p4_joint_word_dissipation_design as WORD
import ou3_p4_source_node_cells as NODES
import ou3_implementation_word_language as LANG

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    p2 = P2.build(path)
    p2f = P2.validate(p2)
    lang = LANG.build(path)
    lf = LANG.validate(lang)
    nodes = NODES.build()
    nf = NODES.validate(nodes)
    failures = (
        [f"P2: {x}" for x in p2f]
        + [f"language: {x}" for x in lf]
        + [f"source-nodes: {x}" for x in nf]
    )

    h = float(lang["word_contract"]["configured_runtime"]["imu_dt_s"])
    samples = int(lang["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"])
    schedule = WORD._schedule(path, samples, h)
    mandatory = sorted(set(schedule["S_steps"] + schedule["vector_steps"]))
    if not mandatory:
        failures.append("terminal mandatory operation schedule is empty")
        first = last = span = None
    else:
        first, last = mandatory[0], mandatory[-1]
        span = last - first

    gaps = list(p2["clock_phase_gap_alphabet_samples"])
    min_gap = min(gaps) if gaps else None
    if span is None or min_gap is None or not span < min_gap:
        failures.append(
            f"terminal cluster span {span} is not shorter than minimum P2 stage gap {min_gap}"
        )

    pair_states = int(p2["stage_boundary_pair_states"])
    physical = int(p2["physical_partition_states"])
    materialized = int(nodes["partition"]["states"])
    if physical != materialized:
        failures.append("P2 and P4 source-node physical partitions differ")

    # A boundary can be at any cluster sample; add one sentinel for no boundary.
    # This is a safe upper count. Exact dynamic enumeration may discard phases
    # inconsistent with the incoming P2 clock state, but can never require a
    # second successor tuple t because span < minimum stage gap.
    boundary_positions_upper = (span + 1) if span is not None else 0
    schedule_shapes_upper = boundary_positions_upper + 1
    reduced_pair_phase_cases_upper = pair_states * schedule_shapes_upper

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_TERMINAL_CLUSTER_SINGLE_P2_STAGE_BOUNDARY_REDUCTION",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "P2_phased_graph_pass": not p2f,
        "P4_source_node_materialization_pass": not nf,
        "physical_source_states": physical,
        "reachable_stage_pair_states": pair_states,
        "minimum_finite_stage_gap_valid_samples": min_gap,
        "terminal_mandatory_operation_samples": mandatory,
        "terminal_cluster_first_sample": first,
        "terminal_cluster_last_sample": last,
        "terminal_cluster_span_intervals": span,
        "cluster_span_strictly_below_minimum_stage_gap": bool(
            span is not None and min_gap is not None and span < min_gap
        ),
        "maximum_stage_boundaries_inside_cluster": 1,
        "maximum_applied_tuple_changes_inside_cluster": 1,
        "second_order_P2_successor_needed_inside_cluster": False,
        "terminal_source_schedule_form": "constant tuple OR one reachable c_to_s step",
        "reachable_pair_correlation_retained": True,
        "clock_phase_retained": True,
        "arbitrary_cartesian_tuner_switching_used": False,
        "boundary_positions_inside_cluster_upper": boundary_positions_upper,
        "including_no_boundary_schedule_shapes_upper": schedule_shapes_upper,
        "reduced_pair_phase_cases_upper": reduced_pair_phase_cases_upper,
        "full_second_order_gap_labelled_edges_avoided": int(
            p2["second_order_gap_labelled_edges_factorized_count"]
        ),
        "ready_for_terminal_nonlinear_source_enumeration": not failures,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_only", "P2_phased_graph_pass", "P4_source_node_materialization_pass",
        "cluster_span_strictly_below_minimum_stage_gap",
        "reachable_pair_correlation_retained", "clock_phase_retained",
        "ready_for_terminal_nonlinear_source_enumeration",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed",
        "second_order_P2_successor_needed_inside_cluster",
        "arbitrary_cartesian_tuner_switching_used",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if int(d.get("maximum_stage_boundaries_inside_cluster", -1)) != 1:
        f.append("terminal cluster does not have one-boundary maximum")
    if int(d.get("maximum_applied_tuple_changes_inside_cluster", -1)) != 1:
        f.append("terminal cluster does not have one applied-tuple-change maximum")
    if not int(d.get("terminal_cluster_span_intervals", 10**9)) < int(
        d.get("minimum_finite_stage_gap_valid_samples", 0)
    ):
        f.append("terminal cluster span does not prove single-boundary reduction")
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
        "mandatory_samples": d["terminal_mandatory_operation_samples"],
        "first": d["terminal_cluster_first_sample"],
        "last": d["terminal_cluster_last_sample"],
        "span": d["terminal_cluster_span_intervals"],
        "min_stage_gap": d["minimum_finite_stage_gap_valid_samples"],
        "reachable_pairs": d["reachable_stage_pair_states"],
        "schedule_shapes_upper": d["including_no_boundary_schedule_shapes_upper"],
        "reduced_pair_phase_cases_upper": d["reduced_pair_phase_cases_upper"],
        "validation_failures": f,
    }, indent=2, sort_keys=True))
    return 0 if not f else 2


if __name__ == "__main__":
    raise SystemExit(main())
