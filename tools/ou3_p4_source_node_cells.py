#!/usr/bin/env python3
"""Exact source-node cell materialization shared by the nonlinear P4 routes.

P2 already certifies a source-dynamic graph over 800 deployed tuner cells, but
its state partition historically remained internal to
:mod:`ou3_p4_source_path_reachability`.  The nonlinear operation-matched and
whole-word routes both need the *same* start/end source coordinates before they
can attach covariance/information metrics to a graph edge.

This module exposes that existing partition without inventing a second graph.
It intentionally calls the P2 backend's source-parsing and cell-partition
primitives, preserving its exact ordering:

    node = ((tau_index * 8) + sigma_raw_index) * 10 + R_S_index.

The raw tuner sigma state remains distinct from the filter-side committed sigma
floor.  Pseudo-update cadence remains coupled to the node's applied tau cell.
No trajectory/replay values, source-path pruning, or nonlinear theorem
promotion are introduced here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_source_path_reachability as PATH
import ou3_source_reachable_matrix_p3 as P3

SCHEMA = 1
TAU_CELLS = 10
SIGMA_RAW_CELLS = 8
RS_CELLS = 10
EXPECTED_STATES = TAU_CELLS * SIGMA_RAW_CELLS * RS_CELLS


def _partition():
    """Reconstruct the exact P2 state partition from the P2 source backend."""
    c = PATH._constants()
    tau_lo = max(
        c["min_tau"],
        PATH.down(c["tau_coeff"] * 0.5 / c["max_freq"]),
    )
    tau = PATH._cells(tau_lo, c["max_tau"], TAU_CELLS)
    sigma_raw = PATH._cells(PATH.RAW_SIGMA_GRAPH_LOWER, c["max_sigma"], SIGMA_RAW_CELLS)
    rs = PATH._cells(c["min_RS"], c["max_RS"], RS_CELLS)
    return c, tau, sigma_raw, rs


def build() -> dict:
    """Return all 800 JSON-safe source nodes in exact P2 indexing order."""
    c, tau, sigma_raw, rs = _partition()
    sched = P3.source_schedule()
    nodes = []
    subfloor_raw = 0
    floor_active = 0

    for ti, t in enumerate(tau):
        tau_iv = Interval(float(t[0]), float(t[1]))
        cadence = P3.cadence_bounds(tau_iv, sched)
        for si, sraw in enumerate(sigma_raw):
            sfilt = PATH._filter_sigma_box(sraw)
            if float(sraw[0]) < PATH.FILTER_SIGMA_FLOOR:
                subfloor_raw += RS_CELLS
            if float(sfilt[0]) <= PATH.FILTER_SIGMA_FLOOR <= float(sfilt[1]):
                floor_active += RS_CELLS
            for ri, r in enumerate(rs):
                index = ((ti * SIGMA_RAW_CELLS) + si) * RS_CELLS + ri
                nodes.append({
                    "index": index,
                    "tau_index": ti,
                    "sigma_raw_index": si,
                    "R_S_index": ri,
                    "tau_s": [float(t[0]), float(t[1])],
                    "sigma_tuner_raw_mps2": [float(sraw[0]), float(sraw[1])],
                    "sigma_filter_committed_mps2": [float(sfilt[0]), float(sfilt[1])],
                    "R_S_filter_std": [float(r[0]), float(r[1])],
                    "pseudo_update_period_s": [float(cadence[0]), float(cadence[1])],
                })

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P2_SOURCE_NODE_CELL_MATERIALIZATION",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "source_graph_rebuilt_or_pruned_here": False,
        "state_order_matches_P2_nested_loops": True,
        "state_index_formula": "((tau_index*8)+sigma_raw_index)*10+R_S_index",
        "partition": {
            "tau": len(tau),
            "sigma_tuner_raw": len(sigma_raw),
            "R_S": len(rs),
            "states": len(nodes),
        },
        "dt_s": float(c["dt"]),
        "raw_tuner_sigma_partition_lower": float(sigma_raw[0][0]),
        "filter_sigma_floor_mps2": float(PATH.FILTER_SIGMA_FLOOR),
        "raw_tuner_sigma_subfloor_node_count": subfloor_raw,
        "filter_sigma_floor_intersecting_node_count": floor_active,
        "nodes": nodes,
        "P4_metric_attached_here": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "attach the actual source-correlated H/A covariance-information metric to each materialized P2 node, then evaluate nonlinear directional pullbacks or whole-word generalized Jacobians on the existing reachable g->h edge family"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    """Fail closed if materialization drifts from the certified P2 partition."""
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P2_SOURCE_NODE_CELL_MATERIALIZATION":
        f.append("wrong qualification")
    for key in ("source_only", "state_order_matches_P2_nested_loops"):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "source_graph_rebuilt_or_pruned_here",
        "P4_metric_attached_here", "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    p = d.get("partition", {})
    expected = {
        "tau": TAU_CELLS,
        "sigma_tuner_raw": SIGMA_RAW_CELLS,
        "R_S": RS_CELLS,
        "states": EXPECTED_STATES,
    }
    if p != expected:
        f.append(f"P2 source partition changed: {p!r}")
    nodes = d.get("nodes", [])
    if len(nodes) != EXPECTED_STATES:
        f.append("source-node list does not contain 800 cells")
        return list(dict.fromkeys(f))
    floor = float(d.get("filter_sigma_floor_mps2", math.nan))
    if floor != PATH.FILTER_SIGMA_FLOOR:
        f.append("filter sigma floor changed")
    if not float(d.get("raw_tuner_sigma_partition_lower", math.inf)) < floor:
        f.append("raw tuner sigma partition no longer extends below filter floor")
    if int(d.get("raw_tuner_sigma_subfloor_node_count", 0)) <= 0:
        f.append("no sub-floor raw tuner sigma nodes materialized")
    if int(d.get("filter_sigma_floor_intersecting_node_count", 0)) <= 0:
        f.append("filter-side sigma floor is not represented")

    for expected_index, node in enumerate(nodes):
        if int(node.get("index", -1)) != expected_index:
            f.append("source-node ordering/index formula drifted")
            break
        ti = int(node["tau_index"])
        si = int(node["sigma_raw_index"])
        ri = int(node["R_S_index"])
        if expected_index != ((ti * SIGMA_RAW_CELLS) + si) * RS_CELLS + ri:
            f.append("source-node index does not match P2 nested-loop formula")
            break
        for key in ("tau_s", "sigma_tuner_raw_mps2", "sigma_filter_committed_mps2", "R_S_filter_std", "pseudo_update_period_s"):
            lo, hi = map(float, node[key])
            if not (math.isfinite(lo) and math.isfinite(hi) and 0.0 < lo <= hi):
                f.append(f"node {expected_index} has invalid {key}")
                break
        if f:
            break
        raw = node["sigma_tuner_raw_mps2"]
        filt = node["sigma_filter_committed_mps2"]
        if float(filt[0]) < floor or float(filt[1]) < float(raw[1]):
            f.append(f"node {expected_index} violates raw/filter sigma relation")
            break
    return list(dict.fromkeys(f))


def node(index: int, payload: dict | None = None) -> dict:
    """Return one JSON-safe P2 source node by exact graph index."""
    d = build() if payload is None else payload
    i = int(index)
    if not 0 <= i < EXPECTED_STATES:
        raise IndexError("P2 source-node index outside [0,800)")
    return d["nodes"][i]


def h18_source_cell(index: int, payload: dict | None = None) -> dict:
    """Convert one P2 node to the interval source-cell schema used by H18 tools."""
    d = build() if payload is None else payload
    n = node(index, d)
    return {
        "source_node_index": int(n["index"]),
        "dt_s": float(d["dt_s"]),
        "tau_s": Interval.outward_bounds(*map(float, n["tau_s"])),
        "sigma_aw_mps2": Interval.outward_bounds(*map(float, n["sigma_filter_committed_mps2"])),
        "R_S_filter_std": Interval.outward_bounds(*map(float, n["R_S_filter_std"])),
        "pseudo_period_s": Interval.outward_bounds(*map(float, n["pseudo_update_period_s"])),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "partition": d["partition"],
        "raw_sigma_lower": d["raw_tuner_sigma_partition_lower"],
        "filter_sigma_floor": d["filter_sigma_floor_mps2"],
        "subfloor_nodes": d["raw_tuner_sigma_subfloor_node_count"],
        "floor_intersecting_nodes": d["filter_sigma_floor_intersecting_node_count"],
        "validation_failures": vf,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
