#!/usr/bin/env python3
"""Recent backward source cones for the OU-III P3/P4 finite-speed tuner graph.

Long-horizon reverse reachability of the coarse 800-cell physical tuner quotient
can eventually become complete.  That does not imply that a *recent* process or
measurement window must forget the finite source slew.  This diagnostic asks a
more local question: for an endpoint physical source cell, which predecessor
cells can occur within at most the number of stage transitions permitted by a
short physical horizon?

The graph is the retained source-only 13..26-sample committed-tuner refinement.
For a requested horizon H, an arbitrary window can begin/end inside a stage, so
we conservatively allow

    ceil(H / (13*dt)) + 2

finite-stage transitions.  Using the minimum stage duration and two extra
boundary crossings can only add predecessor paths.  The frozen-clock branch is
already represented by self-hold edges.

This is a diagnostic only.  It does not turn a short cone into a covariance
certificate, does not assume a tight SpectralMSE powf/sqrtf target, and cannot
promote P3/P4/P5.  Its purpose is to determine whether recent-history source
conditioning is worth carrying into the source-indexed invariant metric.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_sample_clock_source_refinement as CLOCK
import ou3_p4_source_node_cells as NODES

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
DEFAULT_HORIZONS = (0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 3.20)


def _iter_bits(mask: int):
    while mask:
        lsb = mask & -mask
        yield lsb.bit_length() - 1
        mask ^= lsb


def _quantiles(values: list[int]) -> dict:
    if not values:
        raise RuntimeError("empty recent-cone population")
    x = sorted(int(v) for v in values)

    def q(p: float) -> int:
        i = min(len(x) - 1, max(0, int(math.ceil(p * len(x))) - 1))
        return int(x[i])

    return {
        "min": int(x[0]),
        "p50": q(0.50),
        "p90": q(0.90),
        "p99": q(0.99),
        "max": int(x[-1]),
    }


def _cone_snapshot(masks: list[int], nodes: list[dict], n: int) -> dict:
    counts = [m.bit_count() for m in masks]
    tau_spans: list[int] = []
    sigma_spans: list[int] = []
    rs_spans: list[int] = []
    full_tau = full_sigma = full_rs = 0
    worst = None

    for endpoint, mask in enumerate(masks):
        tis: list[int] = []
        sis: list[int] = []
        ris: list[int] = []
        for k in _iter_bits(mask):
            node = nodes[k]
            tis.append(int(node["tau_index"]))
            sis.append(int(node["sigma_raw_index"]))
            ris.append(int(node["R_S_index"]))
        if not tis:
            raise RuntimeError("recent source cone lost endpoint self-state")
        ts = max(tis) - min(tis)
        ss = max(sis) - min(sis)
        rs = max(ris) - min(ris)
        tau_spans.append(ts)
        sigma_spans.append(ss)
        rs_spans.append(rs)
        full_tau += ts == CLOCK.TAU_CELLS - 1
        full_sigma += ss == CLOCK.SIGMA_CELLS - 1
        full_rs += rs == CLOCK.RS_CELLS - 1
        row = {
            "endpoint": endpoint,
            "ancestors": mask.bit_count(),
            "tau_index_range": [min(tis), max(tis)],
            "sigma_raw_index_range": [min(sis), max(sis)],
            "R_S_index_range": [min(ris), max(ris)],
        }
        if worst is None or row["ancestors"] > worst["ancestors"]:
            worst = row

    return {
        "ancestor_count": _quantiles(counts),
        "complete_800_cell_cones": sum(c == n for c in counts),
        "all_endpoints_complete": all(c == n for c in counts),
        "tau_index_span": _quantiles(tau_spans),
        "sigma_raw_index_span": _quantiles(sigma_spans),
        "R_S_index_span": _quantiles(rs_spans),
        "full_tau_span_endpoints": int(full_tau),
        "full_sigma_span_endpoints": int(full_sigma),
        "full_R_S_span_endpoints": int(full_rs),
        "worst_endpoint": worst,
    }


def build(domain_path: Path = DEFAULT_DOMAIN,
          horizons_s: tuple[float, ...] = DEFAULT_HORIZONS) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("recent source cones must not be trajectory fitted")

    nodes_payload = NODES.build()
    nf = NODES.validate(nodes_payload)
    if nf:
        raise RuntimeError(f"source node producer failed: {nf}")
    graph_payload = CLOCK.build(path)
    gf = CLOCK.validate(graph_payload)
    if gf:
        raise RuntimeError(f"finite-speed graph failed: {gf}")

    nodes = nodes_payload["nodes"]
    graph = graph_payload["graph"]
    n = len(nodes)
    if n != NODES.EXPECTED_STATES or len(graph) != n:
        raise RuntimeError("source-node / finite-speed graph cardinality mismatch")

    dt = float(graph_payload["clock"]["dt_binary32_s"])
    min_gap = int(graph_payload["finite_stage_gap_lower_samples"])
    min_stage_s = min_gap * dt
    if not (dt > 0.0 and min_gap > 0 and min_stage_s > 0.0):
        raise RuntimeError("invalid finite-stage timing")

    requested = sorted(set(float(h) for h in horizons_s))
    if not requested or any(not (math.isfinite(h) and h > 0.0) for h in requested):
        raise ValueError("positive finite recent-cone horizons required")
    hops_for = {
        h: int(math.ceil(h / min_stage_s)) + 2
        for h in requested
    }
    max_hops = max(hops_for.values())

    # Reverse adjacency as Python integer bitsets.  The graph already contains
    # the frozen-clock self-hold, and the explicit endpoint bit below guarantees
    # the <=hops cone semantics even if that representation changes later.
    pred = [0] * n
    for i, outs in enumerate(graph):
        bit = 1 << i
        for j in outs:
            pred[int(j)] |= bit
    masks = [1 << j for j in range(n)]

    wanted_by_hop: dict[int, list[float]] = {}
    for h, k in hops_for.items():
        wanted_by_hop.setdefault(k, []).append(h)

    rows: dict[str, dict] = {}
    if 0 in wanted_by_hop:
        snap = _cone_snapshot(masks, nodes, n)
        for h in wanted_by_hop[0]:
            rows[f"{h:.6g}"] = {"max_stage_transitions": 0, **snap}

    saturated_at = None
    for hop in range(1, max_hops + 1):
        new_masks: list[int] = []
        allmask = (1 << n) - 1
        for endpoint in range(n):
            old = masks[endpoint]
            expanded = old
            for k in _iter_bits(old):
                expanded |= pred[k]
                if expanded == allmask:
                    break
            new_masks.append(expanded)
        masks = new_masks
        if saturated_at is None and all(m == allmask for m in masks):
            saturated_at = hop

        if hop in wanted_by_hop:
            snap = _cone_snapshot(masks, nodes, n)
            for h in wanted_by_hop[hop]:
                rows[f"{h:.6g}"] = {"max_stage_transitions": hop, **snap}

        if saturated_at is not None and hop >= max(wanted_by_hop):
            break

    first_all_complete_horizon = None
    for h in requested:
        row = rows[f"{h:.6g}"]
        if row["all_endpoints_complete"]:
            first_all_complete_horizon = h
            break

    one_second_rows = [rows[f"{h:.6g}"] for h in requested if h <= 1.0]
    recent_localization = bool(one_second_rows) and any(
        int(r["complete_800_cell_cones"]) < n for r in one_second_rows
    )

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_RECENT_BACKWARD_SOURCE_CONE_DIAGNOSTIC",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "diagnostic_only": True,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "physical_source_states": n,
        "finite_speed_transition_edges": int(graph_payload["transition_edges"]),
        "base_untimed_transition_edges": int(graph_payload["base_transition_edges"]),
        "clock_dt_binary32_s": dt,
        "minimum_finite_stage_gap_samples": min_gap,
        "minimum_finite_stage_duration_s": min_stage_s,
        "arbitrary_window_phase_extra_stage_transitions": 2,
        "hop_bound_formula": "ceil(H/(min_gap*dt))+2",
        "hop_bound_is_conservative": True,
        "full_deployed_R_S_target_clamp_retained": True,
        "powf_sqrtf_target_tightening_used": False,
        "horizons": rows,
        "all_endpoint_cones_saturate_by_stage_hops": saturated_at,
        "first_requested_horizon_all_endpoints_complete_s": first_all_complete_horizon,
        "recent_le_1s_source_localization_exists": recent_localization,
        "next_obligation": (
            "if recent cones remain localized, consume their source-coordinate extrema in a rigorously propagated source-indexed covariance/information invariant; if they saturate almost immediately, prioritize a correlated target/invariant construction rather than endpoint ancestor conditioning"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_RECENT_BACKWARD_SOURCE_CONE_DIAGNOSTIC":
        f.append("wrong qualification")
    for key in (
        "source_only", "diagnostic_only", "hop_bound_is_conservative",
        "full_deployed_R_S_target_clamp_retained",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "P3_PROMOTED", "P4_PROMOTED", "powf_sqrtf_target_tightening_used",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if int(d.get("physical_source_states", 0)) != NODES.EXPECTED_STATES:
        f.append("physical source-state count changed")
    if int(d.get("finite_speed_transition_edges", 0)) <= 0:
        f.append("finite-speed source graph is empty")
    if int(d.get("finite_speed_transition_edges", 0)) >= int(d.get("base_untimed_transition_edges", 0)):
        f.append("finite-speed graph did not refine untimed source relation")
    if int(d.get("minimum_finite_stage_gap_samples", 0)) != 13:
        f.append("minimum finite stage gap changed")
    if int(d.get("arbitrary_window_phase_extra_stage_transitions", -1)) != 2:
        f.append("arbitrary-window phase cover changed")
    if not d.get("horizons"):
        f.append("no recent-cone horizons emitted")
    for h, row in d.get("horizons", {}).items():
        c = row.get("ancestor_count", {})
        if not (1 <= int(c.get("min", 0)) <= int(c.get("max", 0)) <= NODES.EXPECTED_STATES):
            f.append(f"horizon {h}: invalid ancestor-count range")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--horizons-s", type=float, nargs="*", default=list(DEFAULT_HORIZONS))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, tuple(args.horizons_s))
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "finite_speed_transition_edges": d["finite_speed_transition_edges"],
        "saturate_by_hops": d["all_endpoint_cones_saturate_by_stage_hops"],
        "first_requested_all_complete_s": d["first_requested_horizon_all_endpoints_complete_s"],
        "recent_le_1s_localized": d["recent_le_1s_source_localization_exists"],
        "horizons": d["horizons"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
