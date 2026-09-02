#!/usr/bin/env python3
"""Recent backward-cone diagnostic from the frozen P2 correlation interface.

This diagnostic asks only how quickly the *physical-node projection* of the
pair-state P2 automaton loses endpoint localization over short windows.  It is
not a P3 covariance/information input: the frozen P2 contract forbids replacing
correlated propagation by this projected ancestor hull.

For horizon H we conservatively allow ceil(H/(13*dt))+2 stage transitions to
cover arbitrary start/end phase.  Projection to physical nodes only adds paths,
so the reported cone is an upper reachability diagnostic.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p2_correlation_path_memory as P2
import ou3_p4_source_node_cells as NODES

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
DEFAULT_HORIZONS = (0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 3.20)


def _iter_bits(mask: int):
    while mask:
        lsb = mask & -mask
        yield lsb.bit_length() - 1
        mask ^= lsb


def _quantiles(values: list[int]) -> dict:
    x = sorted(map(int, values))
    if not x:
        raise RuntimeError("empty cone population")
    def q(p):
        return x[min(len(x)-1, max(0, int(math.ceil(p*len(x)))-1))]
    return {"min": x[0], "p50": q(.5), "p90": q(.9), "p99": q(.99), "max": x[-1]}


def _snapshot(masks, nodes, n):
    counts, ts, ss, rs = [], [], [], []
    worst = None
    for endpoint, mask in enumerate(masks):
        idx = list(_iter_bits(mask))
        tis = [int(nodes[k]["tau_index"]) for k in idx]
        sis = [int(nodes[k]["sigma_raw_index"]) for k in idx]
        ris = [int(nodes[k]["R_S_index"]) for k in idx]
        row = {
            "endpoint": endpoint, "ancestors": len(idx),
            "tau_index_range": [min(tis), max(tis)],
            "sigma_raw_index_range": [min(sis), max(sis)],
            "R_S_index_range": [min(ris), max(ris)],
        }
        counts.append(len(idx)); ts.append(max(tis)-min(tis)); ss.append(max(sis)-min(sis)); rs.append(max(ris)-min(ris))
        if worst is None or row["ancestors"] > worst["ancestors"]:
            worst = row
    return {
        "ancestor_count": _quantiles(counts),
        "complete_800_cell_cones": sum(c == n for c in counts),
        "all_endpoints_complete": all(c == n for c in counts),
        "tau_index_span": _quantiles(ts),
        "sigma_raw_index_span": _quantiles(ss),
        "R_S_index_span": _quantiles(rs),
        "full_tau_span_endpoints": sum(x == NODES.TAU_CELLS-1 for x in ts),
        "full_sigma_span_endpoints": sum(x == NODES.SIGMA_RAW_CELLS-1 for x in ss),
        "full_R_S_span_endpoints": sum(x == NODES.RS_CELLS-1 for x in rs),
        "worst_endpoint": worst,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, horizons_s=DEFAULT_HORIZONS):
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text())
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("recent cones must not be trajectory fitted")
    rt = P2.runtime(path)
    nodes = rt["nodes"]
    graph = rt["union_successors"]
    n = len(nodes)
    if n != NODES.EXPECTED_STATES:
        raise RuntimeError("P2 correlation state count changed")
    dt = float(rt["dt"])
    min_gap = min(rt["gaps"])
    requested = sorted(set(map(float, horizons_s)))
    hops_for = {h: int(math.ceil(h/(min_gap*dt)))+2 for h in requested}
    max_hops = max(hops_for.values())

    pred = [0]*n
    for i, outs in enumerate(graph):
        for j in outs:
            pred[int(j)] |= 1 << i
    masks = [1 << i for i in range(n)]
    wanted = {}
    for h,k in hops_for.items(): wanted.setdefault(k, []).append(h)
    rows = {}
    saturated_at = None
    allmask = (1 << n)-1
    for hop in range(1, max_hops+1):
        nxt = []
        for old in masks:
            m = old
            for k in _iter_bits(old):
                m |= pred[k]
                if m == allmask: break
            nxt.append(m)
        masks = nxt
        if saturated_at is None and all(m == allmask for m in masks): saturated_at = hop
        if hop in wanted:
            snap = _snapshot(masks, nodes, n)
            for h in wanted[hop]: rows[f"{h:.6g}"] = {"max_stage_transitions": hop, **snap}

    first_complete = next((h for h in requested if rows[f"{h:.6g}"]["all_endpoints_complete"]), None)
    recent_local = any(not rows[f"{h:.6g}"]["all_endpoints_complete"] for h in requested if h <= 1.0)
    edges = sum(len(x) for x in graph)
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_RECENT_PHYSICAL_PROJECTION_DIAGNOSTIC_FROM_P2_CORRELATION_INTERFACE",
        "source_only": True,
        "diagnostic_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "P2_correlation_interface_version": P2.INTERFACE_VERSION,
        "P2_pair_state_interface_consumed": True,
        "physical_projection_only_adds_paths": True,
        "physical_projection_used_as_P3_covariance_information_input": False,
        "P3_PROMOTED": False, "P4_PROMOTED": False,
        "physical_source_states": n,
        "finite_speed_transition_edges": edges,
        "base_untimed_transition_edges": n*n,
        "clock_dt_binary32_s": dt,
        "minimum_finite_stage_gap_samples": min_gap,
        "arbitrary_window_phase_extra_stage_transitions": 2,
        "hop_bound_is_conservative": True,
        "horizons": rows,
        "all_endpoint_cones_saturate_by_stage_hops": saturated_at,
        "first_requested_horizon_all_endpoints_complete_s": first_complete,
        "recent_le_1s_source_localization_exists": recent_local,
        "next_obligation": "use this only to size P2 path memory; correlated P3 propagation must consume pair-state segment kernels, not these projected cones",
        "failures": [],
    }


def validate(d):
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA: f.append("schema mismatch")
    if d.get("P2_correlation_interface_version") != P2.INTERFACE_VERSION: f.append("P2 interface binding changed")
    for k in ("source_only","diagnostic_only","P2_pair_state_interface_consumed","physical_projection_only_adds_paths","hop_bound_is_conservative"):
        if d.get(k) is not True: f.append(f"{k} is not true")
    for k in ("trajectory_replay_used","filter_changed","physical_projection_used_as_P3_covariance_information_input","P3_PROMOTED","P4_PROMOTED"):
        if d.get(k) is not False: f.append(f"{k} is not false")
    if int(d.get("physical_source_states",0)) != NODES.EXPECTED_STATES: f.append("state count changed")
    if int(d.get("minimum_finite_stage_gap_samples",0)) != 13: f.append("minimum gap changed")
    return list(dict.fromkeys(f))


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); ap.add_argument("--horizons-s",type=float,nargs="*",default=list(DEFAULT_HORIZONS)); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
    d=build(a.domain,tuple(a.horizons_s)); f=validate(d); d["validation_pass"]=not f; d["validation_failures"]=f; a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(d,indent=2,sort_keys=True))
    print(json.dumps({"P2_interface":d["P2_correlation_interface_version"],"saturate_by_hops":d["all_endpoint_cones_saturate_by_stage_hops"],"first_complete_s":d["first_requested_horizon_all_endpoints_complete_s"],"recent_le_1s_localized":d["recent_le_1s_source_localization_exists"],"validation_failures":f},indent=2,sort_keys=True)); return 0 if not f else 2

if __name__ == "__main__": raise SystemExit(main())
