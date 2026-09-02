#!/usr/bin/env python3
"""Diagnostic source-indexed P3 metric scan on the finite-speed tuner graph.

The source-complete P3 route must not collapse the 800 physical tuner cells to
one covariance box merely because long-horizon ancestor reachability eventually
becomes complete.  A switched/multiple-metric certificate may instead attach a
metric to each source cell and prove every *admissible finite-speed transition*
from the start metric into the destination metric.

This producer is deliberately only a design diagnostic.  It builds a static
exact-node P3 comparison at all 800 physical cells, then measures how those
candidate diagonal covariance geometries change across the retained 13..26
sample staged/committed source graph.  The static node comparisons freeze the
source parameters inside the word and therefore are NOT a source-complete P3
certificate.  Likewise, diagonal covariance dominators are metric-design seeds,
not the full Kalman information matrices required for promotion.

The useful outputs are:

* the static exact-node P3 margin distribution, which quantifies the headroom
  available before switched-source costs are paid;
* directed edge ratios between source-indexed covariance seeds;
* the largest tau/sigma/R_S partition jumps actually admitted by the finite
  staging graph.

No replay, operating-domain shrink, filter change, or theorem-gate relaxation is
used here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p3_scaled_process as SCALED
import ou3_p4_sample_clock_source_refinement as CLOCK
import ou3_p4_source_node_cells as NODES
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def _node_seed(mode: str, node: dict, domain: dict, vector: dict, process: dict,
               sched: dict, alpha6: float) -> dict:
    h = float(sched["dt_s"])
    tau_lo, tau_hi = map(float, node["tau_s"])
    x = Interval.outward_bounds(BASE.down(h / tau_hi), BASE.up(h / tau_lo))
    sigma = Interval(*map(float, node["sigma_filter_committed_mps2"]))
    rs = Interval(*map(float, node["R_S_filter_std"]))

    rows = []
    for xcell, rho_x in SCALED.split_x_cell(x):
        rows.append(
            BASE.mode_cell(
                mode,
                xcell,
                rho_x,
                sigma,
                rs,
                domain["normal_live"],
                vector,
                process,
                sched,
                alpha6,
            )
        )
    if not rows:
        raise RuntimeError("static source node produced no x cells")

    dim = len(rows[0]["Sigma_diagonal_upper"])
    upper = [max(float(r["Sigma_diagonal_upper"][k]) for r in rows) for k in range(dim)]
    scales = [min(float(r["comparison_scale_diagonal_squared"][k]) for r in rows) for k in range(dim)]
    if any(not (math.isfinite(x) and x > 0.0) for x in upper + scales):
        raise RuntimeError("static node seed lost finite positivity")
    normalized = [u / s for u, s in zip(upper, scales)]
    margin = min(float(r["relative_Riccati_injection_margin_lower"]) for r in rows)
    return {
        "mode": mode,
        "dimension": dim,
        "static_frozen_source_only": True,
        "P3_PROMOTABLE": False,
        "relative_Riccati_injection_margin_lower": margin,
        "useful_gate": BASE.MIN_USEFUL_DELTA,
        "useful_static_node": margin >= BASE.MIN_USEFUL_DELTA,
        "Sigma_diagonal_upper_seed": upper,
        "comparison_scale_diagonal_squared_lower_seed": scales,
        "Sigma_over_scale_diagonal_upper_seed": normalized,
        "x_subcells": len(rows),
    }


def _directed_ratio(a: list[float], b: list[float]) -> float:
    """max diagonal ratio for W_b/W_a when M_i = diag(1/a_i)."""
    return max(float(x) / float(y) for x, y in zip(a, b))


def _summary(values: list[float]) -> dict:
    if not values:
        raise RuntimeError("empty switching-factor population")
    x = sorted(values)
    def q(p: float) -> float:
        i = min(len(x) - 1, max(0, int(math.ceil(p * len(x))) - 1))
        return float(x[i])
    return {
        "count": len(x),
        "min": float(x[0]),
        "p50": q(0.50),
        "p90": q(0.90),
        "p99": q(0.99),
        "max": float(x[-1]),
        "count_le_1p25": sum(v <= 1.25 for v in x),
        "count_le_2": sum(v <= 2.0 for v in x),
        "count_le_4": sum(v <= 4.0 for v in x),
        "count_le_10": sum(v <= 10.0 for v in x),
        "count_le_100": sum(v <= 100.0 for v in x),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("source-indexed metric scan must not be trajectory fitted")

    nodes_payload = NODES.build()
    nf = NODES.validate(nodes_payload)
    if nf:
        raise RuntimeError(f"source node producer failed: {nf}")
    graph_payload = CLOCK.build(path)
    gf = CLOCK.validate(graph_payload)
    if gf:
        raise RuntimeError(f"finite-speed source graph failed: {gf}")

    vector = BASE.VECTOR.build()
    process = BASE.PROCESS.build()
    sched = BASE.source_schedule()
    alpha6 = BASE.vector_alpha6(domain["normal_live"], vector)

    seeds = {"H": [], "A": []}
    for node in nodes_payload["nodes"]:
        for mode in ("H", "A"):
            seeds[mode].append(
                _node_seed(mode, node, domain, vector, process, sched, alpha6)
            )

    graph = graph_payload["graph"]
    if len(graph) != NODES.EXPECTED_STATES:
        raise RuntimeError("finite-speed graph/source-node cardinality mismatch")

    edge_rows = []
    ratio_pop = {
        "H_covariance_seed": [], "A_covariance_seed": [],
        "H_normalized_seed": [], "A_normalized_seed": [],
    }
    max_partition_jump = {"tau_index": 0, "sigma_raw_index": 0, "R_S_index": 0}
    worst = {k: None for k in ratio_pop}

    nodes = nodes_payload["nodes"]
    for i, outs in enumerate(graph):
        ni = nodes[i]
        for j in outs:
            nj = nodes[j]
            for key in max_partition_jump:
                max_partition_jump[key] = max(
                    max_partition_jump[key], abs(int(ni[key]) - int(nj[key]))
                )
            row = {"start": i, "end": int(j)}
            for mode in ("H", "A"):
                si, sj = seeds[mode][i], seeds[mode][j]
                cov = _directed_ratio(
                    si["Sigma_diagonal_upper_seed"],
                    sj["Sigma_diagonal_upper_seed"],
                )
                norm = _directed_ratio(
                    si["Sigma_over_scale_diagonal_upper_seed"],
                    sj["Sigma_over_scale_diagonal_upper_seed"],
                )
                ck = f"{mode}_covariance_seed"
                nk = f"{mode}_normalized_seed"
                ratio_pop[ck].append(cov)
                ratio_pop[nk].append(norm)
                row[ck] = cov
                row[nk] = norm
                if worst[ck] is None or cov > worst[ck]["ratio"]:
                    worst[ck] = {"ratio": cov, "start": i, "end": int(j)}
                if worst[nk] is None or norm > worst[nk]["ratio"]:
                    worst[nk] = {"ratio": norm, "start": i, "end": int(j)}
            edge_rows.append(row)

    static_summary = {}
    for mode in ("H", "A"):
        margins = [float(s["relative_Riccati_injection_margin_lower"]) for s in seeds[mode]]
        static_summary[mode] = {
            "useful_nodes": sum(m >= BASE.MIN_USEFUL_DELTA for m in margins),
            "worst_margin": min(margins),
            "best_margin": max(margins),
            "margin_distribution": _summary(margins),
        }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_FINITE_SPEED_SOURCE_INDEXED_METRIC_DESIGN_SCAN",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "diagnostic_only": True,
        "static_frozen_source_seed_is_not_source_complete": True,
        "diagonal_seed_is_not_full_information_metric": True,
        "P3_PROMOTED": False,
        "useful_gate": BASE.MIN_USEFUL_DELTA,
        "source_states": len(nodes),
        "finite_speed_transition_edges": len(edge_rows),
        "base_untimed_transition_edges": int(graph_payload["base_transition_edges"]),
        "finite_speed_graph_all_to_all": bool(graph_payload["source_graph_all_to_all"]),
        "clock_gap_samples": [
            int(graph_payload["finite_stage_gap_lower_samples"]),
            int(graph_payload["finite_stage_gap_upper_samples"]),
        ],
        "max_partition_index_jump_per_finite_stage": max_partition_jump,
        "static_node_seed_summary": static_summary,
        "directed_metric_switch_summary": {k: _summary(v) for k, v in ratio_pop.items()},
        "worst_directed_metric_switch_edges": worst,
        "next_obligation": (
            "replace static diagonal seeds by source-indexed invariant covariance/information enclosures and prove every finite-speed source edge maps the start enclosure into the destination enclosure; only then consume the switched metric in P3/P4"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_FINITE_SPEED_SOURCE_INDEXED_METRIC_DESIGN_SCAN":
        f.append("wrong qualification")
    for key in (
        "source_only", "diagnostic_only",
        "static_frozen_source_seed_is_not_source_complete",
        "diagonal_seed_is_not_full_information_metric",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "P3_PROMOTED", "finite_speed_graph_all_to_all",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if int(d.get("source_states", 0)) != NODES.EXPECTED_STATES:
        f.append("source state count changed")
    if int(d.get("finite_speed_transition_edges", 0)) <= 0:
        f.append("finite-speed graph has no edges")
    if int(d.get("finite_speed_transition_edges", 0)) >= int(d.get("base_untimed_transition_edges", 0)):
        f.append("finite-speed graph did not reduce untimed relation")
    if d.get("clock_gap_samples") != [13, 26]:
        f.append("retained finite source clock changed")
    for mode in ("H", "A"):
        s = d.get("static_node_seed_summary", {}).get(mode, {})
        if int(s.get("useful_nodes", -1)) < 0:
            f.append(f"{mode}: missing useful static-node count")
        for key in ("worst_margin", "best_margin"):
            x = s.get(key)
            if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
                f.append(f"{mode}: invalid {key}")
    for key, s in d.get("directed_metric_switch_summary", {}).items():
        x = s.get("max")
        if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
            f.append(f"{key}: invalid switching-factor maximum")
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
        "source_states": d["source_states"],
        "finite_speed_transition_edges": d["finite_speed_transition_edges"],
        "max_partition_index_jump_per_finite_stage": d["max_partition_index_jump_per_finite_stage"],
        "static_node_seed_summary": d["static_node_seed_summary"],
        "directed_metric_switch_summary": d["directed_metric_switch_summary"],
        "worst_directed_metric_switch_edges": d["worst_directed_metric_switch_edges"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
