#!/usr/bin/env python3
"""Fast source-indexed translation-metric scan on the finite-speed tuner graph.

The source-complete P3 route must not collapse the 800 physical tuner cells to
one covariance box merely because long-horizon ancestor reachability eventually
becomes complete.  A switched/multiple-metric construction may attach a metric
to each source cell and prove every admissible finite-speed transition from the
start enclosure into the destination enclosure.

This producer is deliberately only a design diagnostic.  It computes, for each
physical source node, the retained finite-memory translation covariance
dominator in [v,p,S,a_w] and the corresponding dimensionless similarity seed
obtained with D_h=diag(sigma*h,sigma*h^2,sigma*h^3,sigma).  It then measures
how these seeds change across the retained 13..26-sample committed-source graph.

Unlike the earlier version of this diagnostic, it does not call the expensive
scaled-process proof inside every H/A source cell.  That process proof depends
primarily on tau and is not needed to answer the metric-switching question.  The
all-cell e3 scan remains the separate process-headroom diagnostic.

No replay, operating-domain shrink, powf/sqrtf target tightening, filter change,
or theorem-gate relaxation is used.  These diagonal seeds are not the final
18/21-state information matrices and cannot promote P3/P4/P5.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_sample_clock_source_refinement as CLOCK
import ou3_p4_source_node_cells as NODES
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2


def _node_seed(node: dict, domain: dict, sched: dict) -> dict:
    h = float(sched["dt_s"])
    tau = Interval(*map(float, node["tau_s"]))
    sigma = Interval(*map(float, node["sigma_filter_committed_mps2"]))
    rs = Interval(*map(float, node["R_S_filter_std"]))
    Tpe = BASE.pos(domain["normal_live"]["vector_pe_recurrence_window_s"], "PE recurrence")
    upper, timing = BASE.translation_upper(tau, sigma, rs, Tpe, sched)
    upper = [float(x) for x in upper]

    # A covariance upper is divided by the smallest scale square to remain an
    # upper bound in the dimensionless similarity coordinates throughout the
    # source cell.  This is only a diagonal metric-design seed.
    scales = [
        sigma.lo * h,
        sigma.lo * h * h,
        sigma.lo * h * h * h,
        sigma.lo,
    ]
    scales2 = [x * x for x in scales]
    normalized = [BASE.up(u / BASE.down(s)) for u, s in zip(upper, scales2)]
    if any(not (math.isfinite(x) and x > 0.0) for x in upper + normalized):
        raise RuntimeError("translation metric seed lost finite positivity")
    return {
        "Sigma_translation_diagonal_upper": upper,
        "Sigma_over_Dh2_diagonal_upper": normalized,
        "word_horizon_s_lower": float(timing["word_horizon_s_lower"]),
        "word_horizon_s_upper": float(timing["word_horizon_s_upper"]),
    }


def _directed_ratio(start: list[float], end: list[float]) -> float:
    """Metric-switch factor M_end <= r M_start for M_i=diag(1/U_i).

    Since M_end/M_start=U_start/U_end coordinatewise, return max U_start/U_end.
    """
    return max(BASE.up(float(a) / BASE.down(float(b))) for a, b in zip(start, end))


def _summary(values: list[float]) -> dict:
    if not values:
        raise RuntimeError("empty switching-factor population")
    x = sorted(float(v) for v in values)

    def q(p: float) -> float:
        i = min(len(x) - 1, max(0, int(math.ceil(p * len(x))) - 1))
        return float(x[i])

    return {
        "count": len(x),
        "min": x[0],
        "p50": q(0.50),
        "p90": q(0.90),
        "p99": q(0.99),
        "max": x[-1],
        "count_le_1p01": sum(v <= 1.01 for v in x),
        "count_le_1p05": sum(v <= 1.05 for v in x),
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

    sched = BASE.source_schedule()
    nodes = nodes_payload["nodes"]
    seeds = [_node_seed(node, domain, sched) for node in nodes]
    graph = graph_payload["graph"]
    if len(graph) != NODES.EXPECTED_STATES:
        raise RuntimeError("finite-speed graph/source-node cardinality mismatch")

    cov_ratios: list[float] = []
    norm_ratios: list[float] = []
    max_partition_jump = {"tau_index": 0, "sigma_raw_index": 0, "R_S_index": 0}
    worst_cov = None
    worst_norm = None
    edge_count = 0

    for i, outs in enumerate(graph):
        ni = nodes[i]
        for j0 in outs:
            j = int(j0)
            nj = nodes[j]
            edge_count += 1
            for key in max_partition_jump:
                max_partition_jump[key] = max(
                    max_partition_jump[key], abs(int(ni[key]) - int(nj[key]))
                )
            cov = _directed_ratio(
                seeds[i]["Sigma_translation_diagonal_upper"],
                seeds[j]["Sigma_translation_diagonal_upper"],
            )
            norm = _directed_ratio(
                seeds[i]["Sigma_over_Dh2_diagonal_upper"],
                seeds[j]["Sigma_over_Dh2_diagonal_upper"],
            )
            cov_ratios.append(cov)
            norm_ratios.append(norm)
            if worst_cov is None or cov > worst_cov["ratio"]:
                worst_cov = {"ratio": cov, "start": i, "end": j}
            if worst_norm is None or norm > worst_norm["ratio"]:
                worst_norm = {"ratio": norm, "start": i, "end": j}

    # Also expose how broad the node seed family itself is.  This is the loss
    # paid by the current global covariance box before any source-transition
    # correlation is retained.
    coordinate_global_spread = []
    coordinate_normalized_spread = []
    for k in range(4):
        vals = [s["Sigma_translation_diagonal_upper"][k] for s in seeds]
        nvals = [s["Sigma_over_Dh2_diagonal_upper"][k] for s in seeds]
        coordinate_global_spread.append(BASE.up(max(vals) / BASE.down(min(vals))))
        coordinate_normalized_spread.append(BASE.up(max(nvals) / BASE.down(min(nvals))))

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_FINITE_SPEED_SOURCE_INDEXED_TRANSLATION_METRIC_DESIGN_SCAN",
        "source_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "diagnostic_only": True,
        "translation_diagonal_seed_is_not_full_information_metric": True,
        "scaled_process_recomputed_per_source_node": False,
        "powf_sqrtf_target_tightening_used": False,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "source_states": len(nodes),
        "finite_speed_transition_edges": edge_count,
        "base_untimed_transition_edges": int(graph_payload["base_transition_edges"]),
        "finite_speed_graph_all_to_all": bool(graph_payload["source_graph_all_to_all"]),
        "clock_gap_samples": [
            int(graph_payload["finite_stage_gap_lower_samples"]),
            int(graph_payload["finite_stage_gap_upper_samples"]),
        ],
        "state_order": ["v", "p", "S", "a_w"],
        "max_partition_index_jump_per_finite_stage": max_partition_jump,
        "global_node_covariance_spread_by_coordinate": coordinate_global_spread,
        "global_node_normalized_spread_by_coordinate": coordinate_normalized_spread,
        "directed_metric_switch_summary": {
            "physical_covariance_seed": _summary(cov_ratios),
            "dimensionless_Dh_seed": _summary(norm_ratios),
        },
        "worst_directed_metric_switch_edges": {
            "physical_covariance_seed": worst_cov,
            "dimensionless_Dh_seed": worst_norm,
        },
        "next_obligation": (
            "construct source-indexed invariant covariance/information enclosures on the staged/committed automaton and prove every finite-speed source edge maps the start enclosure into the destination enclosure; these diagonal translation seeds only measure how much source-switch geometry must be absorbed"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_FINITE_SPEED_SOURCE_INDEXED_TRANSLATION_METRIC_DESIGN_SCAN":
        f.append("wrong qualification")
    for key in (
        "source_only", "diagnostic_only",
        "translation_diagonal_seed_is_not_full_information_metric",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "scaled_process_recomputed_per_source_node", "powf_sqrtf_target_tightening_used",
        "P3_PROMOTED", "P4_PROMOTED", "finite_speed_graph_all_to_all",
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
    if d.get("state_order") != ["v", "p", "S", "a_w"]:
        f.append("translation state order changed")
    for key, row in d.get("directed_metric_switch_summary", {}).items():
        x = row.get("max")
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
        "global_node_covariance_spread_by_coordinate": d["global_node_covariance_spread_by_coordinate"],
        "global_node_normalized_spread_by_coordinate": d["global_node_normalized_spread_by_coordinate"],
        "directed_metric_switch_summary": d["directed_metric_switch_summary"],
        "worst_directed_metric_switch_edges": d["worst_directed_metric_switch_edges"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
