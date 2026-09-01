#!/usr/bin/env python3
"""Bridge one exact P2 source node into the rigorous direct P3 cell backend.

P2 partitions the shipping tuner state into 800 reachable source cells, while
P3 historically scans its own source-uniform grid.  The nonlinear P4 routes
need a common node identity before they can attach start/end covariance and
information factors.  This producer performs that first source-faithful join.

For one P2 node it preserves exactly:

* the node's applied ``tau`` interval, converted to ``x=h/tau`` with outward
  rounding;
* the filter-side committed ``sigma_aw`` interval (not the raw tuner sigma);
* the node's ``R_S`` interval;
* the deployed source schedule, process and vector bounds.

The x interval is refined only as required by the existing validated P3
``split_x_cell`` backend.  Every resulting subcell is evaluated with the direct
generalized-matrix P3 comparison for H=18 and A=21.  The producer then takes a
componentwise Loewner-safe diagonal covariance upper and the minimum certified
relative injection margin across those subcells.

This is deliberately an intermediate bridge.  A diagonal covariance upper is
not ``Sigma_KF(g)`` itself, and it has no attitude/linear cross terms.  Therefore
this output cannot be used as the final P4 metric or to whiten a nonlinear word.
Its purpose is to bind P3 evidence to the exact P2 node indexing so the next
backend can materialize the full source-correlated covariance/information
factor without inventing a second source partition.
"""
from __future__ import annotations

import argparse
import importlib
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_source_node_cells as NODES
from ou3_proof_module_state import preserve_module_bindings

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def _evaluate(mode: str, xparts, sigma: Interval, rs: Interval, live: dict,
              vector: dict, process: dict, sched: dict, alpha6: float, P3D) -> dict:
    rows = []
    for x, rho in xparts:
        row = P3D.mode_cell(mode, x, rho, sigma, rs, live, vector, process, sched, alpha6)
        rows.append((x, rho, row))
    if not rows:
        raise RuntimeError(f"{mode}: exact-node x interval produced no validated P3 subcells")

    dimension = 18 if mode == "H" else 21
    diag = [0.0] * dimension
    deltas = []
    subcells = []
    for x, rho, row in rows:
        upper = list(row["Sigma_diagonal_upper"])
        if len(upper) != dimension:
            raise RuntimeError(f"{mode}: P3 covariance upper dimension changed")
        for i, value in enumerate(upper):
            diag[i] = max(diag[i], float(value))
        delta = float(row["relative_Riccati_injection_margin_lower"])
        if not (math.isfinite(delta) and delta > 0.0):
            raise RuntimeError(f"{mode}: exact-node P3 subcell lost strict relative margin")
        deltas.append(delta)
        subcells.append({
            "x_h_over_tau": x.as_list(),
            "rho_translation_lower": float(rho),
            "relative_Riccati_injection_margin_lower": delta,
            "limiting_block": row["generalized_matrix_inequality"]["limiting_block"],
            "injected_noise_floor_route": row.get("injected_noise_floor_route"),
        })

    return {
        "dimension": dimension,
        "validated_x_subcell_count": len(rows),
        "relative_Riccati_injection_margin_lower": min(deltas),
        "Sigma_diagonal_upper_componentwise_over_node": diag,
        "full_source_node_covariance_matrix_materialized_here": False,
        "attitude_linear_covariance_cross_terms_materialized_here": False,
        "source_node_information_factor_materialized_here": False,
        "subcells": subcells,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_node_index: int = 0) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P2/P3 node bridge must not be trajectory fitted")

    nodes = NODES.build()
    nf = NODES.validate(nodes)
    if nf:
        raise RuntimeError(f"P2 source-node materialization failed: {nf}")
    node = NODES.node(source_node_index, nodes)

    # The direct P3 module intentionally installs a direct mode_cell backend on
    # import.  Keep that historical side effect inside the same broad scoped
    # module-isolation boundary used by the nonlinear P4 routes.
    with preserve_module_bindings():
        P3D = importlib.import_module("ou3_source_reachable_matrix_p3_direct")
        B = P3D.BASE
        live = domain["normal_live"]
        vector = B.VECTOR.build()
        process = B.PROCESS.build()
        words = B.WORDS.build(path)
        prereq = []
        prereq += [f"vector: {x}" for x in B.VECTOR.validate(vector)]
        prereq += [f"process: {x}" for x in B.PROCESS.validate(process)]
        prereq += [f"word: {x}" for x in B.WORDS.validate(words)]
        if prereq:
            raise RuntimeError(f"P3 node bridge prerequisites failed: {prereq}")

        sched = B.source_schedule()
        h = float(sched["dt_s"])
        tau_lo, tau_hi = map(float, node["tau_s"])
        if not (0.0 < tau_lo <= tau_hi):
            raise RuntimeError("P2 node tau interval is invalid")
        x = Interval.outward_bounds(B.down(h / tau_hi), B.up(h / tau_lo))
        sigma = Interval(float(node["sigma_filter_committed_mps2"][0]),
                         float(node["sigma_filter_committed_mps2"][1]))
        rs = Interval(float(node["R_S_filter_std"][0]), float(node["R_S_filter_std"][1]))
        xparts = B.split_x_cell(x)
        alpha6 = B.vector_alpha6(live, vector)
        modes = {
            mode: _evaluate(mode, xparts, sigma, rs, live, vector, process, sched, alpha6, P3D)
            for mode in ("H", "A")
        }

    failures = []
    for mode, dim in (("H", 18), ("A", 21)):
        m = modes[mode]
        if m["dimension"] != dim:
            failures.append(f"{mode}: wrong dimension")
        if m["validated_x_subcell_count"] <= 0:
            failures.append(f"{mode}: no validated x subcells")
        if not m["relative_Riccati_injection_margin_lower"] > 0.0:
            failures.append(f"{mode}: node-conditioned P3 margin is not strict")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_EXACT_P2_NODE_TO_DIRECT_P3_CELL_BRIDGE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "exact_P2_source_node_partition_used": True,
        "P2_source_node_count": int(nodes["partition"]["states"]),
        "source_node_index": int(source_node_index),
        "source_node": node,
        "x_h_over_tau_from_exact_node": x.as_list(),
        "P3_x_subdivision_is_validation_only_not_domain_shrink": True,
        "modes": modes,
        "actual_source_node_Sigma_KF_materialized_here": False,
        "actual_source_node_information_metric_materialized_here": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "materialize the full source-correlated Sigma_KF(g) enclosure/factor on this exact P2 node, retaining attitude-linear cross terms; then extend the same construction to all 800 nodes and reachable g->h edges before Joseph PSD pullback or whole-word whitening"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_EXACT_P2_NODE_TO_DIRECT_P3_CELL_BRIDGE":
        f.append("wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit", "exact_P2_source_node_partition_used",
        "P3_x_subdivision_is_validation_only_not_domain_shrink",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "declared_domain_changed",
        "actual_source_node_Sigma_KF_materialized_here",
        "actual_source_node_information_metric_materialized_here",
        "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("P2_source_node_count") != 800:
        f.append("P2 source-node count is not 800")
    if not 0 <= int(d.get("source_node_index", -1)) < 800:
        f.append("source-node index outside P2 partition")
    for mode, dim in (("H", 18), ("A", 21)):
        m = d.get("modes", {}).get(mode, {})
        if m.get("dimension") != dim:
            f.append(f"{mode}: wrong mode dimension")
        if int(m.get("validated_x_subcell_count", 0)) <= 0:
            f.append(f"{mode}: missing validated x subcells")
        delta = m.get("relative_Riccati_injection_margin_lower")
        if not isinstance(delta, (int, float)) or not math.isfinite(float(delta)) or float(delta) <= 0.0:
            f.append(f"{mode}: invalid node-conditioned P3 margin")
        diag = m.get("Sigma_diagonal_upper_componentwise_over_node", [])
        if len(diag) != dim or any(not (math.isfinite(float(x)) and float(x) > 0.0) for x in diag):
            f.append(f"{mode}: invalid node-conditioned covariance diagonal upper")
        for key in (
            "full_source_node_covariance_matrix_materialized_here",
            "attitude_linear_covariance_cross_terms_materialized_here",
            "source_node_information_factor_materialized_here",
        ):
            if m.get(key) is not False:
                f.append(f"{mode}: {key} is not false")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve(), source_node_index=args.source_node_index)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "validation_pass": not vf,
        "source_node": d["source_node_index"],
        "x_h_over_tau": d["x_h_over_tau_from_exact_node"],
        "H_subcells": d["modes"]["H"]["validated_x_subcell_count"],
        "H_delta": d["modes"]["H"]["relative_Riccati_injection_margin_lower"],
        "A_subcells": d["modes"]["A"]["validated_x_subcell_count"],
        "A_delta": d["modes"]["A"]["relative_Riccati_injection_margin_lower"],
        "actual_metric_materialized": d["actual_source_node_information_metric_materialized_here"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
