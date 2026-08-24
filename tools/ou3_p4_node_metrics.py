#!/usr/bin/env python3
"""Source-node group-compatible metric construction for OU-III P4.

The nonlinear theorem metric is deliberately not the Kalman inverse covariance.
For each compact (tau,sigma_aw) source node it uses the OU similarity scales

    a_w ~ sigma,
    v   ~ sigma tau,
    p   ~ sigma tau^2,
    S   ~ sigma tau^3,

with isotropic 3x3 blocks.  Gyro-bias and active accelerometer-bias blocks use
the declared theorem error scales.  The attitude scale is balanced against the
accelerometer geometry, theta_scale ~ sigma/f_ref, and a small global multiplier
family is exposed for validated proof-design optimization.  A single multiplier
must be selected for the complete graph; node-by-node arbitrary rescaling is
not permitted because that could manufacture contraction around cycles.

Every resulting local quadratic is exactly of the manuscript form

    Pbar_i = blkdiag((a_R_i/2) I3, P_xi_i),

and therefore lifts to W_i=a_R_i(1-cos(theta))+xi'P_xi_i xi.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_source_reachable_matrix_p3 as P3

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
ATTITUDE_BALANCE_MULTIPLIERS = (0.5, 1.0, 2.0, 4.0)
TAU_NODE_COUNT = 4
SIGMA_NODE_COUNT = 4


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def geom_edges(lo: float, hi: float, n: int) -> list[float]:
    if not (0.0 < lo <= hi) or n < 1:
        raise ValueError("invalid geometric metric partition")
    if lo == hi:
        return [lo, hi]
    r = (hi / lo) ** (1.0 / n)
    out = [lo]
    for _ in range(n - 1):
        out.append(out[-1] * r)
    out.append(hi)
    return out


def gmid(lo: float, hi: float) -> float:
    return math.sqrt(float(lo) * float(hi))


def block3(label: str, scale: float) -> dict:
    if not (math.isfinite(scale) and scale > 0.0):
        raise ValueError(f"invalid {label} metric scale")
    w = 1.0 / (scale * scale)
    return {"label": label, "physical_scale": scale, "diagonal_weight": w, "multiplicity": 3}


def node_metric(mode: str, tau: float, sigma: float, f_ref: float,
                bg_scale: float, ba_scale: float, attitude_mult: float) -> dict:
    theta_scale = attitude_mult * sigma / f_ref
    # Avoid a proof-design singularity at a vanishing sea-amplitude floor.  The
    # floor is not a state-domain assumption; it only prevents the metric weight
    # from diverging.  The value is below one tenth degree and is substantially
    # tighter than any chart eventually promoted by P4.
    theta_scale = max(theta_scale, 1.0e-3)
    a_R = 2.0 / (theta_scale * theta_scale)
    blocks = [
        block3("b_g", bg_scale),
        block3("v", sigma * tau),
        block3("p", sigma * tau * tau),
        block3("S", sigma * tau * tau * tau),
        block3("a_w", sigma),
    ]
    if mode == "A":
        blocks.append(block3("b_a", ba_scale))
    weights = [a_R / 2.0] * 3
    for b in blocks:
        weights.extend([b["diagonal_weight"]] * 3)
    return {
        "kind": "GROUP_COMPATIBLE_NODE_METRIC",
        "node_dependent": True,
        "attitude_block_isotropic": True,
        "attitude_linear_cross_terms": False,
        "common_Euclidean_metric": False,
        "equals_Kalman_inverse_covariance": False,
        "a_R": a_R,
        "theta_physical_scale_rad": theta_scale,
        "P_xi_blocks": blocks,
        "Pbar_diagonal": weights,
        "P_xi_lambda_min": min(weights[3:]),
        "P_xi_lambda_max": max(weights[3:]),
        "Pbar_lambda_min": min(weights),
        "Pbar_lambda_max": max(weights),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P4 metric domain must not be trajectory fitted")
    live = domain["normal_live"]
    startup = domain["startup"]["physical_handoff_coordinate_bounds"]
    f_ref = float(live["specific_force_norm_lower_mps2"])
    bg_scale = float(startup["gyro_bias_error_norm_upper_rad_s"])
    ba_scale = float(startup["accelerometer_bias_error_norm_upper_mps2"])
    sched = P3.source_schedule()
    tau_lo, tau_hi = map(float, sched["tau_applied_invariant_s"])
    sigma_lo, sigma_hi = map(float, sched["sigma_aw_applied_safety"])
    te = geom_edges(tau_lo, tau_hi, TAU_NODE_COUNT)
    se = geom_edges(sigma_lo, sigma_hi, SIGMA_NODE_COUNT)

    candidates = {}
    for mult in ATTITUDE_BALANCE_MULTIPLIERS:
        nodes = {"H": [], "A": []}
        for ti in range(TAU_NODE_COUNT):
            for si in range(SIGMA_NODE_COUNT):
                t0, t1 = te[ti], te[ti+1]
                s0, s1 = se[si], se[si+1]
                tau = gmid(t0, t1)
                sigma = gmid(s0, s1)
                nid = f"t{ti}_s{si}"
                for mode in ("H", "A"):
                    nodes[mode].append({
                        "id": f"{mode}_{nid}",
                        "tau_cell_s": [down(t0), up(t1)],
                        "sigma_aw_cell_mps2": [down(s0), up(s1)],
                        "metric_anchor_tau_s": tau,
                        "metric_anchor_sigma_aw_mps2": sigma,
                        "metric": node_metric(mode, tau, sigma, f_ref, bg_scale, ba_scale, mult),
                    })
        candidates[f"attitude_balance_{mult:g}"] = {
            "global_attitude_balance_multiplier": mult,
            "same_multiplier_on_every_graph_node": True,
            "nodes": nodes,
        }

    return {
        "schema": SCHEMA,
        "qualification": "SOURCE_NODE_GROUP_COMPATIBLE_P4_METRIC_CANDIDATES",
        "source_generated_not_trajectory_fit": True,
        "metric_formula": "OU_SIMILARITY_BLOCKS_WITH_SOURCE_NODE_TAU_SIGMA_AND_GROUP_ATTITUDE_BLOCK",
        "tau_partition": [down(x) if i == 0 else up(x) if i == len(te)-1 else x for i,x in enumerate(te)],
        "sigma_partition": [down(x) if i == 0 else up(x) if i == len(se)-1 else x for i,x in enumerate(se)],
        "candidate_policy": "validated P4 endpoint enclosure selects one global attitude-balance multiplier; no per-node arbitrary scalar rescaling",
        "candidates": candidates,
        "selected_candidate": None,
        "selection_requires_validated_endpoint_word_bound": True,
        "pass": True,
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("metric producer is trajectory fitted")
    candidates = d.get("candidates", {})
    if not candidates:
        failures.append("no metric candidates")
    for cname, c in candidates.items():
        if c.get("same_multiplier_on_every_graph_node") is not True:
            failures.append(f"{cname}: arbitrary node rescaling permitted")
        for mode, dim in (("H",18),("A",21)):
            nodes = c.get("nodes", {}).get(mode, [])
            if len(nodes) != TAU_NODE_COUNT * SIGMA_NODE_COUNT:
                failures.append(f"{cname}.{mode}: source-node partition incomplete")
                continue
            for node in nodes:
                m = node.get("metric", {})
                if m.get("kind") != "GROUP_COMPATIBLE_NODE_METRIC":
                    failures.append(f"{cname}.{node.get('id')}: wrong metric kind")
                if m.get("attitude_linear_cross_terms") is not False:
                    failures.append(f"{cname}.{node.get('id')}: attitude-linear cross term present")
                if m.get("equals_Kalman_inverse_covariance") is not False:
                    failures.append(f"{cname}.{node.get('id')}: metric aliases inverse covariance")
                diag = m.get("Pbar_diagonal", [])
                if len(diag) != dim or any(not (isinstance(x,(int,float)) and math.isfinite(float(x)) and float(x)>0.0) for x in diag):
                    failures.append(f"{cname}.{node.get('id')}: invalid Pbar diagonal")
    if d.get("selection_requires_validated_endpoint_word_bound") is not True:
        failures.append("metric selection can occur without validated endpoint word bound")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve())
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"candidate_count": len(d["candidates"]), "nodes_per_mode": TAU_NODE_COUNT*SIGMA_NODE_COUNT, "failures": failures}, indent=2))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
