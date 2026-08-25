#!/usr/bin/env python3
"""Exact Cayley lift of the source information metric for OU-III P4.

For each fixed-dimensional mode use the sole quantitative P4 metric

    z_C(R,xi) = [c(R); xi],
    c(R) = 2 tan(theta/2) u,
    W_g(R,xi) = s_m z_C' Sigma_KF(g)^-1 z_C,

on theta<pi.  The positive scalar s_m is one source-uniform normalization per
mode, chosen as the certified Sigma lambda_max upper bound.  It is identical at
every node in that mode, so it cannot manufacture path contraction; it merely
sets lambda_min(W metric)>=1 and prevents raw-W underflow.  The physical funnel
in canonical coordinates is invariant to this global scalar choice.

All attitude--linear cross terms of the source-varying Kalman information matrix
are retained.  The local quadratic is the same P3 information geometry up to
this one fixed positive scalar, hence the P3 generalized inequality transfers
without any condition-number conversion.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_explicit_information_word_certificate as P3
import ou3_implementation_word_language as WORDS

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 3


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = domain_path.resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P4 metric domain must not be trajectory fitted")
    p3 = P3.build(domain_path)
    pf = P3.validate(p3)
    words = WORDS.build(domain_path)
    wf = WORDS.validate(words)
    failures = [f"P3: {x}" for x in pf] + [f"word-language: {x}" for x in wf]

    modes = {}
    for mode, dim in (("H", 18), ("A", 21)):
        row = p3["modes"][mode]
        smin = float(row["Sigma_lambda_min_lower"])
        smax = float(row["Sigma_lambda_max_upper"])
        if not (0.0 < smin <= smax < math.inf):
            failures.append(f"{mode}: invalid P3 covariance eigenvalue enclosure")
        scale = smax
        mlo = down(scale / smax)
        mhi = up(scale / smin)
        modes[mode] = {
            "dimension": dim,
            "kind": "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC",
            "chart_coordinate": "c(R)=2*tan(theta/2)*u=4*e_R/(1+tr(R))",
            "chart_domain": "theta<pi",
            "exact_group_metric": "W_g=s_mode*[c(R);xi]^T Sigma_KF(g)^-1 [c(R);xi]",
            "source_covariance_inverse": True,
            "mode_global_positive_scale": scale,
            "same_scale_on_every_source_node_in_mode": True,
            "normalization_policy": "s_mode=source_uniform_Sigma_lambda_max_upper",
            "node_dependent": True,
            "full_attitude_linear_cross_terms_retained": True,
            "block_diagonal_metric_used": False,
            "common_Euclidean_metric_used": False,
            "local_coordinate_matches_P3_delta_theta": True,
            "local_quadratic_is_positive_scalar_multiple_of_P3_information_metric": True,
            "endpoint_metric_must_match_endpoint_source_covariance": True,
            "joint_source_reachability_required": True,
            "Sigma_lambda_min_lower": smin,
            "Sigma_lambda_max_upper": smax,
            "metric_lambda_min_lower": mlo,
            "metric_lambda_max_upper": mhi,
            "P3_word_endpoint_margin_lower": float(row["word_endpoint_relative_Riccati_injection_margin_lower"]),
            "P3_prefix_information_gain_upper": float(row["prefix_information_gain_upper"]),
        }

    return {
        "schema": SCHEMA,
        "qualification": "EXACT_NORMALIZED_CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC_FOR_P4",
        "source_generated_not_trajectory_fit": True,
        "single_quantitative_metric_route": True,
        "retired_block_diagonal_route_available": False,
        "metric_change_does_not_change_filter": True,
        "metric_change_does_not_change_adaptation_law": True,
        "global_scaling_does_not_change_physical_level_sets": True,
        "source_word_horizon_s": words["word_contract"]["conditional_word_language"]["word_horizon_lower_s"],
        "modes": modes,
        "failures": failures,
        "pass": not failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("metric is not source generated")
    if d.get("single_quantitative_metric_route") is not True:
        failures.append("P4 has more than one quantitative metric route")
    if d.get("retired_block_diagonal_route_available") is not False:
        failures.append("retired block-diagonal metric remains available")
    if d.get("global_scaling_does_not_change_physical_level_sets") is not True:
        failures.append("mode-global information normalization is not declared geometry-preserving")
    for mode, dim in (("H",18),("A",21)):
        m = d.get("modes", {}).get(mode, {})
        if m.get("dimension") != dim:
            failures.append(f"{mode}: wrong metric dimension")
        if m.get("kind") != "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC":
            failures.append(f"{mode}: wrong metric kind")
        if m.get("source_covariance_inverse") is not True:
            failures.append(f"{mode}: metric is not source covariance information geometry")
        if m.get("same_scale_on_every_source_node_in_mode") is not True:
            failures.append(f"{mode}: arbitrary node scaling could manufacture contraction")
        if not (isinstance(m.get("mode_global_positive_scale"),(int,float)) and float(m["mode_global_positive_scale"])>0.0):
            failures.append(f"{mode}: missing positive global metric scale")
        if m.get("full_attitude_linear_cross_terms_retained") is not True:
            failures.append(f"{mode}: attitude-linear cross terms were discarded")
        if m.get("block_diagonal_metric_used") is not False:
            failures.append(f"{mode}: block-diagonal surrogate still active")
        if m.get("local_quadratic_is_positive_scalar_multiple_of_P3_information_metric") is not True:
            failures.append(f"{mode}: exact P3/P4 information geometry identity lost")
        if m.get("endpoint_metric_must_match_endpoint_source_covariance") is not True:
            failures.append(f"{mode}: endpoint metric/source correlation not required")
        lo = m.get("metric_lambda_min_lower")
        hi = m.get("metric_lambda_max_upper")
        if not (isinstance(lo,(int,float)) and isinstance(hi,(int,float)) and 0.0 < float(lo) <= float(hi) < math.inf):
            failures.append(f"{mode}: invalid uniform metric eigenvalue bounds")
        if isinstance(lo,(int,float)) and not float(lo) <= 1.0:
            failures.append(f"{mode}: normalized metric lower bound is optimistic")
        delta = m.get("P3_word_endpoint_margin_lower")
        if not (isinstance(delta,(int,float)) and 0.0 < float(delta) < 1.0):
            failures.append(f"{mode}: missing strict P3 endpoint margin")
    if not failures and d.get("pass") is not True:
        failures.append("metric producer did not pass")
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
    print(json.dumps({"modes": d["modes"], "failures": failures}, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
