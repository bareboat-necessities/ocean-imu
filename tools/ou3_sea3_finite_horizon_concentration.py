#!/usr/bin/env python3
"""Finite-horizon Gaussian concentration kernel for the SEA3 -> P1 bridge.

A non-degenerate Gaussian JONSWAP model has unbounded support, so no finite
hard P1 cap can contain *every* realization deterministically.  The declared
operating domain already carries a finite-horizon stochastic failure budget.
This module turns that budget into a rigorous covariance-trace threshold that a
future coupled sea/RAO spectral certificate can target.

The bound is elementary and does not require temporal independence.  For a
centered d-dimensional Gaussian vector X with

    tr Cov[X] <= v,

all component variances are <= v.  If ||X||_2 > A, at least one component has
absolute value > A/sqrt(d).  The one-dimensional Gaussian Chernoff bound and a
coordinate union bound therefore give

    P(||X||_2 > A) <= 2 d exp(-A^2 / (2 d v)).

A union bound across N sampled instants then gives

    P(max_k ||X_k||_2 > A) <= 2 d N exp(-A^2 / (2 d v)),

regardless of correlation between samples.  We choose an integer exponent t
with validated arithmetic such that 2 d N exp(-t) is below the allocated event
budget, then it is sufficient to require

    v <= A^2 / (2 d t).

The exponential is evaluated as repeated interval products of the repository's
validated exp(-1/2), so ordinary libm transcendental calls are not used in the
proof decision.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, down, up
import ou3_sea3_finite_window_response_admission as ADMIT
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
DEFAULT_RESPONSE_DOMAIN = REPO / "tools" / "ou3_sea3_directional_response_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_FINITE_HORIZON_GAUSSIAN_P1_CONCENTRATION_V1"
DIMENSION = 3


def exp_minus_integer(t: int) -> Interval:
    """Validated exp(-t), t a nonnegative integer, from exp(-1/2) products."""
    if not isinstance(t, int) or isinstance(t, bool) or t < 0:
        raise ValueError("t must be a nonnegative integer")
    base = VT.exp_point(-0.5)
    out = Interval.point(1.0)
    for _ in range(2 * t):
        out = out * base
    return out


def select_tail_exponent(samples: int, event_budget: float, dimension: int = DIMENSION) -> dict:
    """Choose the smallest integer t with 2*d*N*exp(-t) <= budget."""
    if not isinstance(samples, int) or isinstance(samples, bool) or samples <= 0:
        raise ValueError("samples must be a positive integer")
    if not isinstance(dimension, int) or isinstance(dimension, bool) or dimension <= 0:
        raise ValueError("dimension must be a positive integer")
    budget = float(event_budget)
    if not (math.isfinite(budget) and 0.0 < budget < 1.0):
        raise ValueError("event budget must lie strictly between zero and one")
    budget_lo = down(budget)

    multiplier = 2 * dimension * samples
    for t in range(1, 257):
        exp_iv = exp_minus_integer(t)
        union_upper = up(float(multiplier) * exp_iv.hi)
        if union_upper <= budget_lo:
            return {
                "integer_tail_exponent": t,
                "validated_exp_minus_t": exp_iv.as_list(),
                "union_failure_probability_upper": union_upper,
                "allocated_event_budget_lower": budget_lo,
                "minimal_integer_exponent": True,
            }
    raise RuntimeError("failed to find finite-horizon tail exponent below 256")


def trace_variance_threshold(norm_cap: float, tail_exponent: int, dimension: int = DIMENSION) -> float:
    """Outward-lower sufficient trace-covariance threshold A^2/(2*d*t)."""
    cap = float(norm_cap)
    if not (math.isfinite(cap) and cap > 0.0):
        raise ValueError("norm cap must be finite and positive")
    if not isinstance(tail_exponent, int) or isinstance(tail_exponent, bool) or tail_exponent <= 0:
        raise ValueError("tail exponent must be a positive integer")
    denom = 2 * dimension * tail_exponent
    cap2_lo = down(cap * cap)
    return down(cap2_lo / float(denom))


def build(
    samples: int,
    domain_path: Path = DEFAULT_DOMAIN,
    response_domain_path: Path = DEFAULT_RESPONSE_DOMAIN,
    repo: Path = REPO,
) -> dict:
    """Build thresholds for one explicitly finite sampled horizon."""
    domain = json.loads(Path(domain_path).read_text(encoding="utf-8"))
    total_budget = float(domain["stochastic"]["finite_horizon_failure_probability_budget"])
    if not (0.0 < total_budget < 1.0):
        raise RuntimeError("invalid finite-horizon failure-probability budget")

    admission = ADMIT.build_contract(domain_path, response_domain_path, repo)
    failures = ADMIT.validate(admission)
    if failures:
        raise RuntimeError(f"finite-window admission contract invalid: {failures}")

    # Equal allocation is deliberately simple and auditable.  It can later be
    # optimized without changing the concentration theorem itself.
    per_event_budget = down(total_budget / 2.0)
    acc_tail = select_tail_exponent(samples, per_event_budget)
    rate_tail = select_tail_exponent(samples, per_event_budget)
    caps = admission["normal_live_caps"]
    acc_v = trace_variance_threshold(
        caps["non_gravitational_cog_acceleration_norm_upper_mps2"],
        acc_tail["integer_tail_exponent"],
    )
    rate_v = trace_variance_threshold(
        caps["body_rate_norm_upper_deg_s"],
        rate_tail["integer_tail_exponent"],
    )
    combined_failure_upper = up(
        acc_tail["union_failure_probability_upper"]
        + rate_tail["union_failure_probability_upper"]
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "samples": samples,
        "dimension": DIMENSION,
        "centered_Gaussian_response_required": True,
        "temporal_independence_required": False,
        "cross_axis_independence_required": False,
        "trajectory_replay_used": False,
        "ordinary_libm_transcendental_used_in_proof_decision": False,
        "finite_horizon_failure_probability_budget": total_budget,
        "budget_allocation": {
            "acceleration": per_event_budget,
            "body_rate": per_event_budget,
            "sum_upper": up(per_event_budget + per_event_budget),
        },
        "acceleration": {
            "norm_cap_mps2": caps["non_gravitational_cog_acceleration_norm_upper_mps2"],
            "tail": acc_tail,
            "required_trace_covariance_upper_m2_s4": acc_v,
            "required_trace_RMS_upper_mps2": down(math.sqrt(acc_v)),
        },
        "body_rate": {
            "norm_cap_deg_s": caps["body_rate_norm_upper_deg_s"],
            "tail": rate_tail,
            "required_trace_covariance_upper_deg2_s2": rate_v,
            "required_trace_RMS_upper_deg_s": down(math.sqrt(rate_v)),
        },
        "combined_failure_probability_upper": combined_failure_upper,
        "combined_failure_within_declared_budget": combined_failure_upper <= down(total_budget),
        "response_parameter_box_sha256": admission["response_parameter_box_sha256"],
        "covariance_trace_producer_attached": False,
        "finite_horizon_good_event_promoted": False,
        "deterministic_left_inclusion_closed": False,
        "infinite_horizon_Gaussian_hard_bound_claimed": False,
        "P2_pruning_promoted": False,
        "P3_promoted": False,
        "P4_promoted": False,
        "P5_promoted": False,
        "next_obligation": (
            "derive validated uniform acceleration and body-rate covariance-trace enclosures on the coupled JONSWAP/RAO domain and compare them to these finite-horizon thresholds"
        ),
    }


def evaluate_covariance_candidate(
    certificate: dict,
    *,
    acceleration_trace_covariance_upper_m2_s4: float,
    body_rate_trace_covariance_upper_deg2_s2: float,
    validated_covariance_trace_enclosures: bool,
    response_parameter_box_sha256: str,
) -> dict:
    """Evaluate validated spectral covariance bounds against the kernel."""
    failures = validate(certificate)
    if failures:
        raise ValueError(f"invalid concentration certificate: {failures}")

    reasons: list[str] = []
    if validated_covariance_trace_enclosures is not True:
        reasons.append("covariance-trace bounds are not validated outward enclosures")
    if response_parameter_box_sha256 != certificate["response_parameter_box_sha256"]:
        reasons.append("covariance candidate is not bound to the certified RAO parameter box")

    acc = float(acceleration_trace_covariance_upper_m2_s4)
    rate = float(body_rate_trace_covariance_upper_deg2_s2)
    if not (math.isfinite(acc) and acc >= 0.0):
        reasons.append("acceleration covariance trace must be finite and nonnegative")
    elif acc > float(certificate["acceleration"]["required_trace_covariance_upper_m2_s4"]):
        reasons.append("acceleration covariance trace exceeds finite-horizon concentration threshold")
    if not (math.isfinite(rate) and rate >= 0.0):
        reasons.append("body-rate covariance trace must be finite and nonnegative")
    elif rate > float(certificate["body_rate"]["required_trace_covariance_upper_deg2_s2"]):
        reasons.append("body-rate covariance trace exceeds finite-horizon concentration threshold")

    passed = not reasons
    return {
        "finite_horizon_good_event_candidate_pass": passed,
        "decision": "PASS_CANDIDATE" if passed else "FAIL_CANDIDATE",
        "validation_failures": reasons,
        "finite_horizon_good_event_probability_lower": (
            down(1.0 - float(certificate["combined_failure_probability_upper"])) if passed else 0.0
        ),
        "deterministic_left_inclusion_promoted": False,
        "infinite_horizon_claim_promoted": False,
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        failures.append("schema/qualification mismatch")
    if not isinstance(d.get("samples"), int) or d.get("samples", 0) <= 0:
        failures.append("invalid finite horizon")
    if d.get("dimension") != DIMENSION:
        failures.append("response dimension changed")
    if d.get("centered_Gaussian_response_required") is not True:
        failures.append("centered-Gaussian hypothesis disappeared")
    if d.get("temporal_independence_required") is not False:
        failures.append("temporal independence was incorrectly introduced")
    if d.get("cross_axis_independence_required") is not False:
        failures.append("cross-axis independence was incorrectly introduced")
    if d.get("trajectory_replay_used") is not False:
        failures.append("trajectory replay entered concentration proof")
    if d.get("ordinary_libm_transcendental_used_in_proof_decision") is not False:
        failures.append("unvalidated libm transcendental entered proof decision")
    if d.get("combined_failure_within_declared_budget") is not True:
        failures.append("combined concentration failure exceeds declared stochastic budget")
    for key in (
        "covariance_trace_producer_attached",
        "finite_horizon_good_event_promoted",
        "deterministic_left_inclusion_closed",
        "infinite_horizon_Gaussian_hard_bound_claimed",
        "P2_pruning_promoted",
        "P3_promoted",
        "P4_promoted",
        "P5_promoted",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} was promoted")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, required=True)
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--response-domain", type=Path, default=DEFAULT_RESPONSE_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    d = build(args.samples, args.domain, args.response_domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(
        json.dumps(
            {
                "samples": d["samples"],
                "failure_budget": d["finite_horizon_failure_probability_budget"],
                "acceleration_trace_covariance_threshold": d["acceleration"][
                    "required_trace_covariance_upper_m2_s4"
                ],
                "body_rate_trace_covariance_threshold": d["body_rate"][
                    "required_trace_covariance_upper_deg2_s2"
                ],
                "combined_failure_probability_upper": d["combined_failure_probability_upper"],
                "validation_failures": failures,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
