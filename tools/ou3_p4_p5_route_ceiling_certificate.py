#!/usr/bin/env python3
"""Route ceiling for the OU-III P4/P5 uniform transported-defect bridge.

P4 certifies an inner level ``W_*`` and P5 must show the P1 handoff family lies
in the P4 strict-decrease domain.  Both use one accounting:

    ||r_word||_M <= B W_0,   B = F N kappa,
    strict decrease requires   sqrt(W_0) <= delta / (2 B),

with ``F`` the prefix bootstrap factor, ``N`` the number of defect-injecting
source operations in the word, ``kappa`` the uniform per-operation metric defect
coefficient and ``delta`` the P3 word endpoint margin.  Sharpening the
*constants* of that accounting has been the whole P5 search so far.  This
certificate asks the prior question: what is the largest attitude capture radius
the accounting itself can ever report, granting every input its most favourable
admissible value?

The certified attitude radius of a level ``W`` is exactly

    theta(W) = sup{ ||c|| : W_g(z) <= W } = sqrt(lambda_max(Sigma_tt) W / s)
             = a_t sqrt(W),      a_t = sqrt(Sigma_tt_upper / s).

Four properties of the accounting are then used, each stated as an explicit
hypothesis and each satisfied by the shipping P4 producer:

  H1  ``delta <= 1``.  The margin is defined by ``Omega_word >= delta Sigma``
      and the covariance contains its own injected noise, so ``Omega <= Sigma``.

  H2  The word must cover the branch in which every attempted vector correction
      is accepted, so ``N`` is at least the number of IMU samples in the word.

  H3  The exact Cayley composition is ``c (+) d = (c+d+0.5 d x c)/(1-d.c/4)``.
      Its homogeneous reference is ``c+d``, so any uniform defect bound is at
      least ``0.5 sup||d|| sup||c||`` -- the cross term is exact, not a slack
      estimate, and is maximal at ``d`` orthogonal to ``c``.

  H4  The uniform accepted-injection bound is at least ``eps a_t sqrt(W)``.  The
      shipping producer has ``eps = 1``: the exact supremum of the attitude part
      of ``K H z`` over the metric ball is ``a_t sqrt(W)``, because
      ``||A^T (A A^T + I)^-1 A|| -> 1`` and ``||E_t^T P^{1/2}|| = sqrt(lambda_max(P_tt))``.

Under H1-H4, an attitude-supported defect costs at least
``sqrt(s lambda_max((Sigma^-1)_tt)) >= sqrt(s/lambda_max(Sigma_tt)) = 1/a_t`` in
the metric, so

    kappa >= (1/a_t) * 0.5 * eps a_t^2 = 0.5 eps a_t,
    B     >= F N 0.5 eps a_t,
    theta_capture = a_t delta/(2B) <= delta / (F N eps).

``a_t`` cancels: the ceiling is independent of how tight the covariance bound
is, of the metric normalization, and of every constant the search has been
sharpening.  With ``delta <= 1`` it depends only on the word's operation count
and on the injection hypothesis.

This is a ceiling on the *proof route*, not a property of the filter.  It says
the uniform transported-defect accounting cannot report the P1 handoff radius,
so P5 cannot close on it and no further sharpening of ``kappa`` can change that.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_nonlinear_word_certificate as P4
import ou3_p5_heading_handoff_contract as HEADING

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1

# Ceiling values of the accounting inputs.  H1: Omega <= Sigma.
MAX_ADMISSIBLE_DELTA = 1.0
# H4 at the shipping producer.  Kept explicit so the ceiling can be re-read at a
# weaker injection hypothesis without editing the derivation.
SHIPPING_INJECTION_FRACTION = 1.0
# Most generous prefix accounting: no overshoot at all (the shipping route pays 4).
MOST_GENEROUS_PREFIX_FACTOR = 1.0
# Exact Cayley cross-term coefficient.  Not adjustable.
CAYLEY_CROSS_COEFF = 0.5


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _ceiling(delta: float, prefix_factor: float, injections: float,
             injection_fraction: float) -> float:
    """theta_capture <= delta / (F N eps).  Rounded upward: it is a ceiling."""
    denom = down(down(prefix_factor * injections) * injection_fraction)
    if not denom > 0.0:
        raise RuntimeError("route ceiling denominator is not positive")
    return up(delta / denom)


def _handoff_rows(heading: dict) -> list[dict]:
    rows = []
    for key, label in (
        ("gauged_quality_handoff", "normal_gauged"),
        ("gauged_timeout_subbranch", "timeout_gauged"),
    ):
        node = heading[key]
        rows.append({
            "handoff_branch": label,
            "full_attitude_cayley_norm_upper": float(node["full_attitude_cayley_norm_upper"]),
        })
    ung = heading.get("ungauged_timeout_subbranch", {})
    rows.append({
        "handoff_branch": "timeout_ungauged_yaw_quotient",
        "full_attitude_cayley_norm_upper": None,
        "full_heading_bound_available": bool(ung.get("full_heading_cayley_bound_available", False)),
    })
    return rows


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    p4 = P4.build(domain_path)
    p4f = P4.validate(p4)
    heading = HEADING.build()
    hf = HEADING.validate(heading)
    failures = [f"P4: {x}" for x in p4f] + [f"heading: {x}" for x in hf]

    handoffs = _handoff_rows(heading)
    bounded = [r for r in handoffs if r["full_attitude_cayley_norm_upper"] is not None]
    largest = max(r["full_attitude_cayley_norm_upper"] for r in bounded)
    smallest = min(r["full_attitude_cayley_norm_upper"] for r in bounded)

    modes = {}
    for mode in ("H", "A"):
        row = p4.get("modes", {}).get(mode)
        if not row:
            failures.append(f"{mode}: P4 mode certificate missing")
            continue
        delta = float(row["word_endpoint_relative_Riccati_injection_margin_lower"])
        B = float(row["transported_word_defect_B_upper"])
        samples = float(row["word_samples_upper"])
        operations = float(row["state_operation_count_upper"])
        a_theta = float(row["metric_consistent_defect_transport"]["attitude_chart_scale"])
        prefix_factor = float(row["prefix_W_factor_upper"])

        # What the shipping certificate reports today.
        sqrt_capture = down(delta / up(2.0 * B))
        theta_now = up(a_theta * sqrt_capture)

        # Ceiling of the same accounting at the shipping word structure and
        # prefix factor, with delta at its maximum admissible value.
        ceiling_shipping = _ceiling(MAX_ADMISSIBLE_DELTA, prefix_factor, samples,
                                    SHIPPING_INJECTION_FRACTION)
        # Ceiling after also granting a perfect prefix accounting.
        ceiling_absolute = _ceiling(MAX_ADMISSIBLE_DELTA, MOST_GENEROUS_PREFIX_FACTOR,
                                    samples, SHIPPING_INJECTION_FRACTION)
        # Ceiling at the certificate's own delta, everything else ideal.
        ceiling_at_source_delta = _ceiling(delta, MOST_GENEROUS_PREFIX_FACTOR, samples,
                                           SHIPPING_INJECTION_FRACTION)

        # What each hypothesis would have to become for the route to reach the
        # largest bounded P1 handoff node.
        breakeven_injections = down(MAX_ADMISSIBLE_DELTA / up(largest * MOST_GENEROUS_PREFIX_FACTOR))
        breakeven_injection_fraction = down(
            MAX_ADMISSIBLE_DELTA / up(largest * MOST_GENEROUS_PREFIX_FACTOR * samples)
        )

        reaches = ceiling_absolute >= largest
        modes[mode] = {
            "mode": mode,
            "word_samples_upper": samples,
            "state_operation_count_upper": operations,
            "prefix_W_factor_upper": prefix_factor,
            "attitude_chart_scale": a_theta,
            "P3_word_endpoint_delta_lower": delta,
            "transported_word_defect_B_upper": B,
            "certified_attitude_capture_radius_now": theta_now,
            "route_ceiling_at_shipping_prefix_factor": ceiling_shipping,
            "route_ceiling_absolute": ceiling_absolute,
            "route_ceiling_at_source_delta": ceiling_at_source_delta,
            "largest_bounded_P1_handoff_cayley_norm": largest,
            "smallest_bounded_P1_handoff_cayley_norm": smallest,
            "shortfall_factor_now_vs_largest_handoff": up(largest / theta_now) if theta_now > 0.0 else math.inf,
            "shortfall_factor_ceiling_vs_largest_handoff": up(largest / ceiling_absolute),
            "shortfall_factor_ceiling_vs_smallest_handoff": up(smallest / ceiling_absolute),
            "breakeven_injecting_operations_per_word": breakeven_injections,
            "breakeven_injection_fraction_at_source_word": breakeven_injection_fraction,
            "route_can_reach_P1_handoff": reaches,
        }
        if reaches:
            failures.append(
                f"{mode}: route ceiling unexpectedly covers the P1 handoff; the "
                "obstruction statement must be re-derived before it is published"
            )

    blocked = bool(modes) and all(m["route_can_reach_P1_handoff"] is False for m in modes.values())
    return {
        "schema": SCHEMA,
        "qualification": "UNIFORM_TRANSPORTED_DEFECT_ROUTE_CEILING_FOR_P4_INNER_LEVEL_AND_P5_CAPTURE",
        "claim": "P5_CANNOT_CLOSE_ON_THE_UNIFORM_TRANSPORTED_DEFECT_ACCOUNTING",
        "source_generated_not_trajectory_fit": True,
        "outward_rounded": True,
        "source_replay_used": False,
        "ceiling_is_about_the_proof_route_not_the_filter": True,
        "ceiling_formula": "theta_capture <= delta / (prefix_factor * injecting_operations * injection_fraction)",
        "attitude_chart_scale_cancels_from_the_ceiling": True,
        "hypotheses": {
            "H1_delta_at_most_one": "Omega_word <= Sigma, so the relative Riccati injection margin is at most 1",
            "H2_injecting_operations_at_least_word_samples": "the word language admits the all-accepted branch, so the uniform bound must cover one accepted vector correction per sample",
            "H3_exact_Cayley_cross_term": "c (+) d = (c+d+0.5 d x c)/(1-d.c/4); 0.5 sup||d|| sup||c|| is exact, not slack",
            "H4_injection_at_least_eps_a_theta_sqrt_W": "sup of the attitude part of K H z over the metric ball is a_t sqrt(W); the shipping producer uses eps=1",
        },
        "P1_handoff_family": handoffs,
        "modes": modes,
        "P4_P5_UNIFORM_TRANSPORT_ROUTE_CEILING": "BELOW_P1_HANDOFF" if blocked else "NOT_ESTABLISHED",
        "required_structural_change": [
            "match each injected defect against the information decrease of the same "
            "operation instead of against the whole-word endpoint margin delta",
            "carry a directional or block margin so the attitude channel is not "
            "rate-limited by the slowest source channel of the word",
            "replace the Cayley quadratic remainder on the attitude channel by an "
            "exact finite-angle sector contraction of the deployed correction",
        ],
        "retired_search_direction": "further sharpening of the uniform per-operation "
                                    "defect coefficient kappa cannot close P5; the ceiling "
                                    "is independent of kappa",
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True or d.get("source_replay_used") is not False:
        failures.append("route ceiling is not source-only")
    if d.get("ceiling_is_about_the_proof_route_not_the_filter") is not True:
        failures.append("route ceiling mis-states its own scope")
    if not d.get("modes"):
        failures.append("route ceiling produced no mode rows")
    for mode, m in d.get("modes", {}).items():
        for key in ("route_ceiling_absolute", "route_ceiling_at_shipping_prefix_factor",
                    "certified_attitude_capture_radius_now"):
            v = m.get(key)
            if not (isinstance(v, (int, float)) and math.isfinite(float(v)) and float(v) > 0.0):
                failures.append(f"{mode}: {key} is not finite positive")
        if m.get("route_can_reach_P1_handoff") is not False:
            failures.append(f"{mode}: route ceiling did not establish the obstruction")
        if not float(m.get("certified_attitude_capture_radius_now", math.inf)) <= float(
            m.get("route_ceiling_absolute", 0.0)
        ):
            failures.append(f"{mode}: reported radius exceeds its own route ceiling")
    if not failures and d.get("P4_P5_UNIFORM_TRANSPORT_ROUTE_CEILING") != "BELOW_P1_HANDOFF":
        failures.append("route ceiling status was not established")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve())
    failures = validate(d)
    out = dict(d)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "P4_P5_UNIFORM_TRANSPORT_ROUTE_CEILING": out["P4_P5_UNIFORM_TRANSPORT_ROUTE_CEILING"],
        "modes": {
            mode: {
                k: out["modes"][mode][k]
                for k in (
                    "certified_attitude_capture_radius_now",
                    "route_ceiling_at_shipping_prefix_factor",
                    "route_ceiling_absolute",
                    "largest_bounded_P1_handoff_cayley_norm",
                    "shortfall_factor_ceiling_vs_largest_handoff",
                    "breakeven_injecting_operations_per_word",
                )
            }
            for mode in out.get("modes", {})
        },
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
