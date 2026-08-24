#!/usr/bin/env python3
"""Independent final composition gate for the OU-III deployment theorem.

This gate deliberately recomputes stochastic concentration and finite capture
from primitive validated constants. It never trusts a supplied final failure
probability, capture time, deployment PASS bit, or a self-asserted source-domain
completeness flag. Hybrid closure is independently normalized against the
current source-domain obligations, with periodic a_w covariance synchronization
discharged by its source-bound PSD/Loewner proof rather than a replay margin.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_hybrid_contract as HYBRID
import ou3_source_domain_contract as SOURCE_DOMAIN

REQUIRED_HYBRID = set(SOURCE_DOMAIN.HYBRID_OBLIGATIONS)
SOURCE_DOMAIN_KEYS = (
    "schema",
    "claim",
    "source_generated_not_trajectory_fit",
    "source_complete_parameter_domain",
    "validated_arithmetic",
    "outward_rounded",
    "implementation_header",
    "continuous_parameters",
    "timing_constants_s",
    "discrete_source_branches",
    "hybrid_obligations",
    "periodic_aw_covariance_sync_proof",
)


def pos(x) -> float:
    y = float(x)
    if not math.isfinite(y) or not y > 0.0:
        raise ValueError(f"expected finite positive value, got {x!r}")
    return y


def nonneg(x) -> float:
    y = float(x)
    if not math.isfinite(y) or y < 0.0:
        raise ValueError(f"expected finite nonnegative value, got {x!r}")
    return y


def t_star_for_radius(radius2: float, m: float, v: float, b: float) -> float:
    if radius2 <= m:
        return 0.0
    # Solve m + 2 sqrt(v t) + 2 b t = radius2 by monotone bisection.
    lo, hi = 0.0, 1.0

    def f(t: float) -> float:
        return m + 2.0 * math.sqrt(max(0.0, v * t)) + 2.0 * b * t

    while f(hi) < radius2:
        hi *= 2.0
        if hi > 1.0e12:
            break
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if f(mid) <= radius2:
            lo = mid
        else:
            hi = mid
    return lo


def validate_source_domain(payload: dict) -> dict:
    """Bind promotion to the contract regenerated from the current source tree."""
    expected = SOURCE_DOMAIN.build(SOURCE_DOMAIN.DEFAULT_HEADER.resolve())
    failures = []
    for key in SOURCE_DOMAIN_KEYS:
        if payload.get(key) != expected.get(key):
            failures.append(f"source-domain field {key!r} does not match current implementation")
    return {
        "pass": not failures,
        "failures": failures,
        "implementation_header": expected["implementation_header"],
        "continuous_parameters": expected["continuous_parameters"],
        "hybrid_obligations": expected["hybrid_obligations"],
    }


def derive_stochastic(s: dict) -> dict:
    lam = pos(s["lambda_W_upper"])
    if not lam < 1.0:
        raise ValueError("lambda_W_upper must be < 1")
    ell = int(s["word_length_samples_upper"])
    if ell < 1:
        raise ValueError("word_length_samples_upper must be >= 1")
    Lx = nonneg(s["L_X_upper"])
    G = nonneg(s["G_bar_upper"])
    czw = nonneg(s["c_zw_upper"])
    rstar = nonneg(s["r_star_upper"])
    cww = nonneg(s["c_ww_upper"])
    s2 = nonneg(s["s2_upper"])
    s4 = nonneg(s["s4_upper"])
    gW = nonneg(s["g_W_upper"])
    hW = nonneg(s["h_W_upper"])
    mSigma = nonneg(s["Sigma_trace_upper"])
    vSigma = nonneg(s["Sigma_trace_square_upper"])
    bSigma = nonneg(s["Sigma_norm_upper"])
    wstar = pos(s["localization_radius_lower"])
    N = int(s["word_horizon_count"])
    if N < 1:
        raise ValueError("word_horizon_count must be >= 1")
    bW = nonneg(s["b_W_upper"])
    vW = nonneg(s["v_W_upper"])
    a = pos(s["funnel_level_a_lower"])
    W0 = nonneg(s["initial_W_upper"])
    budget = pos(s["failure_probability_budget"])
    if not budget < 1.0:
        raise ValueError("failure_probability_budget must be < 1")

    gstar = G + czw * rstar
    nu1 = 2.0 * gstar * gstar * s2 + 2.0 * cww * cww * s4
    if math.isclose(Lx, 1.0, rel_tol=0.0, abs_tol=1.0e-14):
        geom = float(ell)
    else:
        geom = sum(Lx ** (2 * r) for r in range(ell))
    nuW = ell * nu1 * geom
    lambda_s = 0.5 * (1.0 + lam)
    sigma_s2 = (hW + gW * gW * lam / (2.0 * (1.0 - lam))) * nuW

    floor = sigma_s2 / (1.0 - lambda_s)
    envelope_max = max(W0, floor)
    xstar = a - envelope_max
    VW = vW / (1.0 - lambda_s * lambda_s)
    excursion = 1.0
    if xstar > 0.0:
        denom = 2.0 * (VW + bW * xstar / 3.0)
        excursion = 0.0 if denom == 0.0 else min(
            1.0, N * math.exp(-(xstar * xstar) / denom)
        )

    tstar = t_star_for_radius(wstar * wstar, mSigma, vSigma, bSigma)
    localization = min(1.0, N * ell * math.exp(-tstar))
    total = min(1.0, excursion + localization)

    return {
        "lambda_s": lambda_s,
        "nu_W_upper": nuW,
        "sigma_s_squared_upper": sigma_s2,
        "mean_square_floor_upper": floor,
        "drift_envelope_max_upper": envelope_max,
        "excursion_margin_lower": xstar,
        "V_W_upper": VW,
        "gaussian_t_star_lower": tstar,
        "localization_failure_probability_upper": localization,
        "excursion_failure_probability_upper": excursion,
        "finite_horizon_failure_probability_upper": total,
        "failure_probability_budget": budget,
        "pass": bool(xstar > 0.0 and total < 1.0 and total <= budget),
    }


def derive_capture(c: dict) -> dict:
    lam = pos(c["lambda_upper"])
    if not lam < 1.0:
        raise ValueError("capture lambda_upper must be < 1")
    gamma = nonneg(c["gamma_upper"])
    c0 = nonneg(c["initial_level_upper"])
    bstar = gamma / (1.0 - lam)
    beta = pos(c.get("strict_superlevel_factor", 1.001))
    if not beta > 1.0:
        raise ValueError("strict_superlevel_factor must be > 1")
    beta_level = bstar * beta + pos(c.get("strict_superlevel_absolute", 1.0e-12))
    if c0 <= beta_level:
        words = 0
    elif c0 <= bstar:
        words = 0
    else:
        ratio = (beta_level - bstar) / (c0 - bstar)
        words = max(0, int(math.ceil(math.log(ratio) / math.log(lam))))
    word_horizon_s = pos(c["word_horizon_s_upper"])
    return {
        "asymptotic_level_b_star_upper": bstar,
        "strict_capture_level_b_eta": beta_level,
        "capture_words_upper": words,
        "capture_time_s_upper": words * word_horizon_s,
        "pass": math.isfinite(bstar) and beta_level > bstar and words >= 0,
    }


def compose(check: dict, source_domain_payload: dict, primitive: dict) -> dict:
    source_domain = validate_source_domain(source_domain_payload)
    hybrid = HYBRID.validate(check)
    modes = check.get("modes", {})
    continuous_pass = all(
        modes.get(m, {}).get("linear_pass") and modes.get(m, {}).get("nonlinear_pass")
        for m in ("H", "A")
    )
    provenance_pass = bool(check.get("validated_enclosure_provenance_pass"))

    try:
        stochastic = derive_stochastic(primitive["stochastic_primitives"])
        capture = {
            mode: derive_capture(primitive["capture_primitives"][mode])
            for mode in ("H", "A")
        }
        arithmetic_error = None
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        stochastic = {"pass": False}
        capture = {"H": {"pass": False}, "A": {"pass": False}}
        arithmetic_error = str(exc)

    capture_pass = all(capture[m]["pass"] for m in ("H", "A"))
    final_pass = bool(
        source_domain["pass"]
        and provenance_pass
        and continuous_pass
        and hybrid["pass"]
        and stochastic.get("pass")
        and capture_pass
    )
    return {
        "schema": 3,
        "qualification": "INDEPENDENT_DEPLOYMENT_THEOREM_COMPOSITION_GATE",
        "source_domain": source_domain,
        "source_domain_pass": source_domain["pass"],
        "validated_provenance_pass": provenance_pass,
        "continuous_linear_and_nonlinear_pass": continuous_pass,
        "hybrid": hybrid,
        "hybrid_pass": hybrid["pass"],
        "hybrid_required": sorted(REQUIRED_HYBRID),
        "hybrid_seen": hybrid["satisfied"],
        "stochastic": stochastic,
        "capture": capture,
        "finite_capture_pass": capture_pass,
        "arithmetic_error": arithmetic_error,
        "deployment_theorem_certificate": "PASS" if final_pass else "FAIL",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--validated-check", type=Path, required=True)
    ap.add_argument("--source-domain", type=Path, required=True)
    ap.add_argument("--primitive-bounds", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    check = json.loads(args.validated_check.read_text())
    source_domain_payload = json.loads(args.source_domain.read_text())
    primitive = json.loads(args.primitive_bounds.read_text())
    out = compose(check, source_domain_payload, primitive)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True))
    print(json.dumps(out, indent=2, sort_keys=True))
    return 0 if out["deployment_theorem_certificate"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
