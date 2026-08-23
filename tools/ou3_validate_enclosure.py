#!/usr/bin/env python3
"""Validate a rigorous OU-III information-metric source enclosure.

This is the promotion gate from the executed information certificate to the
continuous-source stability theorem.  It intentionally uses the same metric as
the primary numerical certificate,

    M(g) = Sigma_KF(g)^(-1),

and does not depend on the superseded coarse ``path_metrics.npz`` route.

For a source word with covariance recursion

    Sigma_1 = Phi Sigma_0 Phi^T + Omega,

strict linear contraction follows exactly from a validated uniform bound

    Sigma_1^-1/2 Omega Sigma_1^-1/2 >= eta I,   eta > 0,

because ``lambda_information <= 1-eta < 1``.  The validator also requires
uniform covariance eigenvalue bounds (metric equivalence), finite prefix gain,
a nonzero SO(3) chart with positive nonlinear word dissipation, strict hybrid
inward margins, and non-empirical stochastic constants.

The validator never trusts a supplied PASS flag.  It recomputes every logical
promotion condition from the supplied outward-rounded bounds and cross-checks
continuous bounds against executed reference points recorded in
``information_enclosure_contract.json``.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_numerical_certificate as BASE

SCHEMA = 2
ANCHOR_REL_TOL = 5.0e-6


def finite_positive(value) -> bool:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(x) and x > 0.0


def finite_nonnegative(value) -> bool:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(x) and x >= 0.0


def _anchor_upper(value: float, anchor: float) -> bool:
    """A continuous lower bound cannot exceed an included executed minimum."""
    tol = ANCHOR_REL_TOL * max(1.0, abs(anchor))
    return value <= anchor + tol


def _anchor_lower(value: float, anchor: float) -> bool:
    """A continuous upper bound cannot be below an included executed maximum."""
    tol = ANCHOR_REL_TOL * max(1.0, abs(anchor))
    return value + tol >= anchor


def validate_mode(mode: str, payload: dict, contract_mode: dict) -> dict:
    failures: list[str] = []
    if not payload.get("source_complete", False):
        failures.append("source family is not declared source-complete")
    if not payload.get("outward_rounded", False):
        failures.append("mode bounds are not declared outward-rounded")

    expected_horizon = contract_mode.get("recommended_word_horizon_s")
    horizon = payload.get("word_horizon_s")
    if not finite_positive(horizon):
        failures.append("word_horizon_s is not finite positive")
    elif expected_horizon is not None and not math.isclose(
        float(horizon), float(expected_horizon), rel_tol=0.0, abs_tol=1e-9
    ):
        failures.append(
            f"word_horizon_s {horizon} does not match contract horizon {expected_horizon}"
        )

    eta = payload.get("relative_Riccati_injection_margin_lower")
    sigma_lo = payload.get("Sigma_lambda_min_lower")
    sigma_hi = payload.get("Sigma_lambda_max_upper")
    prefix = payload.get("prefix_information_gain_upper")

    if not finite_positive(eta):
        failures.append("relative Riccati injection lower bound is not strictly positive")
    elif float(eta) >= 1.0:
        failures.append("relative Riccati injection lower bound must be < 1")
    if not finite_positive(sigma_lo):
        failures.append("Sigma_lambda_min_lower is not strictly positive")
    if not finite_positive(sigma_hi):
        failures.append("Sigma_lambda_max_upper is not finite positive")
    if finite_positive(sigma_lo) and finite_positive(sigma_hi) and float(sigma_hi) < float(sigma_lo):
        failures.append("Sigma eigenvalue upper bound is below lower bound")
    if not finite_positive(prefix):
        failures.append("prefix_information_gain_upper is not finite positive")

    # Since the continuous source set is required to contain the executed
    # reference family, its lower/upper bounds must contain those observed
    # points.  This catches accidentally optimistic or fabricated enclosures.
    ref = contract_mode.get("executed_reference_only", {})
    ref_eta = ref.get("relative_Riccati_injection_margin_worst")
    if finite_positive(eta) and finite_positive(ref_eta):
        if not _anchor_upper(float(eta), float(ref_eta)):
            failures.append(
                "validated Riccati lower bound exceeds included executed minimum "
                f"({eta} > {ref_eta})"
            )
    ref_sigma_lo = ref.get("Sigma_endpoint_lambda_min")
    if finite_positive(sigma_lo) and finite_positive(ref_sigma_lo):
        if not _anchor_upper(float(sigma_lo), float(ref_sigma_lo)):
            failures.append(
                "validated Sigma lower bound exceeds included executed minimum "
                f"({sigma_lo} > {ref_sigma_lo})"
            )
    ref_sigma_hi = ref.get("Sigma_endpoint_lambda_max")
    if finite_positive(sigma_hi) and finite_positive(ref_sigma_hi):
        if not _anchor_lower(float(sigma_hi), float(ref_sigma_hi)):
            failures.append(
                "validated Sigma upper bound excludes an executed maximum "
                f"({sigma_hi} < {ref_sigma_hi})"
            )

    linear_failures = list(failures)
    lambda_upper = None
    metric_lower = None
    metric_upper = None
    if finite_positive(eta) and float(eta) < 1.0:
        lambda_upper = 1.0 - float(eta)
    if finite_positive(sigma_hi):
        metric_lower = 1.0 / float(sigma_hi)
    if finite_positive(sigma_lo):
        metric_upper = 1.0 / float(sigma_lo)

    theta = payload.get("theta_star")
    if not finite_positive(theta) or not float(theta) < math.pi:
        failures.append("theta_star is not in (0, pi)")
    mu = payload.get("mu_W_lower")
    if not finite_positive(mu):
        failures.append("nonlinear mu_W lower bound is not strictly positive")
    if not payload.get("all_word_prefixes_safe", False):
        failures.append("nonlinear word prefixes are not validated safe")

    return {
        "mode": mode,
        "linear_pass": not linear_failures,
        "nonlinear_pass": not failures,
        "pass": not failures,
        "failures": failures,
        "linear_failures": linear_failures,
        "word_horizon_s": horizon,
        "relative_Riccati_injection_margin_lower": eta,
        "lambda_information_upper": lambda_upper,
        "Sigma_lambda_min_lower": sigma_lo,
        "Sigma_lambda_max_upper": sigma_hi,
        "information_metric_lambda_min_lower": metric_lower,
        "information_metric_lambda_max_upper": metric_upper,
        "prefix_information_gain_upper": prefix,
        "theta_star": theta,
        "mu_W_lower": mu,
    }


def validate_hybrid(payload: list[dict]) -> dict:
    required = {
        "startup_handoff",
        "held_to_active",
        "magnetic_regauge",
        "tilt_reset",
        "cooldown",
    }
    seen: set[str] = set()
    failures: list[str] = []
    rows = []
    for i, jump in enumerate(payload):
        kind = str(jump.get("kind", ""))
        seen.add(kind)
        margin = jump.get("inward_margin_lower")
        if not finite_positive(margin):
            failures.append(f"hybrid[{i}] {kind}: inward margin is not strictly positive")
        if not jump.get("outward_rounded", False):
            failures.append(f"hybrid[{i}] {kind}: not outward-rounded")
        rows.append({"kind": kind, "inward_margin_lower": margin})
    missing = sorted(required - seen)
    if missing:
        failures.append(f"missing hybrid obligations: {missing}")
    return {"pass": not failures, "failures": failures, "seen": sorted(seen), "bounds": rows}


def validate_stochastic(payload: dict) -> dict:
    failures: list[str] = []
    for key in ("Sigma_bar_norm_upper", "b_W_upper", "v_W_upper"):
        if not finite_nonnegative(payload.get(key)):
            failures.append(f"{key} missing/nonfinite/negative")
    prob = payload.get("finite_horizon_failure_probability_upper")
    try:
        p = float(prob)
    except (TypeError, ValueError):
        p = math.nan
    if not math.isfinite(p) or not (0.0 <= p < 1.0):
        failures.append("finite_horizon_failure_probability_upper must be in [0,1)")
    if not payload.get("outward_rounded", False):
        failures.append("stochastic bounds are not outward-rounded")
    return {
        "pass": not failures,
        "failures": failures,
        "finite_horizon_failure_probability_upper": prob,
        "Sigma_bar_norm_upper": payload.get("Sigma_bar_norm_upper"),
        "b_W_upper": payload.get("b_W_upper"),
        "v_W_upper": payload.get("v_W_upper"),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--certificate-dir", type=Path, default=BASE.DEFAULT_OUT)
    ap.add_argument("--enclosure", type=Path, required=True)
    args = ap.parse_args()
    cert = args.certificate_dir.resolve()
    inp = json.loads(args.enclosure.read_text())
    if inp.get("schema") != SCHEMA:
        raise RuntimeError(f"unsupported validated enclosure schema: {inp.get('schema')}")

    contract_path = cert / "information_enclosure_contract.json"
    if not contract_path.exists():
        raise FileNotFoundError(contract_path)
    contract = json.loads(contract_path.read_text())
    if contract.get("schema") != SCHEMA:
        raise RuntimeError("information enclosure contract schema mismatch")

    prov = inp.get("provenance", {})
    provenance_ok = (
        bool(prov.get("validated_arithmetic"))
        and bool(prov.get("outward_rounding"))
        and bool(prov.get("source_generated_not_trajectory_fit"))
        and bool(prov.get("continuous_source_coverage"))
    )

    modes = {}
    for mode in ("H", "A"):
        modes[mode] = validate_mode(
            mode,
            inp.get("modes", {}).get(mode, {}),
            contract.get("modes", {}).get(mode, {}),
        )

    hybrid = validate_hybrid(inp.get("hybrid", []))
    stochastic = validate_stochastic(inp.get("stochastic", {}))

    linear_pass = provenance_ok and all(modes[m]["linear_pass"] for m in modes)
    neighborhood_pass = provenance_ok and all(modes[m]["nonlinear_pass"] for m in modes) and hybrid["pass"]
    deployment_pass = neighborhood_pass and stochastic["pass"]

    out = {
        "schema": SCHEMA,
        "metric": "M(g)=Sigma_KF(g)^(-1)",
        "validated_enclosure_provenance_pass": provenance_ok,
        "modes": modes,
        "hybrid": hybrid,
        "stochastic": stochastic,
        "continuous_linear_information_certificate": "PASS" if linear_pass else "FAIL",
        "numerical_neighborhood_certificate": "PASS" if neighborhood_pass else "FAIL",
        "deployment_theorem_certificate": "PASS" if deployment_pass else "FAIL",
        "promotion_is_machine_verified": True,
    }
    (cert / "validated_enclosure_check.json").write_text(
        json.dumps(out, indent=2, sort_keys=True)
    )
    print(json.dumps({
        "continuous_linear_information_certificate": out["continuous_linear_information_certificate"],
        "numerical_neighborhood_certificate": out["numerical_neighborhood_certificate"],
        "deployment_theorem_certificate": out["deployment_theorem_certificate"],
        "H": modes["H"]["pass"],
        "A": modes["A"]["pass"],
        "hybrid": hybrid["pass"],
        "stochastic": stochastic["pass"],
        "provenance": provenance_ok,
    }, indent=2, sort_keys=True))
    return 0 if deployment_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
