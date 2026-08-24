#!/usr/bin/env python3
"""Validate the rigorous OU-III source-word enclosure and promotion margins.

Schema 3 separates the coordinate-invariant P3 endpoint information inequality
from the P4/P5 Lyapunov geometry.  P3 must provide a strictly positive complete-
word endpoint Riccati injection margin.  Nonlinear promotion must additionally
provide a group-compatible node metric of the paper's form

    Pbar_i = blkdiag((a_R_i/2) I3, P_xi_i),
    W_i = a_R_i (1-cos(theta)) + xi' P_xi_i xi.

No inverse-covariance-as-W alias, repeated one-step contraction, common-
Euclidean fallback, or replay-fitted promotion path is accepted.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_numerical_certificate as BASE

SCHEMA = 3
SOURCE_NOISE_SCHEMA = 1
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
    tol = ANCHOR_REL_TOL * max(1.0, abs(anchor))
    return value <= anchor + tol


def _anchor_lower(value: float, anchor: float) -> bool:
    tol = ANCHOR_REL_TOL * max(1.0, abs(anchor))
    return value + tol >= anchor


def _validate_group_metric(payload: dict) -> tuple[dict, list[str]]:
    failures: list[str] = []
    if payload.get("kind") != "GROUP_COMPATIBLE_NODE_METRIC":
        failures.append("path metric kind is not GROUP_COMPATIBLE_NODE_METRIC")
    if payload.get("node_dependent") is not True:
        failures.append("path metric does not permit/retain node dependence")
    if payload.get("attitude_block_isotropic") is not True:
        failures.append("attitude block is not isotropic (a_R/2) I3")
    if payload.get("attitude_linear_cross_terms") is not False:
        failures.append("attitude-linear cross terms are not excluded")
    if payload.get("common_Euclidean_metric") is not False:
        failures.append("common Euclidean metric fallback is active")
    if payload.get("equals_Kalman_inverse_covariance") is not False:
        failures.append("nonlinear W is incorrectly identified with Kalman inverse covariance")

    aR = payload.get("a_R_lower")
    pmin = payload.get("P_xi_lambda_min_lower")
    pmax = payload.get("P_xi_lambda_max_upper")
    pbar_min = payload.get("Pbar_lambda_min_lower")
    pbar_max = payload.get("Pbar_lambda_max_upper")
    for label, value in (
        ("a_R_lower", aR),
        ("P_xi_lambda_min_lower", pmin),
        ("P_xi_lambda_max_upper", pmax),
        ("Pbar_lambda_min_lower", pbar_min),
        ("Pbar_lambda_max_upper", pbar_max),
    ):
        if not finite_positive(value):
            failures.append(f"{label} is not finite positive")
    if finite_positive(pmin) and finite_positive(pmax) and float(pmax) < float(pmin):
        failures.append("P_xi upper metric bound is below lower bound")
    if finite_positive(pbar_min) and finite_positive(pbar_max) and float(pbar_max) < float(pbar_min):
        failures.append("Pbar upper metric bound is below lower bound")
    return {
        "kind": payload.get("kind"),
        "node_dependent": payload.get("node_dependent"),
        "a_R_lower": aR,
        "P_xi_lambda_min_lower": pmin,
        "P_xi_lambda_max_upper": pmax,
        "Pbar_lambda_min_lower": pbar_min,
        "Pbar_lambda_max_upper": pbar_max,
        "attitude_block_isotropic": payload.get("attitude_block_isotropic"),
        "attitude_linear_cross_terms": payload.get("attitude_linear_cross_terms"),
        "common_Euclidean_metric": payload.get("common_Euclidean_metric"),
        "equals_Kalman_inverse_covariance": payload.get("equals_Kalman_inverse_covariance"),
    }, failures


def validate_mode(mode: str, payload: dict, contract_mode: dict,
                  sampled_mode: dict | None = None) -> dict:
    failures: list[str] = []
    if not payload.get("source_complete", False):
        failures.append("source family is not declared source-complete")
    if not payload.get("outward_rounded", False):
        failures.append("mode bounds are not declared outward-rounded")
    if payload.get("joint_source_reachability") is not True:
        failures.append("mode enclosure is not jointly source-reachable")
    if payload.get("one_sample_decrease_used") is not False:
        failures.append("one-sample decrease/repeated-step shortcut is active")

    expected_horizon = contract_mode.get("recommended_word_horizon_s")
    horizon = payload.get("word_horizon_s")
    if not finite_positive(horizon):
        failures.append("word_horizon_s is not finite positive")
    elif expected_horizon is not None and not math.isclose(
        float(horizon), float(expected_horizon), rel_tol=0.0, abs_tol=1e-9
    ):
        failures.append(f"word_horizon_s {horizon} does not match contract horizon {expected_horizon}")

    # No compatibility alias: schema 3 requires the endpoint field explicitly.
    eta = payload.get("word_endpoint_relative_Riccati_injection_margin_lower")
    sigma_lo = payload.get("Sigma_lambda_min_lower")
    sigma_hi = payload.get("Sigma_lambda_max_upper")
    prefix = payload.get("prefix_information_gain_upper")

    if not finite_positive(eta):
        failures.append("word-endpoint Riccati injection lower bound is not strictly positive")
    elif float(eta) >= 1.0:
        failures.append("word-endpoint Riccati injection lower bound must be < 1")
    if not finite_positive(sigma_lo):
        failures.append("Sigma_lambda_min_lower is not strictly positive")
    if not finite_positive(sigma_hi):
        failures.append("Sigma_lambda_max_upper is not finite positive")
    if finite_positive(sigma_lo) and finite_positive(sigma_hi) and float(sigma_hi) < float(sigma_lo):
        failures.append("Sigma eigenvalue upper bound is below lower bound")
    if not finite_positive(prefix):
        failures.append("prefix_information_gain_upper is not finite positive")

    ref = contract_mode.get("executed_reference_only", {})
    ref_eta = ref.get("relative_Riccati_injection_margin_worst")
    if finite_positive(eta) and finite_positive(ref_eta) and not _anchor_upper(float(eta), float(ref_eta)):
        failures.append(f"validated endpoint Riccati lower bound exceeds included executed minimum ({eta} > {ref_eta})")
    ref_sigma_lo = ref.get("Sigma_endpoint_lambda_min")
    if finite_positive(sigma_lo) and finite_positive(ref_sigma_lo) and not _anchor_upper(float(sigma_lo), float(ref_sigma_lo)):
        failures.append(f"validated Sigma lower bound exceeds included executed minimum ({sigma_lo} > {ref_sigma_lo})")
    ref_sigma_hi = ref.get("Sigma_endpoint_lambda_max")
    if finite_positive(sigma_hi) and finite_positive(ref_sigma_hi) and not _anchor_lower(float(sigma_hi), float(ref_sigma_hi)):
        failures.append(f"validated Sigma upper bound excludes an executed maximum ({sigma_hi} < {ref_sigma_hi})")

    linear_failures = list(failures)
    lambda_upper = 1.0 - float(eta) if finite_positive(eta) and float(eta) < 1.0 else None

    # P3 inverse-covariance conditioning diagnostics remain diagnostics only.
    conditioning_lower = 1.0 / float(sigma_hi) if finite_positive(sigma_hi) else None
    conditioning_upper = 1.0 / float(sigma_lo) if finite_positive(sigma_lo) else None

    metric, metric_failures = _validate_group_metric(payload.get("path_metric", {}))
    failures.extend(metric_failures)

    theta = payload.get("theta_star")
    if not finite_positive(theta) or not float(theta) < math.pi:
        failures.append("theta_star is not in (0, pi)")

    endpoint_ratio = payload.get("endpoint_W_ratio_upper")
    mu = None
    if not finite_nonnegative(endpoint_ratio) or not float(endpoint_ratio) < 1.0:
        failures.append("endpoint_W_ratio_upper is not in [0,1)")
    else:
        mu = 1.0 - float(endpoint_ratio)

    certified_level = payload.get("certified_level_W")
    if not finite_positive(certified_level):
        failures.append("certified_level_W is not finite positive")
    if not payload.get("all_word_prefixes_safe", False):
        failures.append("nonlinear word prefixes are not validated safe")

    if sampled_mode and finite_positive(certified_level):
        first_fail = sampled_mode.get("first_fail_W")
        if finite_positive(first_fail) and float(certified_level) >= float(first_fail):
            failures.append(f"certified nonlinear level reaches/exceeds a sampled failure ({certified_level} >= {first_fail})")

    return {
        "mode": mode,
        "linear_pass": not linear_failures,
        "nonlinear_pass": not failures,
        "pass": not failures,
        "failures": failures,
        "linear_failures": linear_failures,
        "word_horizon_s": horizon,
        "word_endpoint_relative_Riccati_injection_margin_lower": eta,
        "lambda_information_upper": lambda_upper,
        "Sigma_lambda_min_lower": sigma_lo,
        "Sigma_lambda_max_upper": sigma_hi,
        "P3_inverse_covariance_conditioning_lambda_min_lower": conditioning_lower,
        "P3_inverse_covariance_conditioning_lambda_max_upper": conditioning_upper,
        "prefix_information_gain_upper": prefix,
        "path_metric": metric,
        "theta_star": theta,
        "endpoint_W_ratio_upper": endpoint_ratio,
        "mu_W_lower": mu,
        "certified_level_W": certified_level,
    }


def validate_hybrid(payload: list[dict], modes: dict) -> dict:
    required = {"startup_handoff", "held_to_active", "magnetic_regauge", "tilt_reset", "cooldown"}
    seen: set[str] = set()
    failures: list[str] = []
    rows = []

    for i, jump in enumerate(payload):
        kind = str(jump.get("kind", ""))
        seen.add(kind)
        row_failures: list[str] = []
        if not jump.get("source_complete", False):
            row_failures.append("source family is not source-complete")
        if not jump.get("outward_rounded", False):
            row_failures.append("bounds are not outward-rounded")

        source = jump.get("source_level_W_upper")
        gain = jump.get("jump_gain_upper")
        additive = jump.get("additive_W_upper")
        dest = jump.get("destination_level_W")
        new_coord = jump.get("new_coordinate_W_upper", 0.0)
        for label, value in (
            ("source_level_W_upper", source), ("jump_gain_upper", gain),
            ("additive_W_upper", additive), ("new_coordinate_W_upper", new_coord),
        ):
            if not finite_nonnegative(value):
                row_failures.append(f"{label} missing/nonfinite/negative")
        if not finite_positive(dest):
            row_failures.append("destination_level_W missing/nonfinite/nonpositive")

        dest_mode = str(jump.get("destination_mode", ""))
        if dest_mode not in modes:
            row_failures.append("destination_mode must be H or A")
        elif finite_positive(dest):
            mode_level = modes[dest_mode].get("certified_level_W")
            if not finite_positive(mode_level) or float(dest) > float(mode_level):
                row_failures.append("destination level is outside the certified same-mode level")

        if kind == "held_to_active":
            if jump.get("dimension_change_handled_by_embedding") is not True:
                row_failures.append("held_to_active must declare explicit dimension embedding")
            if jump.get("source_dimension") != 18 or jump.get("destination_dimension") != 21:
                row_failures.append("held_to_active dimensions must be 18 -> 21")
            if "new_coordinate_W_upper" not in jump:
                row_failures.append("held_to_active must bound the new bias-coordinate energy")

        if kind == "tilt_reset":
            if jump.get("discarded_pre_reset_tilt_excluded_from_multiplicative_gain") is not True:
                row_failures.append("tilt reset still charges discarded pre-reset tilt in multiplicative gain")
            if jump.get("reset_to_funnel_exact_map") is not True:
                row_failures.append("tilt reset is not certified as exact reset-to-funnel map")

        if kind == "cooldown":
            if jump.get("reachable_word_product_used") is not True:
                row_failures.append("cooldown does not use reachable word products")
            if jump.get("global_worst_word_power_used") is not False:
                row_failures.append("cooldown uses a power of a global worst-word factor")

        post = margin = None
        if (finite_nonnegative(source) and finite_nonnegative(gain)
                and finite_nonnegative(additive) and finite_nonnegative(new_coord)
                and finite_positive(dest)):
            post = float(gain) * float(source) + float(additive) + float(new_coord)
            margin = float(dest) - post
            if not margin > 0.0:
                row_failures.append("recomputed inward margin is not strictly positive")

        if row_failures:
            failures.extend(f"hybrid[{i}] {kind}: {f}" for f in row_failures)
        rows.append({
            "kind": kind, "destination_mode": dest_mode,
            "source_level_W_upper": source, "jump_gain_upper": gain,
            "additive_W_upper": additive, "new_coordinate_W_upper": new_coord,
            "post_jump_W_upper": post, "destination_level_W": dest,
            "inward_margin_lower": margin, "pass": not row_failures,
        })

    missing = sorted(required - seen)
    if missing:
        failures.append(f"missing hybrid obligations: {missing}")
    return {"pass": not failures, "failures": failures, "seen": sorted(seen), "bounds": rows}


def source_noise_moments(payload: dict) -> tuple[dict, list[str]]:
    failures: list[str] = []
    if payload.get("schema") != SOURCE_NOISE_SCHEMA:
        failures.append("source noise certificate schema mismatch")
    if not payload.get("source_generated_not_trajectory_fit", False):
        failures.append("noise certificate is not source-generated")
    std = payload.get("standardized_increment", {})
    d = std.get("dimension")
    if not isinstance(d, int) or d <= 0:
        failures.append("standardized primitive dimension is not positive integer")
        d = 0
    if std.get("covariance_upper_identity") is not True:
        failures.append("standardized covariance is not certified <= identity")
    return {
        "dimension": d,
        "Sigma_bar_norm_upper": 1.0 if d else None,
        "trace_Sigma_bar_upper": float(d) if d else None,
        "trace_Sigma_bar_squared_upper": float(d) if d else None,
        "s2_upper": float(d) if d else None,
        "s4_upper": float(d * d + 2 * d) if d else None,
    }, failures


def _geometric_square_sum(L: float, n: int) -> float:
    if n <= 0:
        return 0.0
    q = L * L
    if abs(q - 1.0) <= 1e-12:
        return float(n)
    try:
        ans = (q ** n - 1.0) / (q - 1.0)
    except OverflowError:
        return math.inf
    return ans if math.isfinite(ans) and ans >= 0.0 else math.inf


def _gaussian_t_star(w_star: float, m: float, v: float, b: float) -> float | None:
    w2 = w_star * w_star
    if not (w2 > m and b > 0.0 and v >= 0.0):
        return None
    root = math.sqrt(v + 2.0 * b * (w2 - m))
    y = (root - math.sqrt(v)) / (2.0 * b)
    return y * y if y > 0.0 and math.isfinite(y) else None


def validate_stochastic(payload: dict, modes: dict, source_noise: dict) -> dict:
    failures: list[str] = []
    if not payload.get("source_complete", False):
        failures.append("stochastic source family is not source-complete")
    if not payload.get("outward_rounded", False):
        failures.append("stochastic sensitivity bounds are not outward-rounded")
    if not payload.get("localization_prefix_safe", False):
        failures.append("localized word prefixes are not validated safe")
    if payload.get("gaussian_localization_used") is not True:
        failures.append("Gaussian localization is not explicitly certified")
    if payload.get("freedman_excursion_used") is not True:
        failures.append("Freedman excursion concentration is not explicitly certified")
    if payload.get("markov_union_fallback_used") is not False:
        failures.append("retired Markov/union stochastic fallback is active")

    moments, noise_failures = source_noise_moments(source_noise)
    failures.extend(noise_failures)
    vals = {}
    for name in ("L_X_upper", "G_bar_upper", "c_zw_upper", "r_star_upper",
                 "c_ww_upper", "g_W_upper", "h_W_upper"):
        value = payload.get(name)
        if not finite_nonnegative(value):
            failures.append(f"{name} missing/nonfinite/negative")
        else:
            vals[name] = float(value)

    w_star = payload.get("localization_radius_standardized")
    if not finite_positive(w_star):
        failures.append("localization_radius_standardized is not finite positive")
    else:
        w_star = float(w_star)
    ell = payload.get("word_samples_upper")
    if not isinstance(ell, int) or ell <= 0:
        failures.append("word_samples_upper must be a positive integer")
        ell = 0
    N = payload.get("finite_horizon_words")
    if not isinstance(N, int) or N <= 0:
        failures.append("finite_horizon_words must be a positive integer")
        N = 0

    a = payload.get("funnel_level_a")
    W0 = payload.get("W0_upper")
    if not finite_positive(a):
        failures.append("funnel_level_a is not finite positive")
    else:
        a = float(a)
        mode_levels = [m.get("certified_level_W") for m in modes.values()]
        if any(not finite_positive(x) for x in mode_levels) or a > min(float(x) for x in mode_levels):
            failures.append("stochastic funnel level exceeds a certified deterministic mode level")
    if not finite_nonnegative(W0):
        failures.append("W0_upper missing/nonfinite/negative")
    elif finite_positive(a) and not float(W0) < float(a):
        failures.append("W0_upper must be strictly inside funnel_level_a")
    else:
        W0 = float(W0)

    ratios = [m.get("endpoint_W_ratio_upper") for m in modes.values()]
    lambda_W = None
    if any(not finite_nonnegative(r) or not float(r) < 1.0 for r in ratios):
        failures.append("deterministic mode contraction factors are unavailable")
    else:
        lambda_W = max(float(r) for r in ratios)

    dt = source_noise.get("physical_scales", {}).get("imu_dt_s")
    horizons = [m.get("word_horizon_s") for m in modes.values()]
    if finite_positive(dt) and ell > 0 and all(finite_positive(h) for h in horizons):
        minimum_samples = math.ceil(max(float(h) for h in horizons) / float(dt) - 1e-12)
        if ell < minimum_samples:
            failures.append(f"word_samples_upper {ell} is below source horizon minimum {minimum_samples}")
    else:
        failures.append("cannot connect word_samples_upper to source sample period")

    nu1 = nu_W = lambda_s = sigma_s2 = None
    b_W = v_W = V_W = drift_max = x_star = None
    t_star = p_localization = p_freedman = p_total = None
    if not failures:
        s2, s4 = moments["s2_upper"], moments["s4_upper"]
        g_star = vals["G_bar_upper"] + vals["c_zw_upper"] * vals["r_star_upper"]
        nu1 = 2.0 * g_star * g_star * s2 + 2.0 * vals["c_ww_upper"] ** 2 * s4
        nu_W = float(ell) * nu1 * _geometric_square_sum(vals["L_X_upper"], ell)
        if not math.isfinite(nu_W):
            failures.append("localized source-word noise moment bound overflowed")
        else:
            lambda_s = 0.5 * (1.0 + lambda_W)
            sigma_s2 = (vals["h_W_upper"] + vals["g_W_upper"] ** 2 * lambda_W /
                         (2.0 * (1.0 - lambda_W))) * nu_W
            fixed = sigma_s2 / (1.0 - lambda_s)
            lamN = lambda_s ** N
            WN = lamN * W0 + (1.0 - lamN) * fixed
            drift_max = max(W0, WN)
            x_star = a - drift_max
            if not x_star > 0.0:
                failures.append("stochastic drift envelope reaches/exceeds the funnel boundary")
            b_W = a
            v_W = a * a / 4.0
            V_W = v_W / (1.0 - lambda_s * lambda_s)
            t_star = _gaussian_t_star(w_star, moments["trace_Sigma_bar_upper"],
                                     moments["trace_Sigma_bar_squared_upper"],
                                     moments["Sigma_bar_norm_upper"])
            if t_star is None:
                failures.append("localization radius does not exceed the source Gaussian RMS threshold")
            if x_star > 0.0 and t_star is not None:
                exponent = -x_star * x_star / (2.0 * (V_W + b_W * x_star / 3.0))
                p_freedman = float(N) * math.exp(exponent)
                p_localization = float(N * ell) * math.exp(-t_star)
                p_total = p_freedman + p_localization
                if not (math.isfinite(p_total) and 0.0 <= p_total < 1.0):
                    failures.append("recomputed finite-horizon stochastic failure probability is not < 1")

    return {
        "pass": not failures, "failures": failures, "source_noise_moments": moments,
        "lambda_W_upper": lambda_W, "lambda_s_upper": lambda_s,
        "nu1_upper": nu1, "nu_W_upper": nu_W, "sigma_s2_upper": sigma_s2,
        "funnel_level_a": a, "W0_upper": W0, "b_W_upper": b_W,
        "v_W_upper": v_W, "V_W_upper": V_W,
        "drift_envelope_max_upper": drift_max,
        "excursion_margin_x_star_lower": x_star, "gaussian_t_star_lower": t_star,
        "localization_failure_probability_upper": p_localization,
        "freedman_failure_probability_upper": p_freedman,
        "finite_horizon_failure_probability_upper": p_total,
        "word_samples_upper": ell, "finite_horizon_words": N,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--certificate-dir", type=Path, default=BASE.DEFAULT_OUT)
    ap.add_argument("--enclosure", type=Path, required=True)
    ap.add_argument("--source-noise", type=Path, default=None)
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
    source_noise_path = args.source_noise or (cert / "source_noise_certificate.json")
    if not source_noise_path.exists():
        raise FileNotFoundError(source_noise_path)
    source_noise = json.loads(source_noise_path.read_text())

    radius_path = cert / "neighborhood_radius_search" / "neighborhood_radius_search.json"
    sampled_modes = {}
    if radius_path.exists():
        sampled_modes = json.loads(radius_path.read_text()).get("modes", {})

    prov = inp.get("provenance", {})
    provenance_ok = (
        bool(prov.get("validated_arithmetic")) and bool(prov.get("outward_rounding"))
        and bool(prov.get("source_generated_not_trajectory_fit"))
        and bool(prov.get("continuous_source_coverage"))
        and bool(prov.get("joint_source_reachability"))
    )
    modes = {
        mode: validate_mode(mode, inp.get("modes", {}).get(mode, {}),
                            contract.get("modes", {}).get(mode, {}), sampled_modes.get(mode))
        for mode in ("H", "A")
    }
    hybrid = validate_hybrid(inp.get("hybrid", []), modes)
    stochastic = validate_stochastic(inp.get("stochastic", {}), modes, source_noise)
    linear_pass = provenance_ok and all(modes[m]["linear_pass"] for m in modes)
    nonlinear_pass = provenance_ok and all(modes[m]["nonlinear_pass"] for m in modes)
    hybrid_pass = provenance_ok and nonlinear_pass and hybrid["pass"]
    stochastic_pass = provenance_ok and nonlinear_pass and hybrid["pass"] and stochastic["pass"]
    deployment_pass = linear_pass and nonlinear_pass and hybrid_pass and stochastic_pass
    out = {
        "schema": SCHEMA,
        "metric_policy": {
            "P3": "word-endpoint generalized information inequality; conditioning congruence only",
            "P4_P5": "node-wise group-compatible W_i=a_R_i(1-cos(theta))+xi^T P_xi_i xi",
            "fallbacks": "NONE",
        },
        "validated_enclosure_provenance_pass": provenance_ok,
        "modes": modes, "hybrid": hybrid, "stochastic": stochastic,
        "continuous_linear_information_certificate": "PASS" if linear_pass else "FAIL",
        "numerical_neighborhood_certificate": "PASS" if nonlinear_pass else "FAIL",
        "hybrid_funnel_certificate": "PASS" if hybrid_pass else "FAIL",
        "stochastic_certificate": "PASS" if stochastic_pass else "FAIL",
        "deployment_theorem_certificate": "PASS" if deployment_pass else "FAIL",
        "promotion_is_machine_verified": True,
    }
    (cert / "validated_enclosure_check.json").write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "continuous_linear_information_certificate": out["continuous_linear_information_certificate"],
        "numerical_neighborhood_certificate": out["numerical_neighborhood_certificate"],
        "hybrid_funnel_certificate": out["hybrid_funnel_certificate"],
        "stochastic_certificate": out["stochastic_certificate"],
        "deployment_theorem_certificate": out["deployment_theorem_certificate"],
        "H": modes["H"]["pass"], "A": modes["A"]["pass"],
        "hybrid": hybrid["pass"], "stochastic": stochastic["pass"],
        "provenance": provenance_ok,
    }, indent=2, sort_keys=True))
    return 0 if deployment_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
