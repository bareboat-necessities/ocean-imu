#!/usr/bin/env python3
"""Validated first outer-H source-word/decrease certificate for OU-III P5.

This is deliberately *not* the local P4 ``B W`` recurrence with a larger
radius. It constructs a finite-angle, anisotropic sufficient inequality on the
actual P1 H-mode handoff nodes.

The error coordinates are

    z = [c(R), b_g, v, p, S, a_w],
    c(R) = 2 tan(theta/2) u,

and the metric is the same normalized Cayley-information geometry used by P4.
The source map is split according to its real nonlinear structure:

* prediction of [v,p,S,a_w] is exactly linear and is charged zero nonlinear
  remainder;
* the accelerometer nonlinear innovation depends on c and a_w, not on v or p;
* the magnetometer nonlinear innovation depends on c only;
* the S=0 innovation r_S=-S is exactly linear. Its complete Kalman correction,
  including S->attitude, remains in the homogeneous map. Only the finite
  quaternion/Cayley injection is nonlinear, and the certificate exposes the
  source-uniform S->attitude prefix bound needed to control it.

For the exact Cayley rotation

    R(c) = I + 4/(4+|c|^2)[c]x + 2/(4+|c|^2)[c]x^2,

one has, for |c|<=q,

    ||R(c)-I-[c]x|| <= q^2 (q+2)/(4+q^2),
    ||R(c)-I||       <= 2 q/sqrt(4+q^2).

These identities give state-proportional (rather than global ``C |z|^2``)
innovation-defect ratios on each P1 node. For a Kalman correction
K=P H'(H P H'+R)^-1=P+ H' R^-1, a defect eta is inserted at the posterior node.
The exact information identity therefore gives

    K' (P+)^-1 K = R^-1 H P+ H' R^-1 <= R^-1,

hence ||K eta||_(P+)^-1 <= ||eta||_R^-1. The left-error covariance reset is an
exact congruence, so this endpoint defect norm is transported through that reset
without an additional condition-number factor.

The resulting outer word test is a valid sufficient test conditional on a
source-safe candidate outer prefix domain. The producer does not silently assume
that bootstrap: full P5 promotion additionally requires a prefix-domain proof.
Failure is reported as a first failing inequality; it never promotes P5 by
replay or by extrapolating the P4 local recurrence.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_explicit_information_word_certificate as P3
import ou3_implementation_proof_manifest as MANIFEST
import ou3_p4_exact_word_map as WORDMAP
import ou3_p4_nonlinear_word_certificate as P4
import ou3_p5_startup_capture_certificate as P5ID
import ou3_startup_stability_certificate as P1

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def mul_up(a: float, b: float) -> float:
    return up(float(a) * float(b))


def add_up(a: float, b: float) -> float:
    return up(float(a) + float(b))


def div_up(a: float, b: float) -> float:
    if not b > 0.0:
        raise ValueError("positive denominator required")
    return up(float(a) / float(b))


def sqrt_up(x: float) -> float:
    if not (math.isfinite(x) and x >= 0.0):
        raise ValueError("sqrt requires finite nonnegative input")
    return up(math.sqrt(x))


def _rotation_bounds(q: float) -> dict:
    if not (math.isfinite(q) and 0.0 <= q < 2.0):
        raise RuntimeError("outer Cayley radius must lie in [0,2)")
    q2 = mul_up(q, q)
    den = down(4.0 + q2)
    rem = div_up(mul_up(q2, add_up(q, 2.0)), den)
    diff = div_up(mul_up(2.0, q), down(math.sqrt(4.0 + q2)))
    rem_over_q = 0.0 if q == 0.0 else div_up(rem, q)
    return {
        "cayley_norm_upper": q,
        "R_minus_I_minus_skew_norm_upper": rem,
        "R_minus_I_norm_upper": diff,
        "R_minus_I_minus_skew_linearized_ratio_upper": rem_over_q,
    }


def _node(name: str, q: float, bounds: dict) -> dict:
    rot = _rotation_bounds(q)
    bg = float(bounds["gyro_bias_error_norm_upper_rad_s"])
    v = float(bounds["velocity_error_norm_upper_mps"])
    p = float(bounds["position_error_norm_upper_m"])
    S = float(bounds["integral_displacement_error_norm_upper_m_s"])
    aw = float(bounds["latent_acceleration_error_norm_upper_mps2"])
    for label, x in (("b_g", bg), ("v", v), ("p", p), ("S", S), ("a_w", aw)):
        if not (math.isfinite(x) and x >= 0.0):
            raise RuntimeError(f"invalid P1 {label} outer radius")
    driver = sqrt_up(q*q + bg*bg + S*S + aw*aw)
    full = sqrt_up(q*q + bg*bg + v*v + p*p + S*S + aw*aw)
    return {
        "name": name,
        "rotation": rot,
        "coordinate_norm_radii": {
            "c": q,
            "b_g": bg,
            "v": v,
            "p": p,
            "S": S,
            "a_w": aw,
        },
        "nonlinear_driver_coordinates": ["c", "b_g", "S", "a_w"],
        "exact_linear_prediction_coordinates": ["v", "p", "S", "a_w"],
        "driver_canonical_norm_upper": driver,
        "full_canonical_norm_upper": full,
    }


def _kalman_posterior_metric_gain(scale: float, measurement_std: float) -> float:
    """sqrt(scale)/sigma from K'(P+)^-1 K <= R^-1."""
    if not (scale > 0.0 and measurement_std > 0.0):
        raise RuntimeError("positive metric scale and measurement std required")
    return div_up(sqrt_up(scale), measurement_std)


def _vector_defect(node: dict, H: dict, domain: dict) -> dict:
    q = float(node["rotation"]["cayley_norm_upper"])
    rot_ratio = float(node["rotation"]["R_minus_I_minus_skew_linearized_ratio_upper"])
    a_ratio = div_up(2.0, down(math.sqrt(4.0 + q*q)))
    live = domain["normal_live"]
    fmax = float(live["specific_force_norm_upper_mps2"])
    mmax = float(live["magnetic_vector_norm_upper_uT"])
    awmax = float(node["coordinate_norm_radii"]["a_w"])

    # Exact finite-angle source structure:
    # eta_acc <= fmax*rem_ratio*|c| + a_ratio*|c|*|a_w|.
    # On the product driver ball |c||a_w| <= min(q,awmax)||z_driver||.
    acc_residual_ratio = add_up(
        mul_up(fmax, rot_ratio),
        mul_up(a_ratio, min(q, awmax)),
    )
    mag_residual_ratio = mul_up(mmax, rot_ratio)

    acc_metric_gain = _kalman_posterior_metric_gain(
        float(H["metric_mode_global_positive_scale"]),
        float(H["measurement_bounds"]["acc_measurement_std_mps2"]),
    )
    mag_metric_gain = _kalman_posterior_metric_gain(
        float(H["metric_mode_global_positive_scale"]),
        float(H["measurement_bounds"]["mag_measurement_std_uT"]),
    )
    acc_info_ratio = mul_up(acc_metric_gain, acc_residual_ratio)
    mag_info_ratio = mul_up(mag_metric_gain, mag_residual_ratio)
    return {
        "finite_angle_accel_residual_defect_per_driver_norm_upper": acc_residual_ratio,
        "finite_angle_mag_residual_defect_per_driver_norm_upper": mag_residual_ratio,
        "kalman_posterior_information_gain_for_accel_defect_upper": acc_metric_gain,
        "kalman_posterior_information_gain_for_mag_defect_upper": mag_metric_gain,
        "accel_nonlinear_information_norm_per_driver_norm_upper": acc_info_ratio,
        "mag_nonlinear_information_norm_per_driver_norm_upper": mag_info_ratio,
        "posterior_metric_lemma": "K^T(P_plus)^-1K <= R^-1",
        "v_or_p_charged_as_vector_nonlinearity": False,
        "S_charged_as_vector_measurement_nonlinearity": False,
    }


def _s_to_attitude_prefix_bound(H: dict, node: dict) -> dict:
    """Full S->attitude gain bound from the source-uniform covariance envelope.

    R_S is isotropic in the proved deployment. For a PSD covariance block,
    P_th,S=P_th,th^(1/2) C P_SS^(1/2), ||C||<=1. With
    R_S=sigma_S^2 I, spectral calculus gives
    ||P_SS^(1/2)(P_SS+sigma_S^2 I)^-1||<=1/(2 sigma_S), hence
    ||K_th,S||<=sqrt(lambda_max(P_th,th))/(2 sigma_S).

    The current calculation deliberately starts with the source-uniform P3
    covariance eigenvalue upper. If it fails, P5 must replace that global box by
    a source-staged post-goLive theta/S cross-covariance enclosure; replay
    covariance is not admissible.
    """
    ptheta_upper = float(H["Sigma_lambda_max_upper"])
    import ou3_source_domain_contract as SOURCE
    text = (REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h").read_text(encoding="utf-8")
    rs_std = float(SOURCE.parse_const(text, "MIN_R_S"))
    Smax = float(node["coordinate_norm_radii"]["S"])
    gain = div_up(sqrt_up(ptheta_upper), mul_up(2.0, rs_std))
    dtheta = mul_up(gain, Smax)
    return {
        "full_S_to_attitude_gain_retained": True,
        "S_innovation_is_exactly_linear": True,
        "source_uniform_covariance_lambda_max_upper": ptheta_upper,
        "R_S_filter_std_lower": rs_std,
        "K_thetaS_operator_norm_upper": gain,
        "S_radius": Smax,
        "S_induced_attitude_correction_norm_upper": dtheta,
        "deployed_axis_angle_helper_certified_correction_limit": 3.0,
        "prefix_correction_inside_current_validated_group_helper": dtheta < 3.0,
        "tightening_required_if_false": (
            "source-staged post-goLive theta/S covariance and cross-covariance enclosure; the global P3 covariance envelope is too broad for outer prefix geometry"
            if not dtheta < 3.0 else None
        ),
    }


def _word_counts(H: dict, domain: dict) -> dict:
    samples = int(H["word_samples_upper"])
    dt = float(domain["configured_runtime"]["imu_dt_s"])
    horizon = float(H["word_horizon_s"])
    # The configured deployment's magnetometer runs at 25 Hz. The source word
    # allows arbitrary accepted/rejected packets, so charge every due packet.
    mag = int(math.ceil(horizon * 25.0)) + 1
    import ou3_source_domain_contract as SOURCE
    text = (REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h").read_text(encoding="utf-8")
    pseudo_min = float(SOURCE.parse_const(text, "PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT"))
    pseudo = int(math.ceil(horizon / pseudo_min)) + 1
    return {
        "word_horizon_s": horizon,
        "imu_dt_s": dt,
        "samples_upper": samples,
        "accepted_accel_corrections_upper": samples,
        "accepted_mag_corrections_upper": mag,
        "S_zero_corrections_upper": pseudo,
        "note": "source-safe maxima; no favorable rejection pattern is selected",
    }


def _node_word_test(node: dict, H: dict, domain: dict) -> dict:
    vec = _vector_defect(node, H, domain)
    sterm = _s_to_attitude_prefix_bound(H, node)
    counts = _word_counts(H, domain)
    delta = float(H["P3_word_endpoint_delta_lower"])

    rho_vec = add_up(
        mul_up(float(counts["accepted_accel_corrections_upper"]),
               float(vec["accel_nonlinear_information_norm_per_driver_norm_upper"])),
        mul_up(float(counts["accepted_mag_corrections_upper"]),
               float(vec["mag_nonlinear_information_norm_per_driver_norm_upper"])),
    )
    homogeneous_sqrt_decrease_lower = down(0.5 * delta)
    vector_margin = down(homogeneous_sqrt_decrease_lower - rho_vec)

    s_prefix_pass = bool(sterm["prefix_correction_inside_current_validated_group_helper"])
    vector_decrease_pass = vector_margin > 0.0
    # The node radii currently come from P1. They are the candidate outer proof
    # domain, but a complete word certificate must also prove every intermediate
    # source image remains in the domain on which these ratios were evaluated.
    # Do not infer that bootstrap from endpoint decrease.
    prefix_domain_bootstrap_pass = False
    pass_ = s_prefix_pass and vector_decrease_pass and prefix_domain_bootstrap_pass
    if not s_prefix_pass:
        first = "S_TO_ATTITUDE_OUTER_PREFIX_BOUND_NOT_CERTIFIED"
    elif not vector_decrease_pass:
        first = "FINITE_ANGLE_VECTOR_PERTURBATION_EXCEEDS_P3_OUTER_WORD_GAP"
    elif not prefix_domain_bootstrap_pass:
        first = "OUTER_SOURCE_PREFIX_DOMAIN_BOOTSTRAP_NOT_CERTIFIED"
    else:
        first = "NONE"
    return {
        "node": node["name"],
        "vector_defect": vec,
        "S_to_attitude_prefix": sterm,
        "word_counts": counts,
        "P3_homogeneous_sqrt_decrease_lower": homogeneous_sqrt_decrease_lower,
        "outer_vector_nonlinear_information_ratio_upper": rho_vec,
        "outer_vector_word_decrease_margin_lower": vector_margin,
        "vector_word_decrease_pass": vector_decrease_pass,
        "S_prefix_group_correction_pass": s_prefix_pass,
        "candidate_outer_prefix_domain_bootstrap_pass": prefix_domain_bootstrap_pass,
        "outer_word_decrease_pass": pass_,
        "first_failure": first,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("outer H domain must not be trajectory fitted")

    p1 = P1.build(domain_path)
    p3 = P3.build(domain_path)
    p4 = P4.build(domain_path)
    manifest = MANIFEST.build()
    wordmap = WORDMAP.build(domain_path)
    prereq = []
    prereq += [f"P1: {x}" for x in P1.validate(p1)]
    prereq += [f"P3: {x}" for x in P3.validate(p3)]
    prereq += [f"P4: {x}" for x in P4.validate(p4)]
    prereq += [f"manifest: {x}" for x in MANIFEST.validate(manifest)]
    prereq += [f"word-map: {x}" for x in WORDMAP.validate(wordmap)]
    if prereq:
        return {
            "schema": SCHEMA,
            "qualification": "OU3_P5_OUTER_H_SOURCE_WORD_DECREASE_CERTIFICATE",
            "source_generated_not_trajectory_fit": True,
            "source_replay_used": False,
            "prerequisite_failures": prereq,
            "P5_OUTER_H_WORD_CERTIFICATE": "NOT_ESTABLISHED",
            "first_failure": "UPSTREAM_PREREQUISITE_FAILURE",
            "failures": prereq,
        }

    H = p4["modes"]["H"]
    bounds = p1["go_live"]["physical_coordinate_bounds"]
    normal_q = P5ID._cayley_norm_upper_from_cos_lower(
        float(p1["normal_handoff"]["true_gravity_cosine_lower"])
    )
    timeout_q = P5ID._cayley_norm_upper_from_cos_lower(
        float(p1["timeout_handoff"]["combined_true_gravity_cosine_lower"])
    )
    nodes = {
        "normal": _node("normal", normal_q, bounds),
        "timeout": _node("timeout", timeout_q, bounds),
    }
    tests = {name: _node_word_test(node, H, domain) for name, node in nodes.items()}
    pass_ = all(t["outer_word_decrease_pass"] for t in tests.values())
    first = next((tests[n]["first_failure"] for n in ("normal", "timeout")
                  if tests[n]["first_failure"] != "NONE"), "NONE")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_OUTER_H_SOURCE_WORD_DECREASE_CERTIFICATE",
        "claim": "FINITE_ANGLE_ANISOTROPIC_OUTER_H_WORD_SUFFICIENT_DECREASE_TEST",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "local_P4_BW_recurrence_reused_as_outer_certificate": False,
        "metric_route": "SAME_NORMALIZED_CAYLEY_SOURCE_INFORMATION_GEOMETRY_AS_P4",
        "shipping_word_map_bound": True,
        "full_S_to_attitude_gain_retained": True,
        "exact_linear_coordinates_not_charged_as_nonlinearity": ["v", "p"],
        "S_innovation_treated_as_exact_linear_selector": True,
        "finite_angle_rotation_identity_used": True,
        "kalman_nonlinear_defect_lemma": "K^T(P_plus)^-1K <= R^-1",
        "candidate_outer_prefix_domain_requires_separate_bootstrap": True,
        "handoff_nodes": nodes,
        "node_word_tests": tests,
        "outer_to_P4_inner_overlap_target_W": float(H["certified_level_W"]),
        "outer_word_decrease_all_nodes": pass_,
        "P5_OUTER_H_WORD_CERTIFICATE": "PASS" if pass_ else "NOT_ESTABLISHED",
        "first_failure": first,
        "next_widening_if_not_established": (
            "replace the failed perturbative term by a validated exact large-angle source correction sector and a source-staged post-goLive covariance/cross-block prefix enclosure; do not enlarge the local P4 radius or drop S-to-attitude"
            if not pass_ else "compute finite outer-to-inner H funnel recursion"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("outer H certificate is not source generated")
    if d.get("source_replay_used") is not False:
        failures.append("outer H certificate uses replay")
    if d.get("filter_changed") is not False:
        failures.append("outer H certificate changes the filter")
    if d.get("local_P4_BW_recurrence_reused_as_outer_certificate") is not False:
        failures.append("outer H certificate reuses the retired local P4 recurrence")
    if d.get("full_S_to_attitude_gain_retained") is not True:
        failures.append("outer H certificate drops S-to-attitude")
    if d.get("S_innovation_treated_as_exact_linear_selector") is not True:
        failures.append("outer H certificate incorrectly charges S innovation as measurement nonlinearity")
    if d.get("finite_angle_rotation_identity_used") is not True:
        failures.append("outer H certificate does not use finite-angle rotation geometry")
    if d.get("kalman_nonlinear_defect_lemma") != "K^T(P_plus)^-1K <= R^-1":
        failures.append("outer H defect is not measured in the posterior information metric")
    if d.get("candidate_outer_prefix_domain_requires_separate_bootstrap") is not True:
        failures.append("outer H certificate silently assumes prefix-domain invariance")
    tests = d.get("node_word_tests", {})
    if set(tests) != {"normal", "timeout"}:
        failures.append("outer H normal/timeout nodes missing")
    for name in ("normal", "timeout"):
        t = tests.get(name, {})
        if not isinstance(t.get("P3_homogeneous_sqrt_decrease_lower"), (int, float)):
            failures.append(f"{name}: missing homogeneous decrease")
        if not isinstance(t.get("outer_vector_word_decrease_margin_lower"), (int, float)):
            failures.append(f"{name}: missing outer word margin")
        s = t.get("S_to_attitude_prefix", {})
        if s.get("full_S_to_attitude_gain_retained") is not True:
            failures.append(f"{name}: S-to-attitude gain omitted")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    failures = validate(out)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "P5_OUTER_H_WORD_CERTIFICATE": out.get("P5_OUTER_H_WORD_CERTIFICATE"),
        "first_failure": out.get("first_failure"),
        "normal": out.get("node_word_tests", {}).get("normal"),
        "timeout": out.get("node_word_tests", {}).get("timeout"),
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
