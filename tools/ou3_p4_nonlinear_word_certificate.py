#!/usr/bin/env python3
"""Validated exact-nonlinear H/A source-word certificate for OU-III P4.

P4 lifts the source-complete P3 homogeneous word certificate to the exact
shipping nonlinear MEKF map.  It uses no replay sampling and no product of
one-sample nonlinear contraction factors.

On theta<pi use the Cayley coordinate

    c(R)=2 tan(theta/2) u = 4 e_R/(1+tr R),
    z=[c;xi].

For each fixed-dimensional mode the sole quantitative metric is

    W_g = s_m z^T Sigma_KF(g)^-1 z,

where s_m>0 is one source-uniform constant shared by every source node in that
mode.  The certificate chooses s_m=Sigma_lambda_max_upper.  This normalization
changes neither physical level sets nor generalized contraction ratios, but it
keeps the very small P3 margin away from binary64 underflow.  Full
attitude--linear information cross terms are retained.

Every source operation is decomposed into its P3 homogeneous tangent map plus
an exact nonlinear defect.  The exact Cayley rotation formula gives a quadratic
vector-residual remainder; the deployed normalized polynomial quaternion branch
has a cubic Cayley defect; and exact Cayley multiplication supplies the finite
correction composition.  P3 proves homogeneous source-prefix information gain
at most one, so defects may be inserted at any admissible accepted/rejected
position and transported to the word endpoint without exponential branch
enumeration.

If C is a uniform Euclidean quadratic defect constant, Nop bounds the number of
state operations in the word, and m_- I <= M_g <= m_+ I, then on the bootstrap
W_s<=4 W_0,

    ||r_word||_{M_end} <= B W_0,
    B = 4 Nop sqrt(m_+) C / m_-.

The explicit choice sqrt(W_*) <= delta/(8B) makes the nonlinear endpoint defect
consume at most delta/8 in sqrt(W), closes every prefix inside 4W_0, and gives

    W_end <= (1-delta/2) W_0

as a direct positive-gap statement.  The implementation never forms 1-delta/2
when delta is below machine epsilon.  Since V_R<=||c||^2 on the certified chart,

    mu_W >= (delta/2) m_- > 0.

This is an inner theorem seed, not a practical-basin claim.  Source subdivision
can later widen W_* without changing this proof route.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import ou3_explicit_information_word_certificate as P3
import ou3_implementation_proof_manifest as MANIFEST
import ou3_p4_exact_word_map as WORDMAP
import ou3_full_process_ucc as PROCESS
import ou3_p4_group_algebra as GROUP
import ou3_p4_metric_defect_transport as TRANSPORT
import ou3_p4_node_metrics as METRIC
import ou3_source_domain_contract as SOURCE
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SCHEMA = 2

ROTATION_REMAINDER_COEFF = 0.75
SOURCE_SERIES_CAYLEY_CUBIC_COEFF = 0.085
PREFIX_BOOTSTRAP_W_FACTOR = 4.0
PROMOTED_CAYLEY_NORM_LIMIT = 1.0
PROMOTED_THETA_STAR_RAD = 1.0
MAX_STATE_OPERATIONS_PER_IMU_SAMPLE = 4  # predict, S, accel, async mag


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def add_up(a: float, b: float) -> float:
    return up(float(a) + float(b))


def mul_up(a: float, b: float) -> float:
    return up(float(a) * float(b))


def div_up(a: float, b: float) -> float:
    if not b > 0.0:
        raise ValueError("positive denominator required")
    return up(float(a) / float(b))


def mul_down(a: float, b: float) -> float:
    if a < 0.0 or b < 0.0:
        raise ValueError("mul_down requires nonnegative inputs")
    return down(float(a) * float(b))


def div_down(a: float, b: float) -> float:
    if a < 0.0 or not b > 0.0:
        raise ValueError("div_down requires nonnegative numerator and positive denominator")
    return down(float(a) / float(b))


def sqrt_up(x: float) -> float:
    return GROUP.sqrt_point(float(x)).hi


def _deployed_member_float(text: str, name: str) -> float:
    """Value of a deployed `float <name> = <literal>f;` member of the wrapper."""
    m = re.search(rf"float\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f\s*;", text)
    if m is None:
        raise RuntimeError(f"cannot extract deployed member {name}")
    return float(m.group(1))


def _source_measurement_bounds() -> dict:
    v = VECTOR.build()
    vf = VECTOR.validate(v)
    if vf:
        raise RuntimeError(f"vector source certificate failed: {vf}")
    c = v["configured_measurement_bounds"]
    acc_std = float(c["acc_measurement_std_mps2"])
    mag_std = float(c["mag_measurement_std_uT"])
    acc_var_lo = down(acc_std * acc_std)
    mag_var_lo = down(mag_std * mag_std)

    wrapper = WRAPPER.read_text(encoding="utf-8")
    rs_min = float(SOURCE.parse_const(wrapper, "MIN_R_S"))
    # The S=0 correction covariance is diag(rho_x r_S, rho_y r_S, r_S)^2, so its
    # smallest eigenvalue is bounded below by (min(rho_x, rho_y, 1) MIN_R_S)^2 --
    # the isotropic case this used to assert is the special case.
    # Read the factors from the deployed members so the bound follows the source
    # rather than a literal that has to be re-asserted on every retune.
    # The two horizontal axes are independent knobs, so the bound has to take
    # the smaller of them rather than assume one horizontal scale.
    rho_x = _deployed_member_float(wrapper, "R_S_x_factor_")
    rho_y = _deployed_member_float(wrapper, "R_S_y_factor_")
    for name_, value in (("R_S_x_factor_", rho_x), ("R_S_y_factor_", rho_y)):
        if not (0.0 < value <= 4.0):
            raise RuntimeError(f"deployed {name_} out of setter range: {value}")
    rho_h = min(rho_x, rho_y)
    rs_min_axis = mul_down(min(rho_h, 1.0), rs_min)
    rs_var_lo = down(rs_min_axis * rs_min_axis)
    rmin = min(acc_var_lo, mag_var_lo, rs_var_lo)
    if not (math.isfinite(rmin) and rmin > 0.0):
        raise RuntimeError("measurement covariance lower bound is not positive")
    return {
        "acc_measurement_std_mps2": acc_std,
        "mag_measurement_std_uT": mag_std,
        "acc_variance_lower": acc_var_lo,
        "mag_variance_lower": mag_var_lo,
        "S_zero_variance_lower": rs_var_lo,
        "all_correction_R_lambda_min_lower": rmin,
    }


def _source_configuration_checks(domain: dict, manifest: dict) -> list[str]:
    failures: list[str] = []
    runtime = domain.get("configured_runtime", {})
    if runtime.get("imu_lever_arm_enabled") is not False:
        failures.append("configured P4 runtime does not explicitly disable optional IMU lever arm")
    mekf = MEKF.read_text(encoding="utf-8")
    if "bool use_imu_lever_arm_ = false;" not in " ".join(mekf.split()):
        failures.append("source default no-lever-arm semantic changed")
    live = domain.get("normal_live", {})
    limit_source = float(manifest["mekf_defaults"]["acc_bias_limit_mps2"])
    limit_domain = float(live.get("active_accelerometer_bias_projection_limit_mps2", math.nan))
    center = float(live.get("active_accelerometer_bias_state_norm_upper_mps2", math.nan))
    if not (math.isfinite(center) and math.isfinite(limit_domain) and 0.0 <= center < limit_domain):
        failures.append("A-mode bias source node is not strictly inside projection ball")
    if not math.isclose(limit_source, limit_domain, rel_tol=0.0, abs_tol=1.0e-7):
        failures.append("declared A-mode projection limit differs from source acc_bias_limit_")
    return failures


def _composition_quadratic_constant(L: float, C_input: float, q_design: float) -> dict:
    """Bound deployed Cayley correction minus its homogeneous additive map."""
    Ld = add_up(L, mul_up(C_input, q_design))
    dmax = mul_up(Ld, q_design)
    if not dmax < 0.005:
        raise RuntimeError(f"design correction not safely inside source series branch: {dmax}")

    dmax2 = mul_up(dmax, dmax)
    cayley_ratio = add_up(1.0, mul_up(SOURCE_SERIES_CAYLEY_CUBIC_COEFF, dmax2))
    Lc = mul_up(Ld, cayley_ratio)
    dot_max = mul_up(Lc, mul_up(q_design, q_design))
    denom_lower = down(1.0 - up(0.25 * dot_max))
    if not denom_lower > 0.99:
        raise RuntimeError("Cayley product denominator is not safely positive on design radius")

    cross_coeff = div_up(mul_up(0.5, Lc), denom_lower)
    denom_effect_coeff = div_up(
        mul_up(0.25, mul_up(Lc, mul_up(add_up(1.0, Lc), q_design))),
        denom_lower,
    )
    cubic_coeff = mul_up(
        SOURCE_SERIES_CAYLEY_CUBIC_COEFF,
        mul_up(mul_up(Ld, mul_up(Ld, Ld)), q_design),
    )
    attitude_defect = add_up(
        C_input, add_up(cubic_coeff, add_up(cross_coeff, denom_effect_coeff))
    )
    full_state_defect = add_up(attitude_defect, C_input)
    return {
        "linear_correction_gain_L": L,
        "input_quadratic_defect_C": C_input,
        "design_error_norm_radius": q_design,
        "corrected_delta_norm_upper_at_design_radius": dmax,
        "cayley_product_denominator_lower": denom_lower,
        "attitude_quadratic_defect_constant_upper": attitude_defect,
        "full_state_quadratic_defect_constant_upper": full_state_defect,
    }


def _mode_certificate(mode: str, p3: dict, metric: dict, wordmap: dict,
                      domain: dict, measurement: dict) -> dict:
    row = p3["modes"][mode]
    m = metric["modes"][mode]
    delta = float(row["word_endpoint_relative_Riccati_injection_margin_lower"])
    smax = float(row["Sigma_lambda_max_upper"])
    smin = float(row["Sigma_lambda_min_lower"])
    mmin = float(m["metric_lambda_min_lower"])
    mmax = float(m["metric_lambda_max_upper"])
    metric_scale = float(m["mode_global_positive_scale"])
    if not (0.0 < mmin <= mmax < math.inf and metric_scale > 0.0):
        raise RuntimeError(f"{mode}: normalized Cayley metric bounds are invalid")
    sqrt_mmax = sqrt_up(mmax)

    rmin = float(measurement["all_correction_R_lambda_min_lower"])
    Kmax = sqrt_up(div_up(smax, rmin))

    live = domain["normal_live"]
    fmax = float(live["specific_force_norm_upper_mps2"])
    magmax = float(live["magnetic_vector_norm_upper_uT"])
    # Broad source-uniform full-H bound.  No attitude or S cross-gain is removed.
    Hmax = add_up(mul_up(2.0, max(fmax, magmax)), 2.0)
    Lcorr = mul_up(Kmax, Hmax)

    # R(c)-I-[c]x <= 3/4 |c|^2.  Acceleration also has the bilinear
    # attitude/a_w term, bounded symmetrically by 0.75(c^2+a_w^2).
    Cvec_acc = add_up(mul_up(ROTATION_REMAINDER_COEFF, fmax), 1.5)
    Cvec_mag = mul_up(ROTATION_REMAINDER_COEFF, magmax)
    Cvec = max(Cvec_acc, Cvec_mag)
    Cinput = mul_up(Kmax, Cvec)

    # Structured inputs for the metric-consistent transport route.  The P3
    # covariance upper is a Loewner diagonal dominator, so a block maximum over
    # its diagonal bounds lambda_max of that marginal block.
    diag_upper = [float(v) for v in row["matrix_comparison"]["Sigma_diagonal_upper"]]
    if len(diag_upper) != int(row["dimension"]):
        raise RuntimeError(f"{mode}: P3 covariance diagonal upper does not match mode dimension")
    sigma_attitude_upper = TRANSPORT._block_upper(diag_upper, range(0, 3))
    sigma_gyro_bias_upper = TRANSPORT._block_upper(diag_upper, range(3, 6))
    sigma_aw_upper = TRANSPORT._block_upper(diag_upper, range(15, 18))
    # The exact word defects are quadratic in the attitude, gyro-bias and a_w
    # coordinates only; the translation block never enters them.
    sigma_defect_input_upper = up(max(sigma_attitude_upper, sigma_gyro_bias_upper, sigma_aw_upper))
    proc = PROCESS.build()
    pf = PROCESS.validate(proc)
    if pf:
        raise RuntimeError(f"{mode}: process certificate failed: {pf}")
    ab = proc["attitude_gyro_bias"]
    q_theta_lower = float(ab["theta_diagonal_lower"])
    q_bias_lower = float(ab["gyro_bias_diagonal_lower"])
    cross_upper = float(ab["cross_norm_upper"])
    rho_attitude = down(1.0 - div_up(cross_upper, down(math.sqrt(q_theta_lower * q_bias_lower))))
    if not rho_attitude > 0.0:
        raise RuntimeError(f"{mode}: scaled attitude/bias process comparison lost positivity")
    H_attitude = up(max(fmax, magmax))

    q_design = min(1.0e-6, div_down(0.002, mul_up(2.0, max(Lcorr, 1.0))))
    if not q_design > 0.0:
        raise RuntimeError(f"{mode}: failed to obtain positive nonlinear design radius")
    corr = _composition_quadratic_constant(Lcorr, Cinput, q_design)

    dt = float(domain["configured_runtime"]["imu_dt_s"])
    omega = float(live["body_rate_norm_upper_deg_s"]) * math.pi / 180.0
    omega_dt = up(omega * dt)
    if not omega_dt < 0.01:
        raise RuntimeError("configured body rotation per sample exceeds P4 prediction envelope")
    Lpred = div_up(dt, down(1.0 - omega_dt))
    pred = _composition_quadratic_constant(Lpred, 0.0, q_design)

    Coperation = max(
        float(corr["full_state_quadratic_defect_constant_upper"]),
        float(pred["full_state_quadratic_defect_constant_upper"]),
    )
    samples = int(p3["source_word_binding"]["word_samples_upper_at_configured_dt"])
    operation_count = MAX_STATE_OPERATIONS_PER_IMU_SAMPLE * samples

    B_isotropic = mul_up(PREFIX_BOOTSTRAP_W_FACTOR, float(operation_count))
    B_isotropic = mul_up(B_isotropic, sqrt_mmax)
    B_isotropic = mul_up(B_isotropic, Coperation)
    B_isotropic = div_up(B_isotropic, mmin)
    if not (math.isfinite(B_isotropic) and B_isotropic > 0.0):
        raise RuntimeError(f"{mode}: nonlinear word defect gain is not finite positive")

    # Metric-consistent structured transport of the same word defect.  Both are
    # upper bounds on ||r_word||_M/W_0, so the certificate keeps the smaller one
    # and can never be widened by the refinement.
    transport = TRANSPORT.build({
        "metric_scale": metric_scale,
        "word_endpoint_delta_lower": delta,
        "correction_R_lambda_min_lower": rmin,
        "Sigma_attitude_upper": sigma_attitude_upper,
        "Sigma_gyro_bias_upper": sigma_gyro_bias_upper,
        "Sigma_defect_input_upper": sigma_defect_input_upper,
        "rho_attitude_scaled_lower": rho_attitude,
        "Q_theta_diagonal_lower": q_theta_lower,
        "H_attitude_norm_upper": H_attitude,
        "vector_residual_quadratic_constant_upper": Cvec,
        "prediction_increment_gain_upper": Lpred,
        "state_operation_count_upper": operation_count,
    })
    transport_failures = TRANSPORT.validate(transport)
    if transport_failures:
        raise RuntimeError(f"{mode}: metric defect transport failed: {transport_failures}")
    B_metric = float(transport["transported_word_defect_B_upper"])
    if B_metric <= B_isotropic:
        B, defect_route = B_metric, "METRIC_CONSISTENT_STRUCTURED_DEFECT_TRANSPORT"
    else:
        B, defect_route = B_isotropic, "ISOTROPIC_EUCLIDEAN_DEFECT_ENVELOPE"

    sqrt_W_star = div_down(delta, mul_up(8.0, B))
    W_star = mul_down(sqrt_W_star, sqrt_W_star)
    if not W_star > 0.0:
        raise RuntimeError(f"{mode}: certified P4 level underflowed to zero")

    q_prefix = mul_up(2.0, sqrt_up(div_up(W_star, mmin)))
    if not q_prefix < q_design:
        raise RuntimeError(f"{mode}: final level does not close nonlinear design-radius bootstrap")
    if not q_prefix < PROMOTED_CAYLEY_NORM_LIMIT:
        raise RuntimeError(f"{mode}: final prefix reaches Cayley chart boundary")

    correction_prefix = add_up(
        mul_up(Lcorr, q_prefix), mul_up(Cinput, mul_up(q_prefix, q_prefix))
    )
    if not correction_prefix < 1.0e-2:
        raise RuntimeError(f"{mode}: source correction branch cannot be fixed to polynomial path")

    nonlinear_sqrt_fraction = div_up(mul_up(B, sqrt_W_star), delta)
    if not nonlinear_sqrt_fraction <= 0.125000000000001:
        raise RuntimeError(f"{mode}: nonlinear endpoint budget exceeds delta/8")

    relative_decrease = mul_down(0.5, delta)
    mu = mul_down(relative_decrease, mmin)
    if not mu > 0.0:
        raise RuntimeError(f"{mode}: exact nonlinear mu_W did not remain positive")

    projection = None
    if mode == "A":
        center = float(live["active_accelerometer_bias_state_norm_upper_mps2"])
        limit = float(live["active_accelerometer_bias_projection_limit_mps2"])
        margin = down(limit - center)
        projection = {
            "source_center_norm_upper_mps2": center,
            "shipping_projection_limit_mps2": limit,
            "interior_margin_lower_mps2": margin,
            "certified_error_norm_prefix_upper": q_prefix,
            "projection_surface_reached_in_certified_funnel": not (q_prefix < margin),
            "exact_projection_branch_in_certified_funnel": "identity_interior_branch",
        }
        if not q_prefix < margin:
            raise RuntimeError("A: certified P4 funnel reaches accelerometer-bias projection surface")

    horizon = float(p3["source_word_binding"]["word_horizon_lower_s"])
    return {
        "mode": mode,
        "dimension": row["dimension"],
        "source_complete": True,
        "outward_rounded": True,
        "joint_source_reachability": True,
        "one_sample_decrease_used": False,
        "source_replay_used": False,
        "word_horizon_s": horizon,
        "word_endpoint_relative_Riccati_injection_margin_lower": delta,
        "Sigma_lambda_min_lower": smin,
        "Sigma_lambda_max_upper": smax,
        "prefix_information_gain_upper": float(row["prefix_information_gain_upper"]),
        "path_metric": m,
        "metric_mode_global_positive_scale": metric_scale,
        "P3_word_endpoint_delta_lower": delta,
        "P3_homogeneous_prefix_information_gain_upper": 1.0,
        "metric_lambda_min_lower": mmin,
        "metric_lambda_max_upper": mmax,
        "measurement_bounds": measurement,
        "full_gain_norm_upper": Kmax,
        "measurement_linear_operator_norm_upper": Hmax,
        "vector_residual_quadratic_constant_upper": Cvec,
        "correction_quadratic_bound": corr,
        "prediction_quadratic_bound": pred,
        "uniform_operation_quadratic_defect_constant_upper": Coperation,
        "word_samples_upper": samples,
        "state_operation_count_upper": operation_count,
        "transported_word_defect_B_upper": B,
        "transported_word_defect_B_isotropic_upper": B_isotropic,
        "transported_word_defect_B_metric_consistent_upper": B_metric,
        "transported_word_defect_route": defect_route,
        "metric_consistent_defect_transport": transport,
        "Sigma_attitude_block_upper": sigma_attitude_upper,
        "Sigma_defect_input_block_upper": sigma_defect_input_upper,
        "certified_attitude_cayley_radius_upper": mul_up(
            transport["attitude_chart_scale"], mul_up(2.0, sqrt_W_star)
        ),
        "nonlinear_sqrt_budget_fraction_of_delta_upper": nonlinear_sqrt_fraction,
        "certified_level_W": W_star,
        "certified_level_sqrt_W": sqrt_W_star,
        "prefix_W_factor_upper": PREFIX_BOOTSTRAP_W_FACTOR,
        "prefix_canonical_error_norm_upper": q_prefix,
        "theta_star": PROMOTED_THETA_STAR_RAD,
        "cayley_norm_limit": PROMOTED_CAYLEY_NORM_LIMIT,
        "theta_bound_reason": "theta=2 atan(||c||/2)<||c||; certified prefixes have ||c||<=||z||<1",
        "accepted_correction_norm_prefix_upper": correction_prefix,
        "accepted_correction_uses_source_series_branch": True,
        "all_word_prefixes_safe": True,
        "endpoint_relative_W_decrease_lower": relative_decrease,
        "endpoint_ratio_not_formed_below_machine_epsilon": True,
        "mu_W_lower": mu,
        "mu_W_denominator": "V_R(R_e)+||xi||^2",
        "mu_W_reason": "W0-W1 >= (delta/2)W0 >= (delta/2)m_minus*(V_R+||xi||^2)",
        "active_bias_projection": projection,
        "exact_nonlinear_word_pass": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P4 theorem domain is trajectory fitted")

    p3 = P3.build(domain_path)
    p3f = P3.validate(p3)
    manifest = MANIFEST.build()
    mf = MANIFEST.validate(manifest)
    wordmap = WORDMAP.build(domain_path)
    wf = WORDMAP.validate(wordmap)
    metric = METRIC.build(domain_path)
    metf = METRIC.validate(metric)
    failures = [f"P3: {x}" for x in p3f] + [f"manifest: {x}" for x in mf]
    failures += [f"word-map: {x}" for x in wf] + [f"metric: {x}" for x in metf]
    failures += _source_configuration_checks(domain, manifest)
    measurement = _source_measurement_bounds()

    modes = {}
    if not failures:
        for mode in ("H", "A"):
            try:
                modes[mode] = _mode_certificate(mode, p3, metric, wordmap, domain, measurement)
            except Exception as exc:
                failures.append(f"{mode}: {exc}")

    passed = not failures and all(
        modes.get(m, {}).get("exact_nonlinear_word_pass") is True for m in ("H", "A")
    )
    return {
        "schema": SCHEMA,
        "qualification": "VALIDATED_EXACT_NORMALIZED_CAYLEY_NONLINEAR_SOURCE_WORD_CERTIFICATE",
        "claim": "P4_EXACT_NONLINEAR_H_A_WORD_DISSIPATION_AND_PREFIX_SAFETY",
        "source_generated_not_trajectory_fit": True,
        "validated_arithmetic": True,
        "outward_rounded": True,
        "source_replay_used": False,
        "same_source_word_language_as_P3": True,
        "exact_deployed_quaternion_injection": True,
        "full_S_to_attitude_cross_gain": True,
        "metric_route": "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC",
        "mode_global_information_normalization": True,
        "block_diagonal_metric_fallback": False,
        "word_branch_coverage": {
            "method": "uniform exact-operation quadratic defect envelope plus P3 unit segment information transport",
            "explicit_exponential_branch_enumeration_required": False,
            "all_admissible_accelerometer_accept_reject_branches": True,
            "all_admissible_magnetometer_not_due_accept_reject_branches": True,
            "all_admissible_S_not_due_due_branches": True,
            "all_admissible_aw_sync_not_due_due_branches": True,
        },
        "source_subdivision": {
            "kind": "ANALYTIC_MONOTONE_SOURCE_UNIFORM_ENVELOPE",
            "reason_no_cartesian_subdivision_needed": "nonlinear defect constants use source-global physical vector, covariance and measurement-noise bounds; P3 supplies joint-source endpoint and segment inequalities",
            "future_widening_allowed": "joint source-node subdivision may reduce C and enlarge W_star but is not a separate theorem route",
        },
        "modes": modes,
        "P4_EXACT_NONLINEAR_WORD_CERTIFICATE": "PASS" if passed else "FAIL",
        "theorem_promotion": "P4_NORMAL_LIVE_EXACT_WORDS" if passed else "NOT_ESTABLISHED",
        "failures": failures,
        "next_obligation": "P5 finite startup-to-inner-funnel capture, then P6 hybrid jumps",
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True or d.get("source_replay_used") is not False:
        failures.append("P4 certificate is not source-only")
    if d.get("exact_deployed_quaternion_injection") is not True:
        failures.append("P4 does not use deployed quaternion injection")
    if d.get("full_S_to_attitude_cross_gain") is not True:
        failures.append("P4 discarded S-to-attitude gain")
    if d.get("metric_route") != "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC":
        failures.append("P4 metric route is not exact Cayley information lift")
    if d.get("mode_global_information_normalization") is not True:
        failures.append("P4 does not use mode-global geometry-preserving information normalization")
    if d.get("block_diagonal_metric_fallback") is not False:
        failures.append("retired block metric remains as fallback")

    coverage = d.get("word_branch_coverage", {})
    for key in (
        "all_admissible_accelerometer_accept_reject_branches",
        "all_admissible_magnetometer_not_due_accept_reject_branches",
        "all_admissible_S_not_due_due_branches",
        "all_admissible_aw_sync_not_due_due_branches",
    ):
        if coverage.get(key) is not True:
            failures.append(f"P4 branch coverage missing {key}")

    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        if m.get("source_complete") is not True or m.get("joint_source_reachability") is not True:
            failures.append(f"{mode}: source word is not complete/joint")
        if m.get("one_sample_decrease_used") is not False:
            failures.append(f"{mode}: one-sample decrease reintroduced")
        pm = m.get("path_metric", {})
        if pm.get("kind") != "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC":
            failures.append(f"{mode}: wrong P4 path metric")
        if pm.get("same_scale_on_every_source_node_in_mode") is not True:
            failures.append(f"{mode}: metric scaling is not source-node uniform")
        if pm.get("full_attitude_linear_cross_terms_retained") is not True:
            failures.append(f"{mode}: information cross terms were discarded")
        if not (isinstance(m.get("theta_star"), (int,float)) and 0.0 < float(m["theta_star"]) < math.pi):
            failures.append(f"{mode}: theta_star is not in (0,pi)")
        if m.get("all_word_prefixes_safe") is not True:
            failures.append(f"{mode}: word prefix safety failed")
        if m.get("accepted_correction_uses_source_series_branch") is not True:
            failures.append(f"{mode}: exact source correction branch not certified")
        if not (isinstance(m.get("certified_level_W"),(int,float)) and float(m["certified_level_W"]) > 0.0):
            failures.append(f"{mode}: certified W level is not positive")
        if not (isinstance(m.get("endpoint_relative_W_decrease_lower"),(int,float)) and float(m["endpoint_relative_W_decrease_lower"]) > 0.0):
            failures.append(f"{mode}: endpoint W decrease is not positive")
        if not (isinstance(m.get("mu_W_lower"),(int,float)) and float(m["mu_W_lower"]) > 0.0):
            failures.append(f"{mode}: mu_W is not positive")
        if not float(m.get("prefix_canonical_error_norm_upper", math.inf)) < float(m.get("cayley_norm_limit", 0.0)):
            failures.append(f"{mode}: prefix chart bootstrap did not close")
        if not float(m.get("accepted_correction_norm_prefix_upper", math.inf)) < 1.0e-2:
            failures.append(f"{mode}: accepted correction can cross source quaternion branch")
        if mode == "A":
            p = m.get("active_bias_projection", {})
            if p.get("projection_surface_reached_in_certified_funnel") is not False:
                failures.append("A: certified funnel reaches nonsmooth bias projection surface")
    if not failures and d.get("P4_EXACT_NONLINEAR_WORD_CERTIFICATE") != "PASS":
        failures.append("P4 status is not PASS")
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
    compact = {
        mode: {
            "W_star": out.get("modes", {}).get(mode, {}).get("certified_level_W"),
            "mu_W_lower": out.get("modes", {}).get(mode, {}).get("mu_W_lower"),
            "endpoint_relative_W_decrease_lower": out.get("modes", {}).get(mode, {}).get("endpoint_relative_W_decrease_lower"),
            "prefix_canonical_error_norm_upper": out.get("modes", {}).get(mode, {}).get("prefix_canonical_error_norm_upper"),
        }
        for mode in ("H", "A")
    }
    print(json.dumps({
        "P4_EXACT_NONLINEAR_WORD_CERTIFICATE": out["P4_EXACT_NONLINEAR_WORD_CERTIFICATE"],
        "numerical": compact,
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
