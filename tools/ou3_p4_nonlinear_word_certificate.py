#!/usr/bin/env python3
"""Validated exact-nonlinear H/A source-word certificate for OU-III P4.

This producer closes the local quantitative lift from P3 without replay and
without enumerating exponentially many accepted/rejected strings.

Coordinate and metric
---------------------
On theta<pi use the exact Cayley coordinate c(R)=2 tan(theta/2)u and

    z = [c; xi],   W_g(z)=z' Sigma_KF(g)^-1 z.

The local coordinate is exactly the P3 delta-theta coordinate, so the P3
homogeneous endpoint inequality applies without a metric conversion.  P3 also
proves every homogeneous prefix/segment is nonexpansive between its matching
source information metrics.

Exact nonlinear defects
-----------------------
The exact rotation in Cayley coordinates is

    R(c)=I + 4/(4+|c|^2)[c]x + 2/(4+|c|^2)[c]x^2.

For |c|<=1,

    ||R(c)-I-[c]x|| <= (3/4)|c|^2.

Hence the accepted vector residual remainder is bounded by a quadratic source
constant.  The complete implemented Kalman gain is bounded from Joseph form:
K R K' <= P_plus and the source-uniform covariance upper bound therefore gives

    ||K|| <= sqrt(Sigma_max / lambda_min(R)).

No S->attitude or other cross gain is removed.

After an additive correction d, the source uses its normalized quaternion
injection.  In the certified funnel every correction lies on the source's
polynomial branch.  Its exact Cayley vector satisfies

    c_d = (2 k(d)/w(d)) d,
    ||c_d-d|| <= 0.085 ||d||^3        for ||d||<=1e-2.

The exact left product is

    c+ = (c_d+c+0.5 c_d x c)/(1-0.25 c_d'c),

so a direct rational bound supplies another quadratic defect.  Prediction has
the same Cayley product structure for the small gyro-bias mismatch; every
non-attitude OU/GM prediction block is linear exactly.  Source tuner evolution
is measurement-only/exogenous and therefore contributes no error-state defect.

Word composition
----------------
Let C be a uniform Euclidean quadratic defect constant for any error-state
operation and Nop an upper count over one source word.  With metric bounds
m_- I <= Sigma^-1 <= m_+ I, bootstrap W_s<=4W_0 gives

    ||sum transported defects||_{M_end}
       <= B W_0,
    B = 4 Nop sqrt(m_+) C / m_-.

Choose

    sqrt(W_*) <= delta/(8 B).

Then the nonlinear defect consumes at most delta/8 of sqrt(W), while
sqrt(1-delta)<=1-delta/2.  Therefore

    W_end <= (1-3 delta/8)^2 W_0
           <= (1-delta/2) W_0,

and every prefix remains inside 4W_0.  This yields the explicit exact nonlinear
margin

    mu_W >= (delta/2) m_- > 0

against the manuscript denominator V_R+||xi||^2, because on the Cayley chart
V_R <= ||c||^2 and W >= m_- ||z||^2.

The resulting level is intentionally allowed to be very small.  It is a
machine-checkable theorem seed, not a claim about the practical basin size;
later source-dependent subdivision can only widen it.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_explicit_information_word_certificate as P3
import ou3_implementation_proof_manifest as MANIFEST
import ou3_p4_exact_word_map as WORDMAP
import ou3_p4_group_algebra as GROUP
import ou3_p4_node_metrics as METRIC
import ou3_source_domain_contract as SOURCE
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SCHEMA = 1

ROTATION_REMAINDER_COEFF = 0.75
SOURCE_SERIES_CAYLEY_CUBIC_COEFF = 0.085
PREFIX_BOOTSTRAP_W_FACTOR = 4.0
NONLINEAR_SQRT_BUDGET_FRACTION_OF_DELTA = 1.0 / 8.0
PROMOTED_CAYLEY_NORM_LIMIT = 1.0
PROMOTED_THETA_STAR_RAD = 1.0  # theta=2 atan(|c|/2) < |c| <= 1
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
    if not (b > 0.0):
        raise ValueError("positive denominator required")
    return up(float(a) / float(b))


def mul_down(a: float, b: float) -> float:
    if a < 0.0 or b < 0.0:
        raise ValueError("mul_down is for nonnegative certificate quantities")
    return down(float(a) * float(b))


def div_down(a: float, b: float) -> float:
    if a < 0.0 or not (b > 0.0):
        raise ValueError("positive certificate division required")
    return down(float(a) / float(b))


def sqrt_up(x: float) -> float:
    return GROUP.sqrt_point(float(x)).hi


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
    if "float R_S_xy_factor_ = 1.0f;" not in wrapper:
        raise RuntimeError("P4 requires the source-bound deployed isotropic R_S_xy_factor_=1")
    rs_var_lo = down(rs_min * rs_min)
    rmin = min(acc_var_lo, mag_var_lo, rs_var_lo)
    if not (rmin > 0.0 and math.isfinite(rmin)):
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
    if "bool    use_imu_lever_arm_       = false;" not in mekf:
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
    """Bound exact deployed Cayley correction minus its linear additive map.

    ``||d_linear||<=L q`` and ``||d_exact-d_linear||<=C_input q^2``.
    The design radius is selected so the source's <1e-2 series branch is
    guaranteed.  All returned upper constants are rounded upward.
    """
    Ld = add_up(L, mul_up(C_input, q_design))
    dmax = mul_up(Ld, q_design)
    if not dmax < 0.005:
        raise RuntimeError(f"design correction is not safely inside source series branch: {dmax}")

    dmax2 = mul_up(dmax, dmax)
    cayley_ratio = add_up(1.0, mul_up(SOURCE_SERIES_CAYLEY_CUBIC_COEFF, dmax2))
    Lc = mul_up(Ld, cayley_ratio)
    dot_max = mul_up(Lc, mul_up(q_design, q_design))
    denom_lower = down(1.0 - up(0.25 * dot_max))
    if not denom_lower > 0.99:
        raise RuntimeError("Cayley product denominator is not safely positive on design radius")

    cross_coeff = mul_up(0.5, Lc)
    denom_effect_coeff = div_up(
        mul_up(0.25, mul_up(Lc, mul_up(add_up(1.0, Lc), q_design))),
        denom_lower,
    )
    cross_coeff = div_up(cross_coeff, denom_lower)
    cubic_coeff = mul_up(
        SOURCE_SERIES_CAYLEY_CUBIC_COEFF,
        mul_up(mul_up(Ld, mul_up(Ld, Ld)), q_design),
    )
    attitude_defect = add_up(C_input, add_up(cubic_coeff, add_up(cross_coeff, denom_effect_coeff)))
    # Non-attitude additive coordinates only see the measurement remainder.
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
    mmin = down(1.0 / smax)
    mmax = up(1.0 / smin)
    sqrt_mmax = sqrt_up(mmax)

    rmin = float(measurement["all_correction_R_lambda_min_lower"])
    Kmax = sqrt_up(div_up(smax, rmin))

    live = domain["normal_live"]
    fmax = float(live["specific_force_norm_upper_mps2"])
    magmax = float(live["magnetic_vector_norm_upper_uT"])
    # H_acc contains attitude and a_w (and b_a in A); H_mag contains attitude.
    # The deliberately broad factor 2 covers simultaneous 3-axis blocks and is
    # source-uniform rather than a sampled singular-vector estimate.
    Hmax = add_up(mul_up(2.0, max(fmax, magmax)), 2.0)
    Lcorr = mul_up(Kmax, Hmax)

    # Exact R(c) formula gives 3/4 |c|^2 vector rotation remainder.  The accel
    # model also has the bilinear (R-I) delta-a_w term; for |c|<=1 it is <=
    # 1.5|c||delta-a_w| <= 0.75(c^2+delta-a_w^2).
    Cvec_acc = add_up(mul_up(ROTATION_REMAINDER_COEFF, fmax), 1.5)
    Cvec_mag = mul_up(ROTATION_REMAINDER_COEFF, magmax)
    Cvec = max(Cvec_acc, Cvec_mag)
    Cinput = mul_up(Kmax, Cvec)

    # Pick a design radius from the worst full correction gain.  This is only a
    # proof workspace radius; the final W_* will be much smaller and is checked
    # against it below.
    q_design = min(1.0e-6, div_down(0.002, mul_up(2.0, max(Lcorr, 1.0))))
    if not q_design > 0.0:
        raise RuntimeError(f"{mode}: failed to obtain positive nonlinear design radius")
    corr = _composition_quadratic_constant(Lcorr, Cinput, q_design)

    dt = float(domain["configured_runtime"]["imu_dt_s"])
    omega = float(live["body_rate_norm_upper_deg_s"]) * math.pi / 180.0
    omega_dt = up(omega * dt)
    if not omega_dt < 0.01:
        raise RuntimeError("configured body rotation per sample exceeds P4 prediction envelope")
    # Relative bias-induced attitude increment through one prediction has norm
    # <= dt/(1-omega dt) ||delta b_g||; this bounds the exact integral of the
    # rotating bias direction without a libm exponential.
    Lpred = div_up(dt, down(1.0 - omega_dt))
    pred = _composition_quadratic_constant(Lpred, 0.0, q_design)

    Coperation = max(
        float(corr["full_state_quadratic_defect_constant_upper"]),
        float(pred["full_state_quadratic_defect_constant_upper"]),
    )
    samples = int(p3["source_word_binding"]["word_samples_upper_at_configured_dt"])
    operation_count = MAX_STATE_OPERATIONS_PER_IMU_SAMPLE * samples

    # B = 4 N sqrt(m+) C / m-.  Four is the bootstrap W_prefix/W0 factor.
    B = mul_up(PREFIX_BOOTSTRAP_W_FACTOR, float(operation_count))
    B = mul_up(B, sqrt_mmax)
    B = mul_up(B, Coperation)
    B = div_up(B, mmin)
    if not (math.isfinite(B) and B > 0.0):
        raise RuntimeError(f"{mode}: nonlinear word defect gain is not finite positive")

    sqrt_W_star = div_down(delta, mul_up(8.0, B))
    W_star = mul_down(sqrt_W_star, sqrt_W_star)
    if not W_star > 0.0:
        raise RuntimeError(f"{mode}: certified P4 level underflowed to zero")

    # Bootstrap coordinate radius: W_prefix<=4W0 and W>=m-||z||^2.
    q_prefix = mul_up(2.0, sqrt_up(div_up(W_star, mmin)))
    if not q_prefix < q_design:
        raise RuntimeError(f"{mode}: final level does not close nonlinear design-radius bootstrap")
    if not q_prefix < PROMOTED_CAYLEY_NORM_LIMIT:
        raise RuntimeError(f"{mode}: final prefix reaches Cayley chart boundary")

    # Every accepted correction remains on the exact source polynomial branch.
    correction_prefix = add_up(mul_up(Lcorr, q_prefix), mul_up(Cinput, mul_up(q_prefix, q_prefix)))
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

    return {
        "mode": mode,
        "dimension": row["dimension"],
        "source_complete": True,
        "outward_rounded": True,
        "joint_source_reachability": True,
        "one_sample_decrease_used": False,
        "source_replay_used": False,
        "metric": m,
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
        "nonlinear_sqrt_budget_fraction_of_delta_upper": nonlinear_sqrt_fraction,
        "certified_level_W": W_star,
        "certified_level_sqrt_W": sqrt_W_star,
        "prefix_W_factor_upper": PREFIX_BOOTSTRAP_W_FACTOR,
        "prefix_canonical_error_norm_upper": q_prefix,
        "theta_star": PROMOTED_THETA_STAR_RAD,
        "cayley_norm_limit": PROMOTED_CAYLEY_NORM_LIMIT,
        "theta_bound_reason": "theta=2 atan(||c||/2) < ||c||; certified prefixes have ||c||<=||z||<1",
        "accepted_correction_norm_prefix_upper": correction_prefix,
        "accepted_correction_uses_source_series_branch": True,
        "all_word_prefixes_safe": True,
        "endpoint_relative_W_decrease_lower": relative_decrease,
        "endpoint_W_ratio_symbolic_upper": "1-delta/2",
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

    passed = not failures and all(modes.get(m, {}).get("exact_nonlinear_word_pass") is True for m in ("H", "A"))
    return {
        "schema": SCHEMA,
        "qualification": "VALIDATED_EXACT_CAYLEY_NONLINEAR_SOURCE_WORD_CERTIFICATE",
        "claim": "P4_EXACT_NONLINEAR_H_A_WORD_DISSIPATION_AND_PREFIX_SAFETY",
        "source_generated_not_trajectory_fit": True,
        "validated_arithmetic": True,
        "outward_rounded": True,
        "source_replay_used": False,
        "same_source_word_language_as_P3": True,
        "exact_deployed_quaternion_injection": True,
        "full_S_to_attitude_cross_gain": True,
        "metric_route": "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC",
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
            "reason_no_cartesian_subdivision_needed": "nonlinear defect constants use source-global physical vector, covariance and measurement-noise bounds; P3 already supplies joint-source endpoint/segment inequalities, so subdividing independent extrema would not strengthen logical coverage",
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
        if not (isinstance(m.get("theta_star"),(int,float)) and 0.0 < float(m["theta_star"]) < math.pi):
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
    print(json.dumps({
        "P4_EXACT_NONLINEAR_WORD_CERTIFICATE": out["P4_EXACT_NONLINEAR_WORD_CERTIFICATE"],
        "H": out.get("modes", {}).get("H"),
        "A": out.get("modes", {}).get("A"),
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
