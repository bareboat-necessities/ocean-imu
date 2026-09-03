#!/usr/bin/env python3
"""Source-complete H=18/A=21 nonlinear whole-word dissipation producer.

This is the numerical complete-word producer required by the frozen canonical
P4 gate.  It consumes only already certified artifacts and evaluates the
retained source-path theorem of ``w3d-stability-widening-source-path``:

* canonical P3 supplies, on the same source-complete word, the information-word
  inequality  Phi_w' Sigma_{i_r}^-1 Phi_w <= (1-delta_m) Sigma_{i_0}^-1  through
  its relative Riccati injection margin ``delta_m``;
* every source-reachable prefix is nonexpansive in its matching source
  information metrics, so intermediate samples need no decrease;
* the exact Cayley lift, the homogeneous vector remainder sector, the exact
  Joseph/quaternion/reset transport and the finite-speed committed sample-clock
  graph bound every implemented operation defect by a source-uniform quadratic
  form; and
* the metric normalization is the retained one, ``M_i = s_m Sigma_i^-1`` with
  ``s_m = lambda_max(Sigma)_m``, hence ``m_- = 1`` exactly and
  ``m_+ = lambda_max(Sigma)_m / lambda_min(Sigma)_m``.

Writing the exact zero-input word map as ``z_{s+1} = A_s z_s + r_s(z_s)`` with
``||r_s(z)|| <= C ||z||^2`` on the declared operation region, ``N_op`` state
operations per word and the ``W_s <= 4 W_0`` prefix bootstrap give

    ||r_w||_{M} <= B_m W_0,     B_m = 4 N_op sqrt(m_+) C / m_-,
    sqrt(W_*)   = delta_m / (8 B_m),

and on ``W_0 <= W_*``

    sqrt(W_out) <= (sqrt(1-delta_m) + delta_m/8) sqrt(W_0)
                <= (1 - 3 delta_m/8) sqrt(W_0),
    W_out       <= (1 - delta_m/2) W_0.

So the certified whole-word contraction factor is ``rho_m = 1 - delta_m/2`` and
the direct margin is ``mu_{W,m} >= (delta_m/2) m_-``.

Representation.  ``delta_m`` is far below binary64 epsilon, so ``1 - delta_m/2``
is not representable strictly below one.  Following the retained verification
contract this producer reports the *positive strict gap* ``delta_m/2`` in
``one_minus_rho_{H,A}_lower`` and never rounds the strict gap back to one.
``rho_{H,A}_upper`` is the correctly rounded-up companion value.

Quadratic defect budget.  ``C`` is the maximum over the four implemented state
operations of one sample (prediction, due S=0, accepted accelerometer, accepted
magnetometer); rejected/not-due branches are identity corrections and are
covered by the same maximum.  Each coefficient is outward rounded from source
envelopes only:

* prediction -- exact Cayley left composition against its linear target,
  ``[ |c_d||c|/2 + |c_d||c|(|c|+|c_d|)/4 ] / (1 - |c_d||c|/4)`` with
  ``|c_d| <= coeff_hi * dt * |delta b_g|``; plus the exact rotation residual
  ``||R(c)-I-[c]x|| <= (3/4)||c||^2`` acting on the specific-force envelope, and
  the bilinear ``dt ||R-I|| |delta b_a|`` term;
* vector measurements -- the homogeneous remainder sector gives
  ``||eta|| <= (||c||^2/2) ||v||`` for the pure rotational residual, and the
  correction charges it through the full gain bound
  ``||K|| <= sqrt(lambda_max(Sigma)/lambda_min(R))``;
* every accepted correction -- the exact quaternion/reset defect of
  ``ou3_p4_exact_reset_transport`` evaluated with ``|d| <= G_v ||z||``.

The cubic-to-quadratic reduction uses the shipping polynomial-quaternion branch
radius ``delta_ref = 1e-2`` as the operation-region scale, and the certificate
is only claimed after verifying that every accepted correction on the certified
funnel is strictly inside that radius, inside the Cayley chart, and inside the
A-mode accelerometer-bias projection interior.

This producer does not promote anything.  The unique canonical P4 gate decides.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_exact_reset_transport as RESET
import ou3_p4_p3_metric_attachment as METRIC
import ou3_p4_sample_clock_source_refinement as CLOCK
import ou3_p4_signed_joseph_feasibility as AUDIT
import ou3_p4_source_word_timing as TIMING
import ou3_p4_vector_remainder_sector as REMAINDER
import ou3_source_reachable_matrix_p3 as BASE
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_NONLINEAR_WORD_DISSIPATION_V1"
P3_QUALIFICATION = "OU3_P3_CANONICAL_THEOREM_INTERFACE"
MODES = ("H", "A")
DIMENSION = {"H": 18, "A": 21}
STATE_GROUPS = METRIC.STATE_GROUPS
# Prediction, due S=0 pseudo-update, accepted accelerometer, accepted
# magnetometer.  Rejected and not-due branches are identity state corrections.
OPERATIONS_PER_SAMPLE = 4
OPERATION_CLASSES = ("prediction", "S_zero_pseudo", "accelerometer", "magnetometer")
# Operation-region correction radius.  Strictly inside the shipping polynomial
# correction-quaternion branch (|d| < 1e-2), used only as the scale of the
# cubic-to-quadratic reduction; the certified funnel is verified to be strictly
# inside it.
CORRECTION_REGION_RADIUS = 1.0e-3
if not CORRECTION_REGION_RADIUS < RESET.SERIES_BRANCH_NORM:
    raise RuntimeError("correction region left the shipping polynomial quaternion branch")


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def add_up(*xs: float) -> float:
    y = 0.0
    for x in xs:
        y = up(y + float(x))
    return y


def mul_up(a: float, b: float) -> float:
    return up(float(a) * float(b))


def div_up(a: float, b: float) -> float:
    if not float(b) > 0.0:
        raise RuntimeError("positive denominator required")
    return up(float(a) / float(b))


def sqrt_up(x: float) -> float:
    if not (math.isfinite(float(x)) and float(x) >= 0.0):
        raise RuntimeError("finite nonnegative square-root input required")
    return up(math.sqrt(float(x)))


def _finite_positive(x) -> bool:
    return (
        not isinstance(x, bool)
        and isinstance(x, (int, float))
        and math.isfinite(float(x))
        and float(x) > 0.0
    )


def _covariance_operator_bounds(metric: dict, mode: str) -> dict:
    """Source-uniform operator bounds on Sigma over the certified metric family.

    For a PSD matrix partitioned into ``n_g`` state groups, block
    Cauchy-Schwarz gives ``Sigma <= n_g blockdiag(Sigma_gg)``, so the retained
    same-history marginal ceilings yield an operator upper.  The retained
    precision-block join already supplies a genuine matrix lower.
    """
    groups = STATE_GROUPS[mode]
    n_g = len(groups)
    upper = {g: 0.0 for g in groups}
    envelopes = 0
    for row in metric.get("endpoint_rows", []):
        for env_key in ("boundary_history_envelope", "positive_phase_history_envelope"):
            env = row[env_key]
            envelopes += 1
            trans = env["translation_covariance_upper_groups"]
            for g in ("v", "p", "S", "a_w"):
                upper[g] = max(upper[g], float(trans[g]))
            bias = env[f"{mode}_bias_covariance_upper"]
            upper["theta"] = max(upper["theta"], float(bias["theta_covariance_upper"]))
            upper["b_g"] = max(upper["b_g"], float(bias["gyro_bias_covariance_upper"]))
            if mode == "A":
                upper["b_a"] = max(upper["b_a"], float(bias["accel_bias_covariance_upper"]))
    if envelopes != 2 * 800 or any(not _finite_positive(x) for x in upper.values()):
        raise RuntimeError(f"{mode}: same-history covariance ceiling family is incomplete")

    lower_row = metric["global_source_phase_covariance_lower_group_diagonal"][mode]
    if list(lower_row.keys()) != groups:
        raise RuntimeError(f"{mode}: global covariance lower group order changed")
    lam_min = down(min(float(x) for x in lower_row.values()))
    lam_max = mul_up(float(n_g), max(upper.values()))
    if not (_finite_positive(lam_min) and _finite_positive(lam_max) and lam_min <= lam_max):
        raise RuntimeError(f"{mode}: covariance operator bounds lost strict positivity")
    return {
        "state_groups": groups,
        "group_covariance_upper_envelope": upper,
        "group_covariance_lower": {g: float(lower_row[g]) for g in groups},
        "block_cauchy_schwarz_group_factor": n_g,
        "covariance_lambda_max_upper": lam_max,
        "covariance_lambda_min_lower": lam_min,
        "same_history_envelopes_scanned": envelopes,
    }


def _measurement_variance_lowers(vector: dict) -> dict:
    """Configured measurement covariance lowers used by the full gain bound."""
    c = vector["configured_measurement_bounds"]
    acc_std = float(c["acc_measurement_std_mps2"])
    mag_std = float(c["mag_measurement_std_uT"])
    if not (_finite_positive(acc_std) and _finite_positive(mag_std)):
        raise RuntimeError("configured vector measurement std lost positivity")
    sched = BASE.source_schedule()
    rs = Interval(*map(float, sched["R_S_applied_invariant"]))
    rs_lower = float(BASE.rs_variance_lower(rs, sched))
    if not _finite_positive(rs_lower):
        raise RuntimeError("source-invariant S pseudo-measurement covariance lower lost positivity")
    return {
        "accelerometer": down(acc_std * acc_std),
        "magnetometer": down(mag_std * mag_std),
        "S_zero_pseudo": rs_lower,
    }


def _reset_series_coefficients() -> dict:
    """Exact correction-quaternion coefficients on the shipping series branch.

    ``ou3_p4_exact_reset_transport`` encloses ``|a| <= coeff_hi |d|`` and
    ``|a-d| <= coeff_err(delta) |d|`` for ``|d| <= delta``.  On the shipping
    branch ``coeff_err(delta)/delta^2 = (1/8 + delta^2/1920)/(1 - delta^2/8)``
    is increasing, so evaluating it once at ``delta_ref`` gives the valid cubic
    coefficient ``|a-d| <= cubic |d|^3`` for every ``|d| <= delta_ref``.
    """
    ref = float(CORRECTION_REGION_RADIUS)
    corr = RESET.correction_cayley_norm_bounds(ref)
    if corr["source_branch_family"] != "source_polynomial_series":
        raise RuntimeError("correction region left the shipping polynomial quaternion branch")
    coeff_hi = float(corr["series_cayley_coefficient_upper"])
    adiff = float(corr["injected_cayley_minus_delta_norm_upper"])
    cubic = div_up(adiff, down(down(ref * ref) * ref))
    if not (_finite_positive(coeff_hi) and _finite_positive(cubic)):
        raise RuntimeError("correction quaternion series coefficients lost positivity")
    return {
        "correction_region_radius": ref,
        "injected_cayley_linear_coefficient_upper": coeff_hi,
        "injected_cayley_cubic_defect_coefficient_upper": cubic,
        "monotone_series_ratio_derivation": "(1/8+delta^2/1920)/(1-delta^2/8) is increasing on [0,1e-2]",
    }


def _reset_defect_coefficient(gain: float, q: float, series: dict) -> dict:
    """Quadratic coefficient of the exact quaternion/reset defect.

    With ``|d| <= gain ||z||``, ``|c| <= min(||z||, q)`` and
    ``gain ||z|| <= delta_ref`` the exact bound of the reset transport
    primitive reduces to a homogeneous quadratic form:

        rho <= [ cubic |d|^3 (1+q/2)
                 + coeff_hi |d| |c| ( |c| + |d| + |d||c|/2 ) / 4 ] / D
            <= C_reset ||z||^2,
        C_reset = [ cubic gain^2 delta_ref (1+q/2)
                    + coeff_hi delta_ref (1 + gain(1+q/2)) / 4 ] / D.
    """
    ref = float(series["correction_region_radius"])
    coeff_hi = float(series["injected_cayley_linear_coefficient_upper"])
    cubic = float(series["injected_cayley_cubic_defect_coefficient_upper"])
    denom = down(1.0 - mul_up(0.25, mul_up(mul_up(coeff_hi, ref), q)))
    if not denom > 0.0:
        raise RuntimeError("correction/state Cayley composition can reach the antipodal denominator")
    tail = add_up(1.0, mul_up(0.5, q))
    cubic_term = mul_up(mul_up(mul_up(cubic, mul_up(gain, gain)), ref), tail)
    bilinear = mul_up(
        mul_up(0.25, mul_up(coeff_hi, ref)),
        add_up(1.0, mul_up(gain, tail)),
    )
    return {
        "correction_gain_upper": gain,
        "cayley_composition_denominator_lower": denom,
        "cubic_quaternion_term_upper": cubic_term,
        "bilinear_composition_term_upper": bilinear,
        "reset_defect_quadratic_coefficient_upper": div_up(add_up(cubic_term, bilinear), denom),
    }


def _operation_defect_budget(mode: str, cov: dict, rvar: dict, cayley: dict,
                             remainder: dict, domain: dict, dt: float,
                             series: dict) -> dict:
    """Per-operation quadratic defect coefficients in the stored SI coordinate."""
    q = float(cayley["cayley_radius_upper"])
    live = domain["normal_live"]
    handoff = domain["startup"]["physical_handoff_coordinate_bounds"]
    f_max = float(live["specific_force_norm_upper_mps2"])
    m_max = float(live["magnetic_vector_norm_upper_uT"])
    bg_max = float(handoff["gyro_bias_error_norm_upper_rad_s"])
    ref = float(series["correction_region_radius"])
    coeff_hi = float(series["injected_cayley_linear_coefficient_upper"])
    lam_max = float(cov["covariance_lambda_max_upper"])

    if not (0.0 < q < 1.0):
        raise RuntimeError("retained sector Cayley radius left the q<1 operation region")
    if float(remainder["acc_eta_aw_quadratic_coefficient_upper"]) != 0.0:
        raise RuntimeError("active remainder still charges a nonlinear a_w eta term")
    if float(remainder["accelerometer_bias_nonlinear_remainder_coefficient"]) != 0.0:
        raise RuntimeError("active remainder still charges a nonlinear accelerometer-bias term")

    gains = {
        "S_zero_pseudo": sqrt_up(div_up(lam_max, rvar["S_zero_pseudo"])),
        "accelerometer": sqrt_up(div_up(lam_max, rvar["accelerometer"])),
        "magnetometer": sqrt_up(div_up(lam_max, rvar["magnetometer"])),
    }
    # Residual operator radii: |y| <= correction_residual_gain * ||z||.
    # S=0 selects the S state exactly; the accelerometer residual adds the two
    # linear a_w and b_a directions; the magnetometer residual is pure rotation.
    residual_gain = {
        "S_zero_pseudo": 1.0,
        "accelerometer": add_up(f_max, sqrt_up(2.0)),
        "magnetometer": m_max,
    }
    correction_gain = {
        name: mul_up(gains[name], residual_gain[name]) for name in gains
    }

    # Prediction: |d| <= dt |delta b_g|, so |c_d| <= coeff_hi dt ||z||.
    g_pred = mul_up(coeff_hi, dt)
    cd_max = mul_up(g_pred, bg_max)
    if not up(dt * bg_max) <= ref:
        raise RuntimeError("prediction rotation increment leaves the shipping series branch")
    denom_pred = down(1.0 - mul_up(0.25, mul_up(cd_max, q)))
    if not denom_pred > 0.0:
        raise RuntimeError("prediction Cayley composition can reach the antipodal denominator")
    prediction_cayley = div_up(
        mul_up(g_pred, add_up(0.5, mul_up(0.25, add_up(q, cd_max)))), denom_pred
    )
    prediction_rotation = mul_up(dt, mul_up(0.75, f_max))
    prediction_bias_coupling = float(dt)
    prediction = add_up(prediction_cayley, prediction_rotation, prediction_bias_coupling)

    resets = {
        name: _reset_defect_coefficient(correction_gain[name], q, series)
        for name in correction_gain
    }
    eta = {
        "S_zero_pseudo": 0.0,
        "accelerometer": mul_up(gains["accelerometer"], mul_up(0.5, f_max)),
        "magnetometer": mul_up(gains["magnetometer"], mul_up(0.5, m_max)),
    }
    coefficients = {"prediction": prediction}
    for name in correction_gain:
        coefficients[name] = add_up(
            eta[name], float(resets[name]["reset_defect_quadratic_coefficient_upper"])
        )
    if sorted(coefficients) != sorted(OPERATION_CLASSES):
        raise RuntimeError("operation class family changed")
    worst = max(coefficients, key=lambda k: coefficients[k])
    return {
        "mode": mode,
        "operation_classes": list(OPERATION_CLASSES),
        "full_gain_bound_formula": "||K||<=sqrt(lambda_max(Sigma)/lambda_min(R))",
        "measurement_covariance_lower": dict(rvar),
        "measurement_gain_upper": gains,
        "residual_operator_radius_upper": residual_gain,
        "accepted_correction_gain_upper": correction_gain,
        "pure_rotation_eta_quadratic_coefficient_upper": eta,
        "prediction_terms": {
            "exact_cayley_left_composition": prediction_cayley,
            "exact_rotation_residual_on_specific_force": prediction_rotation,
            "rotation_times_accel_bias_bilinear": prediction_bias_coupling,
        },
        "reset_terms": resets,
        "prediction_increment_inside_shipping_series_branch": True,
        "operation_quadratic_defect_upper": coefficients,
        "limiting_operation": worst,
        "quadratic_defect_constant_upper": coefficients[worst],
    }


def _mode_word(mode: str, delta: float, cov: dict, budget: dict,
               n_op: int, q: float, active_bias_interior: float) -> dict:
    """Retained source-path word estimate for one fixed-dimension mode."""
    lam_max = float(cov["covariance_lambda_max_upper"])
    lam_min = float(cov["covariance_lambda_min_lower"])
    # M_i = s_m Sigma_i^-1 with s_m = lambda_max(Sigma) gives m_- = 1 exactly.
    s_m = lam_max
    m_minus = 1.0
    m_plus = div_up(lam_max, lam_min)
    c_defect = float(budget["quadratic_defect_constant_upper"])
    b_m = div_up(
        mul_up(mul_up(4.0, float(n_op)), mul_up(sqrt_up(m_plus), c_defect)), m_minus
    )
    sqrt_w_star = down(down(delta / 8.0) / b_m)
    w_star = down(sqrt_w_star * sqrt_w_star)
    theta_star = down(math.sqrt(down(w_star / m_minus)))
    if not (_finite_positive(b_m) and _finite_positive(w_star) and _finite_positive(theta_star)):
        raise RuntimeError(f"{mode}: certified inner level underflowed to zero")

    checks = {
        "inner_level_inside_cayley_sector": theta_star <= q,
        "inner_level_inside_unit_cayley_ball": theta_star < 1.0,
        "accepted_corrections_inside_shipping_series_branch": all(
            up(float(g) * theta_star) <= float(CORRECTION_REGION_RADIUS)
            for g in budget["accepted_correction_gain_upper"].values()
        ),
        "prediction_increment_inside_shipping_series_branch": bool(
            budget["prediction_increment_inside_shipping_series_branch"]
        ),
    }
    if mode == "A":
        checks["inner_level_inside_accel_bias_projection_interior"] = (
            theta_star <= active_bias_interior
        )
    gap = down(delta / 2.0)
    # Both 1.0 and the outward-rounded float difference are valid uppers of the
    # real number 1-gap; take the smaller one and never round the strict gap
    # itself back into the reported companion.
    rho_upper = min(1.0, up(1.0 - gap))
    mu_w = down(gap * m_minus)
    if not _finite_positive(gap):
        raise RuntimeError(f"{mode}: strict word dissipation gap lost positivity")
    return {
        "mode": mode,
        "dimension": DIMENSION[mode],
        "canonical_P3_word_margin_consumed": float(delta),
        "metric_normalization_scalar_s_m": s_m,
        "metric_eigenvalue_lower_m_minus": m_minus,
        "metric_eigenvalue_upper_m_plus": m_plus,
        "covariance_lambda_max_upper": lam_max,
        "covariance_lambda_min_lower": lam_min,
        "word_state_operations": int(n_op),
        "quadratic_defect_constant_upper": c_defect,
        "limiting_operation": budget["limiting_operation"],
        "word_defect_gain_B_m_upper": b_m,
        "inner_level_sqrt_W_star_lower": sqrt_w_star,
        "inner_level_W_star_lower": w_star,
        "inner_attitude_radius_theta_star_lower_rad": theta_star,
        "direct_word_margin_mu_W_lower": mu_w,
        "one_minus_rho_lower": gap,
        "rho_upper": rho_upper,
        "strict_gap_below_binary64_epsilon": rho_upper >= 1.0,
        "funnel_consistency": checks,
        "funnel_consistent": all(checks.values()),
    }


def _flag(reasons: list[str], ok: bool, message: str) -> bool:
    if not ok:
        reasons.append(message)
    return bool(ok)


def build(p3: dict, metric: dict, timing: dict, clock: dict, audit: dict,
          *, domain_path: Path = DEFAULT_DOMAIN, cayley: dict | None = None,
          remainder: dict | None = None) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("complete-word P4 producer must not be trajectory fitted")
    cayley = CAYLEY.build(path) if cayley is None else cayley
    remainder = REMAINDER.build(path) if remainder is None else remainder

    failures: list[str] = []
    for label, errs in (
        ("metric", METRIC.validate(metric)),
        ("timing", TIMING.validate(timing)),
        ("clock", CLOCK.validate(clock)),
        ("audit", AUDIT.validate(audit)),
        ("Cayley", CAYLEY.validate(cayley)),
        ("remainder", REMAINDER.validate(remainder)),
    ):
        failures.extend(f"{label}: {x}" for x in errs)
    if p3.get("validation_pass") is not True or p3.get("validation_failures"):
        failures.append("P3: canonical artifact is not validated")

    reasons: list[str] = []
    p3_consumed = _flag(
        reasons,
        p3.get("qualification") == P3_QUALIFICATION
        and p3.get("P3_CANONICAL_PASS") is True
        and p3.get("P4_MAY_CONSUME_P3") is True
        and p3.get("only_this_module_may_promote_P3_for_P4") is True
        and float(p3.get("useful_gate", 0.0)) == float(BASE.MIN_USEFUL_DELTA),
        "canonical P3 artifact is not a passing unique-authority verdict at the unchanged gate",
    )
    same_history = _flag(
        reasons,
        metric.get("same_history_P3_frontier_consumed") is True
        and metric.get("independent_cartesian_tau_sigma_R_S_extrema_used") is False
        and int(metric.get("finite_source_phase_classes", 0)) == 800 * 26
        and metric.get("canonical_P3_candidate_numeric_pass_observed") is True
        and audit.get("same_history_P3_metric_consumed") is True,
        "metric attachment and signed-Joseph audit are not the same canonical P3 source history",
    )
    joseph_order = _flag(
        reasons,
        metric.get("translation_post_acc_S_at_metric_boundary") is True
        and metric.get("H_A_fresh_process_floor_at_same_endpoint_vector_packet") is True
        and metric.get("magnetometer_translation_jacobian_zero_on_declared_branch") is True,
        "implemented prediction/measurement/Joseph order is not covered by the consumed metric boundary",
    )
    packet_language = _flag(
        reasons,
        timing.get("nonlinear_timing_obligations_reduce_to_vector_measurements") is True
        and timing.get("ready_for_source_complete_nonlinear_remainder_composition") is True
        and timing.get("fixed_minimum_gap_S_schedule_is_source_complete") is False,
        "source-complete vector packet language is not established by the word timing contract",
    )
    s_linear = _flag(
        reasons,
        timing.get("S_residual_exactly_linear_selector") is True
        and timing.get("S_nonlinear_eta_identically_zero") is True
        and timing.get("S_timing_consumed_by_linear_P3_translation_UCO") is True,
        "S=0 linear timing is not discharged by canonical P3",
    )
    cayley_consumed = _flag(
        reasons,
        cayley.get("pass") is True and float(cayley.get("outer_angle_rad", 0.0)) >= 0.80,
        "exact Cayley sector geometry is not available at the required outer angle",
    )
    remainder_consumed = _flag(
        reasons,
        remainder.get("pass") is True
        and remainder.get("penalties_are_homogeneous_quadratic_not_affine_beta") is True
        and float(remainder.get("outer_angle_rad", 0.0)) >= 0.80,
        "homogeneous vector remainder sector is not available at the required outer angle",
    )
    clock_consumed = _flag(
        reasons,
        clock.get("P2_SAMPLE_CLOCK_REFINEMENT_CERTIFICATE") == "PASS"
        and clock.get("source_graph_all_to_all") is False
        and int(clock.get("partition", {}).get("states", 0)) == 800,
        "finite-speed committed sample-clock source graph is not available",
    )
    frozen_branch = _flag(
        reasons,
        clock.get("frozen_clock_self_loop_included") is True
        and clock.get("clock", {}).get("floating_clock_stagnation_verified") is True
        and metric.get("frozen_clock", {}).get("absorbing_hold_arbitrary_duration_covered") is True,
        "absorbing frozen-clock hold branch is not covered by both the source graph and the metric",
    )
    branch_scope = _flag(
        reasons,
        metric.get("zero_lever_arm_branch") is True
        and metric.get("dormant_transparent_vibration_guard_branch") is True,
        "declared zero-lever-arm / dormant vibration-guard branch is not retained",
    )
    outer_angle = float(cayley.get("outer_angle_rad", 0.0))
    angle_bound = _flag(
        reasons,
        float(remainder.get("outer_angle_rad", -1.0)) == outer_angle,
        "Cayley and remainder sector prerequisites disagree on the outer angle",
    )

    dt = float(timing.get("configured_dt_s", math.nan))
    samples = int(timing.get("word_samples_upper", 0) or 0)
    horizon = timing.get("word_horizon_s")
    edges = clock.get("transition_edges")
    if not (_finite_positive(dt) and samples > 0 and _finite_positive(horizon) and isinstance(edges, int)):
        failures.append("word clock, horizon or sample-clock edge family is invalid")
    n_op = OPERATIONS_PER_SAMPLE * samples

    live = domain["normal_live"]
    interior = down(
        float(live["active_accelerometer_bias_projection_limit_mps2"])
        - float(live["active_accelerometer_bias_state_norm_upper_mps2"])
    )
    margins = p3.get("mode_margins", {})
    worst_margin = p3.get("worst_H_A_margin")

    modes: dict[str, dict] = {}
    series: dict | None = None
    if not failures:
        try:
            series = _reset_series_coefficients()
            vector = VECTOR.build()
            vf = VECTOR.validate(vector)
            if vf:
                raise RuntimeError(f"vector UCO certificate invalid: {vf}")
            rvar = _measurement_variance_lowers(vector)
            for mode in MODES:
                delta = margins.get(mode) if isinstance(margins, dict) else None
                if not _finite_positive(delta) or float(delta) < float(BASE.MIN_USEFUL_DELTA):
                    failures.append(f"{mode}: canonical P3 mode margin is missing or below the unchanged gate")
                    continue
                cov = _covariance_operator_bounds(metric, mode)
                budget = _operation_defect_budget(
                    mode, cov, rvar, cayley, remainder, domain, dt, series
                )
                modes[mode] = _mode_word(
                    mode, float(delta), cov, budget, n_op,
                    float(cayley["cayley_radius_upper"]), interior,
                )
                modes[mode]["covariance_operator_bounds"] = cov
                modes[mode]["operation_defect_budget"] = budget
        except Exception as exc:  # fail closed, never silently promote
            failures.append(f"complete-word evaluation: {type(exc).__name__}: {exc}")

    word_evaluated = bool(not failures and set(modes) == set(MODES))
    dissipation = bool(
        not failures
        and not reasons
        and set(modes) == set(MODES)
        and all(modes[m]["funnel_consistent"] for m in MODES)
        and all(_finite_positive(modes[m]["one_minus_rho_lower"]) for m in MODES)
    )
    if worst_margin is None or not _finite_positive(worst_margin):
        failures.append("canonical P3 worst H/A margin is missing")
        dissipation = False

    def _mode_value(mode: str, key: str):
        return modes[mode][key] if mode in modes else None

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "canonical_P3_artifact_consumed": p3_consumed,
        "canonical_P3_qualification": p3.get("qualification"),
        "canonical_P3_worst_H_A_margin_consumed": worst_margin,
        "canonical_P3_mode_margins_consumed": {
            m: (margins.get(m) if isinstance(margins, dict) else None) for m in MODES
        },
        "same_source_history_for_metric_and_nonlinear_word": same_history,
        "implemented_prediction_measurement_Joseph_order_covered": joseph_order,
        "source_complete_vector_packet_language_covered": packet_language,
        "S_linear_timing_discharged_by_canonical_P3": s_linear,
        "Cayley_exact_geometry_consumed": cayley_consumed and angle_bound,
        "homogeneous_vector_remainder_consumed": remainder_consumed and angle_bound,
        "finite_speed_sample_clock_graph_consumed": clock_consumed,
        "frozen_clock_absorbing_hold_branch_covered": frozen_branch,
        "zero_lever_arm_branch": branch_scope,
        "dormant_transparent_vibration_guard_branch": branch_scope,
        "full_H18_state_word_covered": "H" in modes and modes["H"]["dimension"] == 18,
        "full_A21_state_word_covered": "A" in modes and modes["A"]["dimension"] == 21,
        "signed_nonlinear_remainder_charged": word_evaluated,
        "complete_word_generalized_Jacobian_or_equivalent_bound": word_evaluated,
        "instantaneous_one_sample_decrease_required": False,
        "prefix_nonexpansive_in_matched_source_metrics": True,
        "H_dimension": 18,
        "A_dimension": 21,
        "outer_angle_rad": outer_angle,
        "source_word_horizon_s": horizon,
        "source_word_samples_upper": samples,
        "word_state_operations": n_op,
        "operations_per_sample": OPERATIONS_PER_SAMPLE,
        "sample_clock_transition_edges": edges,
        "source_partition_states": clock.get("partition", {}).get("states"),
        "correction_region_radius": None if series is None else series["correction_region_radius"],
        "correction_quaternion_series_coefficients": series,
        "active_accelerometer_bias_projection_interior_margin_mps2": interior,
        "modes": modes,
        "rho_H_upper": _mode_value("H", "rho_upper"),
        "rho_A_upper": _mode_value("A", "rho_upper"),
        "one_minus_rho_H_lower": _mode_value("H", "one_minus_rho_lower"),
        "one_minus_rho_A_lower": _mode_value("A", "one_minus_rho_lower"),
        "strict_dissipation_margin_H_lower": _mode_value("H", "one_minus_rho_lower"),
        "strict_dissipation_margin_A_lower": _mode_value("A", "one_minus_rho_lower"),
        "direct_word_margin_mu_W_lower": {
            m: _mode_value(m, "direct_word_margin_mu_W_lower") for m in MODES
        },
        "inner_level_W_star_lower": {
            m: _mode_value(m, "inner_level_W_star_lower") for m in MODES
        },
        "inner_attitude_radius_theta_star_lower_rad": {
            m: _mode_value(m, "inner_attitude_radius_theta_star_lower_rad") for m in MODES
        },
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED": dissipation,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "P5_FINITE_CAPTURE_ESTABLISHED_HERE": False,
        "next_obligation": (
            "submit this artifact to tools/ou3_p4_canonical_gate.py; P5 must then prove finite "
            "startup/outer-sector capture into the certified inner funnel W_*"
        ),
        "unmet_theorem_obligations": reasons,
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != QUALIFICATION:
        f.append("wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit",
        "prefix_nonexpansive_in_matched_source_metrics",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "instantaneous_one_sample_decrease_required",
        "P4_USABLE_CERTIFICATE_PROMOTED", "P5_FINITE_CAPTURE_ESTABLISHED_HERE",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("H_dimension") != 18 or d.get("A_dimension") != 21:
        f.append("fixed-dimension H/A word sizes changed")
    if d.get("operations_per_sample") != OPERATIONS_PER_SAMPLE:
        f.append("implemented operation family per sample changed")
    established = bool(d.get("P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"))
    unmet = d.get("unmet_theorem_obligations")
    if not isinstance(unmet, list):
        f.append("unmet obligation list is malformed")
        unmet = []
    if established and unmet:
        f.append("complete-word dissipation was declared with unmet theorem obligations")
    modes = d.get("modes", {})
    if established and set(modes) != set(MODES):
        f.append("complete-word dissipation was declared without both fixed-dimension modes")
    for mode in MODES:
        row = modes.get(mode)
        if row is None:
            if established:
                f.append(f"{mode}: complete-word row is missing")
            continue
        if row.get("dimension") != DIMENSION[mode]:
            f.append(f"{mode}: word dimension is not {DIMENSION[mode]}")
        for key in (
            "metric_eigenvalue_upper_m_plus", "quadratic_defect_constant_upper",
            "word_defect_gain_B_m_upper", "inner_level_W_star_lower",
            "inner_attitude_radius_theta_star_lower_rad",
            "direct_word_margin_mu_W_lower", "one_minus_rho_lower",
        ):
            if not _finite_positive(row.get(key)):
                f.append(f"{mode}: {key} is not finite positive")
        if float(row.get("metric_eigenvalue_lower_m_minus", 0.0)) != 1.0:
            f.append(f"{mode}: retained s_m normalization no longer gives m_-=1")
        rho = row.get("rho_upper")
        gap = row.get("one_minus_rho_lower")
        if not (_finite_number_ok(rho) and _finite_positive(gap)):
            f.append(f"{mode}: contraction pair is malformed")
        elif float(rho) > up(1.0 - float(gap)):
            f.append(f"{mode}: reported rho upper is inconsistent with the exact strict gap")
        if established and row.get("funnel_consistent") is not True:
            f.append(f"{mode}: certified funnel left its declared operation region")
        margin = d.get(f"strict_dissipation_margin_{mode}_lower")
        if not _finite_positive(margin) or (
            _finite_positive(gap) and float(margin) > float(gap)
        ):
            f.append(f"{mode}: strict dissipation margin is not inside the certified gap")
    return list(dict.fromkeys(f))


def _finite_number_ok(x) -> bool:
    return (
        not isinstance(x, bool)
        and isinstance(x, (int, float))
        and math.isfinite(float(x))
    )


def _load(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--p3", type=Path, required=True)
    ap.add_argument("--metric", type=Path, required=True)
    ap.add_argument("--timing", type=Path, required=True)
    ap.add_argument("--clock", type=Path, required=True)
    ap.add_argument("--audit", type=Path, required=True)
    ap.add_argument("--cayley", type=Path)
    ap.add_argument("--remainder", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(
        _load(a.p3), _load(a.metric), _load(a.timing), _load(a.clock), _load(a.audit),
        domain_path=a.domain,
        cayley=None if a.cayley is None else _load(a.cayley),
        remainder=None if a.remainder is None else _load(a.remainder),
    )
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "outer_angle_rad": d["outer_angle_rad"],
        "word_state_operations": d["word_state_operations"],
        "quadratic_defect_constant_upper": {
            m: (d["modes"].get(m) or {}).get("quadratic_defect_constant_upper") for m in MODES
        },
        "metric_eigenvalue_upper_m_plus": {
            m: (d["modes"].get(m) or {}).get("metric_eigenvalue_upper_m_plus") for m in MODES
        },
        "word_defect_gain_B_m_upper": {
            m: (d["modes"].get(m) or {}).get("word_defect_gain_B_m_upper") for m in MODES
        },
        "inner_level_W_star_lower": d["inner_level_W_star_lower"],
        "inner_attitude_radius_theta_star_lower_rad": d["inner_attitude_radius_theta_star_lower_rad"],
        "rho_upper": {"H": d["rho_H_upper"], "A": d["rho_A_upper"]},
        "one_minus_rho_lower": {"H": d["one_minus_rho_H_lower"], "A": d["one_minus_rho_A_lower"]},
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED": d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"],
        "unmet_theorem_obligations": d["unmet_theorem_obligations"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
