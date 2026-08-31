#!/usr/bin/env python3
"""Metric-consistent transported-defect gain for the OU-III P4 nonlinear word.

P4 bounds the exact nonlinear source word by its P3 homogeneous tangent map plus
an exact defect, and needs one number: an upper bound ``B`` on

    ||r_word||_{M_end} / W_0,

where ``M_g = s_m Sigma_KF(g)^-1`` is the mode metric and ``W_0`` the entry
level.  The original route obtains ``B`` by leaving the metric entirely:

    ||z||_2 <= sqrt(W/m_-),  ||r||_M <= sqrt(m_+) ||r||_2,

with ``m_- = s/Sigma_max`` and ``m_+ = s/Sigma_min``.  That pays the full
``sqrt(cond(Sigma))`` of an 18/21-state covariance whose entries range over
thirty-four decades, and it bounds the gain by the isotropic
``||K|| <= sqrt(Sigma_max/R_min)``.  Both steps are avoidable: the defects of the
shipping word are *structured*, and the metric that measures them is the source
covariance itself.

This module supplies the structured replacement.  Three exact facts are used,
each stated against the same source-derived P3/process bounds P4 already
consumes.

1.  **Gain transport (the injected correction defect never leaves the metric).**
    For the shipping gain ``K = P H^T S^-1``, ``S = H P H^T + R``, write
    ``A = R^-1/2 H P^1/2``.  Then

        P^-1/2 K = P^1/2 H^T S^-1 = A^T (A A^T + I)^-1 R^-1/2,

    and ``||A^T (A A^T + I)^-1||_2 = max_i sigma_i/(1+sigma_i^2) <= 1/2``.  Hence
    for every residual defect ``q``

        ||K q||_M <= sqrt(s) ||q||_2 / (2 sqrt(lambda_min(R))).

    The node metric uses ``Sigma_KF(g) >= P``, so ``Sigma^-1 <= P^-1`` and the
    same bound holds for the metric actually used.

2.  **Chart transport (the defect inputs are marginal, not global).**
    Marginalising the exact quadratic form,

        min_xi [c;xi]^T Sigma^-1 [c;xi] = c^T (Sigma_cc)^-1 c,

    so for any coordinate block ``c`` of the state,

        ||c||^2 <= lambda_max(Sigma_cc) W / s.

    The P3 covariance upper is a Loewner diagonal dominator, so
    ``lambda_max(Sigma_cc)`` is bounded by the largest diagonal upper on that
    block.  The nonlinear defects of the shipping word are quadratic in the
    attitude, gyro-bias and ``a_w`` coordinates only; the translation block --
    which carries the whole thirty-four-decade spread -- never enters them.

3.  **Attitude-injection cost (an attitude-supported defect is charged on the
    conditional attitude covariance).**  The quaternion injection remainder is
    supported on the attitude coordinates, so its metric cost is
    ``sqrt(s lambda_max((Sigma^-1)_theta,theta))``.  Every source node is either
    a post-prediction node, where ``Sigma >= Q`` and therefore
    ``(Sigma^-1)_tt <= (Q^-1)_tt <= I/(rho_att q_theta)``, or a node reached from
    one by at most ``n_corr`` corrections inside the same sample, each adding at
    most ``||H_theta||^2/lambda_min(R)`` of attitude information in the exact
    information form ``P+^-1 = P^-1 + H^T R^-1 H``.

Everything below is outward rounded.  The module returns an upper bound on the
same ``B`` the parent route bounds, so the consumer takes the minimum of the two
and can never be widened by this refinement.
"""
from __future__ import annotations

import math

SCHEMA = 1

# Exact Cayley composition c (+) d = (c + d + 0.5 d x c)/(1 - d.c/4).
CAYLEY_CROSS_COEFF = 0.5
CAYLEY_DENOMINATOR_COEFF = 0.25
# Deployed normalized polynomial quaternion branch, same constant as the parent.
SOURCE_SERIES_CAYLEY_CUBIC_COEFF = 0.085
# Prefix bootstrap factor shared with the parent route (W_s <= 4 W_0).
PREFIX_BOOTSTRAP_W_FACTOR = 4.0
# Corrections applied inside one IMU sample after its prediction: S, accel, mag.
MAX_CORRECTIONS_PER_SAMPLE = 3
# Largest accepted attitude injection kept inside the deployed series branch.
MAX_DESIGN_INJECTION_NORM = 5.0e-3


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


def sqrt_up(x: float) -> float:
    if x < 0.0:
        raise ValueError("sqrt of negative")
    r = math.sqrt(float(x))
    # One outward step is enough: math.sqrt is correctly rounded.
    return up(r)


def _block_upper(diag_upper: list[float], indices: range | list[int]) -> float:
    vals = [float(diag_upper[i]) for i in indices]
    if not vals or not all(math.isfinite(v) and v > 0.0 for v in vals):
        raise RuntimeError("covariance diagonal upper block is not finite positive")
    return up(max(vals))


def _composition_defect_upper(c_max: float, d_max: float) -> float:
    """Exact Cayley/quaternion composition remainder above the additive map.

    ``c_max`` bounds the current Cayley norm and ``d_max`` the injected
    correction norm.  The homogeneous reference is ``c + d``; everything else is
    defect.
    """
    if not d_max < MAX_DESIGN_INJECTION_NORM:
        raise RuntimeError(f"injection radius leaves the deployed series branch: {d_max}")
    dot_max = mul_up(d_max, c_max)
    denom_lower = down(1.0 - up(CAYLEY_DENOMINATOR_COEFF * dot_max))
    if not denom_lower > 0.99:
        raise RuntimeError("Cayley product denominator is not safely positive")
    cross = div_up(mul_up(CAYLEY_CROSS_COEFF, dot_max), denom_lower)
    numerator_norm = add_up(add_up(c_max, d_max), mul_up(CAYLEY_CROSS_COEFF, dot_max))
    denom_effect = div_up(
        mul_up(numerator_norm, mul_up(CAYLEY_DENOMINATOR_COEFF, dot_max)), denom_lower
    )
    cubic = mul_up(SOURCE_SERIES_CAYLEY_CUBIC_COEFF, mul_up(d_max, mul_up(d_max, d_max)))
    return add_up(cross, add_up(denom_effect, cubic))


def _gain_at_radius(inputs: dict, w: float) -> dict:
    """Transported word-defect gain B_metric valid on the metric ball ||z||_M<=w."""
    s = float(inputs["metric_scale"])
    sqrt_s = sqrt_up(s)
    rmin = float(inputs["correction_R_lambda_min_lower"])
    sqrt_rmin = math.sqrt(rmin)
    if not sqrt_rmin > 0.0:
        raise RuntimeError("measurement noise floor is not positive")

    a_theta = sqrt_up(div_up(inputs["Sigma_attitude_upper"], s))
    a_bias = sqrt_up(div_up(inputs["Sigma_gyro_bias_upper"], s))
    a_input = sqrt_up(div_up(inputs["Sigma_defect_input_upper"], s))

    # Fact 1: metric cost of a defect injected through the shipping gain.
    gain_transport = div_up(sqrt_s, up(2.0 * sqrt_rmin))
    # Fact 3: metric cost of a defect supported on the attitude coordinates.
    attitude_information_upper = add_up(
        div_up(1.0, down(inputs["rho_attitude_scaled_lower"] * inputs["Q_theta_diagonal_lower"])),
        div_up(
            mul_up(float(MAX_CORRECTIONS_PER_SAMPLE), mul_up(inputs["H_attitude_norm_upper"],
                                                             inputs["H_attitude_norm_upper"])),
            rmin,
        ),
    )
    attitude_transport = sqrt_up(mul_up(s, attitude_information_upper))

    # Fact 2: the defect inputs live on the attitude/bias/a_w marginal block.
    c_max = mul_up(a_theta, w)
    input_max = mul_up(a_input, w)

    # Accepted vector correction: residual remainder pushed through the gain.
    residual_defect = mul_up(inputs["vector_residual_quadratic_constant_upper"],
                             mul_up(input_max, input_max))
    corr_state_defect = mul_up(gain_transport, residual_defect)

    # The same correction's exact attitude injection magnitude.  The linear part
    # is P^1/2 A^T (A A^T + I)^-1 A P^-1/2 z, whose attitude rows are bounded by
    # sqrt(lambda_max(Sigma_tt)) since ||A^T (A A^T+I)^-1 A|| <= 1.
    inject_linear = mul_up(a_theta, w)
    inject_nonlinear = mul_up(mul_up(a_theta, sqrt_s),
                              div_up(residual_defect, up(2.0 * sqrt_rmin)))
    inject_max = add_up(inject_linear, inject_nonlinear)
    corr_inject_defect = mul_up(attitude_transport,
                                _composition_defect_upper(c_max, inject_max))
    correction_defect = add_up(corr_state_defect, corr_inject_defect)

    # Prediction: the only nonlinearity is the attitude composition with the
    # propagated increment, whose norm is the gyro-bias error over one sample.
    pred_increment = mul_up(inputs["prediction_increment_gain_upper"], mul_up(a_bias, w))
    prediction_defect = mul_up(attitude_transport,
                               _composition_defect_upper(c_max, pred_increment))

    per_operation = max(correction_defect, prediction_defect)
    if not (math.isfinite(per_operation) and per_operation > 0.0):
        raise RuntimeError("metric-consistent operation defect is not finite positive")

    # ||r||_M <= per_operation, and per_operation is at least quadratic in w, so
    # dividing by w^2 gives a valid coefficient on the whole ball of radius w.
    kappa = div_up(per_operation, down(w * w))
    B = mul_up(PREFIX_BOOTSTRAP_W_FACTOR,
               mul_up(float(inputs["state_operation_count_upper"]), kappa))
    return {
        "design_metric_radius": w,
        "attitude_chart_scale": a_theta,
        "gyro_bias_chart_scale": a_bias,
        "defect_input_chart_scale": a_input,
        "gain_metric_transport_upper": gain_transport,
        "attitude_conditional_information_upper": attitude_information_upper,
        "attitude_metric_transport_upper": attitude_transport,
        "accepted_injection_norm_upper": inject_max,
        "correction_state_defect_upper": corr_state_defect,
        "correction_injection_defect_upper": corr_inject_defect,
        "prediction_defect_upper": prediction_defect,
        "per_operation_metric_defect_upper": per_operation,
        "metric_defect_coefficient_kappa_upper": kappa,
        "transported_word_defect_B_upper": B,
    }


def build(inputs: dict) -> dict:
    """Largest certified metric level for the metric-consistent transport route.

    ``inputs`` carries only source-derived P3/process/domain numbers.  The design
    radius is searched over a fixed decreasing ladder; a radius is admissible
    only when the level it certifies closes its own prefix bootstrap inside it.
    """
    delta = float(inputs["word_endpoint_delta_lower"])
    if not delta > 0.0:
        raise RuntimeError("P3 endpoint margin is not positive")

    best = None
    ladder = []
    exponent = -2
    while exponent >= -160:
        ladder.append(10.0 ** exponent)
        exponent -= 2
    for w in ladder:
        try:
            g = _gain_at_radius(inputs, w)
        except RuntimeError:
            continue
        B = float(g["transported_word_defect_B_upper"])
        sqrt_W = down(delta / up(8.0 * B))
        W = down(sqrt_W * sqrt_W)
        prefix = mul_up(sqrt_up(PREFIX_BOOTSTRAP_W_FACTOR), sqrt_W)
        if not (W > 0.0 and prefix <= w):
            continue
        if best is None or W > best["certified_level_W"]:
            best = dict(g)
            best["certified_level_sqrt_W"] = sqrt_W
            best["certified_level_W"] = W
            best["prefix_metric_norm_upper"] = prefix
    if best is None:
        raise RuntimeError("metric-consistent transport route found no admissible design radius")

    best["schema"] = SCHEMA
    best["route"] = "METRIC_CONSISTENT_STRUCTURED_DEFECT_TRANSPORT"
    best["prefix_attitude_cayley_norm_upper"] = mul_up(
        best["attitude_chart_scale"], best["prefix_metric_norm_upper"]
    )
    best["prefix_defect_input_norm_upper"] = mul_up(
        best["defect_input_chart_scale"], best["prefix_metric_norm_upper"]
    )
    best["exact_facts_used"] = [
        "P^1/2 H^T S^-1 = A^T (A A^T + I)^-1 R^-1/2 with spectral norm <= 1/(2 sqrt(lambda_min R))",
        "min_xi [c;xi]^T Sigma^-1 [c;xi] = c^T (Sigma_cc)^-1 c",
        "Sigma >= Q at every post-prediction node, plus at most three in-sample corrections",
    ]
    return best


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA:
        failures.append("metric defect transport schema mismatch")
    for key in ("certified_level_W", "certified_level_sqrt_W",
                "transported_word_defect_B_upper", "metric_defect_coefficient_kappa_upper"):
        v = d.get(key)
        if not (isinstance(v, (int, float)) and math.isfinite(float(v)) and float(v) > 0.0):
            failures.append(f"metric defect transport {key} is not finite positive")
    if not float(d.get("prefix_metric_norm_upper", math.inf)) <= float(
        d.get("design_metric_radius", 0.0)
    ):
        failures.append("metric defect transport prefix bootstrap did not close")
    if not float(d.get("prefix_attitude_cayley_norm_upper", math.inf)) < 1.0:
        failures.append("metric defect transport prefix leaves the Cayley chart")
    return failures
