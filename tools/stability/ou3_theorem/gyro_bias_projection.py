"""Source-bound gyro mean invariant and one-prediction transport separation.

This lemma enters the historical AG readout prerequisite of the existing
finite-error word inequality. It does not supply the complete signed temporal
margin, an AG covariance ceiling, a contraction factor, or float32 totality.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[3]
COMMON = ROOT / "src/kalman_ou_common/KalmanOUCoreMath.h"


def source_radius():
    source = COMMON.read_text()
    match = re.search(r"inline constexpr double gyro_bias_radius_rad_s = ([0-9.]+);", source)
    if not match:
        raise ValueError("missing fixed shipping gyro-bias radius")
    radius = F(match[1])
    if radius <= 0:
        raise ValueError("positive fixed gyro radius required")
    return radius


def certificate():
    c = json.loads((Path(__file__).with_name("constants.json")).read_text(), parse_float=F)
    sensor = c["sensor_model"]
    h_min, h_max = sensor["sample_period_min_s"], sensor["sample_period_max_s"]
    omega = c["marine_motion"]["Omega_max_rad_s"]
    physical_bias = c["imu_bias"]["B_g_s_rad_s"]
    residual = c["imu_bias"]["B_g_f_rad_s"]
    radius = source_radius()
    angle = h_max * (omega + physical_bias + residual + radius)
    # The shipping quaternion uses a normalized polynomial below .01 rad.
    polynomial_defect = 4*(angle**6/F(46080)+angle**7/F(645120))
    cap = F(7,1000)
    if not (0 < h_min <= h_max and angle < F(1,100) and angle+polynomial_defect < cap):
        raise ValueError("shipping inputs exceed the certified prediction-angle cap")
    # Rodrigues branch: transverse singular values = h sinc(phi/2).
    # sin(x)>=x-x^3/6 on [0,.0035], so this is an exact rational floor.
    relative_floor = 1-cap**2/24
    lower = h_min*relative_floor
    # The small-rate branch is a polynomial, not exact Rodrigues. Its
    # transverse symmetric part is h(1-phi^2/6)I and its axial value is h.
    small_rate = F(1,10_000_000)
    small_floor = h_min*(1-(small_rate*h_max)**2/6)
    if not small_floor > lower > 0:
        raise ArithmeticError("small-rate source branch lacks the claimed floor")
    return {
        "qualification": "OU3_IMPLEMENTED_GYRO_BIAS_INVARIANT_V2_SLOW_FAST",
        "physical_coordinate": "slow gyro bias; fast remains delivered sensor forcing",
        "amplitude_bound_is_not_temporal_qualification": True,
        "estimator_radius_rad_s": str(radius),
        "physical_slow_bias_rad_s": str(physical_bias),
        "fast_measurement_amplitude_rad_s": str(residual),
        "physical_angular_rate_rad_s": str(omega),
        "admitted_step_s": [str(h_min), str(h_max)],
        "corrected_rate_norm_upper_rad_s": str(omega+physical_bias+residual+radius),
        "prediction_argument_upper_rad": str(angle),
        "quaternion_polynomial_angle_defect_upper_rad": str(polynomial_defect),
        "actual_real_source_rotation_angle_upper_rad": str(cap),
        "margin_from_pi_lower_rad": str(F(314159,100000)-cap),
        "margin_from_two_pi_lower_rad": str(2*F(314159,100000)-cap),
        "gyro_transport_relative_singular_floor": str(relative_floor),
        "gyro_transport_singular_floor_s": str(lower),
        "small_rate_polynomial_singular_floor_s": str(small_floor),
        "physical_amplitude_numbers_changed": False,
        "imu_temporal_model_migrated": True,
        "mean_only_covariance_preserved": True,
        "complete_turn_nominal_bias_alias_excluded": True,
        "single_prediction_transport_source_uniform_real_arithmetic": True,
        "signed_temporal_Delta_gyr_closed": False,
        "force_field_collinearity_excluded": False,
        "uniform_historical_AG_action_closed": False,
        "all_positive_device_timesteps_covered": False,
        "float32_word_totality_certified": False,
        "remaining_scope": "The full signed temporal margin includes chronological reset/observation forcing and projection defects; a one-cell singular floor does not bound that word operator.",
    }
