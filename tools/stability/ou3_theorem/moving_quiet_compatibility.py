"""Exact physical compatibility audit; not an all-time shipping-service certificate.

The first MOVING two-epoch separation target is tested BEFORE a norm bound:
R_+ = Rx(theta) Ry(psi), R_- = Rx(-theta) Ry(psi),
p = .02 sin(pi*t/10) e_x, psi = .02 sin(pi*t/10), B = 75 e_x,
b_a,+/- = +/- g sin(theta) e_y, b_g = b_a,f = b_g,f = 0.
These two complete physical histories have identical delivered sensor records.

Their constant BA difference is not a homogeneous BA perturbation: prediction
must retain (1-phi_b)*b_a.  The shipping estimator is never replaced/reseeded.
See docs/ou3-moving-quiet-compatibility.md for the proofs and their exact scope.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import struct

REPO = Path(__file__).resolve().parents[3]
CONSTANTS = REPO / "tools/stability/ou3_theorem/constants.json"
HEADER = REPO / "src/kalman_ou_iii/Kalman3D_Wave_OU_III.h"
ZERO = (0, 0, 0)  # cos(psi), sin(psi), physical a_x


def binary32(value: float) -> F:
    """The exact real value represented by a shipping binary32 literal."""
    return F(struct.unpack("f", struct.pack("f", value))[0])


def _constant(value):
    value = F(value)
    return {ZERO: value} if value else {}


def _add(a, b):
    out = dict(a)
    for monomial, value in b.items():
        out[monomial] = out.get(monomial, F(0)) + value
    return {m: v for m, v in out.items() if v}


def _scale(a, scale):
    return {m: v * F(scale) for m, v in a.items() if v * F(scale)}


def _mul(a, b):
    out = {}
    for m, v in a.items():
        for n, w in b.items():
            key = tuple(x + y for x, y in zip(m, n))
            out[key] = out.get(key, F(0)) + v * w
    return {m: v for m, v in out.items() if v}


def _mv(a, x):
    return [sum_poly(_mul(v, w) for v, w in zip(row, x)) for row in a]


def sum_poly(values):
    out = {}
    for value in values:
        out = _add(out, value)
    return out


def sensor_polynomial_identity(g: F, sine: F, cosine: F):
    """Formal polynomial equality, not a finite sampling of trigonometric phases."""
    c, s, ax = ({(1, 0, 0): F(1)}, {(0, 1, 0): F(1)}, {(0, 0, 1): F(1)})
    u = [[c, {}, _scale(s, -1)], [{}, _constant(1), {}], [s, {}, c]]
    forces, fields = [], []
    for sign in (-1, 1):
        # Ry(-psi) Rx(-sign*theta) (a_x,0,-g) + sign*g*sin(theta)*e_y.
        force = _mv(u, [ax, _constant(-sign * g * sine), _constant(-g * cosine)])
        force[1] = _add(force[1], _constant(sign * g * sine))
        forces.append(force)
        fields.append(_mv(u, [_constant(75), {}, {}]))
    assert forces[0] == forces[1]
    assert fields[0] == fields[1]
    return {
        "formal_polynomial_difference_zero": True,
        "gyro_identity": "R_s^T dR_s/dt = [psi_dot*e_y]_x, independent of sign s",
        "accelerometer": ["a_x*cos(psi)+g*cos(theta)*sin(psi)", "0",
                          "a_x*sin(psi)-g*cos(theta)*cos(psi)"],
        "magnetometer": ["75*cos(psi)", "0", "75*sin(psi)"],
        "gravity_difference": "(R_+^T-R_-^T)g*e_z = 2*g*sin(theta)*e_y",
        "two_epoch_difference_exact": "0",
    }


def centered_slow_charge(cap: F, rate: F, horizon: F, length: F, max_step: F) -> F:
    """C/L + D L/4 + D h_max, for an interior centered constant-record window.

The last term compares held slow samples with the continuous slow history.
It is NOT legitimate to differentiate a held FAST signal. The whole centered
window must be contained in one quiet interval and have L <= H.
"""
    cap, rate, horizon, length, max_step = map(F, (cap, rate, horizon, length, max_step))
    if cap < 0 or rate < 0 or not 0 < length <= horizon or max_step <= 0:
        raise ValueError("nonnegative cap/rate and positive qualified window/mesh required")
    return cap / length + rate * length / 4 + rate * max_step


def continuous_fibre(g: F, cosine: F, sine: F, amplitude: F, nu_hi: F,
                     slow_cap: F, slow_rate: F) -> dict:
    """An exact same-record arc, not a family of perturbed estimator states.

    alpha is constant on each complete physical history. Only the endpoints
    have constant BA. The intermediate slow bias is a smooth body-frame
    history; it is never imposed as an estimator OU law.
    """
    # Formal cancellation with cos(alpha),sin(alpha) left indeterminate.
    # Monomials: cos(psi),sin(psi),a_x,cos(alpha),sin(alpha).
    z = (0,)*5
    def scalar(x):
        return {z: F(x)} if x else {}
    def variable(i):
        m = list(z); m[i] = 1
        return {tuple(m): F(1)}
    c, s, ax, ca, sa = map(variable, range(5))
    U = [[c, {}, _scale(s, -1)], [{}, scalar(1), {}], [s, {}, c]]
    force = _mv(U, [ax, _scale(sa, -g), _scale(ca, -g)])
    normal = [_scale(s, -1), {}, c]
    bias = [_scale(_mul(_add(ca, scalar(-cosine)), n), g) for n in normal]
    bias[1] = _add(bias[1], _scale(sa, g))
    delivered = [_add(x, y) for x, y in zip(force, bias)]
    expected = _mv(U, [ax, {}, scalar(-g*cosine)])
    if delivered != expected:
        raise ValueError('continuous physical sensor identity failed')
    norm2 = (g*sine)**2
    rate = g*(1-cosine)*amplitude*nu_hi
    # On the principal constant-roll family, the exact slow-amplitude cap is
    # cos(alpha)>=(1+c^2-(B/g)^2)/(2c). This is a necessary family bound,
    # not a characterization of every compatible SLOW+FAST history.
    cap_cos = (1+cosine*cosine-(slow_cap/g)**2)/(2*cosine)
    if not 0 < cap_cos < 1:
        raise ValueError('principal small-roll amplitude comparison unavailable')
    cap_tan2 = (1-cap_cos)/(1+cap_cos)
    angle_cap = F(23, 1000)
    if 4*cap_tan2 >= angle_cap**2:
        raise ValueError('principal-roll rational angle cap changed')
    # theta=2 atan(1/200)<1/100. At the central physical chart the
    # attitude increment is exactly -alpha Ry(-psi)e_x. BA adds curvature.
    arc_angle = F(1, 100)
    curvature2 = g*g*(arc_angle**6/36+arc_angle**4/4)
    return {
        'result_type': 'PROVED analytical theorem',
        'parameter': 'constant alpha in [-theta,theta], tan(theta/2)=1/200',
        'rotation': 'R_alpha=Rx(alpha) Ry(psi)',
        'slow_bias': 'g sin(alpha) e_y + g (cos(alpha)-cos(theta)) Ry(-psi) e_z',
        'fast_and_gyro_bias': 'identically zero',
        'formal_sensor_polynomial_identity': delivered == expected,
        'bias_norm_squared_identity': 'g^2 (1+c^2-2c cos(alpha)), c=cos(theta)',
        'slow_bias_norm_squared_upper': str(norm2),
        'slow_bias_rate_upper': str(rate),
        'slow_amplitude_verified': norm2 <= slow_cap**2,
        'slow_rate_verified': rate <= slow_rate,
        'gravity_span_lower_unchanged': True,
        'principal_family_slow_cap_cos_alpha_lower': str(cap_cos),
        'principal_family_tan_half_angle_squared_upper': str(cap_tan2),
        'principal_family_abs_angle_upper_rad': str(angle_cap),
        'arc_abs_angle_upper_rad': str(arc_angle),
        'arc_BA_tangent_remainder_norm_squared_upper': str(curvature2),
        'BA_tangent_remainder': 'g (sin(alpha)-alpha) e_y + g (cos(alpha)-1) Ry(-psi) e_z',
        'physical_angle_bound_is_normalized_gauge_bound': False,
        'same_record_total_estimator_variation': 'zero, by deterministic prefix induction',
        'nominal_root_perturbation_is_physical_fibre_variation': False,
        'C_Q_alpha_may_be_dropped_from_error_word': False,
        'uniform_precision_metric_gauge_bound_verified': False,
        'actually_applied_magnetic_service_all_time_verified': False,
        'structures_preserved': 'one complete physical history per alpha; identical delivered record and complete shipping state/mean/P/K/frontend/clocks',
        'relaxations_introduced': 'none in the physical arc identity; all-time shipping admission remains open',
    }


def certificate() -> dict:
    limits = json.loads(CONSTANTS.read_text())
    m, b, mag = (limits[k] for k in ("marine_motion", "imu_bias", "magnetic_service"))
    q = lambda value: F(str(value))
    g, sine, cosine, amplitude = F(196133, 20000), F(400, 40001), F(39999, 40001), F(1, 50)
    assert sine*sine + cosine*cosine == 1
    # Classical rational pi enclosure; no rounded floating comparison proves membership.
    pi_lo, pi_hi = F(333, 106), F(355, 113)
    nu_lo, nu_hi = pi_lo/10, pi_hi/10
    # Each 30-s window contains a complete 20-s period and its two extrema.
    half_gravity_chord_lo = cosine*(amplitude-amplitude**3/6)
    half_threshold_hi = q(m["attitude_excitation"]["theta_E_rad"])/2
    bounds = {
        "position_max_m": amplitude,
        "velocity_max_mps": amplitude*nu_hi,
        "acceleration_max_mps2": amplitude*nu_hi**2,
        "jerk_max_mps3": amplitude*nu_hi**3,
        "omega_max_rad_s": amplitude*nu_hi,
        "position_primitive_max_m_s": 2*amplitude/nu_lo,
        "position_span_m": 2*amplitude,
        "half_gravity_chord_lower": half_gravity_chord_lo,
        "slow_accel_bias_norm_mps2": g*sine,
        "slow_accel_bias_rate_mps3": F(0),
        "slow_gyro_bias_norm_rad_s": F(0),
        "slow_gyro_bias_rate_rad_s2": F(0),
        "fast_accel_norm_mps2": F(0),
        "fast_gyro_norm_rad_s": F(0),
    }
    checks = {
        "full_period_in_every_attitude_window": q(m["attitude_excitation"]["T_E_s"]) >= 20,
        "full_period_in_every_displacement_window": q(m["displacement_excitation"]["T_P_s"]) >= 20,
        "gravity_span": half_gravity_chord_lo > half_threshold_hi,
        "displacement_span": 2*amplitude >= q(m["displacement_excitation"]["P_E_m"]),
        "position": amplitude <= q(m["P_max_m"]),
        "velocity": amplitude*nu_hi <= q(m["V_max_mps"]),
        "acceleration": amplitude*nu_hi**2 <= q(m["A_max_mps2"]),
        "jerk": amplitude*nu_hi**3 <= q(m["J_max_mps3"]),
        "angular_rate": amplitude*nu_hi <= q(m["Omega_max_rad_s"]),
        "primitive": 2*amplitude/nu_lo <= q(m["P_AC_max_m_s"]),
        "slow_accel": g*sine <= q(b["B_a_s_mps2"]),
        "field_norm": q(mag["field_norm_min_uT"]) <= 75 <= q(mag["field_norm_max_uT"]),
        "horizontal_field": q(mag["horizontal_field_min_uT"]) <= 75,
        "all_placed_FAST_windows": True,  # both FAST histories are identically zero
    }
    # Bind the literal default BA recurrence used for the metric lower comparison.
    source = HEADER.read_text()
    for literal in ("T sigma_bacc0_ = T(0.004)", "T tau_bacc_ = T(5000.0)",
                    "Matrix3 Q_bacc_ = Matrix3::Identity() * T(2.5e-7)"):
        if literal not in source:
            raise ValueError("shipping BA source changed; re-audit before reproducing this certificate")
    p_ceiling = F(1, 1600)
    stationary = binary32(2.5e-7)*binary32(5000.0)/2
    assert stationary <= p_ceiling and binary32(.004)**2 <= p_ceiling
    lower = (g*sine)**2/p_ceiling
    ca, ha = q(b["fast_accel_accumulation_cap_mps"]), q(b["fast_accel_horizon_s"])
    cg, hg = q(b["fast_gyro_accumulation_cap_rad"]), q(b["fast_gyro_horizon_s"])
    da, dg = q(b["D_a_s_mps3"]), q(b["D_g_s_rad_s2"])
    h = q(limits["sensor_model"]["sample_period_max_s"])
    qa = centered_slow_charge(ca, da, ha, F(14), h)
    qg = centered_slow_charge(cg, dg, hg, F(28), h)
    effective_bias = q(b["B_a_s_mps2"]) + ca/ha
    return {
        "qualification": "OU3_MOVING_QUIET_COMPATIBILITY_V1",
        "result_type": "analytical identities with exact-rational bound reproduction",
        "source_constants_sha256": hashlib.sha256(CONSTANTS.read_bytes()).hexdigest(),
        "shipping_BA_header_sha256": hashlib.sha256(HEADER.read_bytes()).hexdigest(),
        "parameters": {"g": str(g), "sin_theta": str(sine), "cos_theta": str(cosine),
                       "pitch_amplitude_rad": str(amplitude), "displacement_amplitude_m": str(amplitude),
                       "period_s": "20", "nu": "pi/10"},
        "exact_bounds": {k: str(v) for k, v in bounds.items()},
        "physical_membership_checks": checks,
        "marine_and_slow_fast_verified": all(checks.values()),
        "sensor_identity": sensor_polynomial_identity(g, sine, cosine),
        "continuous_moving_fibre": continuous_fibre(g, cosine, sine, amplitude, nu_hi,
                                                     q(b['B_a_s_mps2']), da),
        "quiet_constant_record": {
            "scope": "ideal fixed magnetic record; projection onto attitude, not a full estimator boundedness theorem",
            "effective_total_bias_cap_mps2": str(effective_bias),
            "half_angle_squared_upper": str((effective_bias**2-(g-binary32(9.80665))**2)/(4*g*binary32(9.80665))),
            "constant_FAST_leakage_accel_mps2": str(ca/ha),
            "constant_FAST_leakage_gyro_rad_s": str(cg/hg),
            "centered_window_accel_s": "14", "centered_window_gyro_s": "28",
            "instantaneous_slow_accel_to_total_charge_mps2": str(qa),
            "instantaneous_slow_gyro_to_total_charge_rad_s": str(qg),
            "accel_to_canonical_constant_split_charge_mps2": str(qa+ca/ha),
            "gyro_zero_record_to_canonical_split_charge_rad_s": str(qg),
            "uniform_covariance_metric_tube_certified": False,
        },
        "same_estimator_pair_metric": {
            "scope": "real arithmetic, default BA prior, finite SPD shipping covariance at the queried prefix",
            "BA_marginal_ceiling": str(p_ceiling),
            "literal_stationary_BA_variance": str(stationary),
            "max_of_pair_V_lower": str(lower),
            "comparison_target": "9/400",
            "strictly_above_target": lower > F(9, 400),
            "estimator_symmetry_needed_for_pair_lower": False,
            "BA_physical_forcing_dropped": False,
        },
        "classification": {
            "two_epoch_span_plus_persistent_bias_separation": "B_EXACT_PHYSICAL_COMPATIBILITY_OBSTRUCTION",
            "full_shipping_MOVING_absolute_entry": "OPEN_PENDING_ALL_TIME_ADMISSION",
        },
        "actually_applied_magnetic_service_all_time_verified": False,
        "finite_carried_magnetic_audit_is_all_time_certificate": False,
        "shipping_instability_counterexample": False,
        "universal_MOVING_entry_refuted_on_fully_admitted_shipping_history": False,
        "additional_physical_assumptions": [],
        "shipping_behavior_changed": False,
        "theorem_closed": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    text = json.dumps(certificate(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
