"""Regime contracts and exact algebra; no shipping mode switch is enabled.

See docs/ou3-regime-design.md. Quiet evidence is necessary, never sufficient
for physical STILL. All supplied prefix coefficients are conditional bounds,
not certificates manufactured from a finite replay.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import math

from .marine_motion import norm, sub


def _positive(*values):
    if not all(math.isfinite(x) and x > 0 for x in values):
        raise ValueError("positive finite constants required")


def moving_window_requirement(*, episode_start, episode_end, window_start,
                              window_s, theta_e, span=None):
    """Quantify only over complete windows within ONE physical moving episode.

    An absent span or a boundary-crossing window cannot certify MOVING.
    Episode boundaries are truth-side obligations, never inferred here.
    +infinity is permitted only as an episode end.
    """
    _positive(window_s, theta_e)
    if (theta_e > math.pi or not math.isfinite(episode_start)
            or not math.isfinite(window_start) or math.isnan(episode_end)
            or episode_end <= episode_start):
        raise ValueError("valid physical episode and angular threshold required")
    end = window_start + window_s
    if not math.isfinite(end):
        raise ValueError("window endpoint overflow")
    required = episode_start <= window_start and end <= episode_end
    valid_span = (span is not None and math.isfinite(span) and 0 <= span <= math.pi)
    return {"excitation_required": required,
            "complete_excited_window": bool(required and valid_span and span >= theta_e),
            "regime_certificate": False}


def bridge_prefixes(root_radius, operations):
    """Exact conditional sqrt(V) bounds, including each reset/storage change.

    Operation (g,s) means r_next <= g*r+s in the ACTUAL carried storage.
    g need not be <1. Caller must prove every coefficient on the retained set.
    This enters rho via the endpoint gain and enters c_d through the supply.
    """
    root = F(root_radius)
    if root < 0:
        raise ValueError("nonnegative root radius required")
    gain, supply = F(1), F(0)
    out = []
    for g, s in operations:
        g, s = F(g), F(s)
        if min(g, s) < 0:
            raise ValueError("nonnegative operation bounds required")
        gain, supply = g * gain, g * supply + s
        out.append({"gain": gain, "supply": supply,
                    "radius": gain * root + supply})
    return out


def bridge_retained(root_radius, operations, prefix_radii):
    bounds = bridge_prefixes(root_radius, operations)
    radii = tuple(map(F, prefix_radii))
    if not bounds or len(bounds) != len(radii) or any(x <= 0 for x in radii):
        raise ValueError("one positive retained radius per nonempty prefix required")
    return all(row["radius"] <= radius for row, radius in zip(bounds, radii))


def squared_bridge(gain, supply, eta):
    g, s, eta = map(F, (gain, supply, eta))
    if min(g, s) < 0 or eta <= 0:
        raise ValueError("nonnegative bridge and positive Young parameter required")
    return {"storage_gain": (1 + eta) * g * g,
            "storage_supply": (1 + 1 / eta) * s * s}


def stationary_gyro_average_radius(fast_action_bound, bias_rate, sample_ages, weights,
                                  angular_rate_bound=0, *, slow_amplitude):
    """Terminal SLOW-bias error with a supplied same-history FAST action bound.

    `fast_action_bound` bounds |sum weights_i b_g_f,i|, obtained from the
    reachable window set (e.g. imu_temporal.fast_weighted_outer), NOT from
    dividing a noise sigma by sqrt(N). An amplitude-only bound remains a
    conservative inequality, but is not temporal qualification or admission.
    This identity neither certifies STILL nor changes the shipping estimator.
    """
    n, d, w, bs = map(F, (fast_action_bound, bias_rate, angular_rate_bound, slow_amplitude))
    ages, weights = tuple(map(F, sample_ages)), tuple(map(F, weights))
    if (min(n, d, w, bs) < 0 or not ages or len(ages) != len(weights)
            or min(ages) < 0 or min(weights) < 0 or sum(weights) != 1):
        raise ValueError("nonnegative bounds/ages and normalized weights required")
    return n + w + sum(min(2*bs,d*a)*c for a,c in zip(ages,weights))



@dataclass(frozen=True)
class QuietEvidenceLimits:
    """Amplitude-only necessary screen; NEVER fast temporal qualification."""
    gravity: float
    accel_slow_bias: float
    accel_slow_rate: float
    accel_fast_amplitude: float
    gyro_slow_bias: float
    gyro_slow_rate: float
    gyro_fast_amplitude: float
    dwell_s: float
    max_gap_s: float

    def __post_init__(self):
        _positive(self.gravity, self.dwell_s, self.max_gap_s)
        if self.dwell_s <= self.max_gap_s:
            raise ValueError("dwell must exceed one permitted packet gap")
        for x in (self.accel_slow_bias, self.accel_slow_rate, self.accel_fast_amplitude,
                  self.gyro_slow_bias, self.gyro_slow_rate, self.gyro_fast_amplitude):
            if not math.isfinite(x) or x < 0:
                raise ValueError("finite nonnegative sensor bounds required")


class QuietEvidenceMonitor:
    """Constant-memory necessary-evidence monitor used to falsify detectors.

    An embedded equivalent needs only two anchor vectors and elapsed times.
    Anchors are renewed after each completed dwell while preserving compatible
    status. Entry requires a fresh full dwell after any failure; exit is
    immediate. No variance/noise averaging premise is invented. This helper
    has no access to the shipping estimator and cannot modify it.
    """
    def __init__(self, limits: QuietEvidenceLimits):
        self.limits = limits
        self.reset()

    def reset(self):
        self.anchor = None
        self.previous_time = None
        self.compatible = False

    @property
    def certified_fast_history(self):
        return False

    @property
    def certified_still(self):
        return False

    def step(self, time_s, gyro, accel):
        lim = self.limits
        try:
            values_ok = (math.isfinite(time_s)
                         and norm(gyro) <= lim.gyro_slow_bias + lim.gyro_fast_amplitude
                         and abs(norm(accel) - lim.gravity) <= lim.accel_slow_bias + lim.accel_fast_amplitude)
        except (ValueError, TypeError, OverflowError):
            values_ok = False
        if not values_ok:
            self.reset()
            return "TRANSITION"
        if self.previous_time is not None:
            gap = time_s - self.previous_time
            if not 0 < gap <= lim.max_gap_s:
                self.reset()
        self.previous_time = time_s
        if self.anchor is not None:
            t0, g0, a0 = self.anchor
            age = time_s - t0
            acc_bound = min(2 * lim.accel_slow_bias, lim.accel_slow_rate * age) + 2 * lim.accel_fast_amplitude
            gyro_bound = min(2 * lim.gyro_slow_bias, lim.gyro_slow_rate * age) + 2 * lim.gyro_fast_amplitude
            if norm(sub(accel, a0)) > acc_bound or norm(sub(gyro, g0)) > gyro_bound:
                self.compatible = False
                self.anchor = None
            elif age >= lim.dwell_s:
                self.compatible = True
                self.anchor = None
        if self.anchor is None:
            self.anchor = (time_s, tuple(gyro), tuple(accel))
        return "STILL_COMPATIBLE" if self.compatible else "TRANSITION"


def certificate():
    """Exact sin^3 witness bounds plus conditional bridge audit.

    Trigonometric identities are proved in the design note; only their
    rational global majorants are evaluated here, never sampled extrema.
    """
    import json
    from pathlib import Path
    from .sampled_capture_obstruction import service_certificate
    c = json.loads(Path(__file__).with_name("constants.json").read_text(), parse_float=F)
    alpha, nu, gravity = F(1, 1000), F(1, 40), F("9.80665")
    literal_gravity = F("9.8066501617431640625")
    representation_offset = abs(literal_gravity-gravity)
    bounds = {"B_a_s_mps2": gravity * alpha + representation_offset, "D_a_s_mps3": 3 * gravity * alpha * nu,
              "B_g_s_rad_s": 3 * alpha * nu, "D_g_s_rad_s2": 9 * alpha * nu**2}
    margins = {key: F(c["imu_bias"][key]) - value for key, value in bounds.items()}
    assert all(x > 0 for x in margins.values())
    assert bounds["B_g_s_rad_s"] < c["marine_motion"]["Omega_max_rad_s"]
    service = F(service_certificate()["actual_innovation_service_lower"])
    service *= (1 - alpha**2 / 2)**2
    assert service > c["magnetic_service"]["mu_M"]

    return {
        "qualification": "OU3_REGIME_CONTRACT_V2_SLOW_FAST",
        "bias_model": "two-timescale; witness uses slow components only",
        "fast_accel_and_gyro_identically_zero": True,
        "numerical_device_membership_claimed": False,
        "excitation_quantifier": "every complete T_E window contained in one physical moving episode",
        "boundary_crossing_window_requires_excitation": False,
        "witness_alpha_rad": str(alpha), "witness_nu_rad_s": str(nu),
        "witness_T_E_s": "80*pi", "witness_theta_E_rad": str(2 * alpha),
        "witness_bias_global_upper": {k: str(v) for k, v in bounds.items()},
        "witness_bias_margin": {k: str(v) for k, v in margins.items()},
        "witness_translation_jerk_primitive": "identically zero",
        "witness_actual_service_lower": str(service),
        "rest_motion_joins_C2": True,
        "exact_identical_IMU_history": True,
        "constant_slow_gravity_representation_offset_mps2": str(representation_offset),
        "universal_sound_exact_rest_detector_with_entry_and_exit_liveness": False,
        "finite_exit_delay_from_existing_assumptions": False,
        "stationary_gyro_noise_does_not_average_away": True,
        "stationary_attitude_BA_separation": False,
        "stationary_observable_quotient_boundedness": False,
        "finite_bridge_algebra": True,
        "source_uniform_transition_retention": False,
        "finite_bridges_alone_imply_recurring_stability": False,
        "switching_counterexample_norm_cycle_gain": "3/2",
        "shipping_mode_switch_enabled": False,
        "theorem_closed": False,
    }
