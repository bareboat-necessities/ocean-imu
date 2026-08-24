#!/usr/bin/env python3
"""Explicit vector-packet UCO certificate for the configured OU-III deployment.

Full 3-D attitude/gyro-bias stability cannot be unconditional: one vector leaves
rotation about itself unobservable, and no estimator can recover heading from an
indefinitely collinear specific-force/magnetic history.  The theorem therefore
states the weakest missing physical hypothesis explicitly instead of pretending
that a source clamp proves persistent excitation.

The certified family is the configured 200 Hz / 25 Hz-magnetometer deployment
restricted to the following *operating-envelope hypothesis*:

* each proof packet has nonzero specific-force and magnetic vectors;
* their sine separation is at least 0.01 (about 0.57 deg);
* body rate over the two-packet interval is at most 30 deg/s.

The norm floor is only 1e-3 in physical units.  The magnetic value is no stronger
than the wrapper's own ``mag_init_min_mag_norm`` guard, and the acceleration
value merely excludes exact free-fall cancellation.  The angular-separation and
Live rate bounds are theorem hypotheses, not values fitted from the eight wave
replays.  They are intentionally weak; a deployment outside this PE envelope is
outside the full-heading theorem rather than silently declared stable.

For an accepted vector pair f,m with measurement variances r_a,r_m,

  mu_theta = min(f_min^2/r_a, m_min^2/r_m)
             * (1 - sqrt(1-s_fm^2)).

For two packets separated by Delta_g and nominal rate bounded by omega_bar,

  sigma_min(Gamma_g) >= Delta_min * (1 - omega_bar Delta_max/2) = g_min.

With T_bg=1 s, gamma_min=g_min/T_bg and the scaled six-state
attitude/gyro-bias observability Gramian obeys

  alpha_6 >= mu_theta / (1 + 2/gamma_min^2) > 0.

All arithmetic is widened with ``nextafter``.  ``sqrt`` is used only on positive
scalar endpoints and then widened; no trajectory-derived extrema enter.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re

import ou3_source_domain_contract as SOURCE

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
SIM = REPO / "src" / "util" / "W3dSimCommon.h"
CERT_SIM = REPO / "tests" / "kalman_ou_iii" / "ou3-certificate-sim.cpp"
SCHEMA = 1

# The genuinely physical theorem hypothesis.  These are not estimated from a
# replay and are deliberately much weaker than normal marine operation.
PE = {
    "specific_force_norm_lower_mps2": 1.0e-3,
    "vector_sine_separation_lower": 1.0e-2,
    "body_rate_norm_upper_deg_s": 30.0,
    "gyro_bias_time_scale_s": 1.0,
}


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _one(pattern: str, text: str, label: str) -> float:
    m = re.search(pattern, text, flags=re.MULTILINE | re.DOTALL)
    if not m:
        raise RuntimeError(f"cannot extract configured constant {label}")
    return float(m.group(1))


def _configured_constants() -> dict:
    wrapper = WRAPPER.read_text(encoding="utf-8")
    sim = SIM.read_text(encoding="utf-8")
    cert = CERT_SIM.read_text(encoding="utf-8")

    g = 9.80665
    acc_base_mult = _one(
        r"const\s+float\s+acc_sigma\s*=\s*([0-9.eE+-]+)f\s*\*\s*g_std",
        sim, "acc_sigma/g",
    )
    acc_init_mult = _one(
        r"Vector3f\s+sigma_a_init\(\s*([0-9.eE+-]+)f\s*\*\s*acc_sigma",
        sim, "sigma_a_init multiplier",
    )
    acc_rescale = _one(
        r"constexpr\s+float\s+kSigmaARescale\s*=\s*([0-9.eE+-]+)f",
        cert, "certificate accel rescale",
    )
    mag_odr = _one(
        r"constexpr\s+float\s+kMagOdrHz\s*=\s*([0-9.eE+-]+)f",
        cert, "certificate mag ODR",
    )
    mag_hi = _one(
        r"const\s+float\s+mag_sigma_uT\s*=\s*\(mag_odr_hz\s*<=\s*20\.0f\)\s*\?\s*[0-9.eE+-]+f\s*:\s*([0-9.eE+-]+)f",
        sim, "mag high-ODR sigma",
    )
    mag_init_mult = _one(
        r"const\s+float\s+sigma_m_uT\s*=\s*([0-9.eE+-]+)f\s*\*\s*mag_sigma_uT",
        sim, "sigma_m multiplier",
    )
    mag_rescale = _one(
        r"constexpr\s+float\s+kSigmaMRescale\s*=\s*([0-9.eE+-]+)f",
        cert, "certificate mag rescale",
    )
    mag_norm_guard = _one(
        r"float\s+mag_init_min_mag_norm\s*=\s*([0-9.eE+-]+)f",
        wrapper, "mag_init_min_mag_norm",
    )
    configured_motion_guard = _one(
        r"float\s+mag_extreme_gyro_dps\s*=\s*([0-9.eE+-]+)f",
        wrapper, "mag_extreme_gyro_dps",
    )

    return {
        "acc_measurement_std_mps2": acc_base_mult * g * acc_init_mult * acc_rescale,
        "mag_measurement_std_uT": mag_hi * mag_init_mult * mag_rescale,
        "mag_odr_hz": mag_odr,
        "mag_init_norm_guard_uT": mag_norm_guard,
        "startup_extreme_motion_guard_deg_s": configured_motion_guard,
    }


def build() -> dict:
    c = _configured_constants()

    f_min = float(PE["specific_force_norm_lower_mps2"])
    # Do not claim a magnetic norm floor stronger than the actual acquisition
    # guard.  The theorem can use the smaller of its declared floor and source.
    m_min = min(1.0e-3, c["mag_init_norm_guard_uT"])
    s = float(PE["vector_sine_separation_lower"])
    omega_deg = float(PE["body_rate_norm_upper_deg_s"])
    Tbg = float(PE["gyro_bias_time_scale_s"])

    if not (f_min > 0 and m_min > 0 and 0 < s < 1 and omega_deg > 0 and Tbg > 0):
        raise RuntimeError("invalid PE operating envelope")

    ra = up(c["acc_measurement_std_mps2"] ** 2)
    rm = up(c["mag_measurement_std_uT"] ** 2)
    af = down((f_min * f_min) / ra)
    am = down((m_min * m_min) / rm)

    # 1-sqrt(1-s^2), evaluated stably as s^2/(1+sqrt(1-s^2)).
    root = up(math.sqrt(max(0.0, 1.0 - s * s)))
    angular_factor = down((s * s) / up(1.0 + root))
    mu_theta = down(min(af, am) * angular_factor)

    # Every configured 25 Hz magnetic tick is accompanied by the 200 Hz accel
    # stream.  The two-packet analytical bound uses consecutive packets.
    delta = 1.0 / c["mag_odr_hz"]
    delta_lo = down(delta)
    delta_hi = up(delta)
    omega = up(omega_deg * math.pi / 180.0)
    bracket = down(1.0 - up(0.5 * omega * delta_hi))
    g_min = down(delta_lo * bracket)
    gamma = down(g_min / Tbg)
    alpha6 = down(mu_theta / up(1.0 + up(2.0 / down(gamma * gamma))))

    pass_ = all(math.isfinite(x) and x > 0.0 for x in (
        ra, rm, af, am, angular_factor, mu_theta, bracket, g_min, gamma, alpha6
    )) and omega * delta_hi < 2.0

    return {
        "schema": SCHEMA,
        "qualification": "CONDITIONAL_VECTOR_PACKET_UCO_CONFIGURED_DEPLOYMENT",
        "claim": "ATTITUDE_GYRO_BIAS_UCO_UNDER_EXPLICIT_PERSISTENT_EXCITATION",
        "trajectory_fit": False,
        "unconditional_full_heading_claim": False,
        "persistent_excitation_is_theorem_hypothesis": True,
        "operating_envelope": {
            **PE,
            "magnetic_vector_norm_lower_uT": m_min,
            "packet_gap_s": [delta_lo, delta_hi],
            "note": (
                "The sine-separation and Live body-rate bounds are deployment theorem "
                "hypotheses. They are not inferred from the eight validation trajectories."
            ),
        },
        "configured_measurement_bounds": {
            **c,
            "acc_measurement_variance_upper": ra,
            "mag_measurement_variance_upper": rm,
        },
        "vector_pair": {
            "a_f_lower": af,
            "a_m_lower": am,
            "angular_factor_lower": angular_factor,
            "mu_theta_lower": mu_theta,
        },
        "gyro_bias_two_packet": {
            "omega_bar_rad_s_upper": omega,
            "omega_times_gap_upper": up(omega * delta_hi),
            "gamma_bracket_lower": bracket,
            "Gamma_g_sigma_min_lower_s": g_min,
            "scaled_gamma_lower": gamma,
            "alpha_6_information_lower": alpha6,
        },
        "conditional_source_complete": pass_,
        "pass": pass_,
        "theorem_scope": (
            "full-heading Live UES/ISS is conditional on the stated persistent-excitation "
            "operating envelope; exact collinearity/free-fall histories are excluded as "
            "fundamentally unobservable"
        ),
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("trajectory_fit") is not False:
        failures.append("PE certificate must not be trajectory fitted")
    if d.get("unconditional_full_heading_claim") is not False:
        failures.append("certificate must not claim unconditional heading observability")
    if d.get("persistent_excitation_is_theorem_hypothesis") is not True:
        failures.append("PE hypothesis must be explicit")
    if d.get("conditional_source_complete") is not True or d.get("pass") is not True:
        failures.append("conditional vector UCO did not close")
    for section, key in (
        ("vector_pair", "mu_theta_lower"),
        ("gyro_bias_two_packet", "Gamma_g_sigma_min_lower_s"),
        ("gyro_bias_two_packet", "alpha_6_information_lower"),
    ):
        x = d.get(section, {}).get(key)
        if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or not float(x) > 0.0:
            failures.append(f"{section}.{key} is not finite positive")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "pass": d["pass"],
        "mu_theta_lower": d["vector_pair"]["mu_theta_lower"],
        "alpha_6_information_lower": d["gyro_bias_two_packet"]["alpha_6_information_lower"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
