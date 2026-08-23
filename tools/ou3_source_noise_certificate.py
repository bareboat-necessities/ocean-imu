#!/usr/bin/env python3
"""Source-derived stochastic primitive-noise certificate for adaptive OU-III.

This certificate deliberately does not estimate noise from replay trajectories.
It reads the simulator/source constants that define the OU-III validation noise
model and expresses each physical Gaussian draw as a scale times a standardized
primitive coordinate.  With six 3-vectors reserved per IMU sample,

    z_k = [z_a, z_g, z_ba, z_bg, z_m, z_bm] in R^18,

magnetometer coordinates are simply unused (zero gain) on non-mag samples.
Consequently the pre-gate standardized conditional covariance satisfies

    0 <= Cov(z_k | F_k) <= I_18

uniformly.  Physical units stay in the source-map sensitivity, not in this
covariance bound.  This is the normalization used by the stochastic enclosure
backend and avoids mixing m/s^2, rad/s and microtesla in one Euclidean norm.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_CERT = REPO / "reports" / "results" / "ou3_numerical_certificate"
SIM_H = REPO / "src" / "util" / "W3dSimCommon.h"
SIM_CPP = REPO / "src" / "util" / "W3dSimCommon.cpp"
CERT_SIM = REPO / "tests" / "kalman_ou_iii" / "ou3-certificate-sim.cpp"
SCHEMA = 1
PRIMITIVE_DIMENSION = 18


def _one(pattern: str, text: str, label: str) -> float:
    m = re.search(pattern, text, flags=re.MULTILINE | re.DOTALL)
    if not m:
        raise RuntimeError(f"cannot locate source noise constant {label}")
    return float(m.group(1))


def _source_constants() -> dict:
    h = SIM_H.read_text(encoding="utf-8")
    cpp = SIM_CPP.read_text(encoding="utf-8")
    cert = CERT_SIM.read_text(encoding="utf-8")

    g = _one(r"g_std\s*=\s*([0-9.eE+-]+)f", cpp, "g_std")
    dt_den = _one(r"constexpr\s+float\s+kDt\s*=\s*1\.0f\s*/\s*([0-9.eE+-]+)f", cert, "kDt")
    dt = 1.0 / dt_den
    mag_odr = _one(r"constexpr\s+float\s+kMagOdrHz\s*=\s*([0-9.eE+-]+)f", cert, "kMagOdrHz")
    mag_dt = 1.0 / mag_odr

    acc_mult = _one(r"const\s+float\s+acc_sigma\s*=\s*([0-9.eE+-]+)f\s*\*\s*g_std", h, "acc_sigma/g")
    gyro_white = _one(r"const\s+float\s+gyr_sigma\s*=\s*([0-9.eE+-]+)f", h, "gyr_sigma")
    acc_bias_mult = _one(r"const\s+float\s+acc_bias_range\s*=\s*([0-9.eE+-]+)f\s*\*\s*g_std", h, "acc_bias_range/g")
    gyro_bias_deg = _one(r"const\s+float\s+gyr_bias_range\s*=\s*([0-9.eE+-]+)f\s*\*\s*float\(std::numbers::pi_v<float>\s*/\s*180\.0f\)", h, "gyr_bias_range_deg")
    acc_bias_rw = _one(r"const\s+float\s+acc_bias_rw\s*=\s*([0-9.eE+-]+)f", h, "acc_bias_rw")
    gyro_bias_rw = _one(r"const\s+float\s+gyr_bias_rw\s*=\s*([0-9.eE+-]+)f", h, "gyr_bias_rw")

    m = re.search(
        r"const\s+float\s+mag_sigma_uT\s*=\s*\(mag_odr_hz\s*<=\s*20\.0f\)\s*\?\s*([0-9.eE+-]+)f\s*:\s*([0-9.eE+-]+)f",
        h,
    )
    if not m:
        raise RuntimeError("cannot locate mag_sigma_uT source branch")
    mag_sigma_low = float(m.group(1))
    mag_sigma_high = float(m.group(2))
    mag_white = mag_sigma_low if mag_odr <= 20.0 else mag_sigma_high

    m = re.search(
        r"make_mag_noise_model\(\s*mag_sigma_uT\s*,\s*([0-9.eE+-]+)f\s*,\s*([0-9.eE+-]+)f\s*,\s*([0-9.eE+-]+)f\s*,\s*([0-9.eE+-]+)f\s*,\s*([0-9.eE+-]+)f",
        h,
        flags=re.MULTILINE,
    )
    if not m:
        raise RuntimeError("cannot locate make_mag_noise_model source parameters")
    mag_bias_range, mag_bias_rw, mag_scale, mag_cross, mag_misalign_deg = map(float, m.groups())

    return {
        "g_std_mps2": g,
        "imu_dt_s": dt,
        "imu_rate_hz": 1.0 / dt,
        "mag_odr_hz": mag_odr,
        "mag_dt_s": mag_dt,
        "acc_white_std_mps2": acc_mult * g,
        "gyro_white_std_radps": gyro_white,
        "acc_bias_rw_density_mps2_sqrt_s": acc_bias_rw,
        "gyro_bias_rw_density_radps_sqrt_s": gyro_bias_rw,
        "acc_bias_rw_increment_std_mps2": acc_bias_rw * math.sqrt(dt),
        "gyro_bias_rw_increment_std_radps": gyro_bias_rw * math.sqrt(dt),
        "mag_white_std_uT": mag_white,
        "mag_bias_rw_density_uT_sqrt_s": mag_bias_rw,
        "mag_bias_rw_increment_std_uT": mag_bias_rw * math.sqrt(mag_dt),
        "acc_initial_bias_half_range_mps2": acc_bias_mult * g,
        "gyro_initial_bias_half_range_radps": gyro_bias_deg * math.pi / 180.0,
        "mag_initial_bias_half_range_uT": mag_bias_range,
        "mag_scale_error_abs_max": mag_scale,
        "mag_cross_axis_abs_max": mag_cross,
        "mag_misalignment_abs_max_deg": mag_misalign_deg,
    }


def build_certificate() -> dict:
    c = _source_constants()
    d = PRIMITIVE_DIMENSION
    return {
        "schema": SCHEMA,
        "claim": "OU3_SOURCE_GAUSSIAN_PRIMITIVE_NOISE_CERTIFICATE",
        "qualification": "SOURCE_GENERATED_NOT_TRAJECTORY_FIT",
        "source_generated_not_trajectory_fit": True,
        "standardized_increment": {
            "definition": "z=[acc_white(3),gyro_white(3),acc_bias_rw(3),gyro_bias_rw(3),mag_white(3),mag_bias_rw(3)]",
            "dimension": d,
            "covariance_upper_identity": True,
            "Sigma_bar_norm_upper": 1.0,
            "trace_Sigma_bar_upper": float(d),
            "trace_Sigma_bar_squared_upper": float(d),
            "s2_upper": float(d),
            "s4_upper": float(d * d + 2 * d),
            "mag_coordinates_zero_gain_between_mag_ticks": True,
        },
        "physical_scales": {
            k: c[k]
            for k in (
                "imu_dt_s", "imu_rate_hz", "mag_odr_hz", "mag_dt_s",
                "acc_white_std_mps2", "gyro_white_std_radps",
                "acc_bias_rw_density_mps2_sqrt_s", "gyro_bias_rw_density_radps_sqrt_s",
                "acc_bias_rw_increment_std_mps2", "gyro_bias_rw_increment_std_radps",
                "mag_white_std_uT", "mag_bias_rw_density_uT_sqrt_s",
                "mag_bias_rw_increment_std_uT",
            )
        },
        "initialization_support": {
            "distribution": "BOUNDED_INITIALIZATION_SUPPORT_NOT_PER_SAMPLE_GAUSSIAN_NOISE",
            "acc_bias_half_range_mps2": c["acc_initial_bias_half_range_mps2"],
            "gyro_bias_half_range_radps": c["gyro_initial_bias_half_range_radps"],
            "mag_bias_half_range_uT": c["mag_initial_bias_half_range_uT"],
            "mag_scale_error_abs_max": c["mag_scale_error_abs_max"],
            "mag_cross_axis_abs_max": c["mag_cross_axis_abs_max"],
            "mag_misalignment_abs_max_deg": c["mag_misalignment_abs_max_deg"],
        },
        "proof_use": {
            "normalization": "physical increments equal diagonal/source map scales times standardized z coordinates",
            "conditional_covariance": "0 <= Cov(z_k|F_k) <= I_18",
            "gating_note": "bound is pre-gate; accepted innovations need not remain conditionally centered",
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--certificate-dir", type=Path, default=DEFAULT_CERT)
    args = ap.parse_args()
    out = args.certificate_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    report = build_certificate()
    path = out / "source_noise_certificate.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "claim": report["claim"],
        "dimension": report["standardized_increment"]["dimension"],
        "Sigma_bar_norm_upper": report["standardized_increment"]["Sigma_bar_norm_upper"],
        "s2_upper": report["standardized_increment"]["s2_upper"],
        "s4_upper": report["standardized_increment"]["s4_upper"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
