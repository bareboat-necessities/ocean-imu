#!/usr/bin/env python3
"""Strict process-excitation bounds for complete shipping OU-III H/A states.

This producer is source-bound to the same deployed ``SeaStateFusion_OU_III::Config``
that feeds canonical P3.  In particular the attitude process uses the shipping
Config gyro-noise density ``sigma_g=(0.01,0.01,0.01)`` and Config gyro-bias
random-walk density ``b0=1e-11`` passed through ``initialize_ext`` into the MEKF
constructor.  Historical simulation noise and ``ou3-certificate-sim`` rescale
constants are deliberately not inputs to this certificate.

The translational certificate supplies a strict 12-state process covariance
lower.  This producer closes the remaining prediction blocks:

* attitude + gyro bias (6 states), and
* active accelerometer bias (3 states).

For the configured isotropic attitude/bias process covariance the implementation
forms

    Q_AA = [[q_g h I + Q_BB, q_b I_B],
            [q_b I_B^T,       q_b h I]],

where Q_BB is positive semidefinite and ||I_B||<=h^2/2.  Dropping Q_BB and
using the exact block norm gives the rigorous lower

    lambda_min(Q_AA) >= min(q_g h, q_b h) - q_b h^2/2 > 0.

When accelerometer bias is active, shipping uses the exact Gauss--Markov block

    Q_ba,d = Q_ba * [-tau_b/2 expm1(-2h/tau_b)],

which is enclosed with the validated transcendental backend.

The single-sample full-state minimum remains intentionally only a strict UCC
primitive; it is not the canonical P3 contraction ratio.  The canonical word
must retain all process matrices jointly with its measurements.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import struct

from ou3_interval import Interval
import ou3_source_domain_contract as SOURCE
import ou3_translational_uco_ucc as TRANS
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
SCHEMA = 2


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _f32(x: float) -> float:
    return struct.unpack("!f", struct.pack("!f", float(x)))[0]


def _one(pattern: str, text: str, label: str) -> float:
    m = re.search(pattern, text, flags=re.MULTILINE | re.DOTALL)
    if not m:
        raise RuntimeError(f"cannot extract configured constant {label}")
    return float(m.group(1))


def _config_vec3(text: str, name: str) -> list[float]:
    m = re.search(
        rf"Eigen::Vector3f\s+{re.escape(name)}\s*=\s*Eigen::Vector3f\(\s*"
        r"([0-9.eE+-]+)f\s*,\s*([0-9.eE+-]+)f\s*,\s*([0-9.eE+-]+)f\s*\)",
        text,
    )
    if not m:
        raise RuntimeError(f"cannot extract shipping Config {name}")
    out = [_f32(float(m.group(i))) for i in range(1, 4)]
    if any(not (math.isfinite(x) and x > 0.0) for x in out):
        raise RuntimeError(f"shipping Config {name} is not strictly positive")
    return out


def _config_scalar(text: str, name: str) -> float:
    m = re.search(rf"\bfloat\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f\s*;", text)
    if not m:
        raise RuntimeError(f"cannot extract shipping Config {name}")
    x = _f32(float(m.group(1)))
    if not (math.isfinite(x) and x > 0.0):
        raise RuntimeError(f"shipping Config {name} is not strictly positive")
    return x


def _constants() -> dict:
    mekf = MEKF.read_text(encoding="utf-8")
    wrapper = WRAPPER.read_text(encoding="utf-8")

    sigma_g = _config_vec3(wrapper, "sigma_g")
    b0 = _config_scalar(wrapper, "b0")
    q_ba = _f32(_one(
        r"Q_bacc_\s*=\s*Matrix3::Identity\(\)\s*\*\s*T\(([0-9.eE+-]+)\)",
        mekf, "accelerometer bias process density",
    ))
    tau_ba = _f32(_one(
        r"tau_bacc_\s*=\s*T\(([0-9.eE+-]+)\)",
        mekf, "accelerometer bias correlation time",
    ))

    parity = {
        "outer_Config_sigma_g_present": all(x > 0.0 for x in sigma_g),
        "outer_Config_b0_present": b0 > 0.0,
        "outer_wrapper_forwards_sigma_g_and_b0": (
            "impl_.initialize_ext(cfg_.sigma_a, cfg_.sigma_g, cfg_.sigma_m," in wrapper
            and "cfg_.Pq0, cfg_.Pb0, cfg_.b0, cfg_.R_S_noise," in wrapper
        ),
        "inner_initialize_ext_forwards_gyro_density_and_b0": (
            "std::make_unique<Kalman3D_Wave_OU_III<float>>(sigma_a, sigma_g, sigma_m, Pq0, Pb0, b0, R_S_noise, gravity_magnitude)" in wrapper
        ),
        "MEKF_constructor_builds_Qbase_from_forwarded_values": (
            "Qbase(initialize_Q(gyro_noise_density_rad_sqrt_s, b0))" in mekf
        ),
        "MEKF_gyro_density_setter_squares_density": (
            "Qbase.template topLeftCorner<3,3>() = d.array().square().matrix().asDiagonal();" in mekf
        ),
        "active_ba_GM_process_source_present": (
            "Q_bacc_ * qd_scale" in mekf and "std::expm1(-T(2) * Ts / tau_b)" in mekf
        ),
    }
    failures = [k for k, v in parity.items() if not v]
    if failures:
        raise RuntimeError(f"shipping process source parity failed: {failures}")

    return {
        "gyro_noise_density_rad_sqrt_s_per_axis": sigma_g,
        "gyro_noise_density_rad_sqrt_s_lower": min(sigma_g),
        "gyro_bias_rw_variance_density": b0,
        "accel_bias_process_variance_density": q_ba,
        "accel_bias_tau_s": tau_ba,
        "measurement_or_process_simulation_rescale_consumed": False,
        "W3dSimCommon_consumed": False,
        "ou3_certificate_sim_consumed": False,
        "source_parity": parity,
    }


def build() -> dict:
    source = SOURCE.build(SOURCE.DEFAULT_HEADER.resolve())
    runtime = source["configured_runtime_assumption"]
    h = Interval(*runtime["imu_dt_outward_interval_s"])
    c = _constants()
    trans = TRANS.build(TRANS.DEFAULT_HEADER)
    q_trans = float(trans["process_ucc"]["Q_axis_lambda_min_lower"])

    sigma_g_lo = float(c["gyro_noise_density_rad_sqrt_s_lower"])
    qg = down(sigma_g_lo * sigma_g_lo)
    qb = down(c["gyro_bias_rw_variance_density"])
    if not (qg > 0.0 and qb > 0.0 and h.lo > 0.0):
        raise RuntimeError("attitude/bias process source lost positivity")

    diag_theta = down(qg * h.lo)
    diag_bg = down(qb * h.lo)
    cross = up(qb * up(h.hi * h.hi) / 2.0)
    q_att_bg = down(min(diag_theta, diag_bg) - cross)

    tau_ba = Interval.outward_bounds(c["accel_bias_tau_s"], c["accel_bias_tau_s"])
    x = Interval.outward_bounds(2.0 * h.lo, 2.0 * h.hi) / tau_ba
    em1 = VT.expm1_interval(-x)
    qd_scale = Interval.outward_bounds(-0.5, -0.5) * tau_ba * em1
    q_ba = down(c["accel_bias_process_variance_density"] * qd_scale.lo)

    q_H = down(min(q_att_bg, q_trans))
    q_A = down(min(q_att_bg, q_trans, q_ba))
    passed = all(math.isfinite(v) and v > 0.0 for v in (
        qg, qb, diag_theta, diag_bg, q_att_bg, q_trans, q_ba, q_H, q_A
    ))

    return {
        "schema": SCHEMA,
        "qualification": "FULL_STATE_PREDICTION_PROCESS_UCC_SHIPPING_CONFIGURED_RUNTIME",
        "source_generated_not_trajectory_fit": True,
        "validated_arithmetic": True,
        "outward_rounded": True,
        "configured_runtime": runtime,
        "shipping_Config_consumed": True,
        "historical_certificate_simulation_constants_consumed": False,
        "source_constants": c,
        "attitude_gyro_bias": {
            "q_gyro_lower": qg,
            "q_gyro_bias_lower": qb,
            "theta_diagonal_lower": diag_theta,
            "gyro_bias_diagonal_lower": diag_bg,
            "cross_norm_upper": cross,
            "Q_attitude_gyro_bias_lambda_min_lower": q_att_bg,
            "derivation": "min(q_g*h,q_b*h)-q_b*h^2/2; PSD Simpson BQB term dropped",
            "rate_independent": True,
        },
        "translation": {
            "Q_translation_lambda_min_lower": q_trans,
            "upstream_translation_pass": trans["translation_source_complete"],
        },
        "active_accelerometer_bias": {
            "two_h_over_tau_interval": x.as_list(),
            "expm1_minus_interval": em1.as_list(),
            "qd_scale_interval_s": qd_scale.as_list(),
            "Q_accel_bias_lambda_min_lower": q_ba,
        },
        "modes": {
            "H": {
                "dimension": 18,
                "prediction_Q_lambda_min_lower": q_H,
                "pass": q_H > 0.0,
            },
            "A": {
                "dimension": 21,
                "prediction_Q_lambda_min_lower": q_A,
                "pass": q_A > 0.0,
            },
        },
        "full_process_ucc_pass": passed,
        "pass": passed,
        "continuous_word_enclosed": False,
        "theorem_promotion": "NOT_ESTABLISHED",
        "next_obligation": (
            "use these shipping source-bound process matrices inside the same complete D/T/Qc or P/Psi/Omega event word; do not form a one-step scalar contraction ratio"
        ),
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("qualification") != "FULL_STATE_PREDICTION_PROCESS_UCC_SHIPPING_CONFIGURED_RUNTIME":
        failures.append("qualification mismatch")
    for flag in (
        "source_generated_not_trajectory_fit", "validated_arithmetic",
        "outward_rounded", "shipping_Config_consumed", "full_process_ucc_pass", "pass",
    ):
        if d.get(flag) is not True:
            failures.append(f"{flag} is not true")
    if d.get("historical_certificate_simulation_constants_consumed") is not False:
        failures.append("certificate-simulation constants entered shipping process UCC")
    constants = d.get("source_constants", {})
    if constants.get("W3dSimCommon_consumed") is not False:
        failures.append("W3dSimCommon entered shipping process UCC")
    if constants.get("ou3_certificate_sim_consumed") is not False:
        failures.append("ou3-certificate-sim entered shipping process UCC")
    if constants.get("gyro_noise_density_rad_sqrt_s_per_axis") != [
        _f32(0.01), _f32(0.01), _f32(0.01)
    ]:
        failures.append("shipping Config gyro noise changed")
    if not all(constants.get("source_parity", {}).values()):
        failures.append("shipping process source parity failed")
    for mode in ("H", "A"):
        row = d.get("modes", {}).get(mode, {})
        q = row.get("prediction_Q_lambda_min_lower")
        if row.get("pass") is not True or not isinstance(q, (int, float)) or not math.isfinite(float(q)) or not float(q) > 0.0:
            failures.append(f"mode {mode} process lower bound is not strict")
    if d.get("continuous_word_enclosed") is not False:
        failures.append("process stage must not assert full word enclosure")
    if d.get("theorem_promotion") != "NOT_ESTABLISHED":
        failures.append("process stage must not promote theorem")
    return list(dict.fromkeys(failures))


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
        "shipping_sigma_g": d["source_constants"]["gyro_noise_density_rad_sqrt_s_per_axis"],
        "shipping_b0": d["source_constants"]["gyro_bias_rw_variance_density"],
        "H_Q_min": d["modes"]["H"]["prediction_Q_lambda_min_lower"],
        "A_Q_min": d["modes"]["A"]["prediction_Q_lambda_min_lower"],
        "pass": d["pass"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())