#!/usr/bin/env python3
"""Strict process-excitation lower bounds for the complete OU-III H/A states.

The translational certificate supplies a strict 12-state process covariance
lower bound.  This producer closes the remaining prediction blocks using only
the configured estimator/source constants:

* attitude + gyro bias (6 states), and
* active accelerometer bias (3 states).

For the configured isotropic attitude/bias process covariance the implementation
forms

    Q_AA = [[q_g h I + Q_BB, q_b I_B],
            [q_b I_B^T,       q_b h I]],

where Q_BB is a positive-semidefinite Simpson sum and
I_B = integral_0^h B(s) ds.  Rotation matrices have unit norm and therefore
||B(s)|| <= s, so ||I_B|| <= h^2/2.  Dropping Q_BB and applying the scalar
2x2 Gershgorin lower bound gives

    lambda_min(Q_AA) >= min(q_g h, q_b h) - q_b h^2/2 > 0.

This is independent of angular rate; no replay gyro maximum is needed.

When the accelerometer bias is active, the implemented exact Gauss-Markov block
is

    Q_ba,d = Q_ba * [-tau_b/2 expm1(-2h/tau_b)],

which is enclosed with the validated small-argument ``expm1`` backend.

Since prediction process covariance is block diagonal across these groups, the
complete H/A lower bound is the minimum of their block lower bounds.  The
result may be extremely conservative (the integrated-OU S direction is weak on
a single 5 ms sample), but it is strictly positive and source-uniform.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re

from ou3_interval import Interval
import ou3_source_domain_contract as SOURCE
import ou3_translational_uco_ucc as TRANS
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
SIM = REPO / "src" / "util" / "W3dSimCommon.h"
CERT_SIM = REPO / "tests" / "kalman_ou_iii" / "ou3-certificate-sim.cpp"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _one(pattern: str, text: str, label: str) -> float:
    m = re.search(pattern, text, flags=re.MULTILINE | re.DOTALL)
    if not m:
        raise RuntimeError(f"cannot extract configured constant {label}")
    return float(m.group(1))


def _constants() -> dict:
    mekf = MEKF.read_text(encoding="utf-8")
    wrapper = WRAPPER.read_text(encoding="utf-8")
    sim = SIM.read_text(encoding="utf-8")
    cert = CERT_SIM.read_text(encoding="utf-8")

    gyr_sigma = _one(
        r"const\s+float\s+gyr_sigma\s*=\s*([0-9.eE+-]+)f",
        sim, "sim gyro sigma",
    )
    gyr_init_mult = _one(
        r"Vector3f\s+sigma_g\(\s*([0-9.eE+-]+)f\s*\*\s*gyr_sigma",
        sim, "sigma_g multiplier",
    )
    gyr_rescale = _one(
        r"constexpr\s+float\s+kSigmaGRescale\s*=\s*([0-9.eE+-]+)f",
        cert, "certificate gyro rescale",
    )
    b0 = _one(r"float\s+b0\s*=\s*([0-9.eE+-]+)f", wrapper, "gyro bias RW b0")
    q_ba = _one(
        r"Q_bacc_\s*=\s*Matrix3::Identity\(\)\s*\*\s*T\(([0-9.eE+-]+)\)",
        mekf, "accelerometer bias process density",
    )
    tau_ba = _one(
        r"tau_bacc_\s*=\s*T\(([0-9.eE+-]+)\)",
        mekf, "accelerometer bias correlation time",
    )
    return {
        "gyro_noise_density_rad_sqrt_s": gyr_sigma * gyr_init_mult * gyr_rescale,
        "gyro_bias_rw_variance_density": b0,
        "accel_bias_process_variance_density": q_ba,
        "accel_bias_tau_s": tau_ba,
    }


def build() -> dict:
    source = SOURCE.build(SOURCE.DEFAULT_HEADER.resolve())
    runtime = source["configured_runtime_assumption"]
    h = Interval(*runtime["imu_dt_outward_interval_s"])
    c = _constants()
    trans = TRANS.build(TRANS.DEFAULT_HEADER)
    q_trans = float(trans["process_ucc"]["Q_axis_lambda_min_lower"])

    qg = down(c["gyro_noise_density_rad_sqrt_s"] ** 2)
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
        "qualification": "FULL_STATE_PREDICTION_PROCESS_UCC_CONFIGURED_RUNTIME",
        "source_generated_not_trajectory_fit": True,
        "validated_arithmetic": True,
        "outward_rounded": True,
        "configured_runtime": runtime,
        "source_constants": c,
        "attitude_gyro_bias": {
            "q_gyro_lower": qg,
            "q_gyro_bias_lower": qb,
            "theta_diagonal_lower": diag_theta,
            "gyro_bias_diagonal_lower": diag_bg,
            "cross_norm_upper": cross,
            "Q_attitude_gyro_bias_lambda_min_lower": q_att_bg,
            "derivation": "min(q_g*h,q_b*h)-q_b*h^2/2; PSD Simpson BQB term dropped",
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
            "combine full process UCC with translational and conditional vector UCO to derive "
            "uniform covariance bounds and a strict source-word information contraction"
        ),
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for flag in ("source_generated_not_trajectory_fit", "validated_arithmetic",
                 "outward_rounded", "full_process_ucc_pass", "pass"):
        if d.get(flag) is not True:
            failures.append(f"{flag} is not true")
    for mode in ("H", "A"):
        row = d.get("modes", {}).get(mode, {})
        q = row.get("prediction_Q_lambda_min_lower")
        if row.get("pass") is not True or not isinstance(q, (int, float)) or not math.isfinite(float(q)) or not float(q) > 0.0:
            failures.append(f"mode {mode} process lower bound is not strict")
    if d.get("continuous_word_enclosed") is not False:
        failures.append("process stage must not assert full word enclosure")
    if d.get("theorem_promotion") != "NOT_ESTABLISHED":
        failures.append("process stage must not promote theorem")
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
        "H_Q_min": d["modes"]["H"]["prediction_Q_lambda_min_lower"],
        "A_Q_min": d["modes"]["A"]["prediction_Q_lambda_min_lower"],
        "pass": d["pass"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
