#!/usr/bin/env python3
"""Exact joint source/filter-state -> typed executor sample-coordinate map.

A complete Normal-Live sample is not a purely physical BRMM output.  Shipping
uses two different coordinate classes in the same update:

* measurement-only frontend coordinates: de-heeled raw gyro and accelerometer
  specific force;
* nominal estimator geometry: bias-corrected gyro rate plus R_hat and
  f_hat=R_hat*(a_w_hat-g_N), used by the MEKF prediction/Jacobian.

Conflating truth attitude/acceleration with R_hat/f_hat, or equating the raw
gyro with the corrected rate, would break the shipping recurrence.  This module
materializes the exact algebraic map while retaining sensor/model mismatch as
explicit forcing coordinates.  It creates no independent source boxes and does
not itself certify a finite disturbance magnitude.

All vectors here are already expressed in the shipping virtual body frame B'.
The deployed proof branch has the IMU lever arm disabled, so the nominal CoG
specific-force geometry is exactly R_hat*(a_w_hat-g_N).
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_brmm_correlated_window_outer_enclosure as OUTER

REPO = Path(__file__).resolve().parents[2]
DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SCHEMA = 1
QUALIFICATION = "OU3_BRMM_JOINT_SOURCE_FILTER_EXECUTOR_COORDINATE_MAP_V1"
P3_DELTA = 1.0e-18

Vec3 = tuple[Interval, Interval, Interval]
Mat3 = Sequence[Sequence[Interval]]


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _vec3(x: Sequence[Interval], name: str) -> Vec3:
    if len(x) != 3 or any(not isinstance(v, Interval) for v in x):
        raise ValueError(name + " must be Interval[3]")
    return (x[0], x[1], x[2])


def _mat3(A: Mat3, name: str) -> Mat3:
    if len(A) != 3 or any(len(row) != 3 for row in A):
        raise ValueError(name + " must be Interval[3x3]")
    if any(not isinstance(x, Interval) for row in A for x in row):
        raise ValueError(name + " entries must be intervals")
    return A


def add(a: Sequence[Interval], b: Sequence[Interval]) -> Vec3:
    a = _vec3(a, "a"); b = _vec3(b, "b")
    return tuple(a[i] + b[i] for i in range(3))  # type: ignore[return-value]


def sub(a: Sequence[Interval], b: Sequence[Interval]) -> Vec3:
    a = _vec3(a, "a"); b = _vec3(b, "b")
    return tuple(a[i] - b[i] for i in range(3))  # type: ignore[return-value]


def matvec(A: Mat3, x: Sequence[Interval]) -> Vec3:
    A = _mat3(A, "A"); x = _vec3(x, "x")
    out = []
    for i in range(3):
        s = I(0.0)
        for j in range(3):
            s = s + A[i][j] * x[j]
        out.append(s)
    return tuple(out)  # type: ignore[return-value]


def nominal_cog_specific_force_body(
    R_hat_wb: Mat3,
    a_w_hat_world: Sequence[Interval],
    gravity_mps2: Interval,
) -> Vec3:
    """Shipping no-lever-arm nominal f_hat=R_hat*(a_w_hat-g_N)."""
    if gravity_mps2.lo <= 0.0:
        raise ValueError("positive gravity required")
    g = (I(0.0), I(0.0), gravity_mps2)
    return matvec(R_hat_wb, sub(a_w_hat_world, g))


def corrected_gyro_rate_body(
    gyro_measurement_bprime: Sequence[Interval],
    estimated_gyro_bias_bprime: Sequence[Interval],
) -> Vec3:
    """Exact shipping inner relation last_gyr_bias_corrected=gyr-b_g_hat."""
    return sub(gyro_measurement_bprime, estimated_gyro_bias_bprime)


def measured_specific_force_body(
    R_true_wb: Mat3,
    a_true_world: Sequence[Interval],
    gravity_mps2: Interval,
    true_accel_bias_bprime: Sequence[Interval],
    accel_forcing_bprime: Sequence[Interval],
) -> Vec3:
    """Physical accelerometer input retained with explicit additive forcing.

    `true_accel_bias_bprime` includes the admitted centered/deterministic bias
    history. `accel_forcing_bprime` is the declared exogenous sensor/model
    forcing coordinate; this map does not turn configured Racc into a hard
    pathwise noise cap.
    """
    if gravity_mps2.lo <= 0.0:
        raise ValueError("positive gravity required")
    g = (I(0.0), I(0.0), gravity_mps2)
    f_true = matvec(R_true_wb, sub(a_true_world, g))
    return add(add(f_true, true_accel_bias_bprime), accel_forcing_bprime)


def measured_gyro_body(
    omega_true_bprime: Sequence[Interval],
    true_gyro_bias_bprime: Sequence[Interval],
    gyro_forcing_bprime: Sequence[Interval],
) -> Vec3:
    """Raw de-heeled gyro measurement with explicit bias/forcing coordinates."""
    return add(add(omega_true_bprime, true_gyro_bias_bprime), gyro_forcing_bprime)


@dataclass(frozen=True)
class JointExecutorSample:
    gyro_measurement_bprime: Vec3
    omega_body_corrected: Vec3
    specific_force_body: Vec3
    f_cog_body: Vec3
    R_wb: Mat3


def materialize_sample(
    *,
    R_true_wb: Mat3,
    R_hat_wb: Mat3,
    a_true_world: Sequence[Interval],
    a_w_hat_world: Sequence[Interval],
    omega_true_bprime: Sequence[Interval],
    true_gyro_bias_bprime: Sequence[Interval],
    estimated_gyro_bias_bprime: Sequence[Interval],
    true_accel_bias_bprime: Sequence[Interval],
    gyro_forcing_bprime: Sequence[Interval],
    accel_forcing_bprime: Sequence[Interval],
    gravity_mps2: Interval,
) -> JointExecutorSample:
    """Materialize the exact typed-kernel coordinates from one joint witness."""
    R_true_wb = _mat3(R_true_wb, "R_true_wb")
    R_hat_wb = _mat3(R_hat_wb, "R_hat_wb")
    gyro_meas = measured_gyro_body(
        omega_true_bprime, true_gyro_bias_bprime, gyro_forcing_bprime
    )
    corrected = corrected_gyro_rate_body(gyro_meas, estimated_gyro_bias_bprime)
    specific = measured_specific_force_body(
        R_true_wb, a_true_world, gravity_mps2,
        true_accel_bias_bprime, accel_forcing_bprime,
    )
    f_hat = nominal_cog_specific_force_body(
        R_hat_wb, a_w_hat_world, gravity_mps2
    )
    return JointExecutorSample(gyro_meas, corrected, specific, f_hat, R_hat_wb)


def _shipping_parity() -> dict[str, bool]:
    text = MEKF.read_text(encoding="utf-8")
    domain = json.loads(DOMAIN.read_text(encoding="utf-8"))
    return {
        "qref_is_world_to_body_prime": "qref stores WORLD->BODY'." in text,
        "corrected_rate_is_raw_minus_estimated_bg": (
            "last_gyr_bias_corrected = gyr - gyro_bias;" in text
        ),
        "nominal_cog_force_is_Rhat_aw_minus_g": (
            "const Vector3 f_cog_b = R_wb() * (aw - g_world);" in text
        ),
        "measurement_prediction_uses_same_Rhat_aw_minus_g": (
            "const Vector3 f_pred = R_wb() * (aw - g_world) + lever + ba_term;" in text
        ),
        "world_gravity_is_positive_down": (
            "const Vector3 g_world(0,0,+gravity_magnitude_);" in text
        ),
        "attitude_world_to_body_propagates_negative_rate_left": (
            "dq(-omega*dt) * q_wb(k)" in text
            and "(-last_gyr_bias_corrected * Ts).eval()" in text
        ),
        "proof_branch_lever_arm_disabled": (
            domain["configured_runtime"]["imu_lever_arm_enabled"] is False
        ),
        "domain_specific_force_is_pre_bias_noise_cog_model": (
            "before accelerometer bias/noise" in
            domain["normal_live"]["specific_force_norm_definition"]
        ),
    }


def _identity3() -> list[list[Interval]]:
    return [[I(1.0 if i == j else 0.0) for j in range(3)] for i in range(3)]


def _point_smoke() -> dict:
    z = (I(0.0), I(0.0), I(0.0))
    R = _identity3()
    s = materialize_sample(
        R_true_wb=R, R_hat_wb=R,
        a_true_world=(I(1.0), I(-2.0), I(0.5)),
        a_w_hat_world=(I(1.0), I(-2.0), I(0.5)),
        omega_true_bprime=(I(0.1), I(-0.2), I(0.3)),
        true_gyro_bias_bprime=(I(0.01), I(0.02), I(-0.01)),
        estimated_gyro_bias_bprime=(I(0.005), I(0.01), I(-0.02)),
        true_accel_bias_bprime=z,
        gyro_forcing_bprime=z, accel_forcing_bprime=z,
        gravity_mps2=I(9.80665),
    )
    return {
        "gyro_measurement": [x.as_list() for x in s.gyro_measurement_bprime],
        "omega_body_corrected": [x.as_list() for x in s.omega_body_corrected],
        "specific_force_body": [x.as_list() for x in s.specific_force_body],
        "f_cog_body": [x.as_list() for x in s.f_cog_body],
        "R_wb": [[x.as_list() for x in row] for row in s.R_wb],
        "measurement_and_nominal_geometry_roles_distinct": True,
    }


def build() -> dict:
    outer = OUTER.build()
    of = OUTER.validate(outer)
    if of:
        raise RuntimeError("correlated outer-enclosure prerequisite failed: " + repr(of))
    parity = _shipping_parity()
    if not all(parity.values()):
        raise RuntimeError("shipping coordinate-map parity failed: " + repr(parity))
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "correlated_outer_enclosure_consumed": True,
        "source_filter_joint_coordinate_map_materialized": True,
        "measurement_coordinates_and_nominal_geometry_distinct": True,
        "raw_gyro_and_corrected_rate_distinct": True,
        "corrected_rate_depends_on_current_estimated_gyro_bias": True,
        "truth_attitude_and_nominal_R_hat_distinct": True,
        "truth_acceleration_and_nominal_a_w_hat_distinct": True,
        "accelerometer_measurement_retains_true_bias_and_forcing": True,
        "gyro_measurement_retains_true_bias_and_forcing": True,
        "configured_Racc_used_as_hard_pathwise_noise_cap": False,
        "configured_Rgyro_used_as_hard_pathwise_noise_cap": False,
        "lever_arm_disabled_on_certified_branch": True,
        "nominal_accelerometer_geometry_identity": "f_hat_body=R_hat_wb*(a_w_hat_world-[0,0,g])",
        "raw_gyro_identity": "gyro_meas=omega_true+b_g_true+n_g",
        "corrected_gyro_identity": "omega_corrected=gyro_meas-b_g_hat",
        "raw_accelerometer_identity": "f_meas=R_true_wb*(a_true_world-[0,0,g])+b_a_true+n_a",
        "shipping_parity": parity,
        "point_smoke": _point_smoke(),
        "same_joint_witness_required_for_all_coordinates": True,
        "independent_truth_and_estimator_geometry_boxes_allowed": False,
        "independent_raw_and_corrected_gyro_boxes_allowed": False,
        "independent_measurement_and_nominal_force_boxes_allowed": False,
        "joint_source_output_map_closed": True,
        "complete_601_sample_provider_materialized_here": False,
        "sensor_forcing_hard_bound_closed_here": False,
        "BIAS0_assembled_sensor_qualification_closed_here": False,
        "P3_delta": P3_DELTA,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "propagate this exact joint coordinate map through every O^601_BRMM source transition and estimator state, retaining BIAS/forcing coordinates and primitive-prefix ancestry; only then may the canonical provider materialize the 601-sample family"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for k in (
        "correlated_outer_enclosure_consumed",
        "source_filter_joint_coordinate_map_materialized",
        "measurement_coordinates_and_nominal_geometry_distinct",
        "raw_gyro_and_corrected_rate_distinct",
        "corrected_rate_depends_on_current_estimated_gyro_bias",
        "truth_attitude_and_nominal_R_hat_distinct",
        "truth_acceleration_and_nominal_a_w_hat_distinct",
        "accelerometer_measurement_retains_true_bias_and_forcing",
        "gyro_measurement_retains_true_bias_and_forcing",
        "lever_arm_disabled_on_certified_branch",
        "same_joint_witness_required_for_all_coordinates",
        "joint_source_output_map_closed",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "configured_Racc_used_as_hard_pathwise_noise_cap",
        "configured_Rgyro_used_as_hard_pathwise_noise_cap",
        "independent_truth_and_estimator_geometry_boxes_allowed",
        "independent_raw_and_corrected_gyro_boxes_allowed",
        "independent_measurement_and_nominal_force_boxes_allowed",
        "complete_601_sample_provider_materialized_here",
        "sensor_forcing_hard_bound_closed_here",
        "BIAS0_assembled_sensor_qualification_closed_here",
        "P4_MOTION_PASS", "P4_PASS", "P5_MAY_START",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    parity = d.get("shipping_parity", {})
    if not parity or not all(v is True for v in parity.values()):
        f.append("shipping parity not closed")
    if d.get("P3_delta") != P3_DELTA:
        f.append("P3 delta changed")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    d = build()
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "joint_map": d["joint_source_output_map_closed"],
        "provider": d["complete_601_sample_provider_materialized_here"],
        "P4": d["P4_PASS"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
