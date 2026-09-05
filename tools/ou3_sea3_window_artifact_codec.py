#!/usr/bin/env python3
"""Strict representation codec for canonical complete-SEA3 window artifacts.

This module performs representation only.  It does not generate a sea family,
prove SEA0, or establish source reachability of a covariance seed.  The
canonical executor may call :func:`parse_window_artifact` only *after*
``ou3_sea3_hard_finite_window_source.validate_artifact`` has accepted the
provider-owned witness.

The codec exists to remove an otherwise ambiguous implementation seam.  Every
certified interval endpoint is serialized explicitly as ``[lo, hi]``.  Raw gyro
measurement and bias-corrected MEKF body rate remain separate.  The periodic
``a_w`` covariance-floor event serializes only a boolean request; a numerical
floor increment is forbidden because it depends on each mode's current
Riccati covariance after prediction.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any, Sequence

from ou3_interval import Interval, IntervalMatrix, symmetric_positive_definite_ldlt
import ou3_sea3_complete_window_execution_kernel as KERNEL
import ou3_sea3_frontend_state_step as FRONTEND
import ou3_sea3_hard_finite_window_source as SEA0
import ou3_sea3_private_mahony_state_step as MAHONY
import ou3_sea3_tuner_scheduler_step as TUNER
import ou3_sea3_wpe_state_step as WPE

SCHEMA = 1
QUALIFICATION = "OU3_SEA3_COMPLETE_WINDOW_ARTIFACT_CODEC_V1"


@dataclass(frozen=True)
class ParsedWindow:
    frontend_entry: FRONTEND.FrontEndState
    P0_H: IntervalMatrix
    P0_A: IntervalMatrix
    samples: tuple[KERNEL.SampleCoordinates, ...]


def _number(x: Any, label: str) -> float:
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        raise ValueError(f"{label} must be a finite number")
    y = float(x)
    if not math.isfinite(y):
        raise ValueError(f"{label} must be finite")
    return y


def interval_from_json(value: Any, label: str = "interval") -> Interval:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{label} must be [lo,hi]")
    lo = _number(value[0], f"{label}.lo")
    hi = _number(value[1], f"{label}.hi")
    if lo > hi:
        raise ValueError(f"{label} has reversed endpoints")
    # The provider already supplies outward endpoints.  Re-widening during
    # transport would create artificial uncertainty and would corrupt exact
    # algebraic zeros.
    return Interval(lo, hi)


def interval_to_json(x: Interval) -> list[float]:
    return [float(x.lo), float(x.hi)]


def vec3_from_json(value: Any, label: str = "vec3") -> tuple[Interval, Interval, Interval]:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"{label} must contain three intervals")
    return tuple(interval_from_json(value[i], f"{label}[{i}]") for i in range(3))  # type: ignore[return-value]


def vec3_to_json(value: Sequence[Interval]) -> list[list[float]]:
    if len(value) != 3:
        raise ValueError("vec3 serialization requires length three")
    return [interval_to_json(value[i]) for i in range(3)]


def matrix_from_json(value: Any, rows: int, cols: int, label: str) -> IntervalMatrix:
    if not isinstance(value, list) or len(value) != rows:
        raise ValueError(f"{label} must have {rows} rows")
    out: IntervalMatrix = []
    for i, row in enumerate(value):
        if not isinstance(row, list) or len(row) != cols:
            raise ValueError(f"{label}[{i}] must have {cols} columns")
        out.append([
            interval_from_json(row[j], f"{label}[{i}][{j}]")
            for j in range(cols)
        ])
    return out


def matrix_to_json(value: Sequence[Sequence[Interval]]) -> list[list[list[float]]]:
    if not value:
        return []
    cols = len(value[0])
    if any(len(row) != cols for row in value):
        raise ValueError("cannot serialize ragged interval matrix")
    return [[interval_to_json(x) for x in row] for row in value]


def _bool(value: Any, label: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{label} must be boolean")
    return value


def _dict(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    return value


def _witness(payload: dict[str, Any], expected: str, label: str) -> None:
    if not isinstance(expected, str) or not expected:
        raise ValueError(f"{label} expected witness id is empty")
    if payload.get("witness_id") != expected:
        raise ValueError(f"{label} witness id mismatch")


def mahony_from_json(value: Any) -> MAHONY.State:
    d = _dict(value, "front_end_entry.mahony")
    names = (
        "q0", "q1", "q2", "q3",
        "integral_x", "integral_y", "integral_z", "up_ms2",
    )
    vals = [interval_from_json(d.get(name), f"mahony.{name}") for name in names]
    return MAHONY.State(*vals)


def mahony_to_json(state: MAHONY.State) -> dict[str, Any]:
    return {
        name: interval_to_json(getattr(state, name))
        for name in (
            "q0", "q1", "q2", "q3",
            "integral_x", "integral_y", "integral_z", "up_ms2",
        )
    }


def wpe_from_json(value: Any) -> WPE.WPEState:
    d = _dict(value, "front_end_entry.wpe")
    interval_names = (
        "accel_prev", "high_pass_1", "high_pass_1_prev", "high_pass_2",
        "velocity", "elevation", "velocity_mean", "velocity_sq",
        "elevation_mean", "elevation_sq", "weight", "elapsed_s",
        "raw_period_s", "log_period_s", "last_moment_horizon_s",
        "last_log_horizon_s",
    )
    vals = {
        name: interval_from_json(d.get(name), f"wpe.{name}")
        for name in interval_names
    }
    usable = _bool(d.get("usable_period"), "wpe.usable_period")
    if not usable:
        raise ValueError("canonical Normal-Live front-end entry must have usable WPE")
    return WPE.WPEState(
        accel_prev=vals["accel_prev"],
        high_pass_1=vals["high_pass_1"],
        high_pass_1_prev=vals["high_pass_1_prev"],
        high_pass_2=vals["high_pass_2"],
        velocity=vals["velocity"], elevation=vals["elevation"],
        velocity_mean=vals["velocity_mean"], velocity_sq=vals["velocity_sq"],
        elevation_mean=vals["elevation_mean"], elevation_sq=vals["elevation_sq"],
        weight=vals["weight"], elapsed_s=vals["elapsed_s"],
        raw_period_s=vals["raw_period_s"], log_period_s=vals["log_period_s"],
        usable_period=usable,
        last_moment_horizon_s=vals["last_moment_horizon_s"],
        last_log_horizon_s=vals["last_log_horizon_s"],
    )


def wpe_to_json(state: WPE.WPEState) -> dict[str, Any]:
    names = (
        "accel_prev", "high_pass_1", "high_pass_1_prev", "high_pass_2",
        "velocity", "elevation", "velocity_mean", "velocity_sq",
        "elevation_mean", "elevation_sq", "weight", "elapsed_s",
        "raw_period_s", "log_period_s", "last_moment_horizon_s",
        "last_log_horizon_s",
    )
    out = {name: interval_to_json(getattr(state, name)) for name in names}
    out["usable_period"] = bool(state.usable_period)
    return out


def tuner_from_json(value: Any) -> TUNER.TunerState:
    d = _dict(value, "front_end_entry.tuner")
    b = _dict(d.get("band"), "tuner.band")
    band = TUNER.BandState(
        interval_from_json(b.get("lowpass_low"), "tuner.band.lowpass_low"),
        interval_from_json(b.get("band"), "tuner.band.band"),
        interval_from_json(b.get("p00"), "tuner.band.p00"),
        interval_from_json(b.get("p01"), "tuner.band.p01"),
        interval_from_json(b.get("p11"), "tuner.band.p11"),
        _bool(b.get("ready"), "tuner.band.ready"),
    )
    m = _dict(d.get("moments"), "tuner.moments")
    moments = TUNER.MomentState(
        interval_from_json(m.get("mean_value"), "tuner.moments.mean_value"),
        interval_from_json(m.get("mean_weight"), "tuner.moments.mean_weight"),
        interval_from_json(m.get("sq_value"), "tuner.moments.sq_value"),
        interval_from_json(m.get("sq_weight"), "tuner.moments.sq_weight"),
        interval_from_json(m.get("frequency_hz"), "tuner.moments.frequency_hz"),
    )
    c = _dict(d.get("candidate"), "tuner.candidate")
    candidate = TUNER.CandidateState(
        interval_from_json(c.get("tau"), "tuner.candidate.tau"),
        interval_from_json(c.get("sigma"), "tuner.candidate.sigma"),
        interval_from_json(c.get("rs"), "tuner.candidate.rs"),
    )
    a = _dict(d.get("active"), "tuner.active")
    active = TUNER.ActiveSchedule(
        interval_from_json(a.get("tau"), "tuner.active.tau"),
        interval_from_json(a.get("sigma"), "tuner.active.sigma"),
        interval_from_json(a.get("rs_base"), "tuner.active.rs_base"),
        interval_from_json(a.get("pseudo_period"), "tuner.active.pseudo_period"),
    )
    s = _dict(d.get("scheduler"), "tuner.scheduler")
    scheduler = TUNER.SchedulerState(
        interval_from_json(
            s.get("since_last_commit_stage_s"),
            "tuner.scheduler.since_last_commit_stage_s",
        ),
        _bool(s.get("pending_commit"), "tuner.scheduler.pending_commit"),
    )
    return TUNER.TunerState(band, moments, candidate, active, scheduler)


def tuner_to_json(state: TUNER.TunerState) -> dict[str, Any]:
    return {
        "band": {
            "lowpass_low": interval_to_json(state.band.lowpass_low),
            "band": interval_to_json(state.band.band),
            "p00": interval_to_json(state.band.p00),
            "p01": interval_to_json(state.band.p01),
            "p11": interval_to_json(state.band.p11),
            "ready": bool(state.band.ready),
        },
        "moments": {
            "mean_value": interval_to_json(state.moments.mean_value),
            "mean_weight": interval_to_json(state.moments.mean_weight),
            "sq_value": interval_to_json(state.moments.sq_value),
            "sq_weight": interval_to_json(state.moments.sq_weight),
            "frequency_hz": interval_to_json(state.moments.frequency_hz),
        },
        "candidate": {
            "tau": interval_to_json(state.candidate.tau),
            "sigma": interval_to_json(state.candidate.sigma),
            "rs": interval_to_json(state.candidate.rs),
        },
        "active": {
            "tau": interval_to_json(state.active.tau),
            "sigma": interval_to_json(state.active.sigma),
            "rs_base": interval_to_json(state.active.rs_base),
            "pseudo_period": interval_to_json(state.active.pseudo_period),
        },
        "scheduler": {
            "since_last_commit_stage_s": interval_to_json(
                state.scheduler.since_last_commit_stage_s
            ),
            "pending_commit": bool(state.scheduler.pending_commit),
        },
    }


def frontend_from_json(value: Any, expected_witness_id: str) -> FRONTEND.FrontEndState:
    d = _dict(value, "front_end_entry")
    _witness(d, expected_witness_id, "front_end_entry")
    return FRONTEND.FrontEndState(
        mahony_from_json(d.get("mahony")),
        wpe_from_json(d.get("wpe")),
        tuner_from_json(d.get("tuner")),
    )


def frontend_to_json(state: FRONTEND.FrontEndState, witness_id: str) -> dict[str, Any]:
    if not witness_id:
        raise ValueError("front-end witness id required")
    return {
        "witness_id": witness_id,
        "mahony": mahony_to_json(state.mahony),
        "wpe": wpe_to_json(state.wpe),
        "tuner": tuner_to_json(state.tuner),
    }


def _symmetric_exact_box(A: Sequence[Sequence[Interval]], label: str) -> None:
    n = len(A)
    for i in range(n):
        for j in range(i):
            if A[i][j] != A[j][i]:
                raise ValueError(f"{label} must serialize one symmetric interval matrix")


def live_seed_from_json(
    value: Any,
    expected_witness_id: str,
) -> tuple[IntervalMatrix, IntervalMatrix]:
    d = _dict(value, "live_covariance_seed")
    _witness(d, expected_witness_id, "live_covariance_seed")
    PH = matrix_from_json(d.get("P0_H_interval"), 18, 18, "live_covariance_seed.P0_H_interval")
    PA = matrix_from_json(d.get("P0_A_interval"), 21, 21, "live_covariance_seed.P0_A_interval")
    _symmetric_exact_box(PH, "P0_H_interval")
    _symmetric_exact_box(PA, "P0_A_interval")
    h_ok, _ = symmetric_positive_definite_ldlt(PH)
    a_ok, _ = symmetric_positive_definite_ldlt(PA)
    if not h_ok or not a_ok:
        raise ValueError("provider Live covariance seed must certify SPD in both H18 and A21")
    return PH, PA


def live_seed_to_json(
    P0_H: Sequence[Sequence[Interval]],
    P0_A: Sequence[Sequence[Interval]],
    witness_id: str,
) -> dict[str, Any]:
    if not witness_id:
        raise ValueError("Live seed witness id required")
    return {
        "witness_id": witness_id,
        "P0_H_interval": matrix_to_json(P0_H),
        "P0_A_interval": matrix_to_json(P0_A),
    }


def sample_from_transition(value: Any) -> KERNEL.SampleCoordinates:
    sample = _dict(value, "transition")
    transition_id = sample.get("source_transition_witness_id")
    response_id = sample.get("joint_response_witness_id")
    physical = _dict(sample.get("joint_physical_output"), "transition.joint_physical_output")
    events = _dict(sample.get("source_events"), "transition.source_events")
    if physical.get("source_transition_witness_id") != transition_id:
        raise ValueError("physical payload transition witness mismatch")
    if physical.get("joint_response_witness_id") != response_id:
        raise ValueError("physical payload response witness mismatch")
    if events.get("source_transition_witness_id") != transition_id:
        raise ValueError("event payload transition witness mismatch")
    if "aw_covariance_floor_increment" in events:
        raise ValueError("precomputed covariance-floor increment is forbidden")

    raw = vec3_from_json(physical.get("gyro_measurement_interval"), "gyro_measurement_interval")
    corrected = vec3_from_json(
        physical.get("omega_body_corrected_interval"),
        "omega_body_corrected_interval",
    )
    specific = vec3_from_json(
        physical.get("specific_force_body_interval"),
        "specific_force_body_interval",
    )
    fcog = vec3_from_json(physical.get("f_cog_body_interval"), "f_cog_body_interval")
    Rwb = matrix_from_json(physical.get("R_wb_interval"), 3, 3, "R_wb_interval")

    mag_payload = events.get("magnetometer_events_after_imu")
    if not isinstance(mag_payload, list):
        raise ValueError("magnetometer_events_after_imu must be a list")
    mags: list[KERNEL.MagneticEvent] = []
    for i, event in enumerate(mag_payload):
        e = _dict(event, f"magnetometer_events_after_imu[{i}]")
        mags.append(KERNEL.MagneticEvent(vec3_from_json(
            e.get("m_body_interval"), f"magnetometer_events_after_imu[{i}].m_body_interval"
        )))

    return KERNEL.SampleCoordinates(
        gyro_measurement=MAHONY.Vec3(*raw),
        omega_body_corrected=corrected,
        specific_force=MAHONY.Vec3(*specific),
        f_cog_body=fcog,
        R_wb=Rwb,
        due_S=_bool(events.get("S_zero_due"), "S_zero_due"),
        aw_floor_requested=_bool(
            events.get("aw_covariance_floor_requested"),
            "aw_covariance_floor_requested",
        ),
        magnetometer_events_after_imu=tuple(mags),
    )


def sample_to_payload(
    sample: KERNEL.SampleCoordinates,
    *,
    transition_witness_id: str,
    joint_response_witness_id: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    raw = (sample.gyro_measurement.x, sample.gyro_measurement.y, sample.gyro_measurement.z)
    specific = (
        sample.specific_force.x,
        sample.specific_force.y,
        sample.specific_force.z,
    )
    physical = {
        "source_transition_witness_id": transition_witness_id,
        "joint_response_witness_id": joint_response_witness_id,
        "gyro_measurement_interval": vec3_to_json(raw),
        "omega_body_corrected_interval": vec3_to_json(sample.omega_body_corrected),
        "specific_force_body_interval": vec3_to_json(specific),
        "f_cog_body_interval": vec3_to_json(sample.f_cog_body),
        "R_wb_interval": matrix_to_json(sample.R_wb),
    }
    events = {
        "source_transition_witness_id": transition_witness_id,
        "magnetometer_events_after_imu": [
            {"m_body_interval": vec3_to_json(e.m_body)}
            for e in sample.magnetometer_events_after_imu
        ],
        "aw_covariance_floor_requested": bool(sample.aw_floor_requested),
        "S_zero_due": bool(sample.due_S),
    }
    return physical, events


def parse_window_artifact(d: dict[str, Any]) -> ParsedWindow:
    structural = SEA0.validate_candidate_structure(d)
    if structural:
        raise ValueError("SEA3 window structural validation failed: " + "; ".join(structural))
    frontend = frontend_from_json(d.get("front_end_entry"), d["front_end_entry_witness_id"])
    PH, PA = live_seed_from_json(
        d.get("live_covariance_seed"), d["live_covariance_seed_witness_id"]
    )
    samples = tuple(sample_from_transition(x) for x in d["transitions"])
    if len(samples) != SEA0.SAMPLES:
        raise ValueError("canonical parsed window must contain 601 samples")
    return ParsedWindow(frontend, PH, PA, samples)


def _diag(n: int, value: float) -> IntervalMatrix:
    z = Interval.point(0.0)
    d = Interval.point(value)
    return [[d if i == j else z for j in range(n)] for i in range(n)]


def build() -> dict[str, Any]:
    front = FRONTEND._point_state()  # codec fixture only; never source evidence
    encoded_front = frontend_to_json(front, "front-fixture")
    decoded_front = frontend_from_json(encoded_front, "front-fixture")

    PH = _diag(18, 2.0)
    PA = _diag(21, 2.0)
    seed_payload = live_seed_to_json(PH, PA, "seed-fixture")
    PH2, PA2 = live_seed_from_json(seed_payload, "seed-fixture")

    z = Interval.point(0.0)
    sample = KERNEL.SampleCoordinates(
        gyro_measurement=MAHONY.Vec3(z, z, z),
        omega_body_corrected=(z, z, z),
        specific_force=MAHONY.Vec3(z, z, Interval.point(-9.8)),
        f_cog_body=(z, z, Interval.point(-9.8)),
        R_wb=[[Interval.point(1.0 if i == j else 0.0) for j in range(3)] for i in range(3)],
        due_S=True,
        aw_floor_requested=False,
        magnetometer_events_after_imu=(
            KERNEL.MagneticEvent((Interval.point(20.0), z, Interval.point(40.0))),
        ),
    )
    physical, events = sample_to_payload(
        sample,
        transition_witness_id="tr-fixture",
        joint_response_witness_id="resp-fixture",
    )
    decoded_sample = sample_from_transition({
        "source_transition_witness_id": "tr-fixture",
        "joint_response_witness_id": "resp-fixture",
        "joint_physical_output": physical,
        "source_events": events,
    })

    smoke = {
        "frontend_roundtrip_exact": decoded_front == front,
        "H_seed_roundtrip_exact": PH2 == PH,
        "A_seed_roundtrip_exact": PA2 == PA,
        "sample_roundtrip_exact": decoded_sample == sample,
    }
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generator": False,
        "trajectory_replay_used": False,
        "establishes_source_reachability": False,
        "requires_provider_acceptance_before_canonical_use": True,
        "re_widens_provider_interval_endpoints": False,
        "precomputed_aw_floor_increment_accepted": False,
        "raw_gyro_and_corrected_rate_kept_distinct": True,
        "live_seed_requires_full_H18_and_A21_interval_SPD": True,
        "strict_codec_ready": all(smoke.values()),
        "smoke": smoke,
        "P3_promoted": False,
        "next_obligation": (
            "after SEA0 provider acceptance, parse its complete witness with this codec and execute it "
            "through the connected typed kernel; this codec does not establish SEA0 or seed reachability"
        ),
    }


def validate(d: dict[str, Any]) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "requires_provider_acceptance_before_canonical_use",
        "raw_gyro_and_corrected_rate_kept_distinct",
        "live_seed_requires_full_H18_and_A21_interval_SPD",
        "strict_codec_ready",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_generator", "trajectory_replay_used", "establishes_source_reachability",
        "re_widens_provider_interval_endpoints", "precomputed_aw_floor_increment_accepted",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if not all(d.get("smoke", {}).values()):
        f.append("codec roundtrip smoke failed")
    return list(dict.fromkeys(f))


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
        "strict_codec_ready": d["strict_codec_ready"],
        "smoke": d["smoke"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
