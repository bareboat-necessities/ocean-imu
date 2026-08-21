#!/usr/bin/env python3
"""Paired Monte Carlo validation for the OU-II and OU-III filters.

The simulator retains deterministic executable regression gates, but they are
not Monte Carlo inclusion criteria and are not exported as statistical-row
fields. Validation uses a configurable long
window (900 seconds in full mode, matching the window the simulators gate on),
independent wave-realization, IMU-noise, and initialization seeds, and paired
realizations across filters and ablations.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

import ou_evidence_provenance as evidence_provenance


REPO_ROOT = Path(__file__).resolve().parents[1]
DT_SECONDS = 1.0 / 200.0
GRAVITY_MPS2 = 9.80665
SURROGATE_MIN_FREQ_HZ = 0.02
SURROGATE_MAX_FREQ_HZ = 1.60
# Trailing window the OU-II/OU-III simulators score their own quality gates
# over; see RMS_WINDOW_SEC in src/util/W3dSimCommon.cpp.  Recorded in the
# protocol block so a raw row's gate flag can be read without the source.
SIMULATOR_GATE_WINDOW_SEC = 900.0
WAVE_PHASE_METHOD = (
    "common random phase on the 0.02-1.60 Hz world-velocity and Euler "
    "spectra; displacement and acceleration analytically derived from "
    "velocity; body IMU and quaternion reconstructed"
)
TRANSITION_METHOD = (
    "C2 quintic crossfade between two independently phase-randomized "
    "stationary records, with exact first- and second-derivative cross terms; "
    "not a continuously evolving JONSWAP spectrum"
)
PMSTOKES_POOLING = (
    "PM-Stokes is a separate declared ensemble and is not pooled into the "
    "primary JONSWAP aggregate"
)
PRIMARY_ENDPOINT_INFERENCE = (
    "the confirmatory endpoint carries four companion statements on the same "
    "paired seed-level differences: the percentile bootstrap interval, a "
    "Student-t interval, an exact sign test, and an exact paired "
    "randomization (sign-flip) test enumerated over all 2^n sign patterns.  "
    "They are companions, not independent confirmations: all four read the "
    "same n paired differences, and none of them widens the ensemble."
)

# Adaptation modes.  The three primary modes all run the deployed covariance
# policy and vary only whether the OU operating point adapts online.  The two
# *HeldCovariance* modes repeat their partner mode with the periodic
# re-alignment of the posterior a_w marginal switched off, so that "parameters
# adapt" and "part of the covariance is periodically re-aligned" are separated
# instead of confounded.
MODE_SETTINGS: dict[str, tuple[str, str]] = {
    # mode -> (simulator tuning mode, a_w covariance-synchronization policy)
    "Adaptive": ("adaptive", "periodic"),
    "FixedNominal": ("fixed_nominal", "periodic"),
    "FixedOracle": ("fixed_oracle", "periodic"),
    "AdaptiveHeldCovariance": ("adaptive", "reconfigure"),
    "FixedNominalHeldCovariance": ("fixed_nominal", "reconfigure"),
    # Partial-adaptation channels.  The deployed law couples the integral
    # regularization scale to the OU parameters through
    # r_S = clip(0.35 sigma_aw tau^3), so Adaptive-versus-FixedNominal cannot
    # say which channel earns the benefit.  These two modes freeze one channel
    # at the FixedNominal point while the other keeps adapting, completing a
    # 2x2 factorial with Adaptive and FixedNominal.  Only OU-III exposes an
    # explicit r_S channel, so they are OU-III only.
    "AdaptiveRSOnly": ("adaptive_rs_only", "periodic"),
    "AdaptiveOUOnly": ("adaptive_ou_only", "periodic"),
}

PRIMARY_MODES = ("Adaptive", "FixedNominal", "FixedOracle")

# Modes that freeze one adaptation channel at the nominal point.
PARTIAL_MODES = ("AdaptiveRSOnly", "AdaptiveOUOnly")

# Modes each family can run.  OU-II regularizes the second-order chain with
# (R_p0, R_v0) rather than a single integral scale, so the channel ablation is
# not defined for it.
FAMILY_MODES: dict[str, tuple[str, ...]] = {
    "OU_II": tuple(m for m in MODE_SETTINGS if m not in PARTIAL_MODES),
    "OU_III": tuple(MODE_SETTINGS),
}

# Ablation pairs compared at matched tuning mode: (left, right).
COVARIANCE_SYNC_PAIRS = (
    ("Adaptive", "AdaptiveHeldCovariance"),
    ("FixedNominal", "FixedNominalHeldCovariance"),
)

# Crossfade placement of the non-stationary record, as fractions of the replay
# duration.  The blend is centred on the record, so the scored window contains a
# run-in at the start sea, the whole crossfade, and a tail at the endpoint sea.
#
# The crossfade spans a tenth of the record -- 120 s of the 1200 s replay.  It
# used to span 0.35-0.65, i.e. 360 s, and that was too slow to be an instrument:
# the slowest adaptation memory in either filter is the two-stage sigma_a EWMA
# at about 2*K_periods*T_z, some 34 s on the largest reference sea, so a 360 s
# ramp is quasi-static from the schedule's point of view and the score cannot
# resolve tracking lag from steady-state accuracy.  At 120 s the crossfade is
# about three times that memory: still a sea state changing rather than a step,
# but fast enough that a horizon which lags shows up as an error.  Shortening it
# from the front also leaves the settled interval (780-1200 s) exactly where it
# was, so settled numbers stay comparable with the previously published ones.
TRANSITION_START_FRACTION = 0.45
TRANSITION_END_FRACTION = 0.55

# Scoring segments for the non-stationary record, in seconds from the start of
# the run.  A single trailing window mixes the pure start sea, the crossfade,
# and the pure endpoint sea; these split the same window so that each interval
# can be read on its own.
#
# "recover" is the interval immediately after the crossfade, one crossfade
# length long, and "end" is the settled remainder.  A schedule whose averaging
# horizon is too long does not only lag during the blend: it carries the old
# sea into the new one, and that cost lands after the crossfade rather than
# inside it.  Scoring the two apart is what separates "slow to follow" from
# "wrong once it has followed"; pooled into one endpoint interval they partly
# cancel, which is what a single "end" segment reported before.
TRANSITION_SEGMENTS = ("start", "blend", "recover", "end")

SEGMENT_METRIC_NAMES = (
    "disp_z_rms_m",
    "disp_z_pct_hs",
    "disp_z_ref_rms_m",
    "disp_z_pct_refrms",
    "disp_3d_rms_m",
)

DIRECTION_METRIC_NAMES = (
    "dir_axis_error_deg",
    "dir_axis_abs_error_deg",
    "dir_axis_rmse_deg",
    "dir_axis_circ_std_deg",
    "dir_sense_forward_pct",
    "dir_sense_reverse_pct",
    "dir_sense_uncertain_pct",
    "dir_sense_dominant_pct",
    # Travel-sense correctness against the record's physical propagation
    # direction.  The dir_sense_* shares above are relative to the estimator's
    # own axis representative and are stability, not accuracy; these are the
    # accuracy metrics.
    "dir_travel_error_deg",
    "dir_travel_abs_error_deg",
    "dir_travel_rmse_deg",
    "dir_travel_correct_pct",
    "dir_travel_wrong_pct",
    "dir_travel_unresolved_pct",
)

METRIC_NAMES = (
    "disp_x_rms_m",
    "disp_y_rms_m",
    "disp_z_rms_m",
    "disp_3d_rms_m",
    "disp_x_pct_hs",
    "disp_y_pct_hs",
    "disp_z_pct_hs",
    "disp_3d_pct_refmax",
    "disp_z_ref_rms_m",
    "disp_z_pct_refrms",
    "roll_rms_deg",
    "pitch_rms_deg",
    "yaw_rms_deg",
    *DIRECTION_METRIC_NAMES,
    "accel_bias_3d_rms_mps2",
    "gyro_bias_3d_rms_radps",
    *(
        f"seg_{segment}_{metric}"
        for segment in TRANSITION_SEGMENTS
        for metric in SEGMENT_METRIC_NAMES
    ),
)

# Metrics defined for every run.  The segment metrics exist only where extra
# scoring intervals were requested, so studies that do not request them use
# this list instead of carrying empty columns.
NON_SEGMENT_METRIC_NAMES = tuple(
    metric for metric in METRIC_NAMES if not metric.startswith("seg_")
)

DISPLAY_METRICS = (
    "disp_3d_rms_m",
    "disp_z_rms_m",
    "roll_rms_deg",
    "pitch_rms_deg",
    "yaw_rms_deg",
)

FAMILY_BINARY = {
    "OU_II": REPO_ROOT / "tests" / "kalman_ou_ii" / "kalman_ou_ii-sim",
    "OU_III": REPO_ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim",
}

# Time-series filename suffix each simulator appends when
# W3D_WRITE_TIMESERIES is enabled and a magnetometer is present.
FAMILY_TIMESERIES_SUFFIX = {
    "OU_II": "_fusion_ou2",
    "OU_III": "_fusion_ou3",
}

FAMILY_MAKE_DIR = {
    "OU_II": REPO_ROOT / "tests" / "kalman_ou_ii",
    "OU_III": REPO_ROOT / "tests" / "kalman_ou_iii",
}

DEFAULT_FULL_WAVE_SEEDS = (11, 29, 47, 71, 97, 131, 173, 211, 257, 307)
DEFAULT_FULL_IMU_SEEDS = (101, 211, 307, 401, 503, 601, 701, 809, 907, 1009)
DEFAULT_FULL_INIT_SEEDS = (1009, 1103, 1201, 1301, 1409, 1511, 1601, 1709, 1801, 1901)


@dataclass(frozen=True)
class SeedTriplet:
    wave_phase_seed: int
    imu_noise_seed: int
    initialization_seed: int


@dataclass(frozen=True)
class Scenario:
    name: str
    kind: str
    start_input: Path
    end_input: Path | None = None
    end_height_m: float | None = None


@dataclass(frozen=True)
class TuningPoint:
    tau_s: float
    sigma_a_mps2: float
    R_p0_std_m: float | None = None
    R_v0_std_mps: float | None = None
    RS_ms: float | None = None


def parse_int_list(text: str) -> list[int]:
    values = [int(token.strip()) for token in text.split(",") if token.strip()]
    if not values or any(value < 0 for value in values):
        raise argparse.ArgumentTypeError("seed lists require unsigned integers")
    return values


def broadcast_seed_triplets(
    wave_seeds: Sequence[int],
    imu_seeds: Sequence[int],
    init_seeds: Sequence[int],
) -> list[SeedTriplet]:
    count = max(len(wave_seeds), len(imu_seeds), len(init_seeds))
    for label, values in (
        ("wave", wave_seeds),
        ("IMU", imu_seeds),
        ("initialization", init_seeds),
    ):
        if len(values) not in (1, count):
            raise ValueError(
                f"{label} seed count must be one or match the longest seed list ({count})"
            )

    def expanded(values: Sequence[int]) -> Sequence[int]:
        return values if len(values) == count else values * count

    return [
        SeedTriplet(wave, imu, init)
        for wave, imu, init in zip(
            expanded(wave_seeds), expanded(imu_seeds), expanded(init_seeds)
        )
    ]


def read_wave_csv(path: Path, duration_sec: float | None = None) -> tuple[list[str], np.ndarray]:
    with path.open("r", encoding="utf-8") as stream:
        columns = stream.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1, dtype=np.float64)
    if data.ndim != 2 or data.shape[1] != len(columns):
        raise ValueError(f"unexpected wave CSV shape for {path}")
    if duration_sec is not None:
        count = min(data.shape[0], max(2, int(round(duration_sec / DT_SECONDS))))
        data = data[:count].copy()
    return columns, data


def write_wave_csv(path: Path, columns: Sequence[str], data: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        path,
        data,
        delimiter=",",
        header=",".join(columns),
        comments="",
        fmt="%.9g",
    )


def _column_indices(columns: Sequence[str], names: Sequence[str]) -> list[int]:
    lookup = {name: index for index, name in enumerate(columns)}
    missing = [name for name in names if name not in lookup]
    if missing:
        raise ValueError(f"wave CSV is missing columns: {', '.join(missing)}")
    return [lookup[name] for name in names]


def _quaternion_multiply(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    aw, ax, ay, az = np.moveaxis(a, -1, 0)
    bw, bx, by, bz = np.moveaxis(b, -1, 0)
    return np.stack(
        (
            aw * bw - ax * bx - ay * by - az * bz,
            aw * bx + ax * bw + ay * bz - az * by,
            aw * by - ax * bz + ay * bw + az * bx,
            aw * bz + ax * by - ay * bx + az * bw,
        ),
        axis=-1,
    )


def _quaternion_rotate(q: np.ndarray, vectors: np.ndarray) -> np.ndarray:
    q_vector = q[:, 1:]
    first_cross = np.cross(q_vector, vectors)
    return vectors + 2.0 * (
        q[:, :1] * first_cross + np.cross(q_vector, first_cross)
    )


def _world_to_body_quaternion_from_euler_deg(
    roll_deg: np.ndarray,
    pitch_deg: np.ndarray,
    yaw_deg: np.ndarray,
) -> np.ndarray:
    roll = np.deg2rad(roll_deg) * 0.5
    pitch = np.deg2rad(pitch_deg) * 0.5
    yaw = np.deg2rad(yaw_deg) * 0.5
    cr, sr = np.cos(roll), np.sin(roll)
    cp, sp = np.cos(pitch), np.sin(pitch)
    cy, sy = np.cos(yaw), np.sin(yaw)

    # Standard roll-pitch-yaw gives BODY->WORLD. The simulation CSV stores
    # q_wb_zu, so take its conjugate to obtain WORLD->BODY.
    q_body_to_world = np.stack(
        (
            cr * cp * cy + sr * sp * sy,
            sr * cp * cy - cr * sp * sy,
            cr * sp * cy + sr * cp * sy,
            cr * cp * sy - sr * sp * cy,
        ),
        axis=-1,
    )
    q_body_to_world /= np.linalg.norm(q_body_to_world, axis=1, keepdims=True)
    q_world_to_body = q_body_to_world.copy()
    q_world_to_body[:, 1:] *= -1.0
    return q_world_to_body


def _gyro_from_world_to_body_quaternion(q: np.ndarray, dt: float) -> np.ndarray:
    q_next_conjugate = q[1:].copy()
    q_next_conjugate[:, 1:] *= -1.0
    relative = _quaternion_multiply(q[:-1], q_next_conjugate)
    relative /= np.linalg.norm(relative, axis=1, keepdims=True)
    negative = relative[:, 0] < 0.0
    relative[negative] *= -1.0

    vector = relative[:, 1:]
    vector_norm = np.linalg.norm(vector, axis=1)
    angle = 2.0 * np.arctan2(vector_norm, np.clip(relative[:, 0], 0.0, 1.0))
    scale = np.empty_like(vector_norm)
    small = vector_norm < 1e-12
    scale[small] = 2.0 / dt
    scale[~small] = angle[~small] / (vector_norm[~small] * dt)
    gyro = vector * scale[:, None]
    return np.vstack((gyro, gyro[-1]))


def rebuild_body_imu(columns: Sequence[str], data: np.ndarray) -> np.ndarray:
    result = data.copy()
    lookup = {name: index for index, name in enumerate(columns)}
    roll = result[:, lookup["roll_deg"]]
    pitch = result[:, lookup["pitch_deg"]]
    yaw = result[:, lookup["yaw_deg"]]
    q = _world_to_body_quaternion_from_euler_deg(roll, pitch, yaw)

    world_acc = result[:, _column_indices(columns, ("acc_x", "acc_y", "acc_z"))]
    gravity_inclusive = world_acc.copy()
    gravity_inclusive[:, 2] += GRAVITY_MPS2
    body_acc = _quaternion_rotate(q, gravity_inclusive)
    body_gyro = _gyro_from_world_to_body_quaternion(q, DT_SECONDS)

    result[:, _column_indices(columns, ("acc_bx", "acc_by", "acc_bz"))] = body_acc
    result[:, _column_indices(columns, ("gyro_x", "gyro_y", "gyro_z"))] = body_gyro
    result[:, _column_indices(
        columns,
        ("q_wb_zu_w", "q_wb_zu_x", "q_wb_zu_y", "q_wb_zu_z"),
    )] = q
    return result


def phase_randomize_wave(
    columns: Sequence[str], data: np.ndarray, seed: int
) -> np.ndarray:
    """Build a band-limited, kinematically closed phase surrogate.

    The released JONSWAP traces contain components from 0.02 to 0.8 Hz and
    second-order sum harmonics up to 1.6 Hz.  Their finite record boundaries
    are not periodic, so independently rotating the DFTs of displacement,
    velocity, and acceleration redistributes boundary leakage and breaks
    ``v = d/dt(p)`` and ``a = d/dt(v)``.

    Instead, this routine treats world velocity as the translational primitive.
    It applies one common phase rotation to the retained velocity and Euler
    bins, derives displacement and acceleration analytically from that same
    velocity spectrum, and finally rebuilds body-frame IMU fields.  Retained
    velocity/attitude auto- and cross-spectra are preserved while the complete
    translational realization is kinematically consistent by construction.
    """

    velocity_indices = _column_indices(columns, ("vel_x", "vel_y", "vel_z"))
    displacement_indices = _column_indices(
        columns, ("disp_x", "disp_y", "disp_z")
    )
    acceleration_indices = _column_indices(columns, ("acc_x", "acc_y", "acc_z"))
    attitude_indices = _column_indices(
        columns, ("roll_deg", "pitch_deg", "yaw_deg")
    )

    count = data.shape[0]
    frequencies = np.fft.rfftfreq(count, DT_SECONDS)
    omega = 2.0 * np.pi * frequencies
    retained = (
        (frequencies >= SURROGATE_MIN_FREQ_HZ)
        & (frequencies <= SURROGATE_MAX_FREQ_HZ)
    )

    rng = np.random.default_rng(seed)
    phase = rng.uniform(-np.pi, np.pi, frequencies.size)
    phase[0] = 0.0
    if count % 2 == 0:
        phase[-1] = 0.0
    rotation = np.exp(1j * phase)

    source_velocity = np.fft.rfft(data[:, velocity_indices], axis=0)
    velocity_spectrum = np.zeros_like(source_velocity)
    # The state represents oscillatory displacement, so exclude the source
    # model's separate mean Stokes-drift velocity from this closed chain.
    velocity_spectrum[retained] = (
        source_velocity[retained] * rotation[retained, None]
    )

    displacement_spectrum = np.zeros_like(velocity_spectrum)
    source_displacement = np.fft.rfft(data[:, displacement_indices], axis=0)
    displacement_spectrum[0] = source_displacement[0]
    displacement_spectrum[retained] = (
        velocity_spectrum[retained] / (1j * omega[retained, None])
    )
    acceleration_spectrum = np.zeros_like(velocity_spectrum)
    acceleration_spectrum[retained] = (
        1j * omega[retained, None] * velocity_spectrum[retained]
    )

    source_attitude = np.fft.rfft(data[:, attitude_indices], axis=0)
    attitude_spectrum = np.zeros_like(source_attitude)
    attitude_spectrum[0] = source_attitude[0]
    attitude_spectrum[retained] = (
        source_attitude[retained] * rotation[retained, None]
    )

    result = data.copy()
    result[:, displacement_indices] = np.fft.irfft(
        displacement_spectrum, n=count, axis=0
    )
    result[:, velocity_indices] = np.fft.irfft(
        velocity_spectrum, n=count, axis=0
    )
    result[:, acceleration_indices] = np.fft.irfft(
        acceleration_spectrum, n=count, axis=0
    )
    result[:, attitude_indices] = np.fft.irfft(
        attitude_spectrum, n=count, axis=0
    )
    return rebuild_body_imu(columns, result)


def scale_wave_motion(
    columns: Sequence[str], data: np.ndarray, scale: float
) -> np.ndarray:
    if not (math.isfinite(scale) and scale > 0.0):
        raise ValueError("wave scale must be positive")
    result = data.copy()
    linear_names = (
        "disp_x", "disp_y", "disp_z",
        "vel_x", "vel_y", "vel_z",
        "acc_x", "acc_y", "acc_z",
    )
    linear = _column_indices(columns, linear_names)
    result[:, linear] *= scale
    for name in ("roll_deg", "pitch_deg", "yaw_deg"):
        index = columns.index(name)
        mean = float(np.mean(result[:, index]))
        result[:, index] = mean + scale * (result[:, index] - mean)
    return rebuild_body_imu(columns, result)


def smoothstep_profile(
    times: np.ndarray, start_sec: float, end_sec: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a C2 quintic blend and its first two time derivatives."""

    if not end_sec > start_sec:
        raise ValueError("transition end must be later than transition start")
    duration = end_sec - start_sec
    x = np.clip((times - start_sec) / duration, 0.0, 1.0)
    weight = x**3 * (10.0 + x * (-15.0 + 6.0 * x))
    first = 30.0 * x**2 * (1.0 - x) ** 2 / duration
    second = 60.0 * x * (2.0 * x**2 - 3.0 * x + 1.0) / duration**2
    return weight, first, second


def smoothstep_weight(
    times: np.ndarray, start_sec: float, end_sec: float
) -> np.ndarray:
    return smoothstep_profile(times, start_sec, end_sec)[0]


def make_nonstationary_wave(
    columns: Sequence[str],
    start_data: np.ndarray,
    end_data: np.ndarray,
    seed: int,
    end_scale: float,
    transition_start_sec: float,
    transition_end_sec: float,
) -> np.ndarray:
    count = min(start_data.shape[0], end_data.shape[0])
    start = phase_randomize_wave(columns, start_data[:count], seed * 2 + 1)
    end = phase_randomize_wave(columns, end_data[:count], seed * 2 + 2)
    end = scale_wave_motion(columns, end, end_scale)
    times = start[:, columns.index("time")]
    weight, weight_rate, weight_acceleration = smoothstep_profile(
        times, transition_start_sec, transition_end_sec
    )

    displacement = _column_indices(columns, ("disp_x", "disp_y", "disp_z"))
    velocity = _column_indices(columns, ("vel_x", "vel_y", "vel_z"))
    acceleration = _column_indices(columns, ("acc_x", "acc_y", "acc_z"))
    attitude = _column_indices(columns, ("roll_deg", "pitch_deg", "yaw_deg"))
    result = start.copy()
    weight_column = weight[:, None]
    weight_rate_column = weight_rate[:, None]
    weight_acceleration_column = weight_acceleration[:, None]
    displacement_delta = end[:, displacement] - start[:, displacement]
    velocity_delta = end[:, velocity] - start[:, velocity]

    result[:, displacement] = (
        (1.0 - weight_column) * start[:, displacement]
        + weight_column * end[:, displacement]
    )
    result[:, velocity] = (
        (1.0 - weight_column) * start[:, velocity]
        + weight_column * end[:, velocity]
        + weight_rate_column * displacement_delta
    )
    result[:, acceleration] = (
        (1.0 - weight_column) * start[:, acceleration]
        + weight_column * end[:, acceleration]
        + 2.0 * weight_rate_column * velocity_delta
        + weight_acceleration_column * displacement_delta
    )
    result[:, attitude] = (
        (1.0 - weight_column) * start[:, attitude]
        + weight_column * end[:, attitude]
    )
    return rebuild_body_imu(columns, result)


TEXT_METRIC_KEYS = ("family", "tuning_mode", "aw_cov_sync", "input", "segment")


def _parse_metric_line(line: str) -> dict[str, Any]:
    metrics: dict[str, Any] = {}
    for token in line.split()[1:]:
        key, value = token.split("=", 1)
        if key == "samples":
            metrics[key] = int(value)
        elif key in TEXT_METRIC_KEYS:
            metrics[key] = value
        else:
            metrics[key] = float(value)
    return metrics


def parse_validation_metrics(output: str) -> dict[str, Any]:
    lines = [line for line in output.splitlines() if line.startswith("VALIDATION_METRICS ")]
    if not lines:
        raise ValueError("simulator did not emit VALIDATION_METRICS")
    metrics = _parse_metric_line(lines[-1])
    metrics.pop("segment", None)
    metrics.pop("start_s", None)

    # The axis is defined mod 180 deg, so the signed error averages toward zero
    # across scenarios of opposite nominal heading.  The magnitude is what the
    # accuracy claim needs, so it is carried alongside the signed value.
    signed = metrics.get("dir_axis_error_deg")
    if isinstance(signed, float):
        metrics["dir_axis_abs_error_deg"] = abs(signed)

    # Same argument for the directed travel angle: scenarios of opposite
    # nominal heading would cancel in the signed mean.
    signed_travel = metrics.get("dir_travel_error_deg")
    if isinstance(signed_travel, float):
        metrics["dir_travel_abs_error_deg"] = abs(signed_travel)

    # Flatten any extra scoring segments into the same row so the existing
    # summary and paired-effect machinery covers them unchanged.
    for line in output.splitlines():
        if not line.startswith("VALIDATION_SEGMENT "):
            continue
        segment = _parse_metric_line(line)
        name = str(segment.get("segment", "")).strip()
        if not name:
            continue
        for metric in SEGMENT_METRIC_NAMES:
            if metric in segment:
                metrics[f"seg_{name}_{metric}"] = segment[metric]
    return metrics


def _fixed_environment(family: str, point: TuningPoint) -> dict[str, str]:
    env = {
        "W3D_FIXED_TAU_S": f"{point.tau_s:.9g}",
        "W3D_FIXED_SIGMA_A": f"{point.sigma_a_mps2:.9g}",
    }
    if family == "OU_II":
        if point.R_p0_std_m is None or point.R_v0_std_mps is None:
            raise ValueError("OU-II tuning point lacks pseudo-measurement noise")
        env["W3D_FIXED_R_P0_STD"] = f"{point.R_p0_std_m:.9g}"
        env["W3D_FIXED_R_V0_STD"] = f"{point.R_v0_std_mps:.9g}"
    else:
        if point.RS_ms is None:
            raise ValueError("OU-III tuning point lacks RS")
        env["W3D_FIXED_RS"] = f"{point.RS_ms:.9g}"
    return env


def run_simulator(
    family: str,
    input_path: Path,
    window_sec: float,
    imu_seed: int,
    initialization_seed: int,
    tuning_mode: str = "adaptive",
    tuning_point: TuningPoint | None = None,
    no_noise: bool = False,
    aw_cov_sync: str = "reconfigure",
    write_timeseries: bool = False,
    segments: Sequence[tuple[str, float, float]] = (),
) -> tuple[dict[str, Any], bool, int]:
    binary = FAMILY_BINARY[family]
    if not binary.exists():
        raise FileNotFoundError(f"missing simulator binary: {binary}")

    env = os.environ.copy()
    for key in tuple(env):
        if key.startswith("W3D_FIXED_"):
            env.pop(key)
    env.update(
        {
            "W3D_IMU_SEED": str(imu_seed),
            "W3D_INIT_SEED": str(initialization_seed),
            "W3D_WRITE_TIMESERIES": "1" if write_timeseries else "0",
            "W3D_VALIDATION_WINDOW_SEC": f"{window_sec:.9g}",
            "W3D_COLLECT_ALL_GATES": "1",
            "W3D_TUNING_MODE": tuning_mode,
            "W3D_AW_COV_SYNC": aw_cov_sync,
            "W3D_VALIDATION_SEGMENTS": ",".join(
                f"{name}:{start:.9g}:{stop:.9g}" for name, start, stop in segments
            ),
        }
    )
    if tuning_point is not None:
        env.update(_fixed_environment(family, tuning_point))

    command = [str(binary)]
    if no_noise:
        command.append("--no-noise")
    command.extend(("--input", str(input_path.resolve())))
    completed = subprocess.run(
        command,
        cwd=binary.parent,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    try:
        metrics = parse_validation_metrics(completed.stdout)
    except ValueError as error:
        tail = "\n".join(completed.stdout.splitlines()[-30:])
        raise RuntimeError(
            f"{family} simulator failed with exit {completed.returncode}:\n{tail}"
        ) from error

    gate_matches = re.findall(r"QUALITY_GATE: PASS=([01])", completed.stdout)
    quality_gate_pass = bool(gate_matches) and gate_matches[-1] == "1"
    if completed.returncode not in (0, 1):
        tail = "\n".join(completed.stdout.splitlines()[-30:])
        raise RuntimeError(
            f"{family} simulator failed with exit {completed.returncode}:\n{tail}"
        )
    return metrics, quality_gate_pass, completed.returncode


# The regularization laws and their clamps, mirroring
# SeaStateFusionFilter_OU_III.h and SeaStateFusionFilter_OU_II.h.  The
# fixed-tuning modes derive their frozen operating point here rather than
# reading it out of the filter, so a mismatch silently scores every fixed mode
# at a point the deployed filter would never choose.  Keep these in step with
# the C++ defaults.
# The deployed OU-III law is SpectralMSE:
#     r_S = C_J q_eff^(1/14) sigma_a,B^(6/7) tau^(24/7) / sqrt(T_S),
# returned as the filter input, so no cadence renormalization is applied.
# OU_III_RS_COEFF below is the Cubic C_R, kept because the ablation family uses
# it; the fixed-tuning modes derive their frozen r_S from the deployed law.
OU_III_RS_MSE_COEFF = 0.0538
OU_III_RS_QEFF = 2.0 * (0.0148 ** 2) * (1.0 / 200.0)
OU_III_PSEUDO_TAU_RATIO = 0.015 / 1.1
OU_III_PSEUDO_PERIOD_BOUNDS_S = (1.0 / 200.0, 0.25)
OU_III_SIGMA_COEFF = 0.9
OU_III_RS_COEFF = 17.112
# sqrt(R_a): the accelerometer measurement-noise standard deviation the OU-III
# base schedule uses as its acceleration scale.  Mirrors
# R_S_ACCEL_NOISE_DENSITY_DEFAULT / FREQ_SMOOTHER_DT.
OU_III_RS_ACCEL_SIGMA_MPS2 = 0.0148
OU_III_RS_BOUNDS_MS = (0.15, 400.0)
# The deployed OU-II law is PhysicalMSE:
#     r_p = C_P q_eff^(1/10) sigma_a,B^(4/5) tau^(12/5) / sqrt(T_S),
#     r_v = r_p / (C_P/C_V * tau),
# both returned as filter inputs, so no cadence renormalization is applied.
# OU_II_RP0_COEFF/OU_II_RV0_COEFF below are the Empirical c_p and c_v, kept
# because the law ablation uses them; the fixed-tuning modes derive their
# frozen pair from the deployed law.
OU_II_PSEUDO_MSE_COEFF = 0.1116
OU_II_PSEUDO_MSE_RATIO = 0.4611
OU_II_PSEUDO_QEFF = 2.0 * (0.12 ** 2) * (1.0 / 200.0)
OU_II_PSEUDO_TAU_RATIO = 0.015 / 1.1
OU_II_PSEUDO_PERIOD_BOUNDS_S = (1.0 / 200.0, 0.25)
OU_II_SIGMA_COEFF = 0.85
OU_II_RP0_COEFF = 0.65
OU_II_RP0_BOUNDS_M = (0.05, 150.0)
OU_II_RV0_COEFF = 1.3
OU_II_RV0_BOUNDS_MPS = (0.01, 40.0)


def tuning_point_from_pilot(family: str, metrics: Mapping[str, Any]) -> TuningPoint:
    tau = float(metrics["tau_applied_s"])
    sigma = float(metrics["sigma_applied_mps2"])
    if family == "OU_II":
        # PhysicalMSE, the deployed law.  sigma is the OU prior sigma_aw, so
        # divide c_sigma back out to recover the physical band RMS the
        # distortion penalty depends on, exactly as OU-III does below.  See
        # docs/ou-ii-dual-mse-adaptation.md.
        T_S = min(max(OU_II_PSEUDO_TAU_RATIO * tau,
                      OU_II_PSEUDO_PERIOD_BOUNDS_S[0]),
                  OU_II_PSEUDO_PERIOD_BOUNDS_S[1])
        sigma_aB = max(sigma / OU_II_SIGMA_COEFF, 1e-6)
        r_p = (OU_II_PSEUDO_MSE_COEFF
               * OU_II_PSEUDO_QEFF ** 0.1
               * sigma_aB ** 0.8
               * tau ** 2.4
               / math.sqrt(T_S))
        R_p0 = min(max(r_p, OU_II_RP0_BOUNDS_M[0]), OU_II_RP0_BOUNDS_M[1])
        R_v0 = min(max(r_p / (OU_II_PSEUDO_MSE_RATIO * tau),
                       OU_II_RV0_BOUNDS_MPS[0]), OU_II_RV0_BOUNDS_MPS[1])
        return TuningPoint(tau, sigma, R_p0_std_m=R_p0, R_v0_std_mps=R_v0)
    # SpectralMSE, the deployed law.  sigma here is the OU prior sigma_aw, so
    # divide c_sigma back out to recover the physical band RMS the distortion
    # penalty depends on.  See docs/ou-iii-rs-amplitude-retune.md.
    T_S = min(max(OU_III_PSEUDO_TAU_RATIO * tau, OU_III_PSEUDO_PERIOD_BOUNDS_S[0]),
              OU_III_PSEUDO_PERIOD_BOUNDS_S[1])
    sigma_aB = max(sigma / OU_III_SIGMA_COEFF, 1e-6)
    RS = min(max(OU_III_RS_MSE_COEFF
                 * OU_III_RS_QEFF ** (1.0 / 14.0)
                 * sigma_aB ** (6.0 / 7.0)
                 * tau ** (24.0 / 7.0)
                 / math.sqrt(T_S),
                 OU_III_RS_BOUNDS_MS[0]), OU_III_RS_BOUNDS_MS[1])
    return TuningPoint(tau, sigma, RS_ms=RS)


def calibrate_tuning_point(
    family: str,
    input_path: Path,
    window_sec: float,
) -> TuningPoint:
    metrics, _, _ = run_simulator(
        family,
        input_path,
        window_sec,
        imu_seed=0,
        initialization_seed=0,
        no_noise=True,
    )
    return tuning_point_from_pilot(family, metrics)


def _finite_values(rows: Sequence[Mapping[str, Any]], metric: str) -> np.ndarray:
    # Segment metrics only exist for the non-stationary scenario, so a missing
    # key is an expected absence rather than a defect.
    values = np.asarray(
        [float(row.get(metric, math.nan)) for row in rows], dtype=np.float64
    )
    return values[np.isfinite(values)]


def _regularized_incomplete_beta(a: float, b: float, x: float) -> float:
    """I_x(a,b) by the Lentz continued fraction, to double precision.

    Written out rather than imported so that the inference the manuscript
    quotes depends only on NumPy, which is what the validation workflow
    installs.  Only the Student-t tail needs it.
    """

    if not 0.0 <= x <= 1.0:
        raise ValueError("regularized incomplete beta requires 0 <= x <= 1")
    if x in (0.0, 1.0):
        return x
    # The continued fraction converges rapidly only on one side of the
    # symmetry point; reflect onto that side when necessary.
    if x > (a + 1.0) / (a + b + 2.0):
        return 1.0 - _regularized_incomplete_beta(b, a, 1.0 - x)

    log_front = (
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log1p(-x)
    )

    tiny = 1e-300
    c = 1.0
    d = 1.0 - (a + b) * x / (a + 1.0)
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    fraction = d
    for iteration in range(1, 300):
        even = 2 * iteration
        numerator = (
            iteration * (b - iteration) * x / ((a + even - 1.0) * (a + even))
        )
        d = 1.0 + numerator * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + numerator / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        fraction *= d * c

        numerator = (
            -(a + iteration) * (a + b + iteration) * x
            / ((a + even) * (a + even + 1.0))
        )
        d = 1.0 + numerator * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + numerator / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        fraction *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return math.exp(log_front) * fraction / a


def student_t_two_sided_p(t_statistic: float, degrees_of_freedom: int) -> float:
    """Two-sided Student-t tail probability."""

    if degrees_of_freedom < 1 or not math.isfinite(t_statistic):
        return math.nan
    if t_statistic == 0.0:
        return 1.0
    x = degrees_of_freedom / (degrees_of_freedom + t_statistic * t_statistic)
    return float(
        min(1.0, _regularized_incomplete_beta(0.5 * degrees_of_freedom, 0.5, x))
    )


def student_t_quantile(probability: float, degrees_of_freedom: int) -> float:
    """Inverse Student-t CDF by bisection on the two-sided tail."""

    if degrees_of_freedom < 1 or not 0.0 < probability < 1.0:
        return math.nan
    if probability == 0.5:
        return 0.0
    upper_tail = 1.0 - probability if probability > 0.5 else probability
    target = 2.0 * upper_tail
    low, high = 0.0, 1.0
    while student_t_two_sided_p(high, degrees_of_freedom) > target:
        high *= 2.0
        if high > 1e12:
            break
    for _ in range(200):
        middle = 0.5 * (low + high)
        if student_t_two_sided_p(middle, degrees_of_freedom) > target:
            low = middle
        else:
            high = middle
    quantile = 0.5 * (low + high)
    return quantile if probability > 0.5 else -quantile


# The exhaustive sign-flip enumeration materializes a 2^n x n sign matrix, so
# its cost doubles with every added pair: n = 18 is about 37 MB and n = 22 is
# about 740 MB.  Past the bound the test is taken by Monte Carlo instead.  The
# study runs at n = 10, well under it, so its randomization p-value is exact.
MAX_EXACT_SIGN_FLIP_PAIRS = 18
SIGN_FLIP_RESAMPLES = 200_000


def paired_inference(
    differences: np.ndarray, rng: np.random.Generator
) -> dict[str, Any]:
    """Complementary inference for one paired contrast.

    The percentile bootstrap of the paired mean is the study's interval, but
    with ten pairs it is a small-sample interval whose coverage is not
    guaranteed.  Three companions are reported alongside it:

      * a Student-t interval on the paired differences, which is exact under
        normality and is the conventional small-sample choice; and
      * two assumption-light exact tests -- the sign test, which uses only the
        direction of each difference, and the paired randomization
        (sign-flip) test, which is exact under the sharp null of exchangeable
        signs and is enumerated in full at this sample size.

    Agreement across the four is the actual claim: no single one of them turns
    ten seeds into a large sample.
    """

    finite = np.asarray(differences, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    n = int(finite.size)
    result: dict[str, Any] = {"n_pairs": n}
    if n < 2:
        return result

    mean = float(np.mean(finite))
    std = float(np.std(finite, ddof=1))
    standard_error = std / math.sqrt(n)
    degrees_of_freedom = n - 1
    critical = student_t_quantile(0.975, degrees_of_freedom)
    t_statistic = mean / standard_error if standard_error > 0.0 else math.inf
    result.update(
        {
            "t_ci95_low": mean - critical * standard_error,
            "t_ci95_high": mean + critical * standard_error,
            "t_statistic": t_statistic,
            "t_p_value": student_t_two_sided_p(t_statistic, degrees_of_freedom),
            "t_critical_value": critical,
            "standard_error": standard_error,
        }
    )

    # Exact sign test.  Ties are conservatively counted against the observed
    # direction rather than discarded.
    negative = int(np.sum(finite < 0.0))
    positive = int(np.sum(finite > 0.0))
    zero = n - negative - positive
    extreme = min(negative, positive) + zero
    sign_p = min(
        1.0,
        2.0
        * sum(math.comb(n, k) for k in range(0, extreme + 1))
        / float(2**n),
    )
    result.update(
        {
            "sign_negative": negative,
            "sign_positive": positive,
            "sign_zero": zero,
            "sign_p_value": sign_p,
        }
    )

    # Paired randomization test on the mean, under the sharp null that each
    # difference is equally likely to have carried either sign.
    observed = abs(mean)
    if n <= MAX_EXACT_SIGN_FLIP_PAIRS:
        signs = 1.0 - 2.0 * (
            (np.arange(2**n, dtype=np.int64)[:, None] >> np.arange(n)) & 1
        ).astype(np.float64)
        means = np.abs(signs @ finite) / n
        count = int(np.sum(means >= observed - 1e-15))
        result.update(
            {
                "randomization_p_value": count / float(2**n),
                "randomization_exact": True,
                "randomization_patterns": int(2**n),
            }
        )
    else:
        draws = rng.integers(0, 2, size=(SIGN_FLIP_RESAMPLES, n)) * 2.0 - 1.0
        means = np.abs(draws @ finite) / n
        count = int(np.sum(means >= observed - 1e-15))
        result.update(
            {
                # +1/+1 keeps the Monte Carlo p-value strictly positive, which
                # is the standard correction for an estimated permutation tail.
                "randomization_p_value": (count + 1) / float(SIGN_FLIP_RESAMPLES + 1),
                "randomization_exact": False,
                "randomization_patterns": int(SIGN_FLIP_RESAMPLES),
            }
        )
    return result


# Math-mode filler for a quantity the ensemble cannot support.  Used for both
# p-values and the statistics beside them so that "this run has no number
# here" reads the same way wherever it appears.
UNDEFINED_STATISTIC = r"\text{--}"


def _latex_statistic(value: Any, spec: str) -> str:
    """Format one statistic, or the undefined marker when there is none.

    Not every run scores an ensemble the confirmatory statistics are defined
    over.  A smoke run scores a single seed triplet, and a paired interval
    over one pair is not a wide interval but no interval at all, so
    `paired_inference` reports the pair count and stops rather than inventing
    companions.  The macro file still has to define every name the manuscript
    reads, so a statistic the run cannot support is written as the undefined
    marker here instead of failing the write partway through.

    Values are looked up with `.get`, so a key `paired_inference` omitted and
    a key it filled with a non-finite float land in the same place.
    """

    if value is None:
        return UNDEFINED_STATISTIC
    if isinstance(value, float) and not math.isfinite(value):
        return UNDEFINED_STATISTIC
    return format(value, spec)


def _latex_p_value(value: float | None) -> str:
    """Format a p-value as math-mode *content*, without delimiters.

    The manuscript quotes these inside expressions such as `$p=\\macro$`, so a
    macro carrying its own `$` would close that expression and typeset the
    rest as text.  `<` is likewise emitted bare, which is correct in math mode
    and wrong in text mode -- these macros are only ever used in math mode.
    """

    if value is None or not math.isfinite(value):
        return UNDEFINED_STATISTIC
    if value < 1e-4:
        return r"<\!10^{-4}"
    digits = f"{value:.5f}" if value < 1e-3 else f"{value:.4f}"
    return digits.rstrip("0").rstrip(".")


def _bootstrap_mean_ci(
    values: np.ndarray, resamples: int, rng: np.random.Generator
) -> tuple[float, float]:
    if values.size == 0:
        return math.nan, math.nan
    if values.size == 1:
        value = float(values[0])
        return value, value
    indices = rng.integers(0, values.size, size=(resamples, values.size))
    means = values[indices].mean(axis=1)
    low, high = np.quantile(means, (0.025, 0.975))
    return float(low), float(high)


def summarize_rows(
    rows: Sequence[Mapping[str, Any]],
    bootstrap_resamples: int,
    stats_seed: int,
) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str, str], list[Mapping[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[(str(row["scenario"]), str(row["family"]), str(row["mode"]))].append(row)

    rng = np.random.default_rng(stats_seed)
    summary: list[dict[str, Any]] = []
    for (scenario, family, mode), group in sorted(groups.items()):
        for metric in METRIC_NAMES:
            values = _finite_values(group, metric)
            n = int(values.size)
            mean = float(np.mean(values)) if n else math.nan
            std = float(np.std(values, ddof=1)) if n > 1 else 0.0 if n == 1 else math.nan
            half_width = 1.96 * std / math.sqrt(n) if n > 1 else 0.0 if n == 1 else math.nan
            bootstrap_low, bootstrap_high = _bootstrap_mean_ci(
                values, bootstrap_resamples, rng
            )
            summary.append(
                {
                    "scenario": scenario,
                    "family": family,
                    "mode": mode,
                    "metric": metric,
                    "n": n,
                    "mean": mean,
                    "std": std,
                    "ci95_low": mean - half_width,
                    "ci95_high": mean + half_width,
                    "bootstrap_ci95_low": bootstrap_low,
                    "bootstrap_ci95_high": bootstrap_high,
                }
            )
    return summary


def _paired_effect(
    left_rows: Sequence[Mapping[str, Any]],
    right_rows: Sequence[Mapping[str, Any]],
    metric: str,
    pair_keys: Sequence[str],
    bootstrap_resamples: int,
    rng: np.random.Generator,
) -> dict[str, Any] | None:
    def indexed(rows: Sequence[Mapping[str, Any]]) -> dict[tuple[Any, ...], float]:
        return {
            tuple(row[key] for key in pair_keys): float(row[metric])
            for row in rows
            if math.isfinite(float(row.get(metric, math.nan)))
        }

    left = indexed(left_rows)
    right = indexed(right_rows)
    keys = sorted(set(left).intersection(right))
    if not keys:
        return None
    differences = np.asarray([left[key] - right[key] for key in keys], dtype=np.float64)
    n = differences.size
    mean = float(np.mean(differences))
    std = float(np.std(differences, ddof=1)) if n > 1 else 0.0
    cohen_dz = mean / std if n > 1 and std > 0.0 else math.nan
    degrees_of_freedom = n - 1
    correction = (
        math.gamma(degrees_of_freedom / 2.0)
        / (
            math.sqrt(degrees_of_freedom / 2.0)
            * math.gamma((degrees_of_freedom - 1.0) / 2.0)
        )
        if degrees_of_freedom > 1
        else math.nan
    )
    hedges_gz = correction * cohen_dz if math.isfinite(cohen_dz) else math.nan
    bootstrap_low, bootstrap_high = _bootstrap_mean_ci(
        differences, bootstrap_resamples, rng
    )
    return {
        "n_pairs": int(n),
        "mean_paired_difference": mean,
        "std_paired_difference": std,
        "bootstrap_ci95_low": bootstrap_low,
        "bootstrap_ci95_high": bootstrap_high,
        "cohen_dz": cohen_dz,
        "hedges_gz": hedges_gz,
    }


def paired_effect_rows(
    rows: Sequence[Mapping[str, Any]],
    bootstrap_resamples: int,
    stats_seed: int,
) -> list[dict[str, Any]]:
    effects: list[dict[str, Any]] = []
    rng = np.random.default_rng(stats_seed + 1)
    pair_keys = ("scenario", "wave_phase_seed", "imu_noise_seed", "initialization_seed")
    scenarios = sorted({str(row["scenario"]) for row in rows})

    for scenario in scenarios:
        adaptive = [
            row for row in rows
            if row["scenario"] == scenario and row["mode"] == "Adaptive"
        ]
        ou_ii = [row for row in adaptive if row["family"] == "OU_II"]
        ou_iii = [row for row in adaptive if row["family"] == "OU_III"]
        for metric in METRIC_NAMES:
            effect = _paired_effect(
                ou_iii, ou_ii, metric, pair_keys,
                bootstrap_resamples, rng,
            )
            if effect:
                effects.append(
                    {
                        "scenario": scenario,
                        "family": "OU_III_vs_OU_II",
                        "comparison": "OU_III_minus_OU_II",
                        "metric": metric,
                        "left": "OU_III/Adaptive",
                        "right": "OU_II/Adaptive",
                        **effect,
                    }
                )

        for family in sorted({str(row["family"]) for row in adaptive}):
            family_adaptive = [row for row in adaptive if row["family"] == family]

            def mode_rows(mode: str, family: str = family) -> list[Mapping[str, Any]]:
                return [
                    row for row in rows
                    if row["scenario"] == scenario
                    and row["family"] == family
                    and row["mode"] == mode
                ]

            for baseline in ("FixedNominal", "FixedOracle"):
                baseline_rows = mode_rows(baseline)
                for metric in METRIC_NAMES:
                    effect = _paired_effect(
                        family_adaptive, baseline_rows, metric, pair_keys,
                        bootstrap_resamples, rng,
                    )
                    if effect:
                        effects.append(
                            {
                                "scenario": scenario,
                                "family": family,
                                "comparison": f"Adaptive_minus_{baseline}",
                                "metric": metric,
                                "left": f"{family}/Adaptive",
                                "right": f"{family}/{baseline}",
                                **effect,
                            }
                        )

            # Channel ablation: each partially adapting mode against the
            # fully frozen point, so the benefit can be attributed to the
            # channel that was left free rather than to adaptation at large.
            for partial in PARTIAL_MODES:
                partial_rows = mode_rows(partial)
                baseline_rows = mode_rows("FixedNominal")
                for metric in METRIC_NAMES:
                    effect = _paired_effect(
                        partial_rows, baseline_rows, metric, pair_keys,
                        bootstrap_resamples, rng,
                    )
                    if effect:
                        effects.append(
                            {
                                "scenario": scenario,
                                "family": family,
                                "comparison": f"{partial}_minus_FixedNominal",
                                "metric": metric,
                                "left": f"{family}/{partial}",
                                "right": f"{family}/FixedNominal",
                                **effect,
                            }
                        )

            # Covariance-reset ablation at matched tuning mode: isolates the
            # effect of repeatedly overwriting the posterior a_w marginal from
            # the effect of letting the operating point adapt.
            for left_mode, right_mode in COVARIANCE_SYNC_PAIRS:
                left_rows = mode_rows(left_mode)
                right_rows = mode_rows(right_mode)
                for metric in METRIC_NAMES:
                    effect = _paired_effect(
                        left_rows, right_rows, metric, pair_keys,
                        bootstrap_resamples, rng,
                    )
                    if effect:
                        effects.append(
                            {
                                "scenario": scenario,
                                "family": family,
                                "comparison": f"{left_mode}_minus_{right_mode}",
                                "metric": metric,
                                "left": f"{family}/{left_mode}",
                                "right": f"{family}/{right_mode}",
                                **effect,
                            }
                        )
    return effects


# Shard files are the only intermediate between a sharded regeneration and the
# bundle it has to reproduce byte for byte, so they are written with plain json
# rather than write_json: the latter maps non-finite floats to null, and several
# metrics are legitimately NaN.  json's repr-based float formatting and its
# NaN/Infinity literals both round-trip exactly, and CSV would not round-trip
# types at all -- every value would come back a string and land in the published
# JSON as one.
SHARD_ORDER_KEY = "_order"
SHARD_PREFIX = "ou_validation_shard"


def shard_path(shard_dir: Path, prefix: str, index: int) -> Path:
    return shard_dir / f"{prefix}_{index:02d}.json"


def write_shard(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as stream:
        json.dump(list(rows), stream, allow_nan=True)


def read_shards(
    shard_dir: Path, prefix: str, *, ordered: bool
) -> list[dict[str, Any]]:
    """Merge shard files back into the rows a single process would have held.

    A study whose rows are already sorted into a canonical order downstream
    passes ordered=False and is merged by concatenation.  One that relies on
    append order passes ordered=True, and every row then carries the position it
    would have had, so the merge is a sort that does not depend on which shard
    produced what.

    Either way a missing or duplicated shard is a hard error.  The bundle is
    published evidence, and a short one still looks perfectly well formed.
    """
    paths = sorted(shard_dir.glob(f"{prefix}_*.json"))
    if not paths:
        raise FileNotFoundError(f"no shard files under {shard_dir}")

    rows: list[dict[str, Any]] = []
    for path in paths:
        with path.open(encoding="utf-8") as stream:
            rows.extend(json.load(stream))

    if not ordered:
        return rows

    positions = [row[SHARD_ORDER_KEY] for row in rows]
    if len(set(positions)) != len(positions):
        raise ValueError(f"duplicate rows across shards in {shard_dir}")
    if sorted(positions) != list(range(len(positions))):
        missing = sorted(set(range(max(positions) + 1)) - set(positions))
        raise ValueError(
            f"shards in {shard_dir} are not a complete set; "
            f"missing positions {missing[:10]}"
        )

    rows.sort(key=lambda row: row[SHARD_ORDER_KEY])
    for row in rows:
        del row[SHARD_ORDER_KEY]
    return rows


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, restval="")
        writer.writeheader()
        writer.writerows(rows)


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        numeric = float(value)
        return numeric if math.isfinite(numeric) else None
    return value


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(_json_safe(value), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def write_latex_table(path: Path, summary: Sequence[Mapping[str, Any]]) -> None:
    selected = [row for row in summary if row["metric"] in DISPLAY_METRICS]
    indexed = {
        (row["scenario"], row["family"], row["mode"], row["metric"]): row
        for row in selected
    }
    groups = sorted({(row["scenario"], row["family"], row["mode"]) for row in selected})

    def escaped(text: str) -> str:
        return text.replace("_", r"\_")

    def formatted(row: Mapping[str, Any] | None) -> str:
        if row is None:
            return "--"
        return f"{float(row['mean']):.3g} $\\pm$ {float(row['std']):.2g}"

    lines = [
        r"% Generated by tools/ou_validation.py; do not edit by hand.",
        r"\begin{tabular}{lllrrrrr}",
        r"\toprule",
        r"Scenario & Family & Mode & 3D disp. (m) & Z disp. (m) & Roll ($^\circ$) & Pitch ($^\circ$) & Yaw ($^\circ$) \\",
        r"\midrule",
    ]
    for scenario, family, mode in groups:
        values = [
            formatted(indexed.get((scenario, family, mode, metric)))
            for metric in DISPLAY_METRICS
        ]
        lines.append(
            " & ".join((escaped(scenario), escaped(family), escaped(mode), *values))
            + r" \\"
        )
    lines.extend((r"\bottomrule", r"\end{tabular}"))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def scenario_spectrum(scenario: str) -> str:
    if scenario.startswith("nonstationary_"):
        return "jonswap"
    if "pmstokes" in scenario:
        return "pmstokes"
    return "jonswap"


def scenario_display_label(scenario: str) -> str:
    if scenario.startswith("nonstationary_"):
        return "Transition"
    height = re.search(r"_H([0-9]+)_([0-9]{3})_", scenario)
    if height:
        value = float(f"{height.group(1)}.{height.group(2)}")
        return rf"$H_s={value:.2f}$ m"
    return scenario.replace("_", " ")


def scenario_sort_key(scenario: str) -> tuple[int, int, float, str]:
    spectrum_rank = 1 if scenario_spectrum(scenario) == "pmstokes" else 0
    if scenario.startswith("nonstationary_"):
        return spectrum_rank, 1, math.inf, scenario
    height = re.search(r"_H([0-9]+)_([0-9]{3})_", scenario)
    value = (
        float(f"{height.group(1)}.{height.group(2)}")
        if height else math.inf
    )
    return spectrum_rank, 0, value, scenario


# TeX control-sequence names are made of category-11 characters, i.e. letters
# only.  A generated name that contains a digit does not fail: TeX ends the
# name at the digit and typesets the remainder as ordinary text beside the
# value, so `\OUValidationIIHs027NominalDifference` renders as the literal
# "27NominalDifference-13.118".  Every generated name is therefore mapped into
# [A-Za-z] here, and checked again before the file is written.
_DIGIT_WORDS = {
    "0": "Zero",
    "1": "One",
    "2": "Two",
    "3": "Three",
    "4": "Four",
    "5": "Five",
    "6": "Six",
    "7": "Seven",
    "8": "Eight",
    "9": "Nine",
}

_PROVIDED_MACRO_RE = re.compile(r"\\(?:provide|new|renew)command\s*\{\\([^}]*)\}")


def latex_letter_token(text: str) -> str:
    """Return ``text`` as a legal TeX control-sequence fragment."""

    return "".join(
        character if character.isalpha() else _DIGIT_WORDS.get(character, "")
        for character in text
    )


def scenario_macro_token(scenario: str) -> str:
    """Letters-only macro fragment identifying a scenario."""

    if scenario.startswith("nonstationary_"):
        return "Transition"
    prefix = "PMStokes" if scenario_spectrum(scenario) == "pmstokes" else ""
    height = re.search(r"_H([0-9]+)_([0-9]{3})_", scenario)
    if not height:
        return prefix + latex_letter_token(scenario.title())
    value = f"{height.group(1)}.{height.group(2)}".rstrip("0").rstrip(".")
    return f"{prefix}Hs{latex_letter_token(value)}"


def assert_latex_macro_names(lines: Iterable[str]) -> None:
    """Fail generation rather than emit a macro TeX cannot name."""

    illegal = sorted(
        {
            name
            for line in lines
            for name in _PROVIDED_MACRO_RE.findall(line)
            if not name.isalpha()
        }
    )
    if illegal:
        raise ValueError(
            "generated LaTeX macro names must contain letters only; "
            f"offending names: {', '.join(illegal)}"
        )


def stationary_normalized_aggregate(
    rows: Sequence[Mapping[str, Any]],
    bootstrap_resamples: int,
    stats_seed: int,
    spectrum: str = "jonswap",
) -> dict[str, dict[str, Any]]:
    """Seed-level mean normalized vertical error over one stationary family.

    The primary endpoint is the JONSWAP ensemble the study was declared on.
    PM-Stokes is aggregated separately rather than pooled, so that adding it
    does not silently redefine the confirmatory comparison.
    """

    by_seed_family: dict[
        tuple[Any, Any, Any, str], list[float]
    ] = defaultdict(list)
    for row in rows:
        if (
            not str(row["scenario"]).startswith("stationary_")
            or row["mode"] != "Adaptive"
            or scenario_spectrum(str(row["scenario"])) != spectrum
        ):
            continue
        key = (
            row["wave_phase_seed"],
            row["imu_noise_seed"],
            row["initialization_seed"],
            str(row["family"]),
        )
        value = float(row["disp_z_pct_hs"])
        if math.isfinite(value):
            by_seed_family[key].append(value)

    seed_keys = sorted({key[:3] for key in by_seed_family})
    family_values: dict[str, np.ndarray] = {}
    for family in ("OU_II", "OU_III"):
        values = [
            float(np.mean(by_seed_family[(*key, family)]))
            for key in seed_keys
            if by_seed_family.get((*key, family))
        ]
        family_values[family] = np.asarray(values, dtype=np.float64)

    if family_values["OU_II"].size != family_values["OU_III"].size:
        raise ValueError("stationary normalized aggregate is not paired")
    if family_values["OU_II"].size == 0:
        raise ValueError("stationary normalized aggregate has no observations")

    rng = np.random.default_rng(stats_seed + 2)
    result: dict[str, dict[str, Any]] = {}
    for family, values in family_values.items():
        bootstrap_low, bootstrap_high = _bootstrap_mean_ci(
            values, bootstrap_resamples, rng
        )
        result[family] = {
            "n": int(values.size),
            "mean": float(np.mean(values)),
            "std": float(np.std(values, ddof=1)) if values.size > 1 else 0.0,
            "bootstrap_ci95_low": bootstrap_low,
            "bootstrap_ci95_high": bootstrap_high,
        }

    differences = family_values["OU_III"] - family_values["OU_II"]
    difference_std = (
        float(np.std(differences, ddof=1)) if differences.size > 1 else 0.0
    )
    cohen_dz = (
        float(np.mean(differences)) / difference_std
        if difference_std > 0.0 else math.nan
    )
    degrees_of_freedom = differences.size - 1
    correction = (
        math.gamma(degrees_of_freedom / 2.0)
        / (
            math.sqrt(degrees_of_freedom / 2.0)
            * math.gamma((degrees_of_freedom - 1.0) / 2.0)
        )
        if degrees_of_freedom > 1 else math.nan
    )
    bootstrap_low, bootstrap_high = _bootstrap_mean_ci(
        differences, bootstrap_resamples, rng
    )
    result["OU_III_minus_OU_II"] = {
        "n_pairs": int(differences.size),
        "mean_paired_difference": float(np.mean(differences)),
        "std_paired_difference": difference_std,
        "bootstrap_ci95_low": bootstrap_low,
        "bootstrap_ci95_high": bootstrap_high,
        "cohen_dz": cohen_dz,
        "hedges_gz": correction * cohen_dz,
        # Companions to the percentile bootstrap on the confirmatory endpoint
        # only.  Every other contrast in the study is descriptive, and giving
        # each of those four p-values as well would enlarge the family of
        # tests without adding evidence.
        **paired_inference(differences, rng),
    }
    return result


def write_publication_table(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    summary: Sequence[Mapping[str, Any]],
    effects: Sequence[Mapping[str, Any]],
    bootstrap_resamples: int,
    stats_seed: int,
    macros_path: Path | None = None,
) -> None:
    indexed_summary = {
        (
            str(row["scenario"]),
            str(row["family"]),
            str(row["mode"]),
            str(row["metric"]),
        ): row
        for row in summary
    }
    indexed_effects = {
        (
            str(row["scenario"]),
            str(row["family"]),
            str(row["comparison"]),
            str(row["metric"]),
        ): row
        for row in effects
    }
    all_scenarios = sorted(
        {
            str(row["scenario"])
            for row in summary
            if row["mode"] == "Adaptive" and row["metric"] == "disp_z_pct_hs"
        },
        key=scenario_sort_key,
    )
    # The confirmatory comparison is the declared JONSWAP ensemble plus the
    # controlled transition; PM-Stokes is reported separately below.
    scenarios = [s for s in all_scenarios if scenario_spectrum(s) == "jonswap"]
    pmstokes_scenarios = [
        s for s in all_scenarios if scenario_spectrum(s) == "pmstokes"
    ]
    transition_scenarios = [
        s for s in scenarios if s.startswith("nonstationary_")
    ]
    aggregate = stationary_normalized_aggregate(
        rows, bootstrap_resamples, stats_seed
    )
    difference = aggregate["OU_III_minus_OU_II"]
    pmstokes_aggregate = (
        stationary_normalized_aggregate(
            rows, bootstrap_resamples, stats_seed, spectrum="pmstokes"
        )
        if pmstokes_scenarios
        else None
    )

    def mean_std(
        scenario: str, family: str, mode: str, metric: str = "disp_z_pct_hs",
        digits: int = 2,
    ) -> str:
        row = indexed_summary.get((scenario, family, mode, metric))
        if row is None or not math.isfinite(float(row["mean"])):
            return "--"
        return (
            f"{float(row['mean']):.{digits}f} $\\pm$ "
            f"{float(row['std']):.{digits}f}"
        )

    def paired_ci(scenario: str, metric: str, digits: int) -> str:
        row = indexed_effects.get(
            (scenario, "OU_III_vs_OU_II", "OU_III_minus_OU_II", metric)
        )
        if row is None:
            return "--"
        return (
            f"{float(row['mean_paired_difference']):+.{digits}f} "
            f"[{float(row['bootstrap_ci95_low']):+.{digits}f}, "
            f"{float(row['bootstrap_ci95_high']):+.{digits}f}]"
        )

    def value_of(
        scenario: str, family: str, mode: str, metric: str, digits: int = 3
    ) -> str:
        row = indexed_summary.get((scenario, family, mode, metric))
        if row is None or not math.isfinite(float(row["mean"])):
            return "--"
        return f"{float(row['mean']):.{digits}f}"

    def ablation_effect(
        scenario: str, family: str, baseline: str, metric: str = "disp_z_pct_hs"
    ) -> Mapping[str, Any]:
        return indexed_effects[
            (scenario, family, f"Adaptive_minus_{baseline}", metric)
        ]

    def signed(value: float, digits: int = 3) -> str:
        return f"{value:+.{digits}f}"

    scenario_macro = {
        scenario: scenario_macro_token(scenario) for scenario in scenarios
    }

    # Every adaptation-ablation number quoted in the manuscript is emitted
    # here, so the prose cannot drift away from the regenerated rows.
    ablation_macros: list[str] = []
    for scenario, macro_name in scenario_macro.items():
        for family in ("OU_II", "OU_III"):
            short = "II" if family == "OU_II" else "III"
            for baseline, tag in (
                ("FixedNominal", "Nominal"),
                ("FixedOracle", "Oracle"),
            ):
                try:
                    effect = ablation_effect(scenario, family, baseline)
                except KeyError:
                    continue
                stem = f"OUValidation{short}{macro_name}{tag}"
                ablation_macros.extend(
                    (
                        rf"\providecommand{{\{stem}Difference}}"
                        rf"{{{signed(float(effect['mean_paired_difference']))}}}",
                        rf"\providecommand{{\{stem}Low}}"
                        rf"{{{signed(float(effect['bootstrap_ci95_low']))}}}",
                        rf"\providecommand{{\{stem}High}}"
                        rf"{{{signed(float(effect['bootstrap_ci95_high']))}}}",
                    )
                )

        # Share of the FixedNominal-to-FixedOracle mismatch that online tuning
        # actually recovers.  "Adaptation removes most of the mismatch" is only
        # an aggregate statement, and this is what makes the spread visible
        # instead of leaving the reader to infer it from three means.
        for family in ("OU_II", "OU_III"):
            short = "II" if family == "OU_II" else "III"
            means = {}
            for mode in ("Adaptive", "FixedNominal", "FixedOracle"):
                row = indexed_summary.get(
                    (scenario, family, mode, "disp_z_pct_hs")
                )
                if row is not None and math.isfinite(float(row["mean"])):
                    means[mode] = float(row["mean"])
            if len(means) != 3:
                continue
            mismatch = means["FixedNominal"] - means["FixedOracle"]
            if abs(mismatch) < 1e-9:
                continue
            recovered = 100.0 * (means["FixedNominal"] - means["Adaptive"]) / mismatch
            ablation_macros.append(
                rf"\providecommand{{\OUValidation{short}{macro_name}Recovered}}"
                rf"{{{recovered:.0f}}}"
            )

    # LaTeX sets a bare comma in math mode with trailing space, so the
    # thousands separator is braced.  Done here rather than inline because the
    # count is absent altogether when the ensemble is too small to randomize.
    randomization_patterns = _latex_statistic(
        difference.get("randomization_patterns"), ","
    ).replace(",", r"{,}")

    lines = [
        r"% Generated by tools/ou_validation.py from the committed full-study rows.",
        *ablation_macros,
        rf"\providecommand{{\OUValidationStationaryPairs}}{{{difference['n_pairs']}}}",
        rf"\providecommand{{\OUValidationOUIINormalizedMean}}{{{aggregate['OU_II']['mean']:.2f}}}",
        rf"\providecommand{{\OUValidationOUIINormalizedStd}}{{{aggregate['OU_II']['std']:.2f}}}",
        rf"\providecommand{{\OUValidationOUIIINormalizedMean}}{{{aggregate['OU_III']['mean']:.2f}}}",
        rf"\providecommand{{\OUValidationOUIIINormalizedStd}}{{{aggregate['OU_III']['std']:.2f}}}",
        rf"\providecommand{{\OUValidationNormalizedDifference}}{{{difference['mean_paired_difference']:+.3f}}}",
        rf"\providecommand{{\OUValidationNormalizedDifferenceLow}}{{{difference['bootstrap_ci95_low']:+.3f}}}",
        rf"\providecommand{{\OUValidationNormalizedDifferenceHigh}}{{{difference['bootstrap_ci95_high']:+.3f}}}",
        rf"\providecommand{{\OUValidationNormalizedDz}}"
        rf"{{{_latex_statistic(difference.get('cohen_dz'), '+.2f')}}}",
        rf"\providecommand{{\OUValidationNormalizedGz}}"
        rf"{{{_latex_statistic(difference.get('hedges_gz'), '+.2f')}}}",
        # Companion inference on the confirmatory endpoint, so the primary
        # claim does not rest on a single small-sample interval method.  These
        # need at least two pairs to exist at all, which a smoke run does not
        # have; see `_latex_statistic`.
        rf"\providecommand{{\OUValidationNormalizedTLow}}"
        rf"{{{_latex_statistic(difference.get('t_ci95_low'), '+.3f')}}}",
        rf"\providecommand{{\OUValidationNormalizedTHigh}}"
        rf"{{{_latex_statistic(difference.get('t_ci95_high'), '+.3f')}}}",
        rf"\providecommand{{\OUValidationNormalizedTStatistic}}"
        rf"{{{_latex_statistic(difference.get('t_statistic'), '+.2f')}}}",
        rf"\providecommand{{\OUValidationNormalizedTP}}"
        rf"{{{_latex_p_value(difference.get('t_p_value'))}}}",
        rf"\providecommand{{\OUValidationNormalizedSignNegative}}"
        rf"{{{_latex_statistic(difference.get('sign_negative'), 'd')}}}",
        rf"\providecommand{{\OUValidationNormalizedSignP}}"
        rf"{{{_latex_p_value(difference.get('sign_p_value'))}}}",
        rf"\providecommand{{\OUValidationNormalizedRandomizationP}}"
        rf"{{{_latex_p_value(difference.get('randomization_p_value'))}}}",
        rf"\providecommand{{\OUValidationNormalizedRandomizationPatterns}}"
        rf"{{{randomization_patterns}}}",
        # Interval construction, so the manuscript states the method rather
        # than leaving "bootstrap 95% interval" unqualified.
        rf"\providecommand{{\OUValidationBootstrapResamples}}{{{bootstrap_resamples:,}}}".replace(",", r"{,}"),
        rf"\providecommand{{\OUValidationStatsSeed}}{{{stats_seed}}}",
        r"\providecommand{\OUValidationBootstrapMethod}{nonparametric percentile bootstrap of the paired mean, resampling seed triplets with replacement}",
        r"\providecommand{\OUValidationBootstrapRNG}{NumPy PCG64}",
        # Number of distinct paired comparisons, for the multiplicity
        # statement.  Counting effect rows instead would count one comparison
        # once per metric and overstate the family of tests.
        rf"\providecommand{{\OUValidationPairedComparisons}}"
        rf"{{{len({(str(row['scenario']), str(row['family']), str(row['comparison'])) for row in effects})}}}",
        r"",
        r"\begin{table*}[t]",
        r"  \centering",
        r"  \caption{Ten-seed paired OU-family comparison over the final \SI{900}{s}. Values are mean $\pm$ sample standard deviation. $\Delta$ is OU--III Adaptive minus OU--II Adaptive; brackets give the paired bootstrap 95\% interval. Negative $\Delta$Z favors OU--III, whereas positive $\Delta$3D indicates larger OU--III three-dimensional displacement error. Yaw is carried alongside the displacement channels because the ensemble redraws the magnetometer calibration on every seed triplet: its spread is the spread of that draw, not of the estimator, and the two families share it. The single-draw deterministic yaw of Table~\ref{tab:mag-hi-yaw} is one realization of this distribution.}",
        r"  \label{tab:ou_mc_family}",
        r"  \footnotesize",
        r"  \setlength{\tabcolsep}{4.0pt}",
        r"  \begin{tabular}{@{}lrrrrrrr@{}}",
        r"    \toprule",
        r"    Scenario & OU--II Z [\%$H_s$] & OU--III Z [\%$H_s$] & $\Delta$Z [percentage points] & $\Delta$3D [m] & OU--II yaw [\si{\degree}] & OU--III yaw [\si{\degree}] & $\Delta$yaw [\si{\degree}] \\",
        r"    \midrule",
    ]
    for scenario in scenarios:
        lines.append(
            "    "
            + " & ".join(
                (
                    scenario_display_label(scenario),
                    mean_std(scenario, "OU_II", "Adaptive"),
                    mean_std(scenario, "OU_III", "Adaptive"),
                    paired_ci(scenario, "disp_z_pct_hs", 3),
                    paired_ci(scenario, "disp_3d_rms_m", 3),
                    mean_std(scenario, "OU_II", "Adaptive", "yaw_rms_deg"),
                    mean_std(scenario, "OU_III", "Adaptive", "yaw_rms_deg"),
                    paired_ci(scenario, "yaw_rms_deg", 3),
                )
            )
            + r" \\"
        )
    lines.extend(
        (
            r"    \bottomrule",
            r"  \end{tabular}",
            r"\end{table*}",
            r"",
            r"\begin{table*}[t]",
            r"  \centering",
            r"  \caption{Adaptation ablation for vertical-displacement RMS error over the final \SI{900}{s}; entries are mean $\pm$ sample standard deviation in percent of $H_s$ ($n=10$ paired seed triplets). FixedOracle is a noise-free, sea-specific or known-endpoint reference and is not deployable online.}",
            r"  \label{tab:ou_mc_adaptation}",
            r"  \footnotesize",
            r"  \setlength{\tabcolsep}{3.6pt}",
            r"  \begin{tabular}{@{}lrrrrrr@{}}",
            r"    \toprule",
            r"    & \multicolumn{3}{c}{OU--II Z [\%$H_s$]} & \multicolumn{3}{c}{OU--III Z [\%$H_s$]} \\",
            r"    \cmidrule(lr){2-4}\cmidrule(lr){5-7}",
            r"    Scenario & Adaptive & FixedNominal & FixedOracle & Adaptive & FixedNominal & FixedOracle \\",
            r"    \midrule",
        )
    )
    for scenario in scenarios:
        lines.append(
            "    "
            + " & ".join(
                (
                    scenario_display_label(scenario),
                    mean_std(scenario, "OU_II", "Adaptive"),
                    mean_std(scenario, "OU_II", "FixedNominal"),
                    mean_std(scenario, "OU_II", "FixedOracle"),
                    mean_std(scenario, "OU_III", "Adaptive"),
                    mean_std(scenario, "OU_III", "FixedNominal"),
                    mean_std(scenario, "OU_III", "FixedOracle"),
                )
            )
            + r" \\"
        )
    lines.extend((r"    \bottomrule", r"  \end{tabular}", r"\end{table*}"))

    sync_modes = [
        mode
        for pair in COVARIANCE_SYNC_PAIRS
        for mode in pair
        if any(row["mode"] == mode for row in summary)
    ]
    if len(sync_modes) == 4:
        lines.extend(
            (
                r"",
                r"\begin{table*}[t]",
                r"  \centering",
                r"  \caption{Covariance-policy ablation for vertical-displacement RMS error over the final \SI{900}{s}, in percent of $H_s$ (mean $\pm$ sample standard deviation, $n=10$ paired seed triplets). \emph{Held} repeats the paired mode with the periodic re-alignment of the posterior $\mat{P}_{a_w a_w}$ marginal switched off, leaving only the reconfiguration-event alignment. Comparing columns within a tuning mode isolates the covariance policy from online adaptation.}",
                r"  \label{tab:ou_mc_covsync}",
                r"  \footnotesize",
                r"  \setlength{\tabcolsep}{3.2pt}",
                r"  \begin{tabular}{@{}lrrrrrrrr@{}}",
                r"    \toprule",
                r"    & \multicolumn{4}{c}{OU--II Z [\%$H_s$]} & \multicolumn{4}{c}{OU--III Z [\%$H_s$]} \\",
                r"    \cmidrule(lr){2-5}\cmidrule(lr){6-9}",
                r"    Scenario & Adapt. & Adapt.\ held & FixedNom. & FixedNom.\ held"
                r" & Adapt. & Adapt.\ held & FixedNom. & FixedNom.\ held \\",
                r"    \midrule",
            )
        )
        for scenario in scenarios:
            lines.append(
                "    "
                + " & ".join(
                    (
                        scenario_display_label(scenario),
                        *(
                            mean_std(scenario, family, mode)
                            for family in ("OU_II", "OU_III")
                            for mode in (
                                "Adaptive",
                                "AdaptiveHeldCovariance",
                                "FixedNominal",
                                "FixedNominalHeldCovariance",
                            )
                        ),
                    )
                )
                + r" \\"
            )
        lines.extend((r"    \bottomrule", r"  \end{tabular}", r"\end{table*}"))

        # Largest |mean paired difference| attributable to the covariance reset
        # alone, across scenarios, families, and both matched pairs.
        candidates = [
            row
            for row in effects
            if any(
                str(row["comparison"]) == f"{left}_minus_{right}"
                for left, right in COVARIANCE_SYNC_PAIRS
            )
        ]
        vertical = [
            row for row in candidates if row["metric"] == "disp_z_pct_hs"
        ]
        worst = max(
            vertical, key=lambda row: abs(float(row["mean_paired_difference"]))
        )
        transition = [
            row
            for row in vertical
            if str(row["scenario"]).startswith("nonstationary_")
            and str(row["comparison"]).startswith("Adaptive_minus_")
        ]
        worst_transition = max(
            transition, key=lambda row: abs(float(row["mean_paired_difference"]))
        )
        lines.extend(
            (
                r"",
                rf"\providecommand{{\OUValidationCovSyncWorstScenario}}{{{scenario_display_label(str(worst['scenario']))}}}",
                rf"\providecommand{{\OUValidationCovSyncWorstFamily}}{{{str(worst['family']).replace('_', '--')}}}",
                rf"\providecommand{{\OUValidationCovSyncWorstDifference}}{{{float(worst['mean_paired_difference']):+.3f}}}",
                rf"\providecommand{{\OUValidationCovSyncWorstLow}}{{{float(worst['bootstrap_ci95_low']):+.3f}}}",
                rf"\providecommand{{\OUValidationCovSyncWorstHigh}}{{{float(worst['bootstrap_ci95_high']):+.3f}}}",
                rf"\providecommand{{\OUValidationCovSyncTransitionFamily}}{{{str(worst_transition['family']).replace('_', '--')}}}",
                rf"\providecommand{{\OUValidationCovSyncTransitionDifference}}{{{float(worst_transition['mean_paired_difference']):+.3f}}}",
                rf"\providecommand{{\OUValidationCovSyncTransitionLow}}{{{float(worst_transition['bootstrap_ci95_low']):+.3f}}}",
                rf"\providecommand{{\OUValidationCovSyncTransitionHigh}}{{{float(worst_transition['bootstrap_ci95_high']):+.3f}}}",
            )
        )

    lines.extend(
        _three_dimensional_table(scenarios, mean_std, paired_ci, indexed_effects)
    )
    lines.extend(
        _channel_ablation_table(
            scenarios, summary, mean_std, indexed_effects, scenario_macro
        )
    )
    if transition_scenarios:
        lines.extend(
            _transition_segment_table(
                transition_scenarios[0], mean_std, value_of, indexed_effects
            )
        )
    if pmstokes_scenarios and pmstokes_aggregate is not None:
        lines.extend(
            _pmstokes_table(
                pmstokes_scenarios, mean_std, paired_ci, pmstokes_aggregate
            )
        )
    lines.extend(_direction_table(all_scenarios, indexed_summary, mean_std))

    assert_latex_macro_names(lines)
    header = (
        r"% Generated by tools/ou_validation.py from the committed full-study rows."
    )
    if macros_path is None:
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return

    # The manuscript quotes these values in the protocol section, which is
    # typeset before the results tables are input, so the definitions are
    # emitted separately and included early.  Splitting them here keeps a
    # single source of truth for both files.
    definitions = [
        line for line in lines if line.lstrip().startswith(r"\providecommand")
    ]
    tables = [
        line
        for line in lines
        if not line.lstrip().startswith(r"\providecommand")
        and not line.startswith("%")
    ]
    macros_path.write_text(
        "\n".join([header, *definitions]) + "\n", encoding="utf-8"
    )
    path.write_text("\n".join([header, *tables]) + "\n", encoding="utf-8")


def _three_dimensional_table(
    scenarios: Sequence[str],
    mean_std: Any,
    paired_ci: Any,
    indexed_effects: Mapping[tuple[str, str, str, str], Mapping[str, Any]],
) -> list[str]:
    """Absolute per-axis and 3D displacement RMS for both families.

    The paired 3D difference alone does not say whether a positive value is a
    small fraction of a large error or the whole of a small one, and it hides
    which horizontal axis moves.  Both filters are therefore reported in
    absolute metres beside the difference.
    """

    lines = [
        r"",
        r"\begin{table*}[t]",
        r"  \centering",
        r"  \caption{Absolute per-axis and three-dimensional displacement RMS error over the final \SI{900}{s} (Adaptive mode, mean in metres, $n=10$ paired seed triplets). $\Delta$ is OU--III minus OU--II with its paired bootstrap 95\% interval; positive $\Delta$ favors OU--II.}",
        r"  \label{tab:ou_mc_axes}",
        r"  \footnotesize",
        r"  \setlength{\tabcolsep}{3.0pt}",
        r"  \begin{tabular}{@{}lrrrrrrrrr@{}}",
        r"    \toprule",
        r"    & \multicolumn{4}{c}{OU--II RMS [m]} & \multicolumn{4}{c}{OU--III RMS [m]} & \\",
        r"    \cmidrule(lr){2-5}\cmidrule(lr){6-9}",
        r"    Scenario & X & Y & Z & 3D & X & Y & Z & 3D & $\Delta$3D [m] \\",
        r"    \midrule",
    ]
    for scenario in scenarios:
        lines.append(
            "    "
            + " & ".join(
                (
                    scenario_display_label(scenario),
                    *(
                        mean_std(scenario, family, "Adaptive", metric, 3)
                        for family in ("OU_II", "OU_III")
                        for metric in (
                            "disp_x_rms_m",
                            "disp_y_rms_m",
                            "disp_z_rms_m",
                            "disp_3d_rms_m",
                        )
                    ),
                    paired_ci(scenario, "disp_3d_rms_m", 3),
                )
            )
            + r" \\"
        )
    lines.extend((r"    \bottomrule", r"  \end{tabular}", r"\end{table*}"))

    # Sign of the 3D effect per scenario, so the prose can say "four of five"
    # without a hand count that can drift from the table.
    resolved_higher = 0
    unresolved = 0
    for scenario in scenarios:
        row = indexed_effects.get(
            (scenario, "OU_III_vs_OU_II", "OU_III_minus_OU_II", "disp_3d_rms_m")
        )
        if row is None:
            continue
        low = float(row["bootstrap_ci95_low"])
        high = float(row["bootstrap_ci95_high"])
        if low > 0.0:
            resolved_higher += 1
        elif low <= 0.0 <= high:
            unresolved += 1
    lines.extend(
        (
            r"",
            rf"\providecommand{{\OUValidationThreeDScenarios}}{{{len(scenarios)}}}",
            rf"\providecommand{{\OUValidationThreeDHigherCount}}{{{resolved_higher}}}",
            rf"\providecommand{{\OUValidationThreeDUnresolvedCount}}{{{unresolved}}}",
        )
    )
    return lines


def _channel_ablation_table(
    scenarios: Sequence[str],
    summary: Sequence[Mapping[str, Any]],
    mean_std: Any,
    indexed_effects: Mapping[tuple[str, str, str, str], Mapping[str, Any]],
    scenario_macro: Mapping[str, str],
) -> list[str]:
    """Which adaptation channel earns the benefit.

    Adaptive versus FixedNominal moves three parameters at once.  Because the
    deployed law sets r_S from tau and sigma_aw, that comparison cannot say
    whether the vertical benefit comes from adapting the OU process or from
    adapting the integral regularization.  These two extra modes complete the
    2x2 factorial by freezing one channel at a time.
    """

    if not any(row["mode"] in PARTIAL_MODES for row in summary):
        return []

    lines = [
        r"",
        r"\begin{table*}[t]",
        r"  \centering",
        r"  \caption{OU--III adaptation-channel ablation for vertical-displacement RMS error over the final \SI{900}{s}, in percent of $H_s$ (mean $\pm$ sample standard deviation, $n=10$ paired seed triplets). The four columns are a $2\times2$ factorial: each channel is either adapted online or frozen at the FixedNominal operating point. \emph{$r_S$ only} freezes $\tau$ and $\sigma_{aw}$ while the integral pseudo-measurement scale keeps adapting; \emph{OU only} does the reverse. Because the deployed law sets $r_S=\clip(0.35\,\sigma_{aw}\tau^{3},0.15,400)$, $r_S$ is derived from the live $\tau$ and $\sigma_{aw}$ estimates even when those are frozen on the way to the filter.}",
        r"  \label{tab:ou_mc_channels}",
        r"  \footnotesize",
        r"  \setlength{\tabcolsep}{4.0pt}",
        r"  \begin{tabular}{@{}lrrrr@{}}",
        r"    \toprule",
        r"    & \multicolumn{2}{c}{$r_S$ adapted} & \multicolumn{2}{c}{$r_S$ frozen} \\",
        r"    \cmidrule(lr){2-3}\cmidrule(lr){4-5}",
        r"    Scenario & $\tau,\sigma_{aw}$ adapted & $\tau,\sigma_{aw}$ frozen"
        r" & $\tau,\sigma_{aw}$ adapted & $\tau,\sigma_{aw}$ frozen \\",
        r"    & (Adaptive) & ($r_S$ only) & (OU only) & (FixedNominal) \\",
        r"    \midrule",
    ]
    for scenario in scenarios:
        lines.append(
            "    "
            + " & ".join(
                (
                    scenario_display_label(scenario),
                    mean_std(scenario, "OU_III", "Adaptive"),
                    mean_std(scenario, "OU_III", "AdaptiveRSOnly"),
                    mean_std(scenario, "OU_III", "AdaptiveOUOnly"),
                    mean_std(scenario, "OU_III", "FixedNominal"),
                )
            )
            + r" \\"
        )
    lines.extend((r"    \bottomrule", r"  \end{tabular}", r"\end{table*}", r""))

    for scenario in scenarios:
        stem = f"OUValidationChannel{scenario_macro.get(scenario, '')}"
        for mode, tag in (
            ("AdaptiveRSOnly", "RSOnly"),
            ("AdaptiveOUOnly", "OUOnly"),
        ):
            row = indexed_effects.get(
                (
                    scenario,
                    "OU_III",
                    f"{mode}_minus_FixedNominal",
                    "disp_z_pct_hs",
                )
            )
            if row is None:
                continue
            lines.extend(
                (
                    rf"\providecommand{{\{stem}{tag}Difference}}"
                    rf"{{{float(row['mean_paired_difference']):+.3f}}}",
                    rf"\providecommand{{\{stem}{tag}Low}}"
                    rf"{{{float(row['bootstrap_ci95_low']):+.3f}}}",
                    rf"\providecommand{{\{stem}{tag}High}}"
                    rf"{{{float(row['bootstrap_ci95_high']):+.3f}}}",
                )
            )
    return lines


def _transition_segment_table(
    scenario: str,
    mean_std: Any,
    value_of: Any,
    indexed_effects: Mapping[tuple[str, str, str, str], Mapping[str, Any]],
) -> list[str]:
    """Split the transition score into its four sea-state intervals.

    The aggregate normalizes by the final H_s although the window opens in the
    start sea, so its percentage is not comparable with a stationary score.
    Absolute RMS and a reference-RMS normalization are reported per interval.

    The endpoint sea is split at one crossfade length past the blend: the
    run-on interval prices how long the schedule keeps carrying the old sea in
    its averages, and the settled interval is the endpoint sea proper.
    """

    segments = (
        ("start", "Pure start sea"),
        ("blend", "Crossfade"),
        ("recover", "Endpoint sea, run-on"),
        ("end", "Endpoint sea, settled"),
        ("", "Whole window"),
    )
    lines = [
        r"",
        r"\begin{table*}[t]",
        r"  \centering",
        r"  \caption{Controlled transition scored by interval rather than as one window (Adaptive mode, mean over $n=10$ seed triplets). $Z_{\mathrm{ref}}$ is the RMS of the reference vertical displacement in the same interval, so $Z/Z_{\mathrm{ref}}$ is a scale-free normalization that remains meaningful while the sea state changes; $Z/H_s^{\mathrm{end}}$ is the whole-window convention normalized by the final $H_s=\SI{4.0}{m}$ and is reported for continuity only. The window opens in the \SI{1.5}{m} start sea, so the two normalizations disagree by construction. The endpoint sea is split at one crossfade length past the blend: the run-on interval is where a schedule that averages too long is still carrying the start sea, and the settled interval is the endpoint sea proper.}",
        r"  \label{tab:ou_transition_segments}",
        r"  \footnotesize",
        r"  \setlength{\tabcolsep}{3.4pt}",
        r"  \begin{tabular}{@{}lrrrrrrrr@{}}",
        r"    \toprule",
        r"    & \multicolumn{4}{c}{OU--II} & \multicolumn{4}{c}{OU--III} \\",
        r"    \cmidrule(lr){2-5}\cmidrule(lr){6-9}",
        r"    Interval & $Z$ [m] & $Z_{\mathrm{ref}}$ [m] & $Z/Z_{\mathrm{ref}}$ [\%] & $Z/H_s^{\mathrm{end}}$ [\%]"
        r" & $Z$ [m] & $Z_{\mathrm{ref}}$ [m] & $Z/Z_{\mathrm{ref}}$ [\%] & $Z/H_s^{\mathrm{end}}$ [\%] \\",
        r"    \midrule",
    ]
    for key, label in segments:
        prefix = f"seg_{key}_" if key else ""
        lines.append(
            "    "
            + " & ".join(
                (
                    label,
                    *(
                        value_of(
                            scenario,
                            family,
                            "Adaptive",
                            f"{prefix}{metric}",
                            digits,
                        )
                        for family in ("OU_II", "OU_III")
                        for metric, digits in (
                            ("disp_z_rms_m", 3),
                            ("disp_z_ref_rms_m", 3),
                            ("disp_z_pct_refrms", 2),
                            ("disp_z_pct_hs", 2),
                        )
                    ),
                )
            )
            + r" \\"
        )
    lines.extend((r"    \bottomrule", r"  \end{tabular}", r"\end{table*}", r""))

    for key, tag in (("start", "Start"), ("blend", "Blend"), ("end", "End")):
        for family, short in (("OU_II", "II"), ("OU_III", "III")):
            row = indexed_effects.get(
                (
                    scenario,
                    family,
                    "Adaptive_minus_FixedNominal",
                    f"seg_{key}_disp_z_pct_refrms",
                )
            )
            if row is None:
                continue
            stem = f"OUValidationSegment{short}{tag}Nominal"
            lines.extend(
                (
                    rf"\providecommand{{\{stem}Difference}}"
                    rf"{{{float(row['mean_paired_difference']):+.3f}}}",
                    rf"\providecommand{{\{stem}Low}}"
                    rf"{{{float(row['bootstrap_ci95_low']):+.3f}}}",
                    rf"\providecommand{{\{stem}High}}"
                    rf"{{{float(row['bootstrap_ci95_high']):+.3f}}}",
                )
            )
    return lines


def _pmstokes_table(
    scenarios: Sequence[str],
    mean_std: Any,
    paired_ci: Any,
    aggregate: Mapping[str, Mapping[str, Any]],
) -> list[str]:
    """Ten-seed PM-Stokes ensemble, kept out of the primary aggregate."""

    difference = aggregate["OU_III_minus_OU_II"]
    lines = [
        r"",
        r"\begin{table*}[t]",
        r"  \centering",
        r"  \caption{Ten-seed paired OU-family comparison on the PM--Stokes seas, scored over the same final \SI{900}{s} window and the same seed triplets as the JONSWAP ensemble. PM--Stokes carries third-order bound harmonics that JONSWAP does not, so it is reported as a separate declared ensemble and is not pooled into the primary aggregate.}",
        r"  \label{tab:ou_mc_pmstokes}",
        r"  \footnotesize",
        r"  \setlength{\tabcolsep}{4.0pt}",
        r"  \begin{tabular}{@{}lrrrrrrr@{}}",
        r"    \toprule",
        r"    Scenario & OU--II Z [\%$H_s$] & OU--III Z [\%$H_s$] & $\Delta$Z [percentage points] & $\Delta$3D [m] & OU--II yaw [\si{\degree}] & OU--III yaw [\si{\degree}] & $\Delta$yaw [\si{\degree}] \\",
        r"    \midrule",
    ]
    for scenario in scenarios:
        lines.append(
            "    "
            + " & ".join(
                (
                    scenario_display_label(scenario),
                    mean_std(scenario, "OU_II", "Adaptive"),
                    mean_std(scenario, "OU_III", "Adaptive"),
                    paired_ci(scenario, "disp_z_pct_hs", 3),
                    paired_ci(scenario, "disp_3d_rms_m", 3),
                    mean_std(scenario, "OU_II", "Adaptive", "yaw_rms_deg"),
                    mean_std(scenario, "OU_III", "Adaptive", "yaw_rms_deg"),
                    paired_ci(scenario, "yaw_rms_deg", 3),
                )
            )
            + r" \\"
        )
    lines.extend((r"    \bottomrule", r"  \end{tabular}", r"\end{table*}", r""))
    lines.extend(
        (
            rf"\providecommand{{\OUValidationPMStokesOUIINormalizedMean}}{{{aggregate['OU_II']['mean']:.2f}}}",
            rf"\providecommand{{\OUValidationPMStokesOUIINormalizedStd}}{{{aggregate['OU_II']['std']:.2f}}}",
            rf"\providecommand{{\OUValidationPMStokesOUIIINormalizedMean}}{{{aggregate['OU_III']['mean']:.2f}}}",
            rf"\providecommand{{\OUValidationPMStokesOUIIINormalizedStd}}{{{aggregate['OU_III']['std']:.2f}}}",
            rf"\providecommand{{\OUValidationPMStokesDifference}}{{{difference['mean_paired_difference']:+.3f}}}",
            rf"\providecommand{{\OUValidationPMStokesDifferenceLow}}{{{difference['bootstrap_ci95_low']:+.3f}}}",
            rf"\providecommand{{\OUValidationPMStokesDifferenceHigh}}{{{difference['bootstrap_ci95_high']:+.3f}}}",
            rf"\providecommand{{\OUValidationPMStokesPairs}}{{{difference['n_pairs']}}}",
        )
    )
    return lines


def _direction_table(
    scenarios: Sequence[str],
    indexed_summary: Mapping[tuple[str, str, str, str], Mapping[str, Any]],
    mean_std: Any,
) -> list[str]:
    """Ten-seed direction accuracy against the generator azimuth.

    The historical direction report gives an estimated heading and a
    Toward/Away/Uncertain split for one deterministic window, but no error
    against truth and no seed-to-seed spread.  Both are scored here over the
    same \\SI{900}{s} window as the displacement metrics.
    """

    stationary = [s for s in scenarios if not s.startswith("nonstationary_")]
    if not stationary:
        return []
    if not any(
        (s, "OU_III", "Adaptive", "dir_axis_abs_error_deg") in indexed_summary
        for s in stationary
    ):
        return []

    lines = [
        r"",
        r"\begin{table*}[t]",
        r"  \centering",
        r"  \caption{Ten-seed OU--III wave-direction results over the final \SI{900}{s}, against the generator azimuth of each record. The propagation axis is defined modulo \SI{180}{\degree}, so $|\Delta\theta|$ is the absolute axial error of the circular-mean estimate and $\theta_{\mathrm{RMSE}}$ is the sample-wise axial RMS error. \emph{Sense} is a genuine correctness rate: the estimator's directed propagation vector, with the vessel heading removed, scored against the physical propagation direction of the record, which is the generator azimuth plus \SI{180}{\degree}. \emph{Unresolved} is the share below the confidence and amplitude thresholds. The FORWARD/BACKWARD classes the estimator exports are relative to the axis representative it happens to return and invert under a \SI{180}{\degree} heading change, so they are not scored here. Entries are mean $\pm$ sample standard deviation over the seed triplets.}",
        r"  \label{tab:ou_mc_direction}",
        r"  \footnotesize",
        r"  \setlength{\tabcolsep}{3.6pt}",
        r"  \begin{tabular}{@{}llrrrrr@{}}",
        r"    \toprule",
        r"    Spectrum & Scenario & $|\Delta\theta|$ [\si{\degree}] & $\theta_{\mathrm{RMSE}}$ [\si{\degree}]"
        r" & Axial SD [\si{\degree}] & Sense [\%] & Unresolved [\%] \\",
        r"    \midrule",
    ]
    for scenario in stationary:
        spectrum = (
            "PM--Stokes"
            if scenario_spectrum(scenario) == "pmstokes"
            else "JONSWAP"
        )
        lines.append(
            "    "
            + " & ".join(
                (
                    spectrum,
                    scenario_display_label(scenario),
                    *(
                        mean_std(scenario, "OU_III", "Adaptive", metric, 2)
                        for metric in (
                            "dir_axis_abs_error_deg",
                            "dir_axis_rmse_deg",
                            "dir_axis_circ_std_deg",
                            "dir_travel_correct_pct",
                            "dir_travel_unresolved_pct",
                        )
                    ),
                )
            )
            + r" \\"
        )
    lines.extend((r"    \bottomrule", r"  \end{tabular}", r"\end{table*}", r""))

    def scenario_means(metric: str) -> list[float]:
        return [
            float(indexed_summary[(s, "OU_III", "Adaptive", metric)]["mean"])
            for s in stationary
            if (s, "OU_III", "Adaptive", metric) in indexed_summary
        ]

    def pooled(metric: str) -> float:
        values = scenario_means(metric)
        return float(np.mean(values)) if values else math.nan

    def worst(metric: str) -> float:
        values = scenario_means(metric)
        return float(np.max(values)) if values else math.nan

    lines.extend(
        (
            rf"\providecommand{{\OUValidationDirectionAbsError}}{{{pooled('dir_axis_abs_error_deg'):.2f}}}",
            rf"\providecommand{{\OUValidationDirectionWorstAbsError}}{{{worst('dir_axis_abs_error_deg'):.2f}}}",
            rf"\providecommand{{\OUValidationDirectionRMSE}}{{{pooled('dir_axis_rmse_deg'):.2f}}}",
            rf"\providecommand{{\OUValidationDirectionDominant}}{{{pooled('dir_sense_dominant_pct'):.1f}}}",
            rf"\providecommand{{\OUValidationDirectionUncertain}}{{{pooled('dir_sense_uncertain_pct'):.1f}}}",
            rf"\providecommand{{\OUValidationTravelCorrect}}{{{pooled('dir_travel_correct_pct'):.1f}}}",
            rf"\providecommand{{\OUValidationTravelWrong}}{{{pooled('dir_travel_wrong_pct'):.1f}}}",
            rf"\providecommand{{\OUValidationTravelUnresolved}}{{{pooled('dir_travel_unresolved_pct'):.1f}}}",
            rf"\providecommand{{\OUValidationTravelAbsError}}{{{pooled('dir_travel_abs_error_deg'):.2f}}}",
        )
    )
    return lines


def write_tuning_points_table(
    path: Path,
    tuning_points: Mapping[str, Mapping[str, "TuningPoint"]],
    nominal_name: str,
    endpoint_name: str,
) -> None:
    """Emit the exact frozen operating points used by the fixed modes.

    Every fixed-mode run in the study is fully determined by one row of this
    table plus the fixed internal anisotropy factors, so the ablation can be
    reproduced without re-running the calibration pass.
    """

    def label(name: str) -> str:
        if name == endpoint_name:
            return r"transition endpoint $H_s=\SI{4.0}{m}$, $T_p=\SI{11.4}{s}$"
        text = scenario_display_label(name)
        if name == nominal_name:
            return text + r"\,$^{\dagger}$"
        return text

    names = sorted(
        {name for points in tuning_points.values() for name in points},
        key=lambda name: (name == endpoint_name, scenario_sort_key(name)),
    )

    lines = [
        r"% Generated by tools/ou_validation.py; do not edit by hand.",
        r"\begin{table*}[t]",
        r"  \centering",
        r"  \caption{Frozen operating points for the fixed-tuning modes, calibrated once from each noise-free \SI{1200}{s} reference record. $\dagger$ marks the single point used by FixedNominal in \emph{every} scenario; FixedOracle uses the row matching its own scenario, and the transition-endpoint row for the non-stationary case. Values are the vertical/base parameters; the filter scales them by its deployed anisotropy constants before use, $(S_\sigma\sigma_{aw},S_\sigma\sigma_{aw},\sigma_{aw})$ for OU--III acceleration and $\operatorname{diag}(\rho_{xy}r_S,\rho_{xy}r_S,r_S)^2$ for the integral pseudo-measurement, and $1.5\sigma_{aw}$ horizontally with $r_{p0}$ isotropic for OU--II.  The OU--III pair is quoted in Table~\ref{tab:implementation-gates} and is $S_\sigma=\rho_{xy}=1$ for the shipped filter; a bundle produced before that re-gauge was scored at $S_\sigma=1.87$.}",
        r"  \label{tab:ou_fixed_points}",
        r"  \footnotesize",
        r"  \setlength{\tabcolsep}{4pt}",
        r"  \begin{tabular}{@{}lrrrrrr@{}}",
        r"    \toprule",
        r"    & \multicolumn{3}{c}{OU--III} & \multicolumn{3}{c}{OU--II} \\",
        r"    \cmidrule(lr){2-4}\cmidrule(lr){5-7}",
        r"    Calibration record & $\tau$ [s] & $\sigma_{aw}$ [\si{\meter\per\second\squared}] & $r_S$ [\si{\meter\second}]"
        r" & $\tau$ [s] & $\sigma_{aw}$ [\si{\meter\per\second\squared}] & $r_{p0}$ [m] \\",
        r"    \midrule",
    ]
    for name in names:
        third = tuning_points.get("OU_III", {}).get(name)
        second = tuning_points.get("OU_II", {}).get(name)

        def cell(value: float | None, digits: int = 3) -> str:
            return "--" if value is None else f"{value:.{digits}f}"

        lines.append(
            "    "
            + " & ".join(
                (
                    label(name),
                    cell(third.tau_s if third else None),
                    cell(third.sigma_a_mps2 if third else None),
                    cell(third.RS_ms if third else None),
                    cell(second.tau_s if second else None),
                    cell(second.sigma_a_mps2 if second else None),
                    cell(second.R_p0_std_m if second else None),
                )
            )
            + r" \\"
        )
    lines.extend((r"    \bottomrule", r"  \end{tabular}", r"\end{table*}"))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_metric_plot(
    path: Path,
    summary: Sequence[Mapping[str, Any]],
    metric: str,
    ylabel: str,
    modes: Sequence[str] = PRIMARY_MODES,
) -> None:
    cache_dir = Path(tempfile.gettempdir()) / "ocean-imu-matplotlib-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    import matplotlib

    matplotlib.use("Agg")
    # Without these the bundle cannot be reproduced: matplotlib derives svg
    # element ids from a per-process salt and stamps a creation date, so two
    # identical runs disagree on bytes the manifest hashes.
    matplotlib.rcParams["svg.hashsalt"] = "ocean-imu-ou-validation"
    import matplotlib.pyplot as plt

    rows = [
        row for row in summary
        if row["metric"] == metric and str(row["mode"]) in modes
    ]
    scenarios = sorted(
        {str(row["scenario"]) for row in rows}, key=scenario_sort_key
    )
    order = {mode: index for index, mode in enumerate(modes)}
    groups = sorted(
        {(str(row["family"]), str(row["mode"])) for row in rows},
        key=lambda item: (item[0], order[item[1]]),
    )
    figure, axis = plt.subplots(figsize=(max(7.0, 1.7 * len(scenarios)), 4.8))
    x = np.arange(len(scenarios), dtype=float)
    offsets = np.linspace(-0.28, 0.28, max(1, len(groups)))
    for offset, (family, mode) in zip(offsets, groups):
        by_scenario = {
            str(row["scenario"]): row
            for row in rows
            if row["family"] == family and row["mode"] == mode
        }
        means = np.asarray([
            float(by_scenario[scenario]["mean"]) if scenario in by_scenario else np.nan
            for scenario in scenarios
        ])
        lows = np.asarray([
            float(by_scenario[scenario]["bootstrap_ci95_low"])
            if scenario in by_scenario else np.nan
            for scenario in scenarios
        ])
        highs = np.asarray([
            float(by_scenario[scenario]["bootstrap_ci95_high"])
            if scenario in by_scenario else np.nan
            for scenario in scenarios
        ])
        axis.errorbar(
            x + offset,
            means,
            yerr=np.vstack((means - lows, highs - means)),
            marker="o",
            capsize=3,
            linewidth=1.2,
            label=(
                f"{family.replace('OU_II', 'OU–II').replace('OU_III', 'OU–III')} "
                f"{mode.replace('HeldCovariance', ' held-cov').replace('Fixed', 'Fixed ')}"
            ),
        )
    # The two spectral families reuse the same nominal heights, so the tick
    # labels must carry the spectrum whenever both are plotted.
    spectra = {scenario_spectrum(scenario) for scenario in scenarios}
    def tick_label(scenario: str) -> str:
        label = scenario_display_label(scenario)
        if len(spectra) > 1 and not scenario.startswith("nonstationary_"):
            suffix = "PM" if scenario_spectrum(scenario) == "pmstokes" else "JS"
            return f"{label}\n{suffix}"
        return label

    axis.set_xticks(x, [tick_label(scenario) for scenario in scenarios])
    axis.set_xlabel("Scenario")
    axis.set_ylabel(ylabel)
    axis.grid(axis="y", alpha=0.3)
    axis.legend(fontsize=8, ncol=2)
    figure.tight_layout()
    figure.savefig(path, format="svg", metadata={"Date": None})
    plt.close(figure)


def transition_window_composition(
    transition_start_sec: float,
    transition_end_sec: float,
    duration_sec: float,
    window_sec: float,
) -> dict[str, float]:
    """Split the scoring window into pure-start, blended, and pure-end time.

    The transition is a crossfade between two independent stationary records,
    so the score is not a uniform sample of "transitioning" conditions: it
    contains a run-in at the start sea, the blend itself, and a long tail at
    the endpoint sea.  Reporting the split keeps that visible.

    The endpoint tail is reported twice: once whole (``pure_end_sea_sec``) and
    once split into the crossfade-length run-on the schedule needs to shed the
    old sea from its averages (``recovery_sec``) and the settled remainder
    (``settled_end_sea_sec``).  Those two are what the ``recover`` and ``end``
    scoring segments cover.
    """

    window_start = duration_sec - window_sec
    start_end = max(window_start, min(transition_start_sec, duration_sec))
    blend_end = max(start_end, min(transition_end_sec, duration_sec))
    # One crossfade length past the blend, measured on the crossfade the record
    # was built with rather than on the part of it the window happens to see.
    crossfade_sec = max(0.0, transition_end_sec - transition_start_sec)
    recovery_end = max(blend_end, min(blend_end + crossfade_sec, duration_sec))
    return {
        "window_start_sec": window_start,
        "window_end_sec": duration_sec,
        "pure_start_sea_sec": start_end - window_start,
        "blended_sec": blend_end - start_end,
        "pure_end_sea_sec": duration_sec - blend_end,
        "recovery_sec": recovery_end - blend_end,
        "settled_end_sea_sec": duration_sec - recovery_end,
    }


def mixture_significant_height_m(
    start_height_m: float, end_height_m: float, weight: float
) -> float:
    """Effective H_s of the crossfade at blend weight ``weight``.

    The two records are independent, so their variances add under the quintic
    weights.  The blend therefore passes through a lower height than a linear
    H_s ramp between the same endpoints would.
    """

    variance = ((1.0 - weight) * start_height_m) ** 2 + (weight * end_height_m) ** 2
    return math.sqrt(variance)


def _read_timeseries_columns(
    path: Path, names: Sequence[str]
) -> dict[str, np.ndarray]:
    with path.open("r", encoding="utf-8") as stream:
        header = stream.readline().strip().split(",")
    indices = [header.index(name) for name in names]
    data = np.loadtxt(
        path, delimiter=",", skiprows=1, usecols=indices, dtype=np.float64
    )
    return {name: data[:, position] for position, name in enumerate(names)}


def read_diagnostic_timeseries(
    family: str, surrogate_path: Path
) -> dict[str, np.ndarray]:
    """Read and delete the time series a diagnostic simulator run just wrote."""

    # The simulator derives its output name from the input path it was given,
    # so the time series lands next to the surrogate rather than in the cwd.
    output_name = "w3d_" + surrogate_path.name.removeprefix("wave_data_")
    suffix = FAMILY_TIMESERIES_SUFFIX[family]
    output_path = surrogate_path.parent / output_name.replace(
        ".csv", f"{suffix}.csv"
    )
    if not output_path.exists():
        raise FileNotFoundError(f"simulator did not write {output_path}")

    wanted = (
        "time",
        "disp_ref_z",
        "disp_est_z",
        "tau_applied",
        "sigma_a_applied",
        "R_p0_applied",
        "freq_tracker_hz",
    )
    try:
        return _read_timeseries_columns(output_path, wanted)
    finally:
        output_path.unlink(missing_ok=True)


def rolling_significant_height(
    reference: np.ndarray, average_sec: float = 60.0
) -> np.ndarray:
    """Rolling significant height of a reference elevation record.

    This reads the sea state out of the record itself rather than out of the
    blend weight the record was built from, so a figure that shows both can be
    checked against the protocol it claims to follow.
    """

    window_samples = max(2, int(round(average_sec / DT_SECONDS)))
    kernel = np.ones(window_samples) / window_samples
    rolling_mean = np.convolve(reference, kernel, mode="same")
    rolling_mean_square = np.convolve(reference**2, kernel, mode="same")
    return 4.0 * np.sqrt(np.maximum(rolling_mean_square - rolling_mean**2, 0.0))


def reproducible_pyplot():
    """Return pyplot configured so two identical runs agree on SVG bytes."""

    cache_dir = Path(tempfile.gettempdir()) / "ocean-imu-matplotlib-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    import matplotlib

    matplotlib.use("Agg")
    # Without these the bundle cannot be reproduced: matplotlib derives svg
    # element ids from a per-process salt and stamps a creation date, so two
    # identical runs disagree on bytes the manifest hashes.
    matplotlib.rcParams["svg.hashsalt"] = "ocean-imu-ou-validation"
    import matplotlib.pyplot as plt

    return plt


def plot_transition_diagnostic(
    svg_path: Path,
    time: np.ndarray,
    reference: np.ndarray,
    estimate: np.ndarray,
    rolling_hs: np.ndarray,
    mixture_hs: np.ndarray,
    linear_hs: np.ndarray,
    tau_applied: np.ndarray,
    sigma_a_applied: np.ndarray,
    r_s_applied: np.ndarray,
    reference_points: Sequence[tuple[str, TuningPoint, str]],
    shaded_intervals: Sequence[tuple[float, float]],
    window_start_sec: float,
) -> None:
    """Draw the four-panel time-domain view of one transition realization.

    Top to bottom: what the sea does, what the wave and its estimate look
    like, what the resulting vertical error is, and what the tuner does about
    it against the frozen operating points.  Both the one-way validation
    protocol and the round-trip ablation protocol are drawn by this function
    so the two figures stay directly comparable.
    """

    plt = reproducible_pyplot()
    figure, axes = plt.subplots(4, 1, figsize=(9.0, 10.0), sharex=True)

    axes[0].plot(time, rolling_hs, color="#0072B2",
                 label=r"reference rolling $H_s$ (60 s)")
    axes[0].plot(time, mixture_hs, color="#D55E00",
                 linestyle="--", label=r"independent-mixture $H_s$")
    axes[0].plot(
        time, linear_hs,
        color="#999999", linestyle=":", label=r"linear $H_s$ ramp (not simulated)",
    )
    axes[0].set_ylabel(r"$H_s$ (m)")
    axes[0].legend(fontsize=7, ncol=2)

    axes[1].plot(time, reference, color="#000000",
                 linewidth=0.7, label="reference $p_z$")
    axes[1].plot(time, estimate, color="#009E73",
                 linewidth=0.7, label="estimated $p_z$")
    axes[1].set_ylabel("vertical displacement (m)")
    axes[1].legend(fontsize=7)

    axes[2].plot(time, estimate - reference, color="#CC79A7", linewidth=0.7)
    axes[2].set_ylabel("vertical error (m)")

    axes[3].plot(time, tau_applied, color="#0072B2",
                 label=r"applied $\tau$ (s)")
    axes[3].plot(time, sigma_a_applied, color="#D55E00",
                 label=r"applied $\sigma_{aw}$ (m/s$^2$)")
    axes[3].plot(time, r_s_applied, color="#009E73",
                 label=r"applied $r_S$ (m$\cdot$s)")
    for point_label, point, style in reference_points:
        for value, color, name in (
            (point.tau_s, "#0072B2", None),
            (point.sigma_a_mps2, "#D55E00", None),
            (point.RS_ms or point.R_p0_std_m, "#009E73", point_label),
        ):
            axes[3].axhline(
                value, color=color, linestyle=style, linewidth=0.9, label=name
            )
    axes[3].set_ylabel("applied tuning")
    axes[3].set_xlabel("time (s)")
    axes[3].legend(fontsize=7, ncol=2)

    for axis in axes:
        axis.grid(alpha=0.3)
        for shade_start, shade_end in shaded_intervals:
            axis.axvspan(shade_start, shade_end, color="#888888", alpha=0.12)
        axis.axvline(
            window_start_sec, color="#444444", linewidth=0.8, linestyle="-.",
        )
    figure.tight_layout()
    figure.savefig(svg_path, format="svg", metadata={"Date": None})
    plt.close(figure)


def write_transition_diagnostic(
    svg_path: Path,
    csv_path: Path,
    family: str,
    surrogate_path: Path,
    columns: Sequence[str],
    surrogate: np.ndarray,
    seed: SeedTriplet,
    window_sec: float,
    transition_start_sec: float,
    transition_end_sec: float,
    nominal_point: TuningPoint,
    oracle_point: TuningPoint,
    start_height_m: float,
    end_height_m: float,
    decimation: int = 20,
) -> dict[str, Any]:
    """Run one transition realization with time-series output and plot it.

    The aggregate tables cannot separate a genuine adaptation lag from the
    changing composition of the scoring window, so this figure shows a single
    realization directly: what the sea does, what the tuner does about it, and
    what the resulting vertical error looks like.
    """

    metrics, _, _ = run_simulator(
        family,
        surrogate_path,
        window_sec,
        seed.imu_noise_seed,
        seed.initialization_seed,
        tuning_mode="adaptive",
        write_timeseries=True,
    )
    series = read_diagnostic_timeseries(family, surrogate_path)

    time = series["time"]
    weight = smoothstep_weight(time, transition_start_sec, transition_end_sec)
    reference = series["disp_ref_z"]
    error = series["disp_est_z"] - reference

    # Rolling significant height of the reference record, as an independent
    # readout of what the crossfade actually does to the sea state.
    rolling_hs = rolling_significant_height(reference)
    mixture_hs = np.asarray(
        [
            mixture_significant_height_m(start_height_m, end_height_m, value)
            for value in weight
        ]
    )
    linear_hs = (1.0 - weight) * start_height_m + weight * end_height_m

    step = max(1, int(decimation))
    write_csv(
        csv_path,
        [
            {
                "time_s": float(time[index]),
                "blend_weight": float(weight[index]),
                "rolling_hs_m": float(rolling_hs[index]),
                "mixture_hs_m": float(mixture_hs[index]),
                "disp_ref_z_m": float(reference[index]),
                "disp_est_z_m": float(series["disp_est_z"][index]),
                "disp_err_z_m": float(error[index]),
                "tau_applied_s": float(series["tau_applied"][index]),
                "sigma_aw_applied_mps2": float(series["sigma_a_applied"][index]),
                "r_s_applied_ms": float(series["R_p0_applied"][index]),
                "freq_tracker_hz": float(series["freq_tracker_hz"][index]),
            }
            for index in range(0, time.size, step)
        ],
    )

    sliced = slice(None, None, step)
    plot_transition_diagnostic(
        svg_path,
        time[sliced],
        reference[sliced],
        series["disp_est_z"][sliced],
        rolling_hs[sliced],
        mixture_hs[sliced],
        linear_hs[sliced],
        series["tau_applied"][sliced],
        series["sigma_a_applied"][sliced],
        series["R_p0_applied"][sliced],
        (("FixedNominal", nominal_point, ":"), ("FixedOracle", oracle_point, "--")),
        ((transition_start_sec, transition_end_sec),),
        float(time[-1]) - window_sec,
    )

    return {
        "family": family,
        "wave_phase_seed": seed.wave_phase_seed,
        "imu_noise_seed": seed.imu_noise_seed,
        "initialization_seed": seed.initialization_seed,
        "decimation": step,
        "disp_z_pct_hs": float(metrics["disp_z_pct_hs"]),
        "midpoint_mixture_hs_m": mixture_significant_height_m(
            start_height_m, end_height_m, 0.5
        ),
        "midpoint_linear_hs_m": 0.5 * (start_height_m + end_height_m),
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_output(*arguments: str) -> str:
    completed = subprocess.run(
        ("git", *arguments),
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return completed.stdout.strip() if completed.returncode == 0 else "unknown"


def find_default_input(data_dir: Path, height_token: str, length_token: str) -> Path:
    matches = sorted(data_dir.glob(
        f"wave_data_jonswap_H{height_token}_L{length_token}_*.csv"
    ))
    if not matches:
        raise FileNotFoundError(
            f"no JONSWAP H{height_token}/L{length_token} data in {data_dir}; "
            "run `make fetch-sim-data`"
        )
    return matches[0]


def build_simulators(families: Sequence[str], eigen_dir: str | None) -> None:
    for family in families:
        command = ["make", "-B", "build"]
        if eigen_dir:
            command.append(f"EIGEN_DIR={eigen_dir}")
        subprocess.run(command, cwd=FAMILY_MAKE_DIR[family], check=True)


def _scenario_slug(path: Path) -> str:
    stem = path.stem.removeprefix("wave_data_")
    return "stationary_" + re.sub(r"[^A-Za-z0-9]+", "_", stem).strip("_")


def _copy_duration(path: Path, destination: Path, duration_sec: float) -> tuple[list[str], np.ndarray]:
    columns, data = read_wave_csv(path, duration_sec)
    write_wave_csv(destination, columns, data)
    return columns, data


NONSTATIONARY_ENDPOINT_NAME = "nonstationary_endpoint_H4_0_Tp11_4"


def restat_bundle(
    source: Path,
    output_dir: Path,
    bootstrap_resamples: int,
    stats_seed: int,
) -> int:
    """Recompute every statistic and table from an archived bundle's rows.

    The simulator replays are the expensive half of this study and are fully
    determined by their seed triplets, so a change to how the rows are
    *summarized* should not require re-running them.  This path reads
    `raw_runs` back out of a committed bundle and rewrites the derived files
    exactly as the run that produced it would have, with the statistics of the
    current source.  It cannot invent rows: whatever ensemble the source
    bundle scored is the ensemble the restated bundle reports.
    """

    restatement_context = evidence_provenance.begin_restatement(
        "validation", source
    )
    with source.open(encoding="utf-8") as stream:
        bundle = json.load(stream)
    raw_rows = [dict(row) for row in bundle["raw_runs"]]
    if not raw_rows:
        raise ValueError(f"{source} carries no raw runs to restate")

    # The bundle JSON is written with sorted keys, so the rows come back
    # alphabetical.  Column order is not information, but a restated bundle
    # that reshuffles 66 columns is impossible to diff against the run it
    # restates, so the archived CSV header is used to put them back.
    archived_csv = source.parent / "ou_validation_raw.csv"
    if archived_csv.exists():
        with archived_csv.open(encoding="utf-8", newline="") as stream:
            header = next(csv.reader(stream), [])
        if header:
            raw_rows = [
                {
                    **{key: row[key] for key in header if key in row},
                    **{key: value for key, value in row.items() if key not in header},
                }
                for row in raw_rows
            ]

    tuning_points: dict[str, dict[str, TuningPoint]] = {
        family: {
            name: TuningPoint(**{key: point.get(key) for key in (
                "tau_s", "sigma_a_mps2", "R_p0_std_m", "R_v0_std_mps", "RS_ms"
            )})
            for name, point in points.items()
        }
        for family, points in bundle.get("fixed_tuning_points", {}).items()
    }
    calibrated_names = {
        name for points in tuning_points.values() for name in points
    }
    scored_scenarios = {str(row["scenario"]) for row in raw_rows}
    # FixedNominal is calibrated from one record that is also scored as a
    # stationary scenario; the endpoint is calibrated but never scored.
    nominal_candidates = sorted(
        name for name in calibrated_names
        if name in scored_scenarios and "H1_500" in name
    )
    nominal_name = nominal_candidates[0] if nominal_candidates else ""

    output_dir.mkdir(parents=True, exist_ok=True)
    summary = summarize_rows(raw_rows, bootstrap_resamples, stats_seed)
    effects = paired_effect_rows(raw_rows, bootstrap_resamples, stats_seed)
    normalized_aggregate = stationary_normalized_aggregate(
        raw_rows, bootstrap_resamples, stats_seed
    )

    raw_path = output_dir / "ou_validation_raw.csv"
    summary_path = output_dir / "ou_validation_summary.csv"
    effects_path = output_dir / "ou_validation_paired_effects.csv"
    json_path = output_dir / "ou_validation.json"
    tex_path = output_dir / "ou_validation_table.tex"
    publication_tex_path = output_dir / "ou_validation_publication.tex"
    publication_macros_path = output_dir / "ou_validation_macros.tex"
    tuning_tex_path = output_dir / "ou_validation_tuning_points.tex"
    manifest_path = output_dir / "ou_validation_manifest.json"

    evidence_provenance.preserve_raw_rows(restatement_context, raw_path)
    write_csv(summary_path, summary)
    write_csv(effects_path, effects)
    write_latex_table(tex_path, summary)
    write_publication_table(
        publication_tex_path,
        raw_rows,
        summary,
        effects,
        bootstrap_resamples,
        stats_seed,
        macros_path=publication_macros_path,
    )
    if tuning_points:
        write_tuning_points_table(
            tuning_tex_path, tuning_points, nominal_name, NONSTATIONARY_ENDPOINT_NAME
        )

    # Fields that describe what was *run* are carried over untouched -- the
    # restated bundle scored exactly the ensemble the source did.  Fields that
    # describe how the rows are summarized or how the study words its own
    # design come from this file, so that a restated bundle cannot disagree
    # with the manuscript generated beside it.
    protocol = dict(bundle.get("protocol", {}))
    protocol["bootstrap_resamples"] = bootstrap_resamples
    protocol["stats_seed"] = stats_seed
    protocol["primary_endpoint_inference"] = PRIMARY_ENDPOINT_INFERENCE
    protocol["pmstokes_pooling"] = PMSTOKES_POOLING
    protocol["restated_from"] = str(
        source.relative_to(REPO_ROOT) if source.is_relative_to(REPO_ROOT) else source
    )
    write_json(
        json_path,
        {
            "protocol": protocol,
            "fixed_tuning_points": bundle.get("fixed_tuning_points", {}),
            "raw_runs": raw_rows,
            "summary": summary,
            "paired_effects": effects,
            "stationary_normalized_aggregate": normalized_aggregate,
            "transition_diagnostic": bundle.get("transition_diagnostic"),
        },
    )

    result_paths = [
        raw_path, summary_path, effects_path, json_path, tex_path,
        publication_tex_path, publication_macros_path,
    ]
    if tuning_points:
        result_paths.append(tuning_tex_path)
    manifest = dict(json.loads(manifest_path.read_text(encoding="utf-8")))    \
        if manifest_path.exists() else {}
    # A restat rewrites some of the bundle's files and leaves others (the
    # figures, the transition series) exactly as the run produced them.  The
    # manifest has to keep covering all of them, so the carried entries are
    # re-hashed from disk rather than trusted or dropped.
    result_files = {
        name: {"sha256": sha256_file(path), "bytes": path.stat().st_size}
        for name, path in (
            (name, output_dir / name)
            for name in manifest.get("result_files", {})
        )
        if path.is_file()
    }
    result_files.update(
        {
            path.name: {"sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in result_paths
        }
    )
    manifest.update(
        {
            "git_commit": git_output("rev-parse", "HEAD"),
            "git_branch": git_output("branch", "--show-current"),
            "git_diff_stat": git_output("diff", "--stat"),
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "command": [sys.executable, *sys.argv],
            "protocol": protocol,
            "restated_from": protocol["restated_from"],
            "result_files": result_files,
            "stationary_normalized_aggregate": normalized_aggregate,
        }
    )
    manifest = evidence_provenance.finalize_restatement_manifest(
        "validation", manifest, restatement_context
    )
    write_json(manifest_path, manifest)

    for path in (*result_paths, manifest_path):
        print(f"Wrote {path}")
    return 0


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("smoke", "full"), default="smoke")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=REPO_ROOT / "plots" / "kalman_ou_ii",
        help="directory containing oceanography-waves-lib release CSVs",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "reports" / "results" / "ou_validation",
    )
    parser.add_argument(
        "--stationary-input",
        action="append",
        type=Path,
        help="stationary CSV (repeatable); defaults to one smoke or four full JONSWAP seas",
    )
    parser.add_argument("--nonstationary-start", type=Path)
    parser.add_argument("--nonstationary-end", type=Path)
    parser.add_argument("--nonstationary-end-height", type=float, default=4.0)
    parser.add_argument("--skip-nonstationary", action="store_true")
    parser.add_argument("--duration-sec", type=float)
    parser.add_argument("--window-sec", type=float)
    parser.add_argument("--transition-start-sec", type=float)
    parser.add_argument("--transition-end-sec", type=float)
    parser.add_argument("--wave-seeds", type=parse_int_list)
    parser.add_argument("--imu-seeds", type=parse_int_list)
    parser.add_argument("--initialization-seeds", type=parse_int_list)
    parser.add_argument(
        "--families",
        choices=("both", "OU_II", "OU_III"),
        default="both",
    )
    parser.add_argument(
        "--adaptation-modes",
        default=",".join(MODE_SETTINGS),
        help="comma-separated subset of " + ",".join(MODE_SETTINGS),
    )
    parser.add_argument("--skip-pmstokes", action="store_true")
    parser.add_argument(
        "--pmstokes-modes",
        default=",".join(PRIMARY_MODES),
        help="adaptation modes scored on the PM-Stokes seas",
    )
    parser.add_argument("--bootstrap-resamples", type=int, default=10_000)
    parser.add_argument("--stats-seed", type=int, default=20260317)
    parser.add_argument(
        "--jobs",
        type=int,
        default=max(1, min(8, (os.cpu_count() or 1))),
        help="parallel (scenario, seed) work units; results are seed-determined "
        "and independent, so this changes runtime only",
    )
    parser.add_argument(
        "--shard-count",
        type=int,
        default=1,
        help="split the seed repetitions across this many independent runs",
    )
    parser.add_argument(
        "--shard-index",
        type=int,
        default=0,
        help="0-based index of this shard; runs its repetitions and writes a "
        "shard file instead of a bundle",
    )
    parser.add_argument(
        "--shard-dir",
        type=Path,
        help="where shard files are written and read back from "
        "(default: <output-dir>/shards)",
    )
    parser.add_argument(
        "--combine-shards",
        action="store_true",
        help="skip simulation, read every shard file, and write the bundle "
        "the unsharded run would have written",
    )
    parser.add_argument(
        "--restat-from",
        type=Path,
        help="skip simulation entirely and recompute the summaries, paired "
        "effects, tables, macros, and bundle from the raw runs archived in an "
        "existing ou_validation.json",
    )
    parser.add_argument("--eigen-dir")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    args.data_dir = args.data_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.restat_from is not None:
        return restat_bundle(
            args.restat_from.resolve(),
            args.output_dir,
            args.bootstrap_resamples,
            args.stats_seed,
        )

    families = ("OU_II", "OU_III") if args.families == "both" else (args.families,)
    allowed_modes = tuple(MODE_SETTINGS)
    adaptation_modes = tuple(
        token.strip() for token in args.adaptation_modes.split(",") if token.strip()
    )
    unknown_modes = sorted(set(adaptation_modes) - set(allowed_modes))
    if unknown_modes:
        raise ValueError(f"unknown adaptation mode(s): {', '.join(unknown_modes)}")
    if not adaptation_modes:
        raise ValueError("at least one adaptation mode is required")
    pmstokes_modes = tuple(
        token.strip() for token in args.pmstokes_modes.split(",") if token.strip()
    )
    unknown_pmstokes = sorted(set(pmstokes_modes) - set(allowed_modes))
    if unknown_pmstokes:
        raise ValueError(
            f"unknown PM-Stokes adaptation mode(s): {', '.join(unknown_pmstokes)}"
        )

    duration_sec = args.duration_sec or (180.0 if args.mode == "smoke" else 1200.0)
    window_sec = args.window_sec or (60.0 if args.mode == "smoke" else 900.0)
    if not (duration_sec > window_sec > 0.0):
        raise ValueError("duration must be greater than the positive score window")
    transition_start = args.transition_start_sec
    transition_end = args.transition_end_sec
    if transition_start is None:
        transition_start = TRANSITION_START_FRACTION * duration_sec
    if transition_end is None:
        transition_end = TRANSITION_END_FRACTION * duration_sec

    if args.mode == "smoke":
        default_wave, default_imu, default_init = ([11], [101], [1009])
    else:
        default_wave = list(DEFAULT_FULL_WAVE_SEEDS)
        default_imu = list(DEFAULT_FULL_IMU_SEEDS)
        default_init = list(DEFAULT_FULL_INIT_SEEDS)
    seeds = broadcast_seed_triplets(
        args.wave_seeds or default_wave,
        args.imu_seeds or default_imu,
        args.initialization_seeds or default_init,
    )

    nominal_input = args.nonstationary_start or find_default_input(
        args.data_dir, "1.500", "50.710"
    )
    end_input = args.nonstationary_end or find_default_input(
        args.data_dir, "8.500", "202.839"
    )
    nominal_input = nominal_input.resolve()
    end_input = end_input.resolve()

    if args.stationary_input:
        stationary_inputs = [path.resolve() for path in args.stationary_input]
    elif args.mode == "smoke":
        stationary_inputs = [nominal_input]
    else:
        stationary_inputs = sorted(args.data_dir.glob("wave_data_jonswap_*.csv"))
        if not args.skip_pmstokes:
            # PM-Stokes carries third-order bound harmonics that JONSWAP does
            # not, so it is a genuinely different input family rather than
            # another draw from the confirmatory ensemble.  It is scored with
            # the same seeds but kept out of the primary aggregate.
            stationary_inputs += sorted(
                args.data_dir.glob("wave_data_pmstokes_*.csv")
            )
    if not stationary_inputs:
        raise FileNotFoundError(f"no stationary inputs found in {args.data_dir}")

    scenarios = [
        Scenario(_scenario_slug(path), "stationary", path)
        for path in stationary_inputs
    ]
    if not args.skip_nonstationary:
        scenarios.append(
            Scenario(
                "nonstationary_H1_5_to_H4_0_Tp5_7_to_11_4",
                "nonstationary",
                nominal_input,
                end_input,
                args.nonstationary_end_height,
            )
        )

    if not args.skip_build:
        build_simulators(families, args.eigen_dir)

    shard_dir = args.shard_dir or (args.output_dir / "shards")
    shard_count = max(1, args.shard_count)
    if not 0 <= args.shard_index < shard_count:
        raise SystemExit(
            f"--shard-index must be in [0, {shard_count}), got {args.shard_index}"
        )

    raw_rows: list[dict[str, Any]] = []
    tuning_points: dict[str, dict[str, TuningPoint]] = defaultdict(dict)
    # One transition realization is retained after the run loop so the
    # time-domain diagnostic can be produced from the same surrogate that was
    # scored, rather than from a separately regenerated record.
    diagnostic_path: Path | None = None
    diagnostic_columns: list[str] = []
    diagnostic_data: np.ndarray | None = None
    diagnostic_start_height = 0.0
    transition_diagnostic: dict[str, Any] | None = None
    transition_svg_path: Path | None = None
    transition_csv_path: Path | None = None

    with tempfile.TemporaryDirectory(prefix="ocean-imu-ou-validation-") as temporary:
        work_dir = Path(temporary)

        # Calibrate each fixed operating point once from its noise-free true
        # trace. It is then held constant across every stochastic realization.
        calibration_sources: dict[str, Path] = {}
        calibration_data: dict[str, tuple[list[str], np.ndarray]] = {}
        for scenario in scenarios:
            if scenario.kind != "stationary":
                continue
            destination = work_dir / "calibration" / scenario.name / scenario.start_input.name
            columns, data = _copy_duration(scenario.start_input, destination, duration_sec)
            calibration_sources[scenario.name] = destination
            calibration_data[scenario.name] = (columns, data)

        nominal_name = _scenario_slug(nominal_input)
        if nominal_name not in calibration_sources:
            destination = work_dir / "calibration" / nominal_name / nominal_input.name
            columns, data = _copy_duration(nominal_input, destination, duration_sec)
            calibration_sources[nominal_name] = destination
            calibration_data[nominal_name] = (columns, data)

        endpoint_name = "nonstationary_endpoint_H4_0_Tp11_4"
        if any(scenario.kind == "nonstationary" for scenario in scenarios):
            end_columns, end_data = read_wave_csv(end_input, duration_sec)
            height_match = re.search(r"_H([0-9.]+)_", end_input.name)
            if not height_match:
                raise ValueError(f"cannot infer source height from {end_input.name}")
            source_height = float(height_match.group(1))
            end_scaled = scale_wave_motion(
                end_columns, end_data, args.nonstationary_end_height / source_height
            )
            endpoint_filename = (
                "wave_data_jonswap_H4.000_L202.839_A-30.00_P120.00.csv"
            )
            endpoint_path = work_dir / "calibration" / endpoint_name / endpoint_filename
            write_wave_csv(endpoint_path, end_columns, end_scaled)
            calibration_sources[endpoint_name] = endpoint_path
            calibration_data[endpoint_name] = (end_columns, end_scaled)

        for family in families:
            for name, path in calibration_sources.items():
                print(f"Calibrating {family} fixed point: {name}", flush=True)
                tuning_points[family][name] = calibrate_tuning_point(
                    family, path, window_sec
                )

        nominal_points = {
            family: tuning_points[family][nominal_name] for family in families
        }

        # Segments are only meaningful where the sea state changes inside the
        # scored window.  They are keyed to the same window the aggregate uses.
        window_start = duration_sec - window_sec
        # The recovery interval is one crossfade long, so it scales with the
        # transition it follows instead of being a constant somebody picked.
        recovery_end = min(
            duration_sec, transition_end + (transition_end - transition_start)
        )
        transition_segments = tuple(
            (name, lower, upper)
            for name, lower, upper in (
                (
                    name,
                    max(start, window_start),
                    min(stop, duration_sec),
                )
                for name, start, stop in (
                    ("start", window_start, transition_start),
                    ("blend", transition_start, transition_end),
                    ("recover", transition_end, recovery_end),
                    ("end", recovery_end, duration_sec),
                )
            )
            if upper > lower
        )

        def modes_for(scenario: Scenario, family: str) -> tuple[str, ...]:
            allowed = set(FAMILY_MODES[family])
            if scenario_spectrum(scenario.name) == "pmstokes":
                allowed &= set(pmstokes_modes)
            return tuple(m for m in adaptation_modes if m in allowed)

        def run_unit(
            scenario: Scenario,
            repetition: int,
            seed: SeedTriplet,
            start_columns: Sequence[str],
            start_data: np.ndarray,
            end_data: np.ndarray | None,
            end_scale: float,
            generated_path: Path,
            oracle_name: str,
        ) -> list[dict[str, Any]]:
            # The surrogate is built inside the work unit so that only the
            # in-flight realizations are resident; queueing them all up front
            # would hold one full-length record per pending unit.
            if scenario.kind == "stationary":
                generated = phase_randomize_wave(
                    start_columns, start_data, seed.wave_phase_seed
                )
            else:
                assert end_data is not None
                generated = make_nonstationary_wave(
                    start_columns,
                    start_data,
                    end_data,
                    seed.wave_phase_seed,
                    end_scale,
                    transition_start,
                    transition_end,
                )
            write_wave_csv(generated_path, start_columns, generated)
            del generated
            segments = (
                transition_segments if scenario.kind == "nonstationary" else ()
            )
            rows: list[dict[str, Any]] = []
            try:
                for family in families:
                    for mode in modes_for(scenario, family):
                        tuning_mode, aw_cov_sync = MODE_SETTINGS[mode]
                        if tuning_mode == "adaptive":
                            point = None
                        elif tuning_mode == "fixed_oracle":
                            point = tuning_points[family][oracle_name]
                        else:
                            # fixed_nominal and both partial-adaptation modes
                            # freeze their channel at the same nominal point,
                            # which is what makes them a clean factorial.
                            point = nominal_points[family]

                        metrics, gate_pass, return_code = run_simulator(
                            family,
                            generated_path,
                            window_sec,
                            seed.imu_noise_seed,
                            seed.initialization_seed,
                            tuning_mode=tuning_mode,
                            tuning_point=point,
                            aw_cov_sync=aw_cov_sync,
                            segments=segments,
                        )
                        rows.append(
                            {
                                "run_id": (
                                    f"{scenario.name}:{repetition}:{family}:{mode}"
                                ),
                                "scenario": scenario.name,
                                "scenario_kind": scenario.kind,
                                "spectrum": scenario_spectrum(scenario.name),
                                "family": family,
                                "mode": mode,
                                "repetition": repetition,
                                "wave_phase_seed": seed.wave_phase_seed,
                                "imu_noise_seed": seed.imu_noise_seed,
                                "initialization_seed": seed.initialization_seed,
                                "score_window_sec": window_sec,
                                **{
                                    key: value for key, value in metrics.items()
                                    if key not in (
                                        "family", "tuning_mode", "aw_cov_sync",
                                        "input",
                                    )
                                },
                            }
                        )
            finally:
                generated_path.unlink(missing_ok=True)
            return rows

        # Simulator runs are independent given their seeds, so (scenario, seed)
        # work units are dispatched to a small pool.  One scenario is in flight
        # at a time: its source records stay resident for the whole scenario,
        # and only the running units hold a generated surrogate.  Rows are
        # re-sorted afterwards so the committed artifacts do not depend on
        # completion order.
        diagnostic_recipe: tuple[Any, ...] | None = None
        for scenario in scenarios:
            # In combine mode the loop still runs, but only to rebuild the
            # transition diagnostic's recipe and the scenario bookkeeping the
            # bundle needs.  No simulator work is dispatched.
            print(f"Scenario: {scenario.name}", flush=True)
            start_columns, start_data = read_wave_csv(
                scenario.start_input, duration_sec
            )
            if scenario.kind == "nonstationary":
                assert scenario.end_input is not None
                end_columns, end_data = read_wave_csv(
                    scenario.end_input, duration_sec
                )
                if end_columns != start_columns:
                    raise ValueError("non-stationary source columns do not match")
                height_match = re.search(r"_H([0-9.]+)_", scenario.end_input.name)
                if not height_match or scenario.end_height_m is None:
                    raise ValueError("non-stationary endpoint height is unavailable")
                end_scale = scenario.end_height_m / float(height_match.group(1))
                oracle_name = endpoint_name
                generated_name = (
                    "wave_data_jonswap_H4.000_L202.839_A-30.00_P120.00.csv"
                )
            else:
                end_data = None
                end_scale = 1.0
                oracle_name = scenario.name
                generated_name = scenario.start_input.name

            if scenario.kind == "nonstationary":
                # Rebuilt after the loop rather than retained, so the diagnostic
                # does not pin a full record in memory for the whole study.  It
                # is seed-determined, so the rebuild is the scored realization.
                diagnostic_path = work_dir / "runs" / scenario.name / (
                    f"wave_{seeds[0].wave_phase_seed}"
                ) / generated_name
                height_source = re.search(
                    r"_H([0-9.]+)_", scenario.start_input.name
                )
                diagnostic_start_height = (
                    float(height_source.group(1)) if height_source else 0.0
                )
                diagnostic_recipe = (
                    list(start_columns),
                    scenario.start_input,
                    scenario.end_input,
                    end_scale,
                )

            pending = []
            with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
                for repetition, seed in enumerate(seeds, start=1):
                    # `repetition` stays the global index so pairing, sorting
                    # and the committed row all read the same whether or not
                    # the study was sharded.
                    if args.combine_shards:
                        continue
                    if (repetition - 1) % shard_count != args.shard_index:
                        continue
                    generated_path = work_dir / "runs" / scenario.name / (
                        f"wave_{seed.wave_phase_seed}"
                    ) / generated_name
                    pending.append(
                        pool.submit(
                            run_unit,
                            scenario,
                            repetition,
                            seed,
                            list(start_columns),
                            start_data,
                            end_data,
                            end_scale,
                            generated_path,
                            oracle_name,
                        )
                    )
                for index, future in enumerate(pending, start=1):
                    raw_rows.extend(future.result())
                    print(
                        f"  {scenario.name}: {index}/{len(pending)} seeds",
                        flush=True,
                    )
            del start_data, end_data

        if diagnostic_path is not None and diagnostic_recipe is not None:
            diagnostic_columns, start_path, end_path, end_scale = diagnostic_recipe
            _, start_data = read_wave_csv(start_path, duration_sec)
            assert end_path is not None
            _, end_data = read_wave_csv(end_path, duration_sec)
            diagnostic_data = make_nonstationary_wave(
                diagnostic_columns,
                start_data,
                end_data,
                seeds[0].wave_phase_seed,
                end_scale,
                transition_start,
                transition_end,
            )
            write_wave_csv(diagnostic_path, diagnostic_columns, diagnostic_data)
            del start_data, end_data

        if args.combine_shards:
            raw_rows = read_shards(shard_dir, SHARD_PREFIX, ordered=False)
            print(f"Combined {len(raw_rows)} runs from {shard_dir}", flush=True)
        elif shard_count > 1:
            path = shard_path(shard_dir, SHARD_PREFIX, args.shard_index)
            write_shard(path, raw_rows)
            print(f"Wrote {path} ({len(raw_rows)} runs)", flush=True)
            return 0

        # This study already sorts its rows into an order that does not depend
        # on completion order, so a sharded merge needs no position bookkeeping
        # -- the same sort puts the merged rows exactly where one process would
        # have left them.
        raw_rows.sort(
            key=lambda row: (
                scenario_sort_key(str(row["scenario"])),
                int(row["repetition"]),
                str(row["family"]),
                str(row["mode"]),
            )
        )

        if diagnostic_path is not None and not args.no_plots:
            print("Transition diagnostic run", flush=True)
            transition_svg_path = args.output_dir / "ou_validation_transition.svg"
            transition_csv_path = args.output_dir / "ou_validation_transition.csv"
            diagnostic_family = "OU_III" if "OU_III" in families else families[0]
            transition_diagnostic = write_transition_diagnostic(
                transition_svg_path,
                transition_csv_path,
                diagnostic_family,
                diagnostic_path,
                diagnostic_columns,
                diagnostic_data,
                seeds[0],
                window_sec,
                transition_start,
                transition_end,
                nominal_points[diagnostic_family],
                tuning_points[diagnostic_family][endpoint_name],
                diagnostic_start_height,
                args.nonstationary_end_height,
            )
            diagnostic_path.unlink(missing_ok=True)

    summary = summarize_rows(raw_rows, args.bootstrap_resamples, args.stats_seed)
    effects = paired_effect_rows(raw_rows, args.bootstrap_resamples, args.stats_seed)

    raw_path = args.output_dir / "ou_validation_raw.csv"
    summary_path = args.output_dir / "ou_validation_summary.csv"
    effects_path = args.output_dir / "ou_validation_paired_effects.csv"
    json_path = args.output_dir / "ou_validation.json"
    tex_path = args.output_dir / "ou_validation_table.tex"
    publication_tex_path = args.output_dir / "ou_validation_publication.tex"
    publication_macros_path = args.output_dir / "ou_validation_macros.tex"
    tuning_tex_path = args.output_dir / "ou_validation_tuning_points.tex"
    manifest_path = args.output_dir / "ou_validation_manifest.json"
    write_csv(raw_path, raw_rows)
    write_csv(summary_path, summary)
    write_csv(effects_path, effects)
    write_latex_table(tex_path, summary)
    write_publication_table(
        publication_tex_path,
        raw_rows,
        summary,
        effects,
        args.bootstrap_resamples,
        args.stats_seed,
        macros_path=publication_macros_path,
    )
    write_tuning_points_table(
        tuning_tex_path, tuning_points, nominal_name, endpoint_name
    )

    calibration_json = {
        family: {name: asdict(point) for name, point in points.items()}
        for family, points in tuning_points.items()
    }
    normalized_aggregate = stationary_normalized_aggregate(
        raw_rows, args.bootstrap_resamples, args.stats_seed
    )
    protocol = {
        "mode": args.mode,
        "duration_sec": duration_sec,
        "score_window_sec": window_sec,
        "quality_gate_window_sec": SIMULATOR_GATE_WINDOW_SEC,
        "transition_start_sec": transition_start,
        "transition_end_sec": transition_end,
        "wave_phase_method": WAVE_PHASE_METHOD,
        "transition_method": TRANSITION_METHOD,
        "pairing": "identical wave, IMU-noise, and initialization seeds",
        "seed_triplets": [asdict(seed) for seed in seeds],
        "adaptation_modes": list(adaptation_modes),
        "fixed_nominal_definition": (
            "noise-free full-trace operating point from Hs=1.5 m, L=50.710 m"
        ),
        "fixed_oracle_definition": (
            "scenario-calibrated fixed reference: the noise-free full-trace "
            "operating point of the stationary sea being scored, and the known "
            "final Hs=4.0 m/Tp=11.4 s point for the transition; it is selected "
            "from knowledge of the sea, not optimized against displacement error"
        ),
        "aw_covariance_sync_policies": {
            mode: policy for mode, (_, policy) in MODE_SETTINGS.items()
            if mode in adaptation_modes
        },
        "transition_window_composition_sec": transition_window_composition(
            transition_start, transition_end, duration_sec, window_sec
        ),
        "transition_midpoint_mixture_hs_m": mixture_significant_height_m(
            1.5, args.nonstationary_end_height, 0.5
        ),
        "transition_midpoint_linear_hs_m": 0.5 * (
            1.5 + args.nonstationary_end_height
        ),
        "partial_adaptation_definition": (
            "AdaptiveRSOnly freezes tau and sigma_aw at the FixedNominal point "
            "while r_S keeps adapting; AdaptiveOUOnly freezes r_S at the "
            "FixedNominal point while tau and sigma_aw keep adapting.  r_S is "
            "derived from the live tau and sigma_aw estimates in both cases, so "
            "freezing the OU channel does not implicitly freeze r_S through the "
            "coupled law.  OU-III only"
        ),
        "primary_endpoint": (
            "seed-level mean normalized vertical displacement RMS error over "
            "the four stationary JONSWAP seas, OU-III Adaptive minus OU-II "
            "Adaptive"
        ),
        "primary_endpoint_inference": PRIMARY_ENDPOINT_INFERENCE,
        "pmstokes_modes": list(pmstokes_modes),
        "pmstokes_pooling": PMSTOKES_POOLING,
        "transition_segments_sec": {
            name: [start, stop] for name, start, stop in transition_segments
        },
        "bootstrap_method": (
            "nonparametric percentile bootstrap of the mean; paired "
            "comparisons resample the seed-level differences with replacement"
        ),
        "bootstrap_rng": "numpy.random.default_rng (PCG64)",
        "bootstrap_resamples": args.bootstrap_resamples,
        "stats_seed": args.stats_seed,
        "multiplicity": (
            "intervals are per-comparison and are not adjusted for "
            "multiplicity; only the primary endpoint is confirmatory and every "
            "other interval is descriptive"
        ),
    }
    write_json(
        json_path,
        {
            "protocol": protocol,
            "fixed_tuning_points": calibration_json,
            "raw_runs": raw_rows,
            "summary": summary,
            "paired_effects": effects,
            "stationary_normalized_aggregate": normalized_aggregate,
            "transition_diagnostic": transition_diagnostic,
        },
    )

    plot_paths: list[Path] = []
    if not args.no_plots:
        plot_paths = [
            args.output_dir / "ou_validation_displacement.svg",
            args.output_dir / "ou_validation_vertical.svg",
            args.output_dir / "ou_validation_attitude.svg",
        ]
        write_metric_plot(
            plot_paths[0],
            summary,
            "disp_3d_rms_m",
            "3D displacement RMS (m)",
        )
        write_metric_plot(
            plot_paths[1],
            summary,
            "disp_z_pct_hs",
            r"Vertical displacement RMS (% $H_s$)",
        )
        write_metric_plot(
            plot_paths[2],
            summary,
            "pitch_rms_deg",
            "Pitch RMS (deg)",
        )
        if all(
            any(row["mode"] == mode for row in summary)
            for pair in COVARIANCE_SYNC_PAIRS for mode in pair
        ):
            covsync_path = args.output_dir / "ou_validation_covsync.svg"
            write_metric_plot(
                covsync_path,
                summary,
                "disp_z_pct_hs",
                r"Vertical displacement RMS (% $H_s$)",
                modes=[mode for pair in COVARIANCE_SYNC_PAIRS for mode in pair],
            )
            plot_paths.append(covsync_path)
        if transition_svg_path is not None and transition_svg_path.exists():
            plot_paths.append(transition_svg_path)

    source_paths = sorted({scenario.start_input for scenario in scenarios}.union(
        {scenario.end_input for scenario in scenarios if scenario.end_input is not None}
    ))
    result_paths = [
        raw_path,
        summary_path,
        effects_path,
        json_path,
        tex_path,
        publication_tex_path,
        publication_macros_path,
        tuning_tex_path,
        *plot_paths,
    ]
    if transition_csv_path is not None and transition_csv_path.exists():
        result_paths.append(transition_csv_path)
    manifest = {
        "git_commit": git_output("rev-parse", "HEAD"),
        "git_branch": git_output("branch", "--show-current"),
        "git_diff_stat": git_output("diff", "--stat"),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "build_environment": evidence_provenance.environment_metadata(),
        "command": [sys.executable, *sys.argv],
        "protocol": protocol,
        "source_files": {
            str(
                path.relative_to(REPO_ROOT)
                if path.is_relative_to(REPO_ROOT)
                else path
            ): {
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for path in source_paths
        },
        "result_files": {
            path.name: {"sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in result_paths
        },
        "fixed_tuning_points": calibration_json,
        "stationary_normalized_aggregate": normalized_aggregate,
        "transition_diagnostic": transition_diagnostic,
    }
    write_json(manifest_path, manifest)

    print(f"Wrote {raw_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {effects_path}")
    print(f"Wrote {json_path}")
    print(f"Wrote {tex_path}")
    print(f"Wrote {publication_tex_path}")
    print(f"Wrote {publication_macros_path}")
    print(f"Wrote {tuning_tex_path}")
    print(f"Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
