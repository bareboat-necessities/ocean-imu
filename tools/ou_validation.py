#!/usr/bin/env python3
"""Paired Monte Carlo validation for the OU-II and OU-III filters.

The runner deliberately keeps the historical final-60-second executable gates
separate from its inferential score.  Validation uses a configurable long
window (900 seconds in full mode), independent wave-realization, IMU-noise,
and initialization seeds, and paired realizations across filters and ablations.
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
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
DT_SECONDS = 1.0 / 200.0
GRAVITY_MPS2 = 9.80665
SURROGATE_MIN_FREQ_HZ = 0.02
SURROGATE_MAX_FREQ_HZ = 1.60
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
}

PRIMARY_MODES = ("Adaptive", "FixedNominal", "FixedOracle")

# Ablation pairs compared at matched tuning mode: (left, right).
COVARIANCE_SYNC_PAIRS = (
    ("Adaptive", "AdaptiveHeldCovariance"),
    ("FixedNominal", "FixedNominalHeldCovariance"),
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
    "roll_rms_deg",
    "pitch_rms_deg",
    "yaw_rms_deg",
    "accel_bias_3d_rms_mps2",
    "gyro_bias_3d_rms_radps",
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


def parse_validation_metrics(output: str) -> dict[str, Any]:
    lines = [line for line in output.splitlines() if line.startswith("VALIDATION_METRICS ")]
    if not lines:
        raise ValueError("simulator did not emit VALIDATION_METRICS")
    metrics: dict[str, Any] = {}
    for token in lines[-1].split()[1:]:
        key, value = token.split("=", 1)
        if key == "samples":
            metrics[key] = int(value)
        elif key in ("family", "tuning_mode", "aw_cov_sync", "input"):
            metrics[key] = value
        else:
            metrics[key] = float(value)
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


def tuning_point_from_pilot(family: str, metrics: Mapping[str, Any]) -> TuningPoint:
    tau = float(metrics["tau_applied_s"])
    sigma = float(metrics["sigma_applied_mps2"])
    if family == "OU_II":
        R_p0 = min(max(1.6 * sigma * tau * tau, 0.05), 18.0)
        R_v0 = min(max(1.4 * sigma * tau, 0.01), 6.0)
        return TuningPoint(tau, sigma, R_p0_std_m=R_p0, R_v0_std_mps=R_v0)
    RS = min(max(1.2 * sigma * tau * tau * tau, 0.4), 35.0)
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
    values = np.asarray([float(row[metric]) for row in rows], dtype=np.float64)
    return values[np.isfinite(values)]


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
            if math.isfinite(float(row[metric]))
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
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
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


def scenario_display_label(scenario: str) -> str:
    if scenario.startswith("nonstationary_"):
        return "Transition"
    height = re.search(r"_H([0-9]+)_([0-9]{3})_", scenario)
    if height:
        value = float(f"{height.group(1)}.{height.group(2)}")
        return rf"$H_s={value:.2f}$ m"
    return scenario.replace("_", " ")


def scenario_sort_key(scenario: str) -> tuple[int, float, str]:
    if scenario.startswith("nonstationary_"):
        return 1, math.inf, scenario
    height = re.search(r"_H([0-9]+)_([0-9]{3})_", scenario)
    value = (
        float(f"{height.group(1)}.{height.group(2)}")
        if height else math.inf
    )
    return 0, value, scenario


def stationary_normalized_aggregate(
    rows: Sequence[Mapping[str, Any]],
    bootstrap_resamples: int,
    stats_seed: int,
) -> dict[str, dict[str, Any]]:
    by_seed_family: dict[
        tuple[Any, Any, Any, str], list[float]
    ] = defaultdict(list)
    for row in rows:
        if (
            not str(row["scenario"]).startswith("stationary_")
            or row["mode"] != "Adaptive"
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
    }
    return result


def write_publication_table(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    summary: Sequence[Mapping[str, Any]],
    effects: Sequence[Mapping[str, Any]],
    bootstrap_resamples: int,
    stats_seed: int,
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
    scenarios = sorted(
        {
            str(row["scenario"])
            for row in summary
            if row["mode"] == "Adaptive" and row["metric"] == "disp_z_pct_hs"
        },
        key=scenario_sort_key,
    )
    aggregate = stationary_normalized_aggregate(
        rows, bootstrap_resamples, stats_seed
    )
    difference = aggregate["OU_III_minus_OU_II"]
    gate_passes = sum(int(row["historical_60s_gate_pass"]) for row in rows)

    def mean_std(scenario: str, family: str, mode: str) -> str:
        row = indexed_summary[(scenario, family, mode, "disp_z_pct_hs")]
        return f"{float(row['mean']):.2f} $\\pm$ {float(row['std']):.2f}"

    def paired_ci(scenario: str, metric: str, digits: int) -> str:
        row = indexed_effects[
            (scenario, "OU_III_vs_OU_II", "OU_III_minus_OU_II", metric)
        ]
        return (
            f"{float(row['mean_paired_difference']):+.{digits}f} "
            f"[{float(row['bootstrap_ci95_low']):+.{digits}f}, "
            f"{float(row['bootstrap_ci95_high']):+.{digits}f}]"
        )

    def ablation_effect(
        scenario: str, family: str, baseline: str, metric: str = "disp_z_pct_hs"
    ) -> Mapping[str, Any]:
        return indexed_effects[
            (scenario, family, f"Adaptive_minus_{baseline}", metric)
        ]

    def signed(value: float, digits: int = 3) -> str:
        return f"{value:+.{digits}f}"

    scenario_macro = {}
    for scenario in scenarios:
        if scenario.startswith("nonstationary_"):
            scenario_macro[scenario] = "Transition"
            continue
        height = re.search(r"_H([0-9]+)_([0-9]{3})_", scenario)
        if height:
            token = f"{height.group(1)}{height.group(2)}".rstrip("0") or "0"
            scenario_macro[scenario] = f"Hs{token}"

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
        rf"\providecommand{{\OUValidationNormalizedDz}}{{{difference['cohen_dz']:+.2f}}}",
        rf"\providecommand{{\OUValidationNormalizedGz}}{{{difference['hedges_gz']:+.2f}}}",
        rf"\providecommand{{\OUValidationLegacyGatePasses}}{{{gate_passes}}}",
        rf"\providecommand{{\OUValidationLegacyGateFailures}}{{{len(rows) - gate_passes}}}",
        r"",
        r"\begin{table*}[t]",
        r"  \centering",
        r"  \caption{Ten-seed paired OU-family comparison over the final \SI{900}{s}. Values are mean $\pm$ sample standard deviation. $\Delta$ is OU--III Adaptive minus OU--II Adaptive; brackets give the paired bootstrap 95\% interval. Negative $\Delta$Z favors OU--III, whereas positive $\Delta$3D indicates larger OU--III three-dimensional displacement error.}",
        r"  \label{tab:ou_mc_family}",
        r"  \footnotesize",
        r"  \setlength{\tabcolsep}{4.0pt}",
        r"  \begin{tabular}{@{}lrrrr@{}}",
        r"    \toprule",
        r"    Scenario & OU--II Z [\%$H_s$] & OU--III Z [\%$H_s$] & $\Delta$Z [percentage points] & $\Delta$3D [m] \\",
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
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


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
        r"  \caption{Frozen operating points for the fixed-tuning modes, calibrated once from each noise-free \SI{1200}{s} reference record. $\dagger$ marks the single point used by FixedNominal in \emph{every} scenario; FixedOracle uses the row matching its own scenario, and the transition-endpoint row for the non-stationary case. Values are the vertical/base parameters; the filter applies $(1.87\sigma_{aw},1.87\sigma_{aw},\sigma_{aw})$ for OU--III acceleration and $\operatorname{diag}(0.36 r_S,0.36 r_S,r_S)^2$ for the integral pseudo-measurement, and $1.5\sigma_{aw}$ horizontally with $0.31 r_{p0}$ for OU--II.}",
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
    axis.set_xticks(x, [scenario_display_label(scenario) for scenario in scenarios])
    axis.set_xlabel("Scenario")
    axis.set_ylabel(ylabel)
    axis.grid(axis="y", alpha=0.3)
    axis.legend(fontsize=8, ncol=2)
    figure.tight_layout()
    figure.savefig(path, format="svg")
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
    """

    window_start = duration_sec - window_sec
    start_end = max(window_start, min(transition_start_sec, duration_sec))
    blend_end = max(start_end, min(transition_end_sec, duration_sec))
    return {
        "window_start_sec": window_start,
        "window_end_sec": duration_sec,
        "pure_start_sea_sec": start_end - window_start,
        "blended_sec": blend_end - start_end,
        "pure_end_sea_sec": duration_sec - blend_end,
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

    cache_dir = Path(tempfile.gettempdir()) / "ocean-imu-matplotlib-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    metrics, _, _ = run_simulator(
        family,
        surrogate_path,
        window_sec,
        seed.imu_noise_seed,
        seed.initialization_seed,
        tuning_mode="adaptive",
        write_timeseries=True,
    )
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
        series = _read_timeseries_columns(output_path, wanted)
    finally:
        output_path.unlink(missing_ok=True)

    time = series["time"]
    weight = smoothstep_weight(time, transition_start_sec, transition_end_sec)
    error = series["disp_est_z"] - series["disp_ref_z"]

    # Rolling significant height of the reference record, as an independent
    # readout of what the crossfade actually does to the sea state.
    window_samples = max(2, int(round(60.0 / DT_SECONDS)))
    kernel = np.ones(window_samples) / window_samples
    reference = series["disp_ref_z"]
    rolling_mean = np.convolve(reference, kernel, mode="same")
    rolling_mean_square = np.convolve(reference**2, kernel, mode="same")
    rolling_hs = 4.0 * np.sqrt(
        np.maximum(rolling_mean_square - rolling_mean**2, 0.0)
    )
    mixture_hs = np.asarray(
        [
            mixture_significant_height_m(start_height_m, end_height_m, value)
            for value in weight
        ]
    )

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
    figure, axes = plt.subplots(4, 1, figsize=(9.0, 10.0), sharex=True)

    axes[0].plot(time[sliced], rolling_hs[sliced], color="#0072B2",
                 label=r"reference rolling $H_s$ (60 s)")
    axes[0].plot(time[sliced], mixture_hs[sliced], color="#D55E00",
                 linestyle="--", label=r"independent-mixture $H_s$")
    axes[0].plot(
        time[sliced],
        (1.0 - weight[sliced]) * start_height_m + weight[sliced] * end_height_m,
        color="#999999", linestyle=":", label=r"linear $H_s$ ramp (not simulated)",
    )
    axes[0].set_ylabel(r"$H_s$ (m)")
    axes[0].legend(fontsize=7, ncol=2)

    axes[1].plot(time[sliced], reference[sliced], color="#000000",
                 linewidth=0.7, label="reference $p_z$")
    axes[1].plot(time[sliced], series["disp_est_z"][sliced], color="#009E73",
                 linewidth=0.7, label="estimated $p_z$")
    axes[1].set_ylabel("vertical displacement (m)")
    axes[1].legend(fontsize=7)

    axes[2].plot(time[sliced], error[sliced], color="#CC79A7", linewidth=0.7)
    axes[2].set_ylabel("vertical error (m)")

    axes[3].plot(time[sliced], series["tau_applied"][sliced], color="#0072B2",
                 label=r"applied $\tau$ (s)")
    axes[3].plot(time[sliced], series["sigma_a_applied"][sliced], color="#D55E00",
                 label=r"applied $\sigma_{aw}$ (m/s$^2$)")
    axes[3].plot(time[sliced], series["R_p0_applied"][sliced], color="#009E73",
                 label=r"applied $r_S$ (m$\cdot$s)")
    for value, color, name in (
        (nominal_point.tau_s, "#0072B2", None),
        (nominal_point.sigma_a_mps2, "#D55E00", None),
        (nominal_point.RS_ms or nominal_point.R_p0_std_m, "#009E73", "FixedNominal"),
    ):
        axes[3].axhline(value, color=color, linestyle=":", linewidth=0.9, label=name)
    for value, color, name in (
        (oracle_point.tau_s, "#0072B2", None),
        (oracle_point.sigma_a_mps2, "#D55E00", None),
        (oracle_point.RS_ms or oracle_point.R_p0_std_m, "#009E73", "FixedOracle"),
    ):
        axes[3].axhline(value, color=color, linestyle="--", linewidth=0.9, label=name)
    axes[3].set_ylabel("applied tuning")
    axes[3].set_xlabel("time (s)")
    axes[3].legend(fontsize=7, ncol=2)

    for axis in axes:
        axis.grid(alpha=0.3)
        axis.axvspan(
            transition_start_sec, transition_end_sec, color="#888888", alpha=0.12
        )
        axis.axvline(
            float(time[-1]) - window_sec, color="#444444", linewidth=0.8,
            linestyle="-.",
        )
    figure.tight_layout()
    figure.savefig(svg_path, format="svg")
    plt.close(figure)

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
    parser.add_argument("--bootstrap-resamples", type=int, default=10_000)
    parser.add_argument("--stats-seed", type=int, default=20260317)
    parser.add_argument("--eigen-dir")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    args.data_dir = args.data_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

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

    duration_sec = args.duration_sec or (180.0 if args.mode == "smoke" else 1200.0)
    window_sec = args.window_sec or (60.0 if args.mode == "smoke" else 900.0)
    if not (duration_sec > window_sec > 0.0):
        raise ValueError("duration must be greater than the positive score window")
    transition_start = args.transition_start_sec
    transition_end = args.transition_end_sec
    if transition_start is None:
        transition_start = 0.35 * duration_sec
    if transition_end is None:
        transition_end = 0.65 * duration_sec

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

        for scenario in scenarios:
            print(f"Scenario: {scenario.name}", flush=True)
            start_columns, start_data = read_wave_csv(scenario.start_input, duration_sec)
            if scenario.kind == "nonstationary":
                assert scenario.end_input is not None
                end_columns, end_data = read_wave_csv(scenario.end_input, duration_sec)
                if end_columns != start_columns:
                    raise ValueError("non-stationary source columns do not match")
                height_match = re.search(r"_H([0-9.]+)_", scenario.end_input.name)
                if not height_match or scenario.end_height_m is None:
                    raise ValueError("non-stationary endpoint height is unavailable")
                end_scale = scenario.end_height_m / float(height_match.group(1))
                oracle_name = endpoint_name
                generated_name = "wave_data_jonswap_H4.000_L202.839_A-30.00_P120.00.csv"
            else:
                end_data = None
                end_scale = 1.0
                oracle_name = scenario.name
                generated_name = scenario.start_input.name

            for repetition, seed in enumerate(seeds, start=1):
                generated_path = work_dir / "runs" / scenario.name / (
                    f"wave_{seed.wave_phase_seed}"
                ) / generated_name
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
                if scenario.kind == "nonstationary" and repetition == 1:
                    diagnostic_path = generated_path
                    diagnostic_columns = list(start_columns)
                    diagnostic_data = generated
                    height_source = re.search(
                        r"_H([0-9.]+)_", scenario.start_input.name
                    )
                    diagnostic_start_height = (
                        float(height_source.group(1)) if height_source else 0.0
                    )

                for family in families:
                    for mode in adaptation_modes:
                        tuning_mode, aw_cov_sync = MODE_SETTINGS[mode]
                        if tuning_mode == "adaptive":
                            point = None
                        elif tuning_mode == "fixed_nominal":
                            point = nominal_points[family]
                        else:
                            point = tuning_points[family][oracle_name]

                        print(
                            f"  rep={repetition}/{len(seeds)} {family} {mode}",
                            flush=True,
                        )
                        metrics, gate_pass, return_code = run_simulator(
                            family,
                            generated_path,
                            window_sec,
                            seed.imu_noise_seed,
                            seed.initialization_seed,
                            tuning_mode=tuning_mode,
                            tuning_point=point,
                            aw_cov_sync=aw_cov_sync,
                        )
                        raw_rows.append(
                            {
                                "run_id": (
                                    f"{scenario.name}:{repetition}:{family}:{mode}"
                                ),
                                "scenario": scenario.name,
                                "scenario_kind": scenario.kind,
                                "family": family,
                                "mode": mode,
                                "repetition": repetition,
                                "wave_phase_seed": seed.wave_phase_seed,
                                "imu_noise_seed": seed.imu_noise_seed,
                                "initialization_seed": seed.initialization_seed,
                                "score_window_sec": window_sec,
                                "historical_60s_gate_pass": int(gate_pass),
                                "simulator_return_code": return_code,
                                **{
                                    key: value for key, value in metrics.items()
                                    if key not in (
                                        "family", "tuning_mode", "aw_cov_sync",
                                        "input",
                                    )
                                },
                            }
                        )
                if generated_path == diagnostic_path:
                    continue
                generated_path.unlink()

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
        "historical_gate_window_sec": 60.0,
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
        "bootstrap_resamples": args.bootstrap_resamples,
        "stats_seed": args.stats_seed,
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
    print(f"Wrote {tuning_tex_path}")
    print(f"Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
