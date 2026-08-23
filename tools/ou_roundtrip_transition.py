#!/usr/bin/env python3
"""Generate the OU-III deployed-law bidirectional transition evidence.

The instrument is strictly the low--high--low transition used by the OU-III
article. It always runs the deployed SpectralMSE configuration through the
normal adaptive simulator path; it does not sweep or compare regularizer laws.
"""

from __future__ import annotations

import argparse
import csv
import sys
import tempfile
from pathlib import Path
from typing import Any, Sequence

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
OU_III_DATA = REPO_ROOT / "tests" / "kalman_ou_iii"

TRANSITION_DURATION_SEC = 120.0
DEFAULT_TRANSITION_START_SEC = 400.0
DEFAULT_TRANSITION_RETURN_START_SEC = 800.0
ROUNDTRIP_DURATION_SEC = 1200.0
ROUNDTRIP_WINDOW_SEC = 900.0
ROUNDTRIP_WAVE_NAME = "wave_data_jonswap_H4.000_L202.839_A30.00_P120.00.csv"
DIAGNOSTIC_LOW_HEIGHT_M = 1.5
DIAGNOSTIC_HIGH_HEIGHT_M = 4.0
DIAGNOSTIC_SOURCE_HEIGHT_M = 8.5
SEGMENT_METRICS = (
    "disp_z_rms_m",
    "disp_z_ref_rms_m",
    "disp_z_pct_refrms",
    "disp_z_pct_hs",
    "disp_3d_rms_m",
)
SEGMENT_LABELS = {
    "low_start": "Low start",
    "rise": "Rise crossfade",
    "high_recover": "High recovery",
    "high": "High settled",
    "fall": "Fall crossfade",
    "low_recover": "Low recovery",
    "low_return": "Low return",
}


def _validation_module():
    tools = str(REPO_ROOT / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    import ou_validation as ov
    return ov


def transition_bounds(
    start_sec: float = DEFAULT_TRANSITION_START_SEC,
    return_start_sec: float = DEFAULT_TRANSITION_RETURN_START_SEC,
) -> tuple[float, float, float, float]:
    end_sec = start_sec + TRANSITION_DURATION_SEC
    return_end_sec = return_start_sec + TRANSITION_DURATION_SEC
    if not (
        0.0 <= start_sec < end_sec <= return_start_sec
        < return_end_sec <= ROUNDTRIP_DURATION_SEC
    ):
        raise ValueError(
            "round-trip transition requires 0 <= rise < high hold < fall <= record"
        )
    return start_sec, end_sec, return_start_sec, return_end_sec


def roundtrip_profile(
    ov,
    times: np.ndarray,
    start_sec: float = DEFAULT_TRANSITION_START_SEC,
    return_start_sec: float = DEFAULT_TRANSITION_RETURN_START_SEC,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    start, end, return_start, return_end = transition_bounds(
        start_sec, return_start_sec
    )
    up, up_rate, up_accel = ov.smoothstep_profile(times, start, end)
    down, down_rate, down_accel = ov.smoothstep_profile(
        times, return_start, return_end
    )
    return up - down, up_rate - down_rate, up_accel - down_accel


def make_roundtrip_wave(
    ov,
    columns: Sequence[str],
    low_data: np.ndarray,
    high_data: np.ndarray,
    seed: int,
    high_scale: float,
    start_sec: float = DEFAULT_TRANSITION_START_SEC,
    return_start_sec: float = DEFAULT_TRANSITION_RETURN_START_SEC,
) -> np.ndarray:
    """Build a kinematically closed low->high->same-low realization."""
    count = min(low_data.shape[0], high_data.shape[0])
    low = ov.phase_randomize_wave(columns, low_data[:count], seed * 2 + 1)
    high = ov.phase_randomize_wave(columns, high_data[:count], seed * 2 + 2)
    high = ov.scale_wave_motion(columns, high, high_scale)
    times = low[:, columns.index("time")]
    weight, weight_rate, weight_acceleration = roundtrip_profile(
        ov, times, start_sec, return_start_sec
    )

    def indices(names: Sequence[str]) -> list[int]:
        return [columns.index(name) for name in names]

    displacement = indices(("disp_x", "disp_y", "disp_z"))
    velocity = indices(("vel_x", "vel_y", "vel_z"))
    acceleration = indices(("acc_x", "acc_y", "acc_z"))
    attitude = indices(("roll_deg", "pitch_deg", "yaw_deg"))

    result = low.copy()
    w = weight[:, None]
    dw = weight_rate[:, None]
    ddw = weight_acceleration[:, None]
    displacement_delta = high[:, displacement] - low[:, displacement]
    velocity_delta = high[:, velocity] - low[:, velocity]

    result[:, displacement] = (
        (1.0 - w) * low[:, displacement] + w * high[:, displacement]
    )
    result[:, velocity] = (
        (1.0 - w) * low[:, velocity]
        + w * high[:, velocity]
        + dw * displacement_delta
    )
    result[:, acceleration] = (
        (1.0 - w) * low[:, acceleration]
        + w * high[:, acceleration]
        + 2.0 * dw * velocity_delta
        + ddw * displacement_delta
    )
    result[:, attitude] = (
        (1.0 - w) * low[:, attitude] + w * high[:, attitude]
    )
    return ov.rebuild_body_imu(columns, result)


def roundtrip_segments(
    start_sec: float = DEFAULT_TRANSITION_START_SEC,
    return_start_sec: float = DEFAULT_TRANSITION_RETURN_START_SEC,
) -> tuple[tuple[str, float, float], ...]:
    """Expose both crossfades, both recovery intervals, and settled seas."""
    start, end, return_start, return_end = transition_bounds(
        start_sec, return_start_sec
    )
    window_start = ROUNDTRIP_DURATION_SEC - ROUNDTRIP_WINDOW_SEC
    return (
        ("low_start", window_start, start),
        ("rise", start, end),
        ("high_recover", end, end + TRANSITION_DURATION_SEC),
        ("high", end + TRANSITION_DURATION_SEC, return_start),
        ("fall", return_start, return_end),
        ("low_recover", return_end, return_end + TRANSITION_DURATION_SEC),
        ("low_return", return_end + TRANSITION_DURATION_SEC, ROUNDTRIP_DURATION_SEC),
    )


def _write_rows(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def roundtrip_score_rows(
    metrics: dict[str, Any],
    segments: Sequence[tuple[str, float, float]],
) -> list[dict[str, Any]]:
    """Extract the simulator's full-rate segment metrics without re-scoring CSV."""
    rows: list[dict[str, Any]] = []
    for name, start, stop in segments:
        row: dict[str, Any] = {
            "segment": name,
            "start_s": start,
            "stop_s": stop,
        }
        for metric in SEGMENT_METRICS:
            key = f"seg_{name}_{metric}"
            if key not in metrics:
                raise RuntimeError(f"simulator did not emit round-trip score {key}")
            row[metric] = float(metrics[key])
        rows.append(row)
    return rows


def write_roundtrip_score_table(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    """Write the complete seven-segment publication table, including both crossfades."""
    lines = [
        r"\begin{table}[t]",
        r"  \centering",
        r"  \caption{Bidirectional low--high--low transition scores for the deployed OU--III SpectralMSE configuration. Rise and fall are the two \SI{120}{s} crossfade scores; recovery and settled intervals are reported separately.}",
        r"  \label{tab:ou-roundtrip-scores}",
        r"  \footnotesize",
        r"  \setlength{\tabcolsep}{2.7pt}",
        r"  \begin{tabular}{@{}lrrrrr@{}}",
        r"    \toprule",
        r"    Segment & $Z$ [m] & $Z_{\rm ref}$ [m] & $Z/Z_{\rm ref}$ [\%] & $Z/H_s$ [\%] & 3-D [m] \\",
        r"    \midrule",
    ]
    for row in rows:
        label = SEGMENT_LABELS[str(row["segment"])]
        lines.append(
            f"    {label} & {row['disp_z_rms_m']:.3f} & "
            f"{row['disp_z_ref_rms_m']:.3f} & {row['disp_z_pct_refrms']:.2f} & "
            f"{row['disp_z_pct_hs']:.2f} & {row['disp_3d_rms_m']:.3f} \\\\"
        )
    lines.extend((r"    \bottomrule", r"  \end{tabular}", r"\end{table}", ""))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def write_roundtrip_diagnostic(
    svg_path: Path,
    csv_path: Path,
    wave_seed: int,
    imu_seed: int,
    init_seed: int,
    start_sec: float = DEFAULT_TRANSITION_START_SEC,
    return_start_sec: float = DEFAULT_TRANSITION_RETURN_START_SEC,
    decimation: int = 20,
) -> dict[str, Any]:
    """Replay one deployed OU-III adaptive round trip and write all evidence."""
    ov = _validation_module()
    low_path = ov.find_default_input(OU_III_DATA, "1.500", "50.710")
    high_path = ov.find_default_input(OU_III_DATA, "8.500", "202.839")
    columns, low_data = ov.read_wave_csv(low_path, ROUNDTRIP_DURATION_SEC)
    _, high_data = ov.read_wave_csv(high_path, ROUNDTRIP_DURATION_SEC)
    high_scale = DIAGNOSTIC_HIGH_HEIGHT_M / DIAGNOSTIC_SOURCE_HEIGHT_M
    segments = roundtrip_segments(start_sec, return_start_sec)

    _, _, _, return_end_sec = transition_bounds(start_sec, return_start_sec)
    window_start_sec = ROUNDTRIP_DURATION_SEC - ROUNDTRIP_WINDOW_SEC

    with tempfile.TemporaryDirectory(prefix="ou_roundtrip_") as tmp:
        root = Path(tmp)
        low_reference = root / "calibration" / "low" / low_path.name
        ov.write_wave_csv(low_reference, columns, low_data)
        high_reference = root / "calibration" / "high" / ROUNDTRIP_WAVE_NAME
        ov.write_wave_csv(
            high_reference,
            columns,
            ov.scale_wave_motion(columns, high_data, high_scale),
        )
        low_point = ov.calibrate_tuning_point(
            "OU_III", low_reference, ROUNDTRIP_WINDOW_SEC
        )
        high_point = ov.calibrate_tuning_point(
            "OU_III", high_reference, ROUNDTRIP_WINDOW_SEC
        )

        generated = make_roundtrip_wave(
            ov, columns, low_data, high_data, wave_seed, high_scale,
            start_sec, return_start_sec,
        )
        surrogate = root / f"wave_{wave_seed}" / ROUNDTRIP_WAVE_NAME
        ov.write_wave_csv(surrogate, columns, generated)
        metrics, _, _ = ov.run_simulator(
            "OU_III",
            surrogate,
            ROUNDTRIP_WINDOW_SEC,
            imu_seed=imu_seed,
            initialization_seed=init_seed,
            tuning_mode="adaptive",
            aw_cov_sync="periodic",
            segments=segments,
            write_timeseries=True,
        )
        series = ov.read_diagnostic_timeseries("OU_III", surrogate)

    score_rows = roundtrip_score_rows(metrics, segments)
    _write_rows(svg_path.parent / "ou_rs_roundtrip_scores.csv", score_rows)
    write_roundtrip_score_table(
        svg_path.parent / "ou_rs_roundtrip_scores.tex", score_rows
    )

    time = series["time"]
    weight, _, _ = roundtrip_profile(ov, time, start_sec, return_start_sec)
    reference = series["disp_ref_z"]
    estimate = series["disp_est_z"]
    error = estimate - reference
    rolling_hs = ov.rolling_significant_height(reference)
    mixture_hs = np.asarray([
        ov.mixture_significant_height_m(
            DIAGNOSTIC_LOW_HEIGHT_M, DIAGNOSTIC_HIGH_HEIGHT_M, value
        )
        for value in weight
    ])
    linear_hs = (
        (1.0 - weight) * DIAGNOSTIC_LOW_HEIGHT_M
        + weight * DIAGNOSTIC_HIGH_HEIGHT_M
    )

    step = max(1, int(decimation))
    _write_rows(
        csv_path,
        [
            {
                "time_s": float(time[index]),
                "blend_weight": float(weight[index]),
                "rolling_hs_m": float(rolling_hs[index]),
                "mixture_hs_m": float(mixture_hs[index]),
                "disp_ref_z_m": float(reference[index]),
                "disp_est_z_m": float(estimate[index]),
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
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    ov.plot_transition_diagnostic(
        svg_path,
        time[sliced],
        reference[sliced],
        estimate[sliced],
        rolling_hs[sliced],
        mixture_hs[sliced],
        linear_hs[sliced],
        series["tau_applied"][sliced],
        series["sigma_a_applied"][sliced],
        series["R_p0_applied"][sliced],
        (
            (r"low-sea fixed point ($H_s=1.5$ m)", low_point, ":"),
            (r"high-sea fixed point ($H_s=4.0$ m)", high_point, "--"),
        ),
        (
            (start_sec, start_sec + TRANSITION_DURATION_SEC),
            (return_start_sec, return_end_sec),
        ),
        window_start_sec,
    )

    return {
        "family": "OU_III",
        "wave_phase_seed": wave_seed,
        "imu_noise_seed": imu_seed,
        "initialization_seed": init_seed,
        "decimation": step,
        "disp_z_pct_hs": float(metrics["disp_z_pct_hs"]),
        "segment_scores": score_rows,
        "low_tau_s": low_point.tau_s,
        "low_sigma_aw_mps2": low_point.sigma_a_mps2,
        "low_r_s_ms": low_point.RS_ms,
        "high_tau_s": high_point.tau_s,
        "high_sigma_aw_mps2": high_point.sigma_a_mps2,
        "high_r_s_ms": high_point.RS_ms,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--transition-start-sec", type=float, default=DEFAULT_TRANSITION_START_SEC
    )
    parser.add_argument(
        "--transition-return-start-sec",
        type=float,
        default=DEFAULT_TRANSITION_RETURN_START_SEC,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "reports" / "results" / "ou_rs_law",
    )
    parser.add_argument("--diagnostic-seeds", default="11,101,1009")
    parser.add_argument("--decimation", type=int, default=20)
    args = parser.parse_args(argv)

    transition_bounds(args.transition_start_sec, args.transition_return_start_sec)
    seeds = [int(value) for value in args.diagnostic_seeds.split(",") if value.strip()]
    if len(seeds) != 3:
        raise SystemExit("--diagnostic-seeds takes wave,imu,init")

    write_roundtrip_diagnostic(
        args.output_dir / "ou_rs_roundtrip_transition.svg",
        args.output_dir / "ou_rs_roundtrip_transition.csv",
        *seeds,
        start_sec=args.transition_start_sec,
        return_start_sec=args.transition_return_start_sec,
        decimation=args.decimation,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
