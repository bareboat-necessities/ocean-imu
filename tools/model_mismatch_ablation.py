#!/usr/bin/env python3
"""Noise-free model-mismatch ablation for OU-II, OU-III, and TFG.

Runs the eight versioned JONSWAP / PM-Stokes stationary seas through each
filter family with the simulator's ``--no-noise`` switch.  The physical wave
motion, attitude motion, ideal accelerometer/gyro signals, and ideal magnetic
field remain present.  What is removed is the simulator-side sensor corruption:
white noise, initial sensor bias, bias random walk, magnetic residual bias,
magnetic scale/cross-axis error, and magnetic misalignment.

The filter itself is deliberately *not* retuned for perfect sensors.  Its normal
Q/R assumptions, adaptive operating point, pseudo-measurements, startup logic,
and regularization stay exactly as deployed.  The resulting RMS error is
therefore a practical noise-free estimator/model-mismatch floor.  It includes
model mismatch, regularization bias, finite adaptation/startup residue, and
numerical effects; it should not be described as a pure plant-model error.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import ou_validation as ouv

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "reports" / "results" / "model_mismatch_ablation"
DEFAULT_WINDOW_SEC = 900.0


@dataclass(frozen=True)
class Record:
    spectrum: str
    hs_m: float
    filename: str


RECORDS = (
    Record("JONSWAP", 0.27, "wave_data_jonswap_H0.270_L14.047_A30.00_P60.00.csv"),
    Record("JONSWAP", 1.50, "wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv"),
    Record("JONSWAP", 4.00, "wave_data_jonswap_H4.000_L112.766_A30.00_P30.00.csv"),
    Record("JONSWAP", 8.50, "wave_data_jonswap_H8.500_L202.839_A-30.00_P72.00.csv"),
    Record("PM-Stokes", 0.27, "wave_data_pmstokes_H0.270_L14.047_A30.00_P60.00.csv"),
    Record("PM-Stokes", 1.50, "wave_data_pmstokes_H1.500_L50.710_A-30.00_P120.00.csv"),
    Record("PM-Stokes", 4.00, "wave_data_pmstokes_H4.000_L112.766_A30.00_P30.00.csv"),
    Record("PM-Stokes", 8.50, "wave_data_pmstokes_H8.500_L202.839_A-30.00_P72.00.csv"),
)

FAMILIES = {
    "OU-II": ROOT / "tests" / "kalman_ou_ii" / "kalman_ou_ii-sim",
    "OU-III": ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim",
    "TFG": ROOT / "tests" / "kalman_tfg" / "kalman_tfg-sim",
}

RMS_FIELDS = (
    "disp_x_rms_m",
    "disp_y_rms_m",
    "disp_z_rms_m",
    "disp_3d_rms_m",
    "roll_rms_deg",
    "pitch_rms_deg",
    "yaw_rms_deg",
    "accel_bias_3d_rms_mps2",
    "gyro_bias_3d_rms_radps",
)

ROW_FIELDS = (
    "family",
    "spectrum",
    "hs_m",
    "input",
    *RMS_FIELDS,
    "disp_z_ref_rms_m",
    "disp_z_pct_refrms",
    "disp_z_pct_hs",
    "wave_period_s",
    "tau_applied_s",
    "sigma_applied_mps2",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="directory containing the versioned wave CSVs (searched recursively)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"result directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--window-sec",
        type=float,
        default=DEFAULT_WINDOW_SEC,
        help="trailing scoring window in seconds (default: 900)",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=3,
        help="parallel simulator processes (default: 3)",
    )
    return parser.parse_args()


def find_record(data_dir: Path, record: Record) -> Path:
    matches = list(data_dir.rglob(record.filename))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"expected exactly one {record.filename} under {data_dir}, found {len(matches)}"
        )
    return matches[0]


def source_commit() -> str | None:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        return None
    value = completed.stdout.strip()
    return value or None


def run_one(family: str, binary: Path, record: Record, input_path: Path,
            window_sec: float) -> dict[str, Any]:
    if not binary.exists():
        raise FileNotFoundError(f"missing simulator binary: {binary}")

    env = os.environ.copy()
    env.update(
        {
            "W3D_WRITE_TIMESERIES": "0",
            "W3D_VALIDATION_WINDOW_SEC": f"{window_sec:.9g}",
            # Ablations can cross regression sentinels; collect the metrics
            # rather than letting a fitted noisy-input gate stop the study.
            "W3D_COLLECT_ALL_GATES": "1",
        }
    )
    command = [str(binary), "--no-noise", "--input", str(input_path.resolve())]
    completed = subprocess.run(
        command,
        cwd=binary.parent,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if "noise=false" not in completed.stdout:
        tail = "\n".join(completed.stdout.splitlines()[-30:])
        raise RuntimeError(
            f"{family} {record.filename}: simulator did not confirm noise=false:\n{tail}"
        )
    try:
        metrics = ouv.parse_validation_metrics(completed.stdout)
    except ValueError as error:
        tail = "\n".join(completed.stdout.splitlines()[-30:])
        raise RuntimeError(
            f"{family} {record.filename}: no validation metrics (exit "
            f"{completed.returncode}):\n{tail}"
        ) from error
    # A quality-gate failure returns 1 under W3D_COLLECT_ALL_GATES.  Those gates
    # were fitted to the noisy deployed validation and are not acceptance
    # criteria for this ablation.
    if completed.returncode not in (0, 1):
        tail = "\n".join(completed.stdout.splitlines()[-30:])
        raise RuntimeError(
            f"{family} {record.filename}: simulator exit {completed.returncode}:\n{tail}"
        )

    row: dict[str, Any] = {
        "family": family,
        "spectrum": record.spectrum,
        "hs_m": record.hs_m,
        "input": record.filename,
    }
    for field in ROW_FIELDS:
        if field in row:
            continue
        value = metrics.get(field, float("nan"))
        row[field] = float(value) if isinstance(value, (int, float)) else value
    return row


def pooled_rms(rows: Iterable[dict[str, Any]], field: str) -> float:
    values = [float(row[field]) for row in rows]
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return float("nan")
    return math.sqrt(sum(value * value for value in finite) / len(finite))


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for family in FAMILIES:
        family_rows = [row for row in rows if row["family"] == family]
        if len(family_rows) != len(RECORDS):
            raise ValueError(f"{family}: expected {len(RECORDS)} rows, found {len(family_rows)}")
        summary: dict[str, Any] = {"family": family, "records": len(family_rows)}
        for field in RMS_FIELDS:
            summary[field] = pooled_rms(family_rows, field)
        z_ref = pooled_rms(family_rows, "disp_z_ref_rms_m")
        summary["disp_z_ref_rms_m"] = z_ref
        summary["disp_z_pct_refrms"] = (
            100.0 * summary["disp_z_rms_m"] / z_ref if z_ref > 0.0 else float("nan")
        )
        result.append(summary)
    return result


def fmt(value: Any, digits: int = 4) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        return "n/a"
    return f"{number:.{digits}f}"


def markdown_report(rows: list[dict[str, Any]], summaries: list[dict[str, Any]],
                    window_sec: float, commit: str | None) -> str:
    lines = [
        "# Noise-free model-mismatch ablation",
        "",
        "This study runs the standard eight stationary JONSWAP / PM-Stokes wave",
        "records through OU-II, OU-III, and TFG with simulator-side sensor",
        "corruption disabled (`--no-noise`).  The filters themselves keep their",
        "normal deployed covariance assumptions, adaptation laws, pseudo-measurements,",
        "startup, and regularization.",
        "",
        f"Scoring uses the trailing **{window_sec:.0f} s** of each 1200 s record.",
        "Magnetometer updates remain enabled, but the magnetic measurements are ideal.",
        "The reported floor therefore contains model/estimator mismatch, intentional",
        "regularization bias, residual adaptation/startup effects, and numerical error;",
        "it is not a claim of pure plant-model mismatch in isolation.",
        "",
    ]
    if commit:
        lines.extend((f"Source commit used for the replay: `{commit}`.", ""))

    lines.extend(
        [
            "## Pooled RMS across the eight equal-duration records",
            "",
            "Because every record contributes the same 900 s window at the same sample",
            "rate, `sqrt(mean(record_RMS^2))` is the exact pooled RMS over their concatenation.",
            "",
            "| Family | X disp [m] | Y disp [m] | Z disp [m] | 3D disp [m] | Z / ref RMS [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for summary in summaries:
        lines.append(
            "| {family} | {x} | {y} | {z} | {d3} | {zp} | {r} | {p} | {yw} | {ab} | {gb} |".format(
                family=summary["family"],
                x=fmt(summary["disp_x_rms_m"]),
                y=fmt(summary["disp_y_rms_m"]),
                z=fmt(summary["disp_z_rms_m"]),
                d3=fmt(summary["disp_3d_rms_m"]),
                zp=fmt(summary["disp_z_pct_refrms"], 3),
                r=fmt(summary["roll_rms_deg"]),
                p=fmt(summary["pitch_rms_deg"]),
                yw=fmt(summary["yaw_rms_deg"]),
                ab=fmt(summary["accel_bias_3d_rms_mps2"], 6),
                gb=fmt(summary["gyro_bias_3d_rms_radps"], 7),
            )
        )

    lines.extend(
        [
            "",
            "## Per-record RMS",
            "",
            "| Family | Sea | Hs [m] | X [m] | Y [m] | Z [m] | 3D [m] | Z / Hs [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    order = {family: index for index, family in enumerate(FAMILIES)}
    sorted_rows = sorted(rows, key=lambda row: (order[row["family"]], row["spectrum"], row["hs_m"]))
    for row in sorted_rows:
        lines.append(
            "| {family} | {sea} | {hs:.2f} | {x} | {y} | {z} | {d3} | {zh} | {r} | {p} | {yw} | {ab} | {gb} |".format(
                family=row["family"],
                sea=row["spectrum"],
                hs=float(row["hs_m"]),
                x=fmt(row["disp_x_rms_m"]),
                y=fmt(row["disp_y_rms_m"]),
                z=fmt(row["disp_z_rms_m"]),
                d3=fmt(row["disp_3d_rms_m"]),
                zh=fmt(row["disp_z_pct_hs"], 3),
                r=fmt(row["roll_rms_deg"]),
                p=fmt(row["pitch_rms_deg"]),
                yw=fmt(row["yaw_rms_deg"]),
                ab=fmt(row["accel_bias_3d_rms_mps2"], 6),
                gb=fmt(row["gyro_bias_3d_rms_radps"], 7),
            )
        )

    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "`--no-noise` removes the simulator's stochastic and calibration-error",
            "sensor terms before they reach the filters.  It does **not** remove wave",
            "nonlinearity, attitude/translation coupling, OU/TFG prior mismatch, the",
            "integral pseudo-measurements, finite adaptation bandwidth, startup residue,",
            "or discretization.  Those effects are intentionally what this ablation",
            "leaves visible.",
            "",
        ]
    )
    return "\n".join(lines)


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in writer.fieldnames})


def main() -> int:
    args = parse_args()
    if not (math.isfinite(args.window_sec) and args.window_sec > 0.0):
        raise SystemExit("--window-sec must be positive")
    if args.jobs < 1:
        raise SystemExit("--jobs must be >= 1")

    data_dir = args.data_dir.resolve()
    inputs = {record: find_record(data_dir, record) for record in RECORDS}
    for family, binary in FAMILIES.items():
        if not binary.exists():
            raise SystemExit(f"missing {family} simulator: {binary}; build it first")

    tasks = []
    rows: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        for family, binary in FAMILIES.items():
            for record in RECORDS:
                tasks.append(
                    executor.submit(
                        run_one, family, binary, record, inputs[record], args.window_sec
                    )
                )
        for future in as_completed(tasks):
            row = future.result()
            rows.append(row)
            print(
                f"{row['family']:6s} {row['spectrum']:9s} Hs={row['hs_m']:4.2f} "
                f"Z={row['disp_z_rms_m']:.5f} m 3D={row['disp_3d_rms_m']:.5f} m "
                f"RPY=({row['roll_rms_deg']:.4f},{row['pitch_rms_deg']:.4f},"
                f"{row['yaw_rms_deg']:.4f}) deg"
            )

    family_order = {family: index for index, family in enumerate(FAMILIES)}
    spectrum_order = {"JONSWAP": 0, "PM-Stokes": 1}
    rows.sort(key=lambda row: (family_order[row["family"]], spectrum_order[row["spectrum"]], row["hs_m"]))
    summaries = summarize(rows)
    commit = source_commit()

    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "model_mismatch_runs.csv", rows, ROW_FIELDS)
    summary_fields = ("family", "records", *RMS_FIELDS, "disp_z_ref_rms_m", "disp_z_pct_refrms")
    write_csv(out / "model_mismatch_summary.csv", summaries, summary_fields)
    (out / "model_mismatch_report.md").write_text(
        markdown_report(rows, summaries, args.window_sec, commit), encoding="utf-8"
    )
    manifest = {
        "study": "noise-free model-mismatch ablation",
        "source_commit": commit,
        "simulation_data": "oceanography-waves-lib v1.1.3",
        "families": list(FAMILIES),
        "records": [record.__dict__ for record in RECORDS],
        "record_count_per_family": len(RECORDS),
        "total_replays": len(RECORDS) * len(FAMILIES),
        "window_sec": args.window_sec,
        "magnetometer_enabled": True,
        "sensor_noise_injected": False,
        "simulator_switch": "--no-noise",
        "filter_covariance_and_tuning_unchanged": True,
        "outputs": [
            "model_mismatch_runs.csv",
            "model_mismatch_summary.csv",
            "model_mismatch_report.md",
        ],
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Wrote noise-free mismatch study to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
