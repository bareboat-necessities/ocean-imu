#!/usr/bin/env python3
"""Engine-noise degradation study for OU-II, OU-III, and TFG.

Replays the eight versioned stationary JONSWAP / PM-Stokes seas through each
filter family while the simulator injects the vibration a hull-mounted IMU
records on a mid-size recreational cruising sailboat motoring under its
inboard auxiliary diesel.  The archetype modeled in ``W3dSimCommon`` is a
naturally aspirated three-cylinder four-stroke on flexible mounts, driving a
three-blade fixed propeller through a 2.6:1 reduction gear.

Everything else is the deployed configuration.  The ordinary sensor noise
models stay on, the filters keep their normal covariances, adaptation laws,
pseudo-measurements, startup logic, and regularization, and the vessel's
rigid-body wave response is untouched: only the sensor path changes.

Four arms run against a common engine-off baseline.

speed
    Engine speed from idle to near maximum at the nominal vibration level.

level
    Vibration level at cruise rpm, from a quiet, well-isolated installation to
    a sensor bolted near the engine bed.

bandwidth
    The sensor's anti-alias bandwidth ahead of the 200 Hz sample rate, at
    cruise rpm and the nominal level.  A narrower filter both removes power
    and stops the high crank orders from folding.

matched
    The control that separates those two effects: the same bandwidth sweep
    with the hull level rescaled so every cell records the *same* vibration
    RMS as the nominal cruise cell.  If degradation tracked where the folded
    lines land, these cells would still differ; if it tracks recorded power,
    they collapse.

path
    Attribution between the two sensors.  The cruise cell rerun with the
    model's gyroscope terms -- hull angular vibration and the gyroscope's own
    linear-acceleration sensitivity -- switched off, so only the accelerometer
    is perturbed.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import ou_validation as ouv

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "reports" / "results" / "engine_noise_degradation"
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

# The nominal operating point.  0.60 m/s^2 of hull broadband RMS is about
# 0.061 g; at the 30-60 Hz where this engine puts most of its energy that is
# roughly 2 mm/s RMS velocity, i.e. inside the ISO 6954 range normally called
# comfortable for small-craft accommodation.  It is a routine cruise level for
# a sensor on a cabin bulkhead, not a worst case.
CRUISE_RPM = 2400.0
NOMINAL_LEVEL_MPS2 = 0.60
NOMINAL_BANDWIDTH_HZ = 80.0

SPEED_RPM = (800.0, 1200.0, 1600.0, 2000.0, 2400.0, 2800.0, 3200.0)
LEVEL_MPS2 = (0.15, 0.30, 0.60, 1.20, 2.40)
BANDWIDTH_HZ = (20.0, 40.0, 80.0, 160.0)


@dataclass(frozen=True)
class Setting:
    """One engine configuration, shared by every family and record."""

    arm: str
    label: str
    rpm: float = 0.0
    level_mps2: float = NOMINAL_LEVEL_MPS2
    bandwidth_hz: float = NOMINAL_BANDWIDTH_HZ
    # Set for the matched arm: the recorded RMS this cell was rescaled onto.
    matched_to_mps2: float = float("nan")
    # Set for the path arm: drop the model's gyroscope coupling terms.
    gyro_path: bool = True

    @property
    def engine_on(self) -> bool:
        return self.rpm > 0.0

    def env(self) -> dict[str, str]:
        if not self.engine_on:
            return {}
        env = {
            "W3D_ENGINE_RPM": f"{self.rpm:.9g}",
            "W3D_ENGINE_LEVEL_MPS2": f"{self.level_mps2:.9g}",
            "W3D_ENGINE_BANDWIDTH_HZ": f"{self.bandwidth_hz:.9g}",
        }
        if not self.gyro_path:
            env["W3D_ENGINE_GYRO_LEVER_M"] = "0"
            env["W3D_ENGINE_GYRO_G_SENS"] = "0"
        return env


# Mean error over the scored window.  The RMS alone cannot say whether the
# vibration cost is a fluctuation or a constant offset, and the answer turns
# out to be the whole mechanism, so both are carried.
MEAN_FIELDS = (
    "disp_x_mean_m",
    "disp_y_mean_m",
    "disp_z_mean_m",
    "roll_mean_deg",
    "pitch_mean_deg",
    "yaw_mean_deg",
)

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

SETTING_FIELDS = (
    "arm",
    "label",
    "rpm",
    "level_mps2",
    "bandwidth_hz",
    "hull_rms_mps2",
    "recorded_rms_mps2",
    "firing_hz",
    "blade_hz",
    "vre_offset_mps2",
)

ROW_FIELDS = (
    "family",
    *SETTING_FIELDS,
    "spectrum",
    "hs_m",
    "input",
    *RMS_FIELDS,
    *MEAN_FIELDS,
    "disp_z_ref_rms_m",
    "disp_z_pct_refrms",
    "disp_z_pct_hs",
    "tau_applied_s",
    "sigma_applied_mps2",
    "accel_variance_m2ps4",
)

SUMMARY_FIELDS = (
    "family",
    *SETTING_FIELDS,
    "records",
    *RMS_FIELDS,
    *MEAN_FIELDS,
    "disp_z_ref_rms_m",
    "disp_z_pct_refrms",
    "disp_3d_ratio_to_baseline",
    "disp_z_ratio_to_baseline",
    "pitch_ratio_to_baseline",
    "yaw_ratio_to_baseline",
)

BANNER_RE = re.compile(r"^ENGINE_VIBRATION (.*)$", re.MULTILINE)
MEANS_RE = re.compile(r"^VALIDATION_METRICS_MEANS (.*)$", re.MULTILINE)


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
        "--no-plots",
        action="store_true",
        help="skip the SVG charts (for environments without Matplotlib)",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=4,
        help="parallel simulator processes (default: 4)",
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
    return completed.stdout.strip() or None


def parse_banner(stdout: str) -> dict[str, float]:
    """Pull the simulator's ENGINE_VIBRATION line into a dict of floats."""

    match = BANNER_RE.search(stdout)
    if not match:
        raise ValueError("simulator printed no ENGINE_VIBRATION banner")
    values: dict[str, float] = {}
    for token in match.group(1).split():
        key, _, raw = token.partition("=")
        try:
            values[key] = float(raw)
        except ValueError:
            continue
    return values


def parse_means(stdout: str) -> dict[str, float]:
    """Pull the simulator's VALIDATION_METRICS_MEANS line into a dict."""

    match = MEANS_RE.search(stdout)
    if not match:
        raise ValueError("simulator printed no VALIDATION_METRICS_MEANS line")
    values: dict[str, float] = {}
    for token in match.group(1).split():
        key, _, raw = token.partition("=")
        try:
            values[key] = float(raw)
        except ValueError:
            continue
    return values


def invoke(binary: Path, input_path: Path, setting: Setting,
           window_sec: float) -> str:
    env = os.environ.copy()
    env.update(
        {
            "W3D_WRITE_TIMESERIES": "0",
            "W3D_VALIDATION_WINDOW_SEC": f"{window_sec:.9g}",
            # Every gate in this study was fitted to a vibration-free input, so
            # collect the metrics instead of stopping on the first failure.
            "W3D_COLLECT_ALL_GATES": "1",
            # All three families now arm the front-end vibration guard by
            # default.  This study is the controlled comparison that motivates
            # it, so it pins the guard off in every family and keeps measuring
            # the unconditioned measurement path.  What the default actually
            # does is tools/ou3_engine_noise_mitigation.py.  Each simulator
            # ignores the prefixes that are not its own, so setting all three
            # here keeps one env block for one protocol.
            "OU_II_ACC_GUARD_HZ": "0",
            "OU_III_ACC_GUARD_HZ": "0",
            "TFG_ACC_GUARD_HZ": "0",
        }
    )
    env.update(setting.env())
    completed = subprocess.run(
        [str(binary), "--input", str(input_path.resolve())],
        cwd=binary.parent,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if completed.returncode not in (0, 1):
        tail = "\n".join(completed.stdout.splitlines()[-30:])
        raise RuntimeError(
            f"{binary.name} exit {completed.returncode} for {setting.label}:\n{tail}"
        )
    return completed.stdout


def probe_setting(binary: Path, input_path: Path, setting: Setting) -> dict[str, float]:
    """Read back the vibration the simulator will actually inject.

    The banner is emitted before the replay loop, so a one-second scoring
    window is enough to recover it; the level and bandwidth it reports are
    record-independent.
    """

    stdout = invoke(binary, input_path, setting, window_sec=1.0)
    return parse_banner(stdout)


def run_one(family: str, binary: Path, record: Record, input_path: Path,
            setting: Setting, banner: dict[str, float],
            window_sec: float) -> dict[str, Any]:
    stdout = invoke(binary, input_path, setting, window_sec)
    if setting.engine_on and "ENGINE_VIBRATION" not in stdout:
        raise RuntimeError(
            f"{family} {setting.label} {record.filename}: engine model did not engage"
        )
    if not setting.engine_on and "ENGINE_VIBRATION" in stdout:
        raise RuntimeError(
            f"{family} baseline {record.filename}: engine model engaged unexpectedly"
        )
    try:
        metrics = ouv.parse_validation_metrics(stdout)
    except ValueError as error:
        tail = "\n".join(stdout.splitlines()[-30:])
        raise RuntimeError(
            f"{family} {setting.label} {record.filename}: no validation metrics:\n{tail}"
        ) from error

    row: dict[str, Any] = {
        "family": family,
        "arm": setting.arm,
        "label": setting.label,
        "rpm": setting.rpm,
        "level_mps2": setting.level_mps2 if setting.engine_on else 0.0,
        "bandwidth_hz": setting.bandwidth_hz if setting.engine_on else float("nan"),
        "hull_rms_mps2": banner.get("hull_rms_mps2", 0.0),
        "recorded_rms_mps2": banner.get("recorded_rms_mps2", 0.0),
        "firing_hz": banner.get("firing_hz", 0.0),
        "blade_hz": banner.get("blade_hz", 0.0),
        "vre_offset_mps2": banner.get("vre_offset_mps2", 0.0),
        "spectrum": record.spectrum,
        "hs_m": record.hs_m,
        "input": record.filename,
    }
    means = parse_means(stdout)
    for name in ROW_FIELDS:
        if name in row:
            continue
        source = means if name in MEAN_FIELDS else metrics
        value = source.get(name, float("nan"))
        row[name] = float(value) if isinstance(value, (int, float)) else value
    return row


def pooled_rms(rows: Iterable[dict[str, Any]], name: str) -> float:
    values = [float(row[name]) for row in rows]
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return float("nan")
    return math.sqrt(sum(value * value for value in finite) / len(finite))


def summarize(rows: list[dict[str, Any]], settings: list[Setting]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    baseline: dict[str, dict[str, Any]] = {}

    for family in FAMILIES:
        for setting in settings:
            cells = [
                row for row in rows
                if row["family"] == family and row["label"] == setting.label
            ]
            if len(cells) != len(RECORDS):
                raise ValueError(
                    f"{family} {setting.label}: expected {len(RECORDS)} rows, "
                    f"found {len(cells)}"
                )
            summary: dict[str, Any] = {"family": family, "records": len(cells)}
            for name in SETTING_FIELDS:
                summary[name] = cells[0][name]
            for name in RMS_FIELDS:
                summary[name] = pooled_rms(cells, name)
            # Pooled as the RMS across records of each record's mean error:
            # the typical size of the systematic offset, without records of
            # opposite nominal heading cancelling in a signed average.
            for name in MEAN_FIELDS:
                summary[name] = pooled_rms(cells, name)
            z_ref = pooled_rms(cells, "disp_z_ref_rms_m")
            summary["disp_z_ref_rms_m"] = z_ref
            summary["disp_z_pct_refrms"] = (
                100.0 * summary["disp_z_rms_m"] / z_ref if z_ref > 0.0 else float("nan")
            )
            if setting.arm == "baseline":
                baseline[family] = summary
            result.append(summary)

    for summary in result:
        reference = baseline[summary["family"]]
        for name, source in (
            ("disp_3d_ratio_to_baseline", "disp_3d_rms_m"),
            ("disp_z_ratio_to_baseline", "disp_z_rms_m"),
            ("pitch_ratio_to_baseline", "pitch_rms_deg"),
            ("yaw_ratio_to_baseline", "yaw_rms_deg"),
        ):
            denominator = float(reference[source])
            summary[name] = (
                float(summary[source]) / denominator
                if denominator > 0.0 else float("nan")
            )
    return result


def build_settings(probe_binary: Path, probe_input: Path) -> list[Setting]:
    """Enumerate the study arms, resolving the matched arm against a probe."""

    settings: list[Setting] = [Setting(arm="baseline", label="engine off")]

    for rpm in SPEED_RPM:
        settings.append(
            Setting(
                arm="speed",
                label=f"speed {rpm:.0f} rpm",
                rpm=rpm,
            )
        )
    for level in LEVEL_MPS2:
        if level == NOMINAL_LEVEL_MPS2:
            # Already covered by the cruise cell of the speed arm.
            continue
        settings.append(
            Setting(
                arm="level",
                label=f"level {level:.2f} m/s^2",
                rpm=CRUISE_RPM,
                level_mps2=level,
            )
        )
    for bandwidth in BANDWIDTH_HZ:
        if bandwidth == NOMINAL_BANDWIDTH_HZ:
            continue
        settings.append(
            Setting(
                arm="bandwidth",
                label=f"bandwidth {bandwidth:.0f} Hz",
                rpm=CRUISE_RPM,
                bandwidth_hz=bandwidth,
            )
        )

    settings.append(
        Setting(
            arm="path",
            label="accelerometer only",
            rpm=CRUISE_RPM,
            gyro_path=False,
        )
    )

    # The matched arm needs the recorded RMS each bandwidth actually produces,
    # so probe the simulator rather than reimplementing the model here.
    nominal = Setting(arm="speed", label="probe", rpm=CRUISE_RPM)
    target = probe_setting(probe_binary, probe_input, nominal)["recorded_rms_mps2"]
    for bandwidth in BANDWIDTH_HZ:
        if bandwidth == NOMINAL_BANDWIDTH_HZ:
            continue
        probe = Setting(
            arm="matched", label="probe", rpm=CRUISE_RPM, bandwidth_hz=bandwidth
        )
        recorded = probe_setting(probe_binary, probe_input, probe)["recorded_rms_mps2"]
        if not (recorded > 0.0):
            raise RuntimeError(f"probe at {bandwidth} Hz recorded no vibration")
        settings.append(
            Setting(
                arm="matched",
                label=f"matched {bandwidth:.0f} Hz",
                rpm=CRUISE_RPM,
                level_mps2=NOMINAL_LEVEL_MPS2 * target / recorded,
                bandwidth_hz=bandwidth,
                matched_to_mps2=target,
            )
        )
    return settings


def fmt(value: Any, digits: int = 4) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        return "n/a"
    return f"{number:.{digits}f}"


def arm_table(summaries: list[dict[str, Any]], arm: str, axis_header: str,
              axis_key: str, axis_digits: int) -> list[str]:
    lines = [
        "",
        f"| Family | {axis_header} | Recorded vib. [m/s²] | Z [m] | 3-D [m] "
        "| Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] "
        "| 3-D / baseline |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    selected = [s for s in summaries if s["arm"] in (arm, "baseline")]
    order = {family: index for index, family in enumerate(FAMILIES)}
    selected.sort(key=lambda s: (order[s["family"]], float(s[axis_key] or 0.0)))
    for summary in selected:
        axis = "off" if summary["arm"] == "baseline" else fmt(summary[axis_key], axis_digits)
        lines.append(
            "| {family} | {axis} | {vib} | {z} | {d3} | {p} | {pm} | {dm} | {yw} "
            "| {ratio} |".format(
                family=summary["family"],
                axis=axis,
                vib=fmt(summary["recorded_rms_mps2"], 4),
                z=fmt(summary["disp_z_rms_m"]),
                d3=fmt(summary["disp_3d_rms_m"]),
                p=fmt(summary["pitch_rms_deg"]),
                pm=fmt(summary["pitch_mean_deg"], 3),
                dm=fmt(
                    math.sqrt(
                        float(summary["disp_x_mean_m"]) ** 2
                        + float(summary["disp_y_mean_m"]) ** 2
                        + float(summary["disp_z_mean_m"]) ** 2
                    ),
                    3,
                ),
                yw=fmt(summary["yaw_rms_deg"]),
                ratio=fmt(summary["disp_3d_ratio_to_baseline"], 2),
            )
        )
    return lines


def markdown_report(summaries: list[dict[str, Any]], window_sec: float,
                    commit: str | None, plots: bool) -> str:
    lines = [
        "# Engine-noise degradation study",
        "",
        "How far the deployed estimators degrade when the IMU also records the",
        "vibration of an inboard auxiliary diesel.  The eight stationary JONSWAP /",
        "PM-Stokes records are replayed through OU-II, OU-III, and TFG with the",
        "ordinary sensor noise models on and the engine vibration model added on",
        "top; the filters keep their deployed covariances, adaptation, startup, and",
        "regularization throughout.",
        "",
        "The vessel modeled is a mid-size recreational cruising sailboat: a",
        "naturally aspirated three-cylinder four-stroke diesel on flexible mounts,",
        "a 2.6:1 reduction gear, and a three-blade fixed propeller.  The model is",
        "a sensor-path model.  It adds crank orders, driveline shaft- and",
        "blade-rate lines, a broadband structural floor, governor hunting, mount",
        "transmissibility, the sensor's finite anti-alias bandwidth, accelerometer",
        "vibration rectification, and gyroscope g-sensitivity.  It does not change",
        "the vessel's rigid-body response to the sea.",
        "",
        f"Scoring uses the trailing **{window_sec:.0f} s** of each 1200 s record, and",
        "every reported number pools the eight equal-duration records as",
        "`sqrt(mean(record_RMS^2))`, which is the exact RMS over their concatenation.",
        "",
    ]
    if commit:
        lines.extend((f"Source commit used for the replay: `{commit}`.", ""))

    lines.extend(
        [
            "## Engine speed",
            "",
            f"Vibration level fixed at {NOMINAL_LEVEL_MPS2:.2f} m/s² of hull broadband",
            f"RMS at {CRUISE_RPM:.0f} rpm, sensor bandwidth"
            f" {NOMINAL_BANDWIDTH_HZ:.0f} Hz.",
        ]
    )
    lines.extend(arm_table(summaries, "speed", "Engine speed [rpm]", "rpm", 0))

    lines.extend(
        [
            "",
            "## Vibration level",
            "",
            f"Engine speed fixed at {CRUISE_RPM:.0f} rpm, sensor bandwidth"
            f" {NOMINAL_BANDWIDTH_HZ:.0f} Hz.  The level is the hull broadband RMS",
            "before the sensor's anti-alias filter, so a quiet, well-isolated",
            "installation sits at the low end and a sensor near the engine bed at the",
            "high end.",
        ]
    )
    lines.extend(arm_table(summaries, "level", "Hull level [m/s²]", "level_mps2", 2))

    lines.extend(
        [
            "",
            "## Sensor anti-alias bandwidth",
            "",
            f"Engine speed fixed at {CRUISE_RPM:.0f} rpm and the hull level at",
            f"{NOMINAL_LEVEL_MPS2:.2f} m/s².  A narrower filter removes power *and*",
            "stops the high crank orders from folding below Nyquist.",
        ]
    )
    lines.extend(arm_table(summaries, "bandwidth", "Bandwidth [Hz]", "bandwidth_hz", 0))

    lines.extend(
        [
            "",
            "## Matched-power control",
            "",
            "The same bandwidth sweep with the hull level rescaled so every cell",
            "records the same vibration RMS as the nominal cruise cell.  The folded",
            "line frequencies still differ from cell to cell; the recorded power no",
            "longer does.",
        ]
    )
    lines.extend(arm_table(summaries, "matched", "Bandwidth [Hz]", "bandwidth_hz", 0))

    lines.extend(
        [
            "",
            "## Sensor attribution",
            "",
            "The nominal cruise cell rerun with the model's gyroscope terms switched",
            "off, so the accelerometer is the only perturbed sensor.  Compare against",
            f"the {CRUISE_RPM:.0f} rpm row of the engine-speed table.",
        ]
    )
    lines.extend(arm_table(summaries, "path", "Engine speed [rpm]", "rpm", 0))

    if plots:
        lines.extend(
            [
                "",
                "## Figures",
                "",
                f"- `{SPEED_PLOT_NAME}`: recorded vibration and the displacement,",
                "  pitch, and yaw response against engine speed.",
                f"- `{MECHANISM_PLOT_NAME}`: the level sweep, the bandwidth sweep with",
                "  its matched-power control, every cell of the study collapsed onto",
                "  recorded vibration RMS, and the rectified tilt offset that drives",
                "  all of it.",
                "",
                "Both are mirrored byte-for-byte into `doc/kalman_ou_iii/` for the article.",
            ]
        )

    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "The engine model perturbs the accelerometer and gyroscope only.  The",
            "magnetometer, the wave records, and the vessel's rigid-body motion are",
            "unchanged, so this study bounds the sensor-path cost of motoring and not",
            "the full difference between sailing and motoring.  A real passage under",
            "power also changes the encounter spectrum, adds propeller-induced surge,",
            "and runs the engine at a speed that itself varies with the sea; none of",
            "that is modeled here.",
            "",
        ]
    )
    return "\n".join(lines)


# Charts published with the article.
FAMILY_COLORS = {
    "OU-II": "#1b6ca8",
    "OU-III": "#c0392b",
    "TFG": "#2e8b57",
}
FAMILY_MARKERS = {
    "OU-II": "o",
    "OU-III": "s",
    "TFG": "^",
}
ARM_MARKERS = {
    "baseline": "*",
    "speed": "s",
    "level": "o",
    "bandwidth": "^",
    "matched": "D",
}

SPEED_PLOT_NAME = "ou_engine_noise_speed.svg"
MECHANISM_PLOT_NAME = "ou_engine_noise_mechanism.svg"


def reproducible_pyplot():
    """Matplotlib configured for byte-reproducible published SVGs."""

    cache_dir = Path(tempfile.gettempdir()) / "ocean-imu-matplotlib-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    import matplotlib

    matplotlib.use("Agg")
    # Without these the figure cannot be reproduced: matplotlib derives svg
    # element ids from a per-process salt and stamps a creation date, so two
    # identical runs disagree on bytes the evidence fingerprint hashes.
    matplotlib.rcParams["svg.hashsalt"] = "ocean-imu-engine-noise"
    import matplotlib.pyplot as plt

    return plt


def _save_svg(figure: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, format="svg", metadata={"Date": None})


def _arm(summaries: list[dict[str, Any]], family: str, arm: str,
         axis_key: str) -> tuple[list[float], list[dict[str, Any]]]:
    selected = sorted(
        (s for s in summaries if s["family"] == family and s["arm"] == arm),
        key=lambda s: float(s[axis_key]),
    )
    return [float(s[axis_key]) for s in selected], selected


def write_speed_plot(path: Path, summaries: list[dict[str, Any]]) -> None:
    """Degradation against engine speed, with the vibration it comes from."""

    plt = reproducible_pyplot()
    families = list(FAMILIES)
    figure, axes = plt.subplots(2, 2, figsize=(7.8, 5.6))

    # The vibration itself is a property of the engine model, not the filter,
    # so read it off any one family.
    rpm, cells = _arm(summaries, families[0], "speed", "rpm")
    axis = axes[0][0]
    axis.plot(
        rpm,
        [float(cell["hull_rms_mps2"]) for cell in cells],
        marker="o",
        markersize=4.5,
        linewidth=1.4,
        color="#444444",
        label="hull",
    )
    axis.plot(
        rpm,
        [float(cell["recorded_rms_mps2"]) for cell in cells],
        marker="s",
        markersize=4.5,
        linewidth=1.4,
        color="#c0392b",
        label="recorded",
    )
    axis.set_xlabel("engine speed (rpm)")
    axis.set_ylabel(r"vibration RMS (m s$^{-2}$)")
    axis.set_title("Injected vibration", fontsize=10)
    axis.set_ylim(bottom=0.0)
    axis.grid(True, alpha=0.25)
    axis.legend(frameon=False, fontsize=8)

    def response(target, name, ylabel, title, log=False) -> None:
        for family in families:
            speeds, selected = _arm(summaries, family, "speed", "rpm")
            base = [s for s in summaries if s["family"] == family and s["arm"] == "baseline"]
            target.plot(
                speeds,
                [float(cell[name]) for cell in selected],
                marker=FAMILY_MARKERS[family],
                markersize=4.5,
                linewidth=1.4,
                color=FAMILY_COLORS[family],
                label=family,
            )
            if base:
                target.axhline(
                    float(base[0][name]),
                    color=FAMILY_COLORS[family],
                    linewidth=0.9,
                    linestyle=":",
                )
        target.set_xlabel("engine speed (rpm)")
        target.set_ylabel(ylabel)
        target.set_title(title, fontsize=10)
        target.grid(True, alpha=0.25)
        if log:
            target.set_yscale("log")
        else:
            target.set_ylim(bottom=0.0)

    # Dotted horizontals are each family's engine-off baseline.
    response(axes[0][1], "disp_3d_rms_m", "RMS error (m)", "3-D displacement", log=True)
    response(axes[1][0], "pitch_rms_deg", "RMS error (deg)", "Pitch")
    response(axes[1][1], "yaw_rms_deg", "RMS error (deg)", "Yaw")

    handles, labels = axes[0][1].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        frameon=False,
        ncol=len(families),
        fontsize=9,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.0),
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    _save_svg(figure, path)
    plt.close(figure)


def write_mechanism_plot(path: Path, summaries: list[dict[str, Any]]) -> None:
    """Level, bandwidth with its matched control, and the collapse."""

    plt = reproducible_pyplot()
    families = list(FAMILIES)
    figure, axes_grid = plt.subplots(2, 2, figsize=(8.4, 6.2))
    axes = [axes_grid[0][0], axes_grid[0][1], axes_grid[1][0], axes_grid[1][1]]

    axis = axes[0]
    for family in families:
        levels, cells = _arm(summaries, family, "level", "level_mps2")
        cruise = [
            s for s in summaries
            if s["family"] == family and s["arm"] == "speed"
            and float(s["rpm"]) == CRUISE_RPM
        ]
        merged = sorted(cells + cruise, key=lambda s: float(s["level_mps2"]))
        axis.plot(
            [float(cell["level_mps2"]) for cell in merged],
            [float(cell["disp_3d_rms_m"]) for cell in merged],
            marker=FAMILY_MARKERS[family],
            markersize=5.0,
            linewidth=1.4,
            color=FAMILY_COLORS[family],
            label=family,
        )
        base = [s for s in summaries if s["family"] == family and s["arm"] == "baseline"]
        if base:
            axis.axhline(
                float(base[0]["disp_3d_rms_m"]),
                color=FAMILY_COLORS[family],
                linewidth=0.9,
                linestyle=":",
            )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xticks(list(LEVEL_MPS2))
    axis.set_xticklabels([f"{level:g}" for level in LEVEL_MPS2])
    # Explicit major ticks do not suppress the log minor labels, which would
    # otherwise be drawn straight through them.
    axis.minorticks_off()
    axis.set_xlabel(r"hull vibration level (m s$^{-2}$)")
    axis.set_ylabel("3-D RMS error (m)")
    axis.set_title("Level sweep at cruise rpm", fontsize=10)
    axis.grid(True, which="major", alpha=0.25)

    axis = axes[1]
    for family in families:
        for arm, style in (("bandwidth", "-"), ("matched", "--")):
            bandwidths, cells = _arm(summaries, family, arm, "bandwidth_hz")
            cruise = [
                s for s in summaries
                if s["family"] == family and s["arm"] == "speed"
                and float(s["rpm"]) == CRUISE_RPM
            ]
            merged = sorted(cells + cruise, key=lambda s: float(s["bandwidth_hz"]))
            axis.plot(
                [float(cell["bandwidth_hz"]) for cell in merged],
                [float(cell["disp_3d_rms_m"]) for cell in merged],
                marker=ARM_MARKERS[arm],
                markersize=5.0,
                linewidth=1.4,
                linestyle=style,
                color=FAMILY_COLORS[family],
                markerfacecolor=FAMILY_COLORS[family] if arm == "bandwidth" else "none",
                label=f"{family} {'as-is' if arm == 'bandwidth' else 'matched'}",
            )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xticks(list(BANDWIDTH_HZ))
    axis.set_xticklabels([f"{bandwidth:g}" for bandwidth in BANDWIDTH_HZ])
    axis.minorticks_off()
    axis.set_xlabel("sensor anti-alias bandwidth (Hz)")
    axis.set_ylabel("3-D RMS error (m)")
    axis.set_title("Bandwidth, and the same at matched power", fontsize=10)
    axis.grid(True, which="major", alpha=0.25)
    axis.legend(frameon=False, fontsize=6.5, ncol=1)

    axis = axes[2]
    for family in families:
        for arm in ("speed", "level", "bandwidth", "matched"):
            cells = [s for s in summaries if s["family"] == family and s["arm"] == arm]
            if not cells:
                continue
            axis.scatter(
                [float(cell["recorded_rms_mps2"]) for cell in cells],
                [float(cell["disp_3d_rms_m"]) for cell in cells],
                marker=ARM_MARKERS[arm],
                s=26.0,
                color=FAMILY_COLORS[family],
                edgecolors="#222222",
                linewidths=0.4,
                label=arm if family == families[0] else None,
            )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel(r"recorded vibration RMS (m s$^{-2}$)")
    axis.set_ylabel("3-D RMS error (m)")
    axis.set_title("Every cell, against recorded power", fontsize=10)
    axis.grid(True, which="major", alpha=0.25)
    axis.legend(frameon=False, fontsize=7, title="arm", title_fontsize=7)

    # The mechanism itself.  What vibration buys is a *static* tilt offset, and
    # the striking feature is that it appears almost fully formed at the
    # quietest level tested and then grows very slowly: there is no gentle
    # onset to trade against.
    axis = axes[3]
    for family in families:
        for arm in ("speed", "level", "bandwidth", "matched"):
            cells = [s for s in summaries if s["family"] == family and s["arm"] == arm]
            if not cells:
                continue
            axis.scatter(
                [float(cell["recorded_rms_mps2"]) for cell in cells],
                [abs(float(cell["pitch_mean_deg"])) for cell in cells],
                marker=ARM_MARKERS[arm],
                s=26.0,
                color=FAMILY_COLORS[family],
                edgecolors="#222222",
                linewidths=0.4,
            )
    # Dotted horizontals are each family's engine-off offset.
    for family in families:
        base = [s for s in summaries if s["family"] == family and s["arm"] == "baseline"]
        if base:
            axis.axhline(
                abs(float(base[0]["pitch_mean_deg"])),
                color=FAMILY_COLORS[family],
                linewidth=0.9,
                linestyle=":",
            )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel(r"recorded vibration RMS (m s$^{-2}$)")
    axis.set_ylabel("pitch offset (deg)")
    axis.set_title("Static tilt offset", fontsize=10)
    axis.grid(True, which="major", alpha=0.25)

    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        frameon=False,
        ncol=len(families),
        fontsize=9,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.0),
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    _save_svg(figure, path)
    plt.close(figure)


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

    probe_binary = FAMILIES["OU-III"]
    probe_input = inputs[RECORDS[0]]
    settings = build_settings(probe_binary, probe_input)

    # The banner is record- and family-independent, so probe each setting once.
    banners: dict[str, dict[str, float]] = {}
    for setting in settings:
        banners[setting.label] = (
            probe_setting(probe_binary, probe_input, setting)
            if setting.engine_on else {}
        )

    rows: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        tasks = [
            executor.submit(
                run_one, family, binary, record, inputs[record], setting,
                banners[setting.label], args.window_sec,
            )
            for family, binary in FAMILIES.items()
            for setting in settings
            for record in RECORDS
        ]
        for index, future in enumerate(as_completed(tasks), start=1):
            row = future.result()
            rows.append(row)
            print(
                f"[{index:4d}/{len(tasks)}] {row['family']:6s} {row['label']:22s} "
                f"{row['spectrum']:9s} Hs={row['hs_m']:4.2f} "
                f"Z={row['disp_z_rms_m']:.5f} 3D={row['disp_3d_rms_m']:.5f} "
                f"pitch={row['pitch_rms_deg']:.4f} yaw={row['yaw_rms_deg']:.4f}"
            )

    family_order = {family: index for index, family in enumerate(FAMILIES)}
    label_order = {setting.label: index for index, setting in enumerate(settings)}
    spectrum_order = {"JONSWAP": 0, "PM-Stokes": 1}
    rows.sort(
        key=lambda row: (
            family_order[row["family"]],
            label_order[row["label"]],
            spectrum_order[row["spectrum"]],
            row["hs_m"],
        )
    )
    summaries = summarize(rows, settings)
    commit = source_commit()

    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "engine_noise_runs.csv", rows, ROW_FIELDS)
    write_csv(out / "engine_noise_summary.csv", summaries, SUMMARY_FIELDS)
    (out / "engine_noise_report.md").write_text(
        markdown_report(summaries, args.window_sec, commit, not args.no_plots),
        encoding="utf-8",
    )
    outputs = [
        "engine_noise_runs.csv",
        "engine_noise_summary.csv",
        "engine_noise_report.md",
    ]
    if not args.no_plots:
        write_speed_plot(out / SPEED_PLOT_NAME, summaries)
        write_mechanism_plot(out / MECHANISM_PLOT_NAME, summaries)
        outputs.extend((SPEED_PLOT_NAME, MECHANISM_PLOT_NAME))

    manifest = {
        "study": "engine-noise degradation",
        "source_commit": commit,
        "simulation_data": "oceanography-waves-lib v1.1.3",
        "families": list(FAMILIES),
        "records": [record.__dict__ for record in RECORDS],
        "record_count_per_family": len(RECORDS),
        "settings": [
            {
                "arm": setting.arm,
                "label": setting.label,
                "rpm": setting.rpm,
                "level_mps2": setting.level_mps2 if setting.engine_on else 0.0,
                "bandwidth_hz": (
                    setting.bandwidth_hz if setting.engine_on else None
                ),
                "recorded_rms_mps2": banners[setting.label].get(
                    "recorded_rms_mps2", 0.0
                ),
                "hull_rms_mps2": banners[setting.label].get("hull_rms_mps2", 0.0),
            }
            for setting in settings
        ],
        "total_replays": len(rows),
        "window_sec": args.window_sec,
        "magnetometer_enabled": True,
        "sensor_noise_injected": True,
        "engine": {
            "cylinders": 3,
            "stroke": 4,
            "gear_ratio": 2.6,
            "blades": 3,
            "cruise_rpm": CRUISE_RPM,
            "nominal_level_mps2": NOMINAL_LEVEL_MPS2,
            "nominal_bandwidth_hz": NOMINAL_BANDWIDTH_HZ,
        },
        "filter_covariance_and_tuning_unchanged": True,
        "outputs": outputs,
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Wrote engine-noise degradation study to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
