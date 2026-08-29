#!/usr/bin/env python3
"""OU-III IMU lever-arm installation study.

The versioned wave records describe rigid-body specific force at the vessel
centre of gravity (CG).  An IMU rigidly mounted at body-frame offset r measures
that force plus the rotational terms

    a_r = alpha x r + omega x (omega x r),

so a few decimetres of installation offset inject a deterministic,
attitude-correlated error into the one channel OU-III uses for attitude and
wave acceleration at the same time.

The installation and the filter's model of it live on opposite sides of the
sensor, and the simulator now implements them there rather than by rewriting
records:

  * ``W3D_IMU_LEVER_ARM_M`` moves the accelerometer truth from the CG to the
    sensor location, before sensor corruption;
  * ``W3D_IMU_LEVER_ARM_MODEL`` selects the lever-arm model the filter's input
    stage applies, after sensor corruption and immediately before fusion.

Four arms are replayed over the same eight seas, with the sensor-noise
realization, OU-III configuration, startup logic, adaptation, and scoring
identical in every one:

  baseline   IMU at the CG;
  unmodeled  IMU off the CG, no filter-side model;
  gyro       IMU off the CG, compensated from the measured (noisy, biased)
             rate through a band-limited derivative -- the deployable model;
  exact      IMU off the CG, compensated from the record's own angular
             kinematics -- the oracle bound on what any model can recover.

The ``exact`` arm is deliberately an upper bound: it answers how much of the
installation penalty is deterministically recoverable when r and the rotational
kinematics are known exactly.  The ``gyro`` arm answers the deployable
question, and the cutoff sweep shows why its one design parameter matters.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import ou_validation as ouv

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "reports" / "results" / "ou3_lever_arm_study"
DOC_DIR = ROOT / "doc" / "kalman_ou_iii"
BINARY = ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim"
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

AXES: dict[str, tuple[float, float, float]] = {
    "x-athwartships": (1.0, 0.0, 0.0),
    "y-fore-aft": (0.0, 1.0, 0.0),
    "z-vertical": (0.0, 0.0, 1.0),
}
DISTANCES_M = (0.10, 0.20, 0.30)

# Modelling arms, in the order they are reported.  The value is what the
# simulator's W3D_IMU_LEVER_ARM_MODEL accepts.
MODES: dict[str, str] = {
    "unmodeled": "none",
    "gyro": "gyro",
    "exact": "exact",
}
MODE_LABELS = {
    "baseline": "IMU at CG",
    "unmodeled": "off-CG, unmodeled",
    "gyro": "off-CG, gyro-derived model",
    "exact": "off-CG, exact model",
}

# The derivative band the deployable model runs ahead of its difference.  The
# sweep is the study's answer to "why 15 Hz": below it the low-pass phase lag
# dominates, above it the differentiated gyro noise does.
SWEEP_CUTOFFS_HZ = (1.0, 2.0, 5.0, 10.0, 15.0, 25.0, 50.0, 100.0)
SWEEP_AXIS = "y-fore-aft"
SWEEP_DISTANCE_M = 0.30
DEFAULT_CUTOFF_HZ = 15.0

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
    "mode",
    "axis",
    "distance_m",
    "spectrum",
    "hs_m",
    "input",
    *RMS_FIELDS,
    "disp_z_ref_rms_m",
    "disp_z_pct_refrms",
    "disp_z_pct_hs",
    "tau_applied_s",
    "sigma_applied_mps2",
)

FIGURES = (
    "ou3_lever_arm_penalty.svg",
    "ou3_lever_arm_tilt.svg",
    "ou3_lever_arm_sea_state.svg",
    "ou3_lever_arm_mechanism.svg",
    "ou3_lever_arm_cutoff.svg",
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    p.add_argument("--window-sec", type=float, default=DEFAULT_WINDOW_SEC)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--mode", choices=("smoke", "full"), default="full")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument(
        "--mirror-doc",
        action="store_true",
        help="copy the article figures into doc/kalman_ou_iii/ byte for byte",
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# Rigid-body kinematics, kept here so the study's own arithmetic is testable
# without the simulator.  The simulator owns the deployed copy.
# ---------------------------------------------------------------------------


def cross(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def add(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(s: float, a: Sequence[float]) -> tuple[float, float, float]:
    return (s * a[0], s * a[1], s * a[2])


def lever_acceleration(
    omega_radps: Sequence[float],
    alpha_radps2: Sequence[float],
    r_body_m: Sequence[float],
) -> tuple[float, float, float]:
    """Rigid-body acceleration of a sensor at r relative to the CG."""
    return add(
        cross(alpha_radps2, r_body_m),
        cross(omega_radps, cross(omega_radps, r_body_m)),
    )


def find_record(data_dir: Path, record: Record) -> Path:
    matches = list(data_dir.rglob(record.filename))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"expected exactly one {record.filename} under {data_dir}, found {len(matches)}"
        )
    return matches[0]


def source_commit() -> str | None:
    p = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    return (p.stdout.strip() or None) if p.returncode == 0 else None


def invoke(
    input_path: Path,
    window_sec: float,
    lever_env: dict[str, str] | None,
) -> str:
    if not BINARY.exists():
        raise FileNotFoundError(f"missing simulator binary: {BINARY}")
    env = os.environ.copy()
    # Drop any inherited lever-arm settings so a run is defined entirely by
    # the case it is meant to be.
    for name in (
        "W3D_IMU_LEVER_ARM_M",
        "W3D_IMU_LEVER_ARM_MODEL",
        "W3D_IMU_LEVER_ARM_CUTOFF_HZ",
    ):
        env.pop(name, None)
    env.update(
        {
            "W3D_WRITE_TIMESERIES": "0",
            "W3D_VALIDATION_WINDOW_SEC": f"{window_sec:.9g}",
            "W3D_COLLECT_ALL_GATES": "1",
        }
    )
    if lever_env:
        env.update(lever_env)
    p = subprocess.run(
        [str(BINARY), "--input", str(input_path.resolve())],
        cwd=BINARY.parent,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if p.returncode not in (0, 1):
        raise RuntimeError(
            f"simulator exit {p.returncode}:\n" + "\n".join(p.stdout.splitlines()[-30:])
        )
    return p.stdout


def parse_lever_arm_result(stdout: str) -> dict[str, float]:
    """Reads the simulator's IMU_LEVER_ARM_RESULT diagnostic line.

    Returns the installed and residual specific-force RMS, i.e. what the
    installation added and what the filter's model failed to remove.  A
    baseline run emits no such line and reports zeros.
    """
    out = {"installed_rms_mps2": 0.0, "residual_rms_mps2": 0.0}
    for line in stdout.splitlines():
        if not line.startswith("IMU_LEVER_ARM_RESULT "):
            continue
        for token in line.split()[1:]:
            key, _, value = token.partition("=")
            if key in out:
                out[key] = float(value)
    return out


def lever_env(
    unit: Sequence[float],
    distance_m: float,
    model: str,
    cutoff_hz: float,
) -> dict[str, str]:
    r = scale(distance_m, unit)
    return {
        "W3D_IMU_LEVER_ARM_M": ",".join(f"{v:.9g}" for v in r),
        "W3D_IMU_LEVER_ARM_MODEL": model,
        "W3D_IMU_LEVER_ARM_CUTOFF_HZ": f"{cutoff_hz:.9g}",
    }


def run_case(
    record: Record,
    source: Path,
    mode: str,
    axis: str,
    distance_m: float,
    unit: Sequence[float],
    window_sec: float,
    cutoff_hz: float = DEFAULT_CUTOFF_HZ,
) -> dict[str, Any]:
    env = (
        None
        if mode == "baseline"
        else lever_env(unit, distance_m, MODES[mode], cutoff_hz)
    )
    stdout = invoke(source, window_sec, env)
    metrics = ouv.parse_validation_metrics(stdout)
    row: dict[str, Any] = {
        "mode": mode,
        "axis": axis,
        "distance_m": distance_m,
        "spectrum": record.spectrum,
        "hs_m": record.hs_m,
        "input": record.filename,
    }
    for field in ROW_FIELDS:
        if field in row:
            continue
        value = metrics.get(field, float("nan"))
        row[field] = float(value) if isinstance(value, (float, int)) else value
    row.update(parse_lever_arm_result(stdout))
    row["cutoff_hz"] = cutoff_hz if mode == "gyro" else float("nan")
    return row


def pooled_rms(rows: Iterable[dict[str, Any]], field: str) -> float:
    vals = [float(r[field]) for r in rows if math.isfinite(float(r[field]))]
    return math.sqrt(sum(v * v for v in vals) / len(vals)) if vals else float("nan")


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str, float], list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(
            (row["mode"], row["axis"], float(row["distance_m"])), []
        ).append(row)
    baseline = groups[("baseline", "cg", 0.0)]
    base3d = pooled_rms(baseline, "disp_3d_rms_m")
    basez = pooled_rms(baseline, "disp_z_rms_m")
    basetilt = max(
        pooled_rms(baseline, "roll_rms_deg"),
        pooled_rms(baseline, "pitch_rms_deg"),
    )
    out: list[dict[str, Any]] = []
    for key in sorted(groups, key=lambda k: (k[0], k[1], k[2])):
        mode, axis, distance = key
        rs = groups[key]
        item: dict[str, Any] = {
            "mode": mode,
            "axis": axis,
            "distance_m": distance,
            "records": len(rs),
        }
        for field in RMS_FIELDS:
            item[field] = pooled_rms(rs, field)
        item["disp_3d_ratio_to_baseline"] = item["disp_3d_rms_m"] / base3d
        item["disp_z_ratio_to_baseline"] = item["disp_z_rms_m"] / basez
        tilt = max(item["roll_rms_deg"], item["pitch_rms_deg"])
        item["max_tilt_rms_deg"] = tilt
        item["tilt_ratio_to_baseline"] = tilt / basetilt
        item["installed_rms_mps2"] = pooled_rms(rs, "installed_rms_mps2")
        item["residual_rms_mps2"] = pooled_rms(rs, "residual_rms_mps2")
        # Fraction of the unmodeled excess this arm removes.  1.0 is a full
        # return to the CG baseline; 0.0 is no model at all.
        item["excess_removed_fraction"] = float("nan")
        item["tilt_excess_removed_fraction"] = float("nan")
        out.append(item)

    unmodeled = {
        (s["axis"], s["distance_m"]): s for s in out if s["mode"] == "unmodeled"
    }
    for item in out:
        if item["mode"] not in ("gyro", "exact"):
            continue
        reference = unmodeled.get((item["axis"], item["distance_m"]))
        if reference is None:
            continue
        for field, base, name in (
            ("disp_3d_rms_m", base3d, "excess_removed_fraction"),
            ("max_tilt_rms_deg", basetilt, "tilt_excess_removed_fraction"),
        ):
            excess = reference[field] - base
            # A fraction of an excess that is itself run-to-run scatter says
            # nothing, and printing 280% of nothing is worse than saying so.
            if excess > 0.01 * base:
                item[name] = (reference[field] - item[field]) / excess
    return out


def summarize_sweep(rows: list[dict[str, Any]], base3d: float) -> list[dict[str, Any]]:
    groups: dict[float, list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(float(row["cutoff_hz"]), []).append(row)
    out: list[dict[str, Any]] = []
    for cutoff in sorted(groups):
        rs = groups[cutoff]
        item: dict[str, Any] = {
            "cutoff_hz": cutoff,
            "axis": SWEEP_AXIS,
            "distance_m": SWEEP_DISTANCE_M,
            "records": len(rs),
            "disp_3d_rms_m": pooled_rms(rs, "disp_3d_rms_m"),
            "installed_rms_mps2": pooled_rms(rs, "installed_rms_mps2"),
            "residual_rms_mps2": pooled_rms(rs, "residual_rms_mps2"),
        }
        item["disp_3d_ratio_to_baseline"] = item["disp_3d_rms_m"] / base3d
        item["residual_fraction"] = (
            item["residual_rms_mps2"] / item["installed_rms_mps2"]
            if item["installed_rms_mps2"] > 0.0
            else float("nan")
        )
        out.append(item)
    return out


def fmt(v: Any, digits: int = 4) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    return "n/a" if not math.isfinite(f) else f"{f:.{digits}f}"


def markdown_report(
    summaries: list[dict[str, Any]],
    sweep: list[dict[str, Any]],
    window_sec: float,
    commit: str | None,
    full: bool,
) -> str:
    lines = [
        "# OU-III IMU lever-arm installation study",
        "",
        "The wave records define motion at the vessel centre of gravity (CG).  A rigidly",
        "mounted IMU displaced by body-frame vector $r$ additionally measures",
        "`alpha x r + omega x (omega x r)`.  The simulator applies that term to the",
        "accelerometer truth before sensor corruption, and applies the filter's own",
        "lever-arm model after corruption, immediately before fusion.  Nothing else",
        "changes: the noise realization, OU-III configuration, adaptation, startup",
        "logic, pseudo-measurements, vibration guard, and scoring are identical in",
        "every arm.",
        "",
        "| Arm | What the filter receives |",
        "| --- | --- |",
        "| baseline | IMU at the CG |",
        "| unmodeled | off-CG specific force, no filter-side model |",
        "| gyro | off-CG specific force, compensated from the measured rate |",
        "| exact | off-CG specific force, compensated from truth kinematics |",
        "",
        "The canonical body directions are x = athwartships, y = fore-aft, and z = vertical.",
        f"Scoring uses the trailing **{window_sec:.0f} s** of each 1200 s record.",
        "",
        "The `exact` arm is an oracle bound on what any lever-arm model can recover.",
        "The `gyro` arm is the deployable one: it sees only the noisy, biased rate and",
        f"reconstructs `alpha` through a two-pole low-pass at {DEFAULT_CUTOFF_HZ:.0f} Hz",
        "followed by a causal second-order difference.",
        "",
    ]
    if commit:
        lines += [f"Source commit: `{commit}`.", ""]
    lines += [
        "## Pooled results",
        "",
        "| Arm | Axis | Offset [cm] | 3D disp [m] | 3D / CG | Max roll/pitch RMS [deg] "
        "| Tilt / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed "
        "| Tilt excess removed |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]

    def pct(value: Any) -> str:
        f = float(value)
        return "n/a" if not math.isfinite(f) else f"{100.0 * f:.1f}%"

    for s in summaries:
        lines.append(
            f"| {s['mode']} | {s['axis']} | {100*float(s['distance_m']):.0f} | "
            f"{fmt(s['disp_3d_rms_m'])} | {fmt(s['disp_3d_ratio_to_baseline'],3)}x | "
            f"{fmt(s['max_tilt_rms_deg'])} | {fmt(s['tilt_ratio_to_baseline'],3)}x | "
            f"{fmt(s['installed_rms_mps2'],4)} | {fmt(s['residual_rms_mps2'],4)} | "
            f"{pct(s['excess_removed_fraction'])} | "
            f"{pct(s['tilt_excess_removed_fraction'])} |"
        )
    if sweep:
        lines += [
            "",
            "## Derivative band of the deployable model",
            "",
            f"Fore-aft arm at {100*SWEEP_DISTANCE_M:.0f} cm, pooled over the same eight seas.",
            "",
            "| Cutoff [Hz] | 3D disp [m] | 3D / CG | Residual [m/s^2] | Residual / installed |",
            "| ---: | ---: | ---: | ---: | ---: |",
        ]
        for s in sweep:
            lines.append(
                f"| {s['cutoff_hz']:.0f} | {fmt(s['disp_3d_rms_m'])} | "
                f"{fmt(s['disp_3d_ratio_to_baseline'],3)}x | "
                f"{fmt(s['residual_rms_mps2'],4)} | {fmt(s['residual_fraction'],3)} |"
            )
    lines += [
        "",
        "## Interpretation",
        "",
        "The comparison isolates one installation effect: rigid-body rotational",
        "acceleration at the sensor location.  No filter covariance, OU schedule,",
        "pseudo-measurement, vibration guard, startup rule, or quality threshold is",
        "retuned for the off-CG cases.  The exact-model arm returns to the CG baseline",
        "to numerical precision, so the whole unmodeled penalty is deterministic and",
        "recoverable rather than an intrinsic OU-III limit.  The gyro-derived arm shows",
        "how much of that is available to firmware that has only the measured rate, and",
        "the cutoff sweep shows that its single design parameter is two-sided: too",
        "narrow a derivative band and the low-pass phase lag misaligns a correction of",
        "the right size, too wide and differentiated gyro noise dominates.",
        "",
        "Read the pooled ratios above with the per-sea figure beside them.  An RMS over",
        "all eight seas is dominated by the largest, where the injected term is smallest",
        "relative to the wave signal, so pooling understates a penalty that is severe in",
        "the mildest seas for displacement and in the steepest seas for attitude.",
        "",
        f"Study matrix: {'8 records, 3 axes, 3 offsets, 3 modelling arms' if full else 'smoke subset'}.",
        "",
        "Figures are written here and mirrored byte-for-byte into `doc/kalman_ou_iii/`",
        "for the article.",
    ]
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

MODE_STYLE = {
    "unmodeled": {"color": "#c1272d", "marker": "o", "linestyle": "-"},
    "gyro": {"color": "#0072b2", "marker": "s", "linestyle": "--"},
    "exact": {"color": "#1a7f37", "marker": "^", "linestyle": ":"},
}
AXIS_LABEL = {
    "x-athwartships": "athwartships",
    "y-fore-aft": "fore-aft",
    "z-vertical": "vertical",
}
# Eight sea labels have to share one axis; the full names collide.
SPECTRUM_SHORT = {"JONSWAP": "JON", "PM-Stokes": "PM"}


def reproducible_pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matplotlib.rcParams["svg.hashsalt"] = "ocean-imu-ou3-lever-arm"
    return plt


def _save(fig: Any, path: Path) -> None:
    fig.savefig(path, format="svg", metadata={"Date": None})


def pick(
    summaries: list[dict[str, Any]], mode: str, axis: str, distance: float
) -> dict[str, Any] | None:
    for s in summaries:
        if (
            s["mode"] == mode
            and s["axis"] == axis
            and abs(float(s["distance_m"]) - distance) < 1e-9
        ):
            return s
    return None


def write_ratio_plot(
    path: Path,
    summaries: list[dict[str, Any]],
    field: str,
    ylabel: str,
    title: str,
) -> None:
    """One panel per canonical direction: penalty, then its two removals."""
    plt = reproducible_pyplot()
    fig, axarr = plt.subplots(1, 3, figsize=(9.6, 3.6), sharey=True)
    for ax, axis in zip(axarr, AXES):
        for mode, style in MODE_STYLE.items():
            xs = [0.0]
            ys = [1.0]
            for d in DISTANCES_M:
                match = pick(summaries, mode, axis, d)
                if match is None:
                    continue
                xs.append(100.0 * d)
                ys.append(float(match[field]))
            ax.plot(xs, ys, label=MODE_LABELS[mode], markersize=5, **style)
        ax.axhline(1.0, color="#666666", linewidth=0.8)
        ax.set_title(AXIS_LABEL[axis], fontsize=10)
        ax.set_xlabel("IMU offset from CG [cm]")
        ax.grid(True, alpha=0.3)
    axarr[0].set_ylabel(ylabel)
    handles, labels = axarr[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=3,
        fontsize=8,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0.0, 0.07, 1.0, 1.0))
    _save(fig, path)
    plt.close(fig)


def _per_sea_ratio(
    rows: list[dict[str, Any]],
    mode: str,
    axis: str,
    distance: float,
    field: str,
) -> tuple[list[tuple[float, str]], list[float]]:
    """Per-sea error of one arm, divided by the CG baseline for that sea.

    Pooling hides the shape of this effect: an RMS over all eight seas is
    dominated by the largest, where the rotational term is smallest relative
    to the wave signal.  The ratio has to be formed sea by sea.
    """
    seas = sorted({(float(r["hs_m"]), r["spectrum"]) for r in rows})

    def value(row: dict[str, Any]) -> float:
        if field == "max_tilt_rms_deg":
            return max(float(row["roll_rms_deg"]), float(row["pitch_rms_deg"]))
        return float(row[field])

    def find(want_mode: str, want_axis: str, want_distance: float, sea) -> float:
        hs, spectrum = sea
        for r in rows:
            if (
                r["mode"] == want_mode
                and r["axis"] == want_axis
                and abs(float(r["distance_m"]) - want_distance) < 1e-9
                and abs(float(r["hs_m"]) - hs) < 1e-9
                and r["spectrum"] == spectrum
            ):
                return value(r)
        return float("nan")

    ratios = []
    for sea in seas:
        base = find("baseline", "cg", 0.0, sea)
        ratios.append(find(mode, axis, distance, sea) / base if base > 0.0 else float("nan"))
    return seas, ratios


def write_sea_state_plot(
    path: Path,
    rows: list[dict[str, Any]],
    disp_axis: str,
    tilt_axis: str,
    distance: float,
) -> None:
    """Per-sea penalty and its removal, for both channels the effect reaches."""
    plt = reproducible_pyplot()
    fig, panels = plt.subplots(1, 2, figsize=(9.6, 3.8))
    width = 0.26
    for panel, (field, axis, title) in zip(
        panels,
        (
            (
                "disp_3d_rms_m",
                disp_axis,
                f"3-D displacement, {AXIS_LABEL[disp_axis]} arm",
            ),
            (
                "max_tilt_rms_deg",
                tilt_axis,
                f"Max roll/pitch, {AXIS_LABEL[tilt_axis]} arm",
            ),
        ),
    ):
        seas: list[tuple[float, str]] = []
        lowest = 1.0
        for index, mode in enumerate(("unmodeled", "gyro", "exact")):
            seas, ratios = _per_sea_ratio(rows, mode, axis, distance, field)
            lowest = min([lowest] + [v for v in ratios if math.isfinite(v)])
            positions = [i + (index - 1) * width for i in range(len(seas))]
            panel.bar(
                positions,
                ratios,
                width=width,
                label=MODE_LABELS[mode],
                color=MODE_STYLE[mode]["color"],
            )
        panel.axhline(1.0, color="#666666", linewidth=0.8)
        panel.set_xticks(range(len(seas)))
        panel.set_xticklabels(
            [f"{SPECTRUM_SHORT[spectrum]}\n$H_s$ {hs:g}" for hs, spectrum in seas],
            fontsize=7,
        )
        panel.set_ylabel("RMS / CG baseline")
        panel.set_title(title, fontsize=10)
        panel.grid(True, axis="y", alpha=0.3)
        # A ratio plot anchored at zero hides the very deviations it is for.
        panel.set_ylim(bottom=min(0.92, lowest - 0.02))
    handles, labels = panels[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=3,
        fontsize=8,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.suptitle(
        f"Per-sea penalty with the IMU {100*distance:.0f} cm off the CG",
        fontsize=11,
    )
    fig.tight_layout(rect=(0.0, 0.07, 1.0, 0.94))
    _save(fig, path)
    plt.close(fig)


def write_mechanism_plot(path: Path, summaries: list[dict[str, Any]]) -> None:
    """What the installation injects, and what each model leaves behind."""
    plt = reproducible_pyplot()
    fig, (left, right) = plt.subplots(1, 2, figsize=(9.6, 3.8))

    for axis in AXES:
        xs = [0.0]
        ys = [0.0]
        for d in DISTANCES_M:
            match = pick(summaries, "unmodeled", axis, d)
            if match is None:
                continue
            xs.append(100.0 * d)
            ys.append(float(match["installed_rms_mps2"]))
        left.plot(xs, ys, marker="o", markersize=5, label=AXIS_LABEL[axis])
    left.set_xlabel("IMU offset from CG [cm]")
    left.set_ylabel("Injected specific force RMS [m/s$^2$]")
    left.set_title("What the installation adds", fontsize=10)
    left.grid(True, alpha=0.3)
    left.legend(fontsize=8)

    # As a fraction of what was injected, not in absolute units: the exact
    # model's residual is identically zero, which no logarithmic axis can
    # draw, and the fraction is the quantity the comparison is about anyway.
    axis_names = list(AXES)
    positions = range(len(axis_names))
    width = 0.27
    for index, mode in enumerate(("unmodeled", "gyro", "exact")):
        heights = []
        for axis in axis_names:
            match = pick(summaries, mode, axis, DISTANCES_M[-1])
            if match is None or not (float(match["installed_rms_mps2"]) > 0.0):
                heights.append(float("nan"))
                continue
            heights.append(
                float(match["residual_rms_mps2"]) / float(match["installed_rms_mps2"])
            )
        offsets = [p + (index - 1) * width for p in positions]
        right.bar(
            offsets,
            heights,
            width=width,
            label=MODE_LABELS[mode],
            color=MODE_STYLE[mode]["color"],
        )
        for x, height in zip(offsets, heights):
            if not math.isfinite(height):
                continue
            right.annotate(
                f"{height:.2f}" if height >= 0.005 else "0",
                xy=(x, height),
                xytext=(0, 2),
                textcoords="offset points",
                ha="center",
                fontsize=7,
            )
    right.set_xticks(list(positions))
    right.set_xticklabels([AXIS_LABEL[a] for a in axis_names], fontsize=9)
    right.set_ylabel("Residual / injected specific force")
    right.set_ylim(0.0, 1.5)
    right.set_title(
        f"What each model leaves at {100*DISTANCES_M[-1]:.0f} cm", fontsize=10
    )
    right.grid(True, axis="y", alpha=0.3)
    right.legend(fontsize=8, loc="upper center", ncol=1, framealpha=0.9)

    fig.tight_layout()
    _save(fig, path)
    plt.close(fig)


def write_cutoff_plot(
    path: Path,
    sweep: list[dict[str, Any]],
    unmodeled_ratio: float,
    exact_ratio: float,
) -> None:
    """The deployable model's one design parameter, and both ways to lose."""
    plt = reproducible_pyplot()
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    xs = [float(s["cutoff_hz"]) for s in sweep]
    ys = [float(s["disp_3d_ratio_to_baseline"]) for s in sweep]
    ax.plot(
        xs,
        ys,
        color=MODE_STYLE["gyro"]["color"],
        marker="s",
        markersize=5,
        label="gyro-derived model",
    )
    ax.axhline(
        unmodeled_ratio,
        color=MODE_STYLE["unmodeled"]["color"],
        linestyle="-",
        linewidth=1.2,
        label=MODE_LABELS["unmodeled"],
    )
    ax.axhline(
        exact_ratio,
        color=MODE_STYLE["exact"]["color"],
        linestyle=":",
        linewidth=1.4,
        label=MODE_LABELS["exact"],
    )
    ax.axvline(DEFAULT_CUTOFF_HZ, color="#666666", linewidth=0.8, linestyle="--")
    # Low on the axis: the legend owns the top-right corner.
    ax.annotate(
        f"deployed {DEFAULT_CUTOFF_HZ:.0f} Hz",
        xy=(DEFAULT_CUTOFF_HZ, min(ys)),
        xytext=(5, 6),
        textcoords="offset points",
        fontsize=8,
        color="#666666",
    )
    ax.set_xscale("log")
    ax.set_xlabel("Derivative low-pass corner [Hz]")
    ax.set_ylabel("Pooled 3-D RMS / CG baseline")
    ax.set_title(
        f"Derivative band of the deployable model "
        f"({AXIS_LABEL[SWEEP_AXIS]} arm, {100*SWEEP_DISTANCE_M:.0f} cm)",
        fontsize=11,
    )
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    _save(fig, path)
    plt.close(fig)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def worst_unmodeled_axis(
    summaries: list[dict[str, Any]], distance: float, field: str
) -> str:
    candidates = [
        s
        for s in summaries
        if s["mode"] == "unmodeled" and abs(float(s["distance_m"]) - distance) < 1e-9
    ]
    return max(candidates, key=lambda s: float(s[field]))["axis"]


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs must be >= 1")
    if not (math.isfinite(args.window_sec) and args.window_sec > 0):
        raise SystemExit("--window-sec must be positive")

    full = args.mode == "full"
    records = RECORDS if full else (RECORDS[1],)
    axes = (
        AXES
        if full
        else {
            "x-athwartships": AXES["x-athwartships"],
            "z-vertical": AXES["z-vertical"],
        }
    )
    distances = DISTANCES_M if full else (0.10, 0.30)
    cutoffs = SWEEP_CUTOFFS_HZ if full else (2.0, 15.0)
    sources = {r: find_record(args.data_dir, r) for r in records}

    out = args.output_dir
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    cases: list[tuple[Record, str, str, float, Sequence[float], float]] = []
    for record in records:
        cases.append((record, "baseline", "cg", 0.0, (0.0, 0.0, 0.0), DEFAULT_CUTOFF_HZ))
        for axis, unit in axes.items():
            for distance in distances:
                for mode in MODES:
                    cases.append(
                        (record, mode, axis, distance, unit, DEFAULT_CUTOFF_HZ)
                    )

    sweep_cases: list[tuple[Record, str, str, float, Sequence[float], float]] = []
    sweep_unit = AXES[SWEEP_AXIS]
    for record in records:
        for cutoff in cutoffs:
            sweep_cases.append(
                (record, "gyro", SWEEP_AXIS, SWEEP_DISTANCE_M, sweep_unit, cutoff)
            )

    def launch(
        pool: ThreadPoolExecutor,
        batch: list[tuple[Record, str, str, float, Sequence[float], float]],
    ) -> list[dict[str, Any]]:
        pending = {
            pool.submit(
                run_case, rec, sources[rec], mode, axis, distance, unit,
                args.window_sec, cutoff,
            ): (rec, mode, axis, distance, cutoff)
            for rec, mode, axis, distance, unit, cutoff in batch
        }
        collected: list[dict[str, Any]] = []
        for future in as_completed(pending):
            rec, mode, axis, distance, cutoff = pending[future]
            row = future.result()
            collected.append(row)
            print(
                f"DONE {rec.spectrum} Hs={rec.hs_m:.2f} {mode} {axis} "
                f"{100*distance:.0f}cm fc={cutoff:.0f} "
                f"3d={row['disp_3d_rms_m']:.6g}",
                flush=True,
            )
        return collected

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows = launch(pool, cases)
        sweep_rows = launch(pool, sweep_cases)

    rows.sort(
        key=lambda r: (
            r["mode"],
            r["axis"],
            float(r["distance_m"]),
            r["spectrum"],
            float(r["hs_m"]),
        )
    )
    sweep_rows.sort(
        key=lambda r: (float(r["cutoff_hz"]), r["spectrum"], float(r["hs_m"]))
    )
    summaries = summarize(rows)
    base3d = pooled_rms(
        [r for r in rows if r["mode"] == "baseline"], "disp_3d_rms_m"
    )
    sweep = summarize_sweep(sweep_rows, base3d)

    runs_path = out / "lever_arm_runs.csv"
    summary_path = out / "lever_arm_summary.csv"
    sweep_runs_path = out / "lever_arm_cutoff_runs.csv"
    sweep_path = out / "lever_arm_cutoff_summary.csv"
    report_path = out / "lever_arm_report.md"
    write_csv(runs_path, rows)
    write_csv(summary_path, summaries)
    write_csv(sweep_runs_path, sweep_rows)
    write_csv(sweep_path, sweep)
    report_path.write_text(
        markdown_report(summaries, sweep, args.window_sec, source_commit(), full)
    )

    generated = [runs_path, summary_path, sweep_runs_path, sweep_path, report_path]
    if not args.no_plots and full:
        # Displacement and attitude do not peak on the same installation
        # direction, so each panel gets the direction that maximizes it.
        worst_disp_axis = worst_unmodeled_axis(
            summaries, DISTANCES_M[-1], "disp_3d_ratio_to_baseline"
        )
        worst_tilt_axis = worst_unmodeled_axis(
            summaries, DISTANCES_M[-1], "tilt_ratio_to_baseline"
        )
        penalty = out / "ou3_lever_arm_penalty.svg"
        tilt = out / "ou3_lever_arm_tilt.svg"
        sea_state = out / "ou3_lever_arm_sea_state.svg"
        mechanism = out / "ou3_lever_arm_mechanism.svg"
        cutoff_plot = out / "ou3_lever_arm_cutoff.svg"
        write_ratio_plot(
            penalty,
            summaries,
            "disp_3d_ratio_to_baseline",
            "Pooled 3-D RMS / CG baseline",
            "Displacement penalty of an off-CG IMU, and its removal",
        )
        write_ratio_plot(
            tilt,
            summaries,
            "tilt_ratio_to_baseline",
            "Pooled max tilt RMS / CG baseline",
            "Attitude penalty of an off-CG IMU, and its removal",
        )
        write_sea_state_plot(
            sea_state, rows, worst_disp_axis, worst_tilt_axis, DISTANCES_M[-1]
        )
        write_mechanism_plot(mechanism, summaries)
        sweep_reference_unmodeled = pick(
            summaries, "unmodeled", SWEEP_AXIS, SWEEP_DISTANCE_M
        )
        sweep_reference_exact = pick(summaries, "exact", SWEEP_AXIS, SWEEP_DISTANCE_M)
        write_cutoff_plot(
            cutoff_plot,
            sweep,
            float(sweep_reference_unmodeled["disp_3d_ratio_to_baseline"]),
            float(sweep_reference_exact["disp_3d_ratio_to_baseline"]),
        )
        generated += [penalty, tilt, sea_state, mechanism, cutoff_plot]

        if args.mirror_doc:
            for name in FIGURES:
                shutil.copyfile(out / name, DOC_DIR / name)

    manifest = {
        "study": "OU-III IMU lever-arm installation",
        "source_commit": source_commit(),
        "simulation_data": "oceanography-waves-lib v1.1.3",
        "mode": args.mode,
        "scoring_window_sec": args.window_sec,
        "axes": axes,
        "distances_m": distances,
        "modes": MODES,
        "derivative_cutoff_hz": DEFAULT_CUTOFF_HZ,
        "cutoff_sweep_hz": cutoffs,
        "model": "a_imu = a_cg + alpha x r + omega x (omega x r)",
        "installation_stage": (
            "simulator applies the rigid-body term to the accelerometer truth "
            "before sensor corruption (W3D_IMU_LEVER_ARM_M)"
        ),
        "model_stage": (
            "simulator removes the filter's modelled term after sensor "
            "corruption, immediately before fusion (W3D_IMU_LEVER_ARM_MODEL)"
        ),
        "exact_arm": "oracle: truth angular kinematics; an upper bound, not a claim",
        "gyro_arm": (
            "deployable: two-pole low-pass on the measured rate followed by a "
            "causal second-order difference"
        ),
        "files": {p.name: sha256(p) for p in generated},
    }
    manifest_path = out / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(report_path.read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
