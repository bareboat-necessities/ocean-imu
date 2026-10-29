#!/usr/bin/env python3
"""Off-CoG IMU installation study of the experimental CARMA filters.

The versioned wave records give rigid-body specific force at the vessel centre
of gravity (CG).  An IMU bolted at body offset r measures that force plus

    alpha x r + omega x (omega x r),

which the simulator applies to the accelerometer truth before sensor
corruption (``W3D_IMU_LEVER_ARM_M``), exactly as in the OU-III study
(tools/ou3_lever_arm_study.py, docs/ou-iii-imu-lever-arm.md).  The filters are
CARMA3 (carma3-sim's A41) and TFG-CARMA3 (tfg_carma3-sim's T41), replayed by
tests/kalman_carma/carma_lever_arm-sim with the configuration of their gated
replays.

Arms, over the same eight seas, noise realization and scoring:

  baseline   IMU at the CG;
  unmodeled  IMU off the CG, no correction;
  exact      the simulator removes the installed term from truth kinematics
             (oracle: what any model could recover);
  wrapper    SeaStateFusionFilter_CARMA3::setImuLeverArmBody, the correction
             the CARMA sketches deploy: the measured rate through a two-pole
             low-pass and a second-order difference, removed at the wrapper's
             input so every consumer sees CoG specific force;
  mekf       ablation: the same model in the inner filter's accelerometer
             update only (OU-III's placement);
  mekf_ou3   ablation: OU-III's own model in OU-III's placement (one-pole
             first difference of the rate, alpha term in the bias Jacobian).

Offsets are given in the sketches' body frame (x forward, y starboard, z down)
and passed to the simulator through its z-up map.  A sweep varies the
wrapper's derivative corner on the vertical arm, where the installed term is
largest.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Iterable, Sequence

import ou3_lever_arm_study as ou3
import ou_validation as ouv
from sim_dataset import input_provenance

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "reports" / "results" / "carma_lever_arm_study"
# A smoke run must not replace the committed full-study evidence.
DEFAULT_SMOKE_OUTPUT_DIR = Path(tempfile.gettempdir()) / "carma_lever_arm_study_smoke"
BINARY = ROOT / "tests" / "kalman_carma" / "carma_lever_arm-sim"
DEFAULT_WINDOW_SEC = 900.0

RECORDS = ou3.RECORDS

# Filter -> carma_lever_arm-sim variant.
FILTERS: dict[str, str] = {"carma3": "A41", "tfg_carma3": "T41"}
FILTER_LABELS = {"carma3": "CARMA3", "tfg_carma3": "TFG-CARMA3"}

# Body NED (x forward, y starboard, z down), the frame of the sketches'
# IMU_LEVER_ARM_*_M.  The vertical arm points up: an IMU above the CG.
AXES: dict[str, tuple[float, float, float]] = {
    "x-forward": (1.0, 0.0, 0.0),
    "y-starboard": (0.0, 1.0, 0.0),
    "z-up": (0.0, 0.0, -1.0),
}
AXIS_LABEL = {"x-forward": "forward", "y-starboard": "starboard", "z-up": "above"}
AXIS_WHERE = {"x-forward": "forward of", "y-starboard": "starboard of", "z-up": "above"}
DISTANCES_M = (0.10, 0.20, 0.30, 0.50, 1.00, 2.00)

# Arm -> (simulator input-stage model, estimator placement or None).
ARMS: dict[str, tuple[str, str | None]] = {
    "unmodeled": ("none", None),
    "exact": ("exact", None),
    "wrapper": ("none", "wrapper"),
    "mekf": ("none", "mekf"),
    "mekf_ou3": ("none", "ou3"),
}
ARM_LABELS = {
    "baseline": "IMU at CG",
    "unmodeled": "unmodeled",
    "exact": "exact (oracle)",
    "wrapper": "wrapper (deployed)",
    "mekf": "inner filter",
    "mekf_ou3": "inner filter, OU-III model",
}

DEFAULT_CUTOFF_HZ = 15.0
SWEEP_CUTOFFS_HZ = (1.0, 2.0, 5.0, 10.0, 15.0, 25.0, 50.0, 100.0)
SWEEP_AXIS = "z-up"
SWEEP_DISTANCE_M = 1.00

RMS_FIELDS = ou3.RMS_FIELDS
ROW_FIELDS = ou3.ROW_FIELDS
RATIO_FIELDS = ("disp_3d_ratio_to_baseline", "tilt_ratio_to_baseline", "yaw_ratio_to_baseline")
# The headline also publishes the worst single sea; the cutoff sweep does not.
SUMMARY_COMPARE_FIELDS = RATIO_FIELDS + ("disp_3d_ratio_worst_sea", "tilt_ratio_worst_sea")
COMPARE_TOLERANCE = 0.01
# A removed fraction is reported only where the uncorrected excess clears the
# run-to-run scatter of a correction; below that it is a ratio of noise.
EXCESS_FLOOR = 0.05

FIGURES = (
    "carma_lever_arm_penalty.svg",
    "carma_lever_arm_tilt.svg",
    "carma_lever_arm_yaw.svg",
    "carma_lever_arm_sea_state.svg",
    "carma_lever_arm_cutoff.svg",
)
# The per-sea figure: the direction each channel suffers most.
SEA_STATE_CASES = (("disp_3d_rms_m", "z-up", "3-D displacement"), ("max_tilt_rms_deg", "x-forward", "max roll/pitch"))
SEA_STATE_DISTANCE_M = 1.00


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, default=None,
                   help=f"default {DEFAULT_OUTPUT_DIR} for --mode full, {DEFAULT_SMOKE_OUTPUT_DIR} for smoke")
    p.add_argument("--window-sec", type=float, default=DEFAULT_WINDOW_SEC)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--mode", choices=("smoke", "full"), default="full")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument(
        "--compare-to",
        type=Path,
        default=None,
        help="directory holding a previous lever_arm_summary.csv and lever_arm_cutoff_summary.csv; "
        "fail if any published ratio moved by more than --compare-tolerance",
    )
    p.add_argument("--compare-tolerance", type=float, default=COMPARE_TOLERANCE)
    return p.parse_args()


def ned_to_zu(v: Sequence[float]) -> tuple[float, float, float]:
    """The simulator's body z-up frame from body NED (src/ahrs/FrameConversions.h)."""
    return (v[1], v[0], -v[2])


def arm_env(arm: str, unit_ned: Sequence[float], distance_m: float, cutoff_hz: float) -> dict[str, str]:
    model, placement = ARMS[arm]
    r_zu = ned_to_zu(ou3.scale(distance_m, unit_ned))
    env = {
        "W3D_IMU_LEVER_ARM_M": ",".join(f"{v:.9g}" for v in r_zu),
        "W3D_IMU_LEVER_ARM_MODEL": model,
        "W3D_IMU_LEVER_ARM_CUTOFF_HZ": f"{cutoff_hz:.9g}",
    }
    if placement is not None:
        env["W3D_IMU_LEVER_ARM_FILTER"] = placement
    return env


def invoke(variant: str, input_path: Path, window_sec: float, env_extra: dict[str, str] | None) -> str:
    import subprocess

    if not BINARY.exists():
        raise FileNotFoundError(f"missing simulator binary: {BINARY}")
    env = os.environ.copy()
    for name in (
        "W3D_IMU_LEVER_ARM_M",
        "W3D_IMU_LEVER_ARM_MODEL",
        "W3D_IMU_LEVER_ARM_CUTOFF_HZ",
        "W3D_IMU_LEVER_ARM_FILTER",
        "W3D_IMU_LEVER_ARM_ALPHA_TAU_S",
    ):
        env.pop(name, None)
    env.update(
        {
            "W3D_CARMA_LEVER_VARIANT": variant,
            "W3D_WRITE_TIMESERIES": "0",
            "W3D_VALIDATION_WINDOW_SEC": f"{window_sec:.9g}",
        }
    )
    if env_extra:
        env.update(env_extra)
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
        raise RuntimeError(f"simulator exit {p.returncode}:\n" + "\n".join(p.stdout.splitlines()[-30:]))
    return p.stdout


def run_case(
    filt: str,
    record: ou3.Record,
    source: Path,
    arm: str,
    axis: str,
    distance_m: float,
    window_sec: float,
    cutoff_hz: float = DEFAULT_CUTOFF_HZ,
) -> dict[str, Any]:
    env = None if arm == "baseline" else arm_env(arm, AXES[axis], distance_m, cutoff_hz)
    stdout = invoke(FILTERS[filt], source, window_sec, env)
    if arm != "baseline" and ARMS[arm][1] is not None and "IMU_LEVER_ARM_FILTER " not in stdout:
        raise RuntimeError(f"{filt} {arm}: the estimator's lever arm was not installed")
    metrics = ouv.parse_validation_metrics(stdout)
    row: dict[str, Any] = {
        "filter": filt,
        "mode": arm,
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
    row.update(ou3.parse_lever_arm_result(stdout))
    # The simulator only sees its own input stage; the estimator's correction
    # happens after it, so its residual there is not the model's residual.
    if arm != "baseline" and ARMS[arm][1] is not None:
        row["residual_rms_mps2"] = float("nan")
    row["cutoff_hz"] = cutoff_hz if arm in ("wrapper", "mekf") else float("nan")
    return row


def max_tilt(item: dict[str, Any]) -> float:
    return max(float(item["roll_rms_deg"]), float(item["pitch_rms_deg"]))


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for filt in FILTERS:
        frows = [r for r in rows if r["filter"] == filt]
        if not frows:
            continue
        groups: dict[tuple[str, str, float], list[dict[str, Any]]] = {}
        for row in frows:
            groups.setdefault((row["mode"], row["axis"], float(row["distance_m"])), []).append(row)
        baseline = groups[("baseline", "cg", 0.0)]
        base = {f: ou3.pooled_rms(baseline, f) for f in RMS_FIELDS}
        base_by_sea = {(r["spectrum"], float(r["hs_m"])): r for r in baseline}
        base_tilt = max(base["roll_rms_deg"], base["pitch_rms_deg"])
        items: list[dict[str, Any]] = []
        for key in sorted(groups):
            arm, axis, distance = key
            rs = groups[key]
            item: dict[str, Any] = {"filter": filt, "mode": arm, "axis": axis, "distance_m": distance,
                                    "records": len(rs)}
            for field in RMS_FIELDS:
                item[field] = ou3.pooled_rms(rs, field)
            item["max_tilt_rms_deg"] = max_tilt(item)
            item["disp_3d_ratio_to_baseline"] = item["disp_3d_rms_m"] / base["disp_3d_rms_m"]
            item["disp_z_ratio_to_baseline"] = item["disp_z_rms_m"] / base["disp_z_rms_m"]
            item["tilt_ratio_to_baseline"] = item["max_tilt_rms_deg"] / base_tilt
            item["yaw_ratio_to_baseline"] = item["yaw_rms_deg"] / base["yaw_rms_deg"]
            # Pooling weights the largest seas, where the installed term is
            # smallest against the wave signal; the worst single sea is
            # carried alongside.
            per_sea = [(r, base_by_sea[(r["spectrum"], float(r["hs_m"]))]) for r in rs]
            item["disp_3d_ratio_worst_sea"] = max(
                float(r["disp_3d_rms_m"]) / float(b["disp_3d_rms_m"]) for r, b in per_sea)
            item["tilt_ratio_worst_sea"] = max(max_tilt(r) / max_tilt(b) for r, b in per_sea)
            item["installed_rms_mps2"] = ou3.pooled_rms(rs, "installed_rms_mps2")
            item["residual_rms_mps2"] = ou3.pooled_rms(rs, "residual_rms_mps2")
            item["excess_removed_fraction"] = float("nan")
            item["tilt_excess_removed_fraction"] = float("nan")
            items.append(item)
        unmodeled = {(s["axis"], s["distance_m"]): s for s in items if s["mode"] == "unmodeled"}
        for item in items:
            if item["mode"] in ("baseline", "unmodeled"):
                continue
            reference = unmodeled.get((item["axis"], item["distance_m"]))
            if reference is None:
                continue
            for field, b, name in (
                ("disp_3d_rms_m", base["disp_3d_rms_m"], "excess_removed_fraction"),
                ("max_tilt_rms_deg", base_tilt, "tilt_excess_removed_fraction"),
            ):
                excess = reference[field] - b
                if excess > EXCESS_FLOOR * b:
                    item[name] = (reference[field] - item[field]) / excess
        out += items
    return out


def summarize_sweep(rows: list[dict[str, Any]], summaries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for filt in FILTERS:
        base = ou3.pick([s for s in summaries if s["filter"] == filt], "baseline", "cg", 0.0)
        if base is None:
            continue
        groups: dict[float, list[dict[str, Any]]] = {}
        for row in rows:
            if row["filter"] == filt:
                groups.setdefault(float(row["cutoff_hz"]), []).append(row)
        for cutoff in sorted(groups):
            rs = groups[cutoff]
            item: dict[str, Any] = {"filter": filt, "cutoff_hz": cutoff, "axis": SWEEP_AXIS,
                                    "distance_m": SWEEP_DISTANCE_M, "records": len(rs)}
            for field in RMS_FIELDS:
                item[field] = ou3.pooled_rms(rs, field)
            item["max_tilt_rms_deg"] = max_tilt(item)
            item["disp_3d_ratio_to_baseline"] = item["disp_3d_rms_m"] / base["disp_3d_rms_m"]
            item["tilt_ratio_to_baseline"] = item["max_tilt_rms_deg"] / base["max_tilt_rms_deg"]
            item["yaw_ratio_to_baseline"] = item["yaw_rms_deg"] / base["yaw_rms_deg"]
            out.append(item)
    return out


def pick(summaries: Iterable[dict[str, Any]], filt: str, arm: str, axis: str, distance: float):
    return ou3.pick([s for s in summaries if s["filter"] == filt], arm, axis, distance)


def markdown_report(
    summaries: list[dict[str, Any]],
    sweep: list[dict[str, Any]],
    window_sec: float,
    commit: str | None,
    full: bool,
) -> str:
    fmt = ou3.fmt

    def pct(value: Any) -> str:
        f = float(value)
        return "n/a" if not math.isfinite(f) else f"{100.0 * f:.0f}%"

    lines = [
        "# CARMA IMU lever-arm installation study",
        "",
        "The wave records define motion at the vessel centre of gravity (CG).  An IMU",
        "rigidly mounted at body offset $r$ additionally measures",
        "`alpha x r + omega x (omega x r)`; the simulator applies that term to the",
        "accelerometer truth before sensor corruption.  Nothing else changes between",
        "arms: noise realization, CARMA configuration, startup, adaptation and scoring",
        "are those of the gated `carma3-sim` (A41) and `tfg_carma3-sim` (T41) replays.",
        "",
        "| Arm | Correction |",
        "| --- | --- |",
        "| baseline | none needed: IMU at the CG |",
        "| unmodeled | none |",
        "| exact | simulator removes the installed term from truth kinematics (oracle) |",
        "| wrapper | `setImuLeverArmBody`, the deployed correction: measured rate, "
        f"{DEFAULT_CUTOFF_HZ:.0f} Hz two-pole low-pass, second-order difference, at the wrapper input |",
        "| mekf | ablation: the same model inside the inner filter's accelerometer update |",
        "| mekf_ou3 | ablation: OU-III's model (one-pole first difference, alpha bias Jacobian) "
        "inside the inner filter |",
        "",
        "Offsets are in the sketches' body frame: forward, starboard, and above the CG.",
        f"Scoring uses the trailing **{window_sec:.0f} s** of each 1200 s record; every",
        "entry is the RMS pooled over the eight seas; ratios are to the filter's own",
        "CG baseline; the worst-sea columns divide each sea by its own baseline and",
        "keep the largest.  The inner-filter arms report no residual: the simulator cannot",
        "see a correction made inside the estimator.  The removed fraction of the",
        f"uncorrected 3-D excess is shown where that excess exceeds {100 * EXCESS_FLOOR:.0f} % of the baseline.",
        "",
    ]
    if commit:
        lines += [f"Source commit: `{commit}`.", ""]
    for filt in FILTERS:
        fs = [s for s in summaries if s["filter"] == filt]
        if not fs:
            continue
        base = pick(fs, filt, "baseline", "cg", 0.0)
        lines += [
            f"## {FILTER_LABELS[filt]}",
            "",
            f"CG baseline: 3-D displacement {fmt(base['disp_3d_rms_m'])} m, max roll/pitch "
            f"{fmt(base['max_tilt_rms_deg'])} deg, yaw {fmt(base['yaw_rms_deg'])} deg.",
            "",
            "| Arm | Axis | Offset [cm] | 3D / CG | 3D / CG, worst sea | Tilt / CG | Tilt / CG, worst sea "
            "| Yaw / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
        order = {a: i for i, a in enumerate(("unmodeled", *[a for a in ARMS if a != "unmodeled"]))}
        rows = sorted((s for s in fs if s["mode"] != "baseline"),
                      key=lambda s: (list(AXES).index(s["axis"]), float(s["distance_m"]), order[s["mode"]]))
        for s in rows:
            lines.append(
                f"| {s['mode']} | {s['axis']} | {100 * float(s['distance_m']):.0f} | "
                f"{fmt(s['disp_3d_ratio_to_baseline'], 3)}x | {fmt(s['disp_3d_ratio_worst_sea'], 3)}x | "
                f"{fmt(s['tilt_ratio_to_baseline'], 3)}x | {fmt(s['tilt_ratio_worst_sea'], 3)}x | "
                f"{fmt(s['yaw_ratio_to_baseline'], 3)}x | {fmt(s['installed_rms_mps2'], 4)} | "
                f"{fmt(s['residual_rms_mps2'], 4)} | {pct(s['excess_removed_fraction'])} |"
            )
        lines.append("")
    if sweep:
        lines += [
            "## Derivative corner of the deployed correction",
            "",
            f"Wrapper arm, IMU {100 * SWEEP_DISTANCE_M:.0f} cm above the CG, pooled over the same seas.",
            "",
            "| Filter | Corner [Hz] | 3D / CG | Tilt / CG | Yaw / CG |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
        for s in sweep:
            lines.append(
                f"| {FILTER_LABELS[s['filter']]} | {float(s['cutoff_hz']):.0f} | {fmt(s['disp_3d_ratio_to_baseline'], 3)}x | "
                f"{fmt(s['tilt_ratio_to_baseline'], 3)}x | {fmt(s['yaw_ratio_to_baseline'], 3)}x |"
            )
        lines.append("")
    lines += [
        "## Reading the tables",
        "",
        "The installed force is deterministic, but it combines with each record's",
        "existing estimator residual, so an arm can score slightly below its CG",
        "baseline; that is not evidence that the installation helps.  The exact arm",
        "bounds what any model can recover.  docs/carma-imu-lever-arm.md interprets",
        "the results.",
        "",
        f"Study matrix: {'8 records, 2 filters, 3 axes, 6 offsets, 5 arms' if full else 'smoke subset'}.",
    ]
    return "\n".join(lines) + "\n"


HEADLINE_DISTANCES_M = (0.30, 1.00, 2.00)


def headline_markdown(summaries: list[dict[str, Any]], sweep: list[dict[str, Any]]) -> str:
    """The table docs/carma-imu-lever-arm.md quotes, rendered from the summaries."""
    fmt = ou3.fmt
    lines: list[str] = []
    for filt in FILTERS:
        lines += [
            f"**{FILTER_LABELS[filt]}** (RMS / CG baseline: eight seas pooled, and the worst single sea)",
            "",
            "| IMU | Uncorrected 3-D | worst sea | Corrected 3-D | worst sea | Corrected tilt, worst sea "
            "| Corrected yaw | Inner-filter tilt, worst sea | OU-III model tilt |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
        for d in HEADLINE_DISTANCES_M:
            for axis in AXES:
                cell = {arm: pick(summaries, filt, arm, axis, d) for arm in ("unmodeled", "wrapper", "mekf", "mekf_ou3")}
                if any(v is None for v in cell.values()):
                    continue
                lines.append(
                    f"| {100 * d:.0f} cm {AXIS_LABEL[axis]} "
                    f"| {fmt(cell['unmodeled']['disp_3d_ratio_to_baseline'], 2)} "
                    f"| {fmt(cell['unmodeled']['disp_3d_ratio_worst_sea'], 2)} "
                    f"| {fmt(cell['wrapper']['disp_3d_ratio_to_baseline'], 3)} "
                    f"| {fmt(cell['wrapper']['disp_3d_ratio_worst_sea'], 3)} "
                    f"| {fmt(cell['wrapper']['tilt_ratio_worst_sea'], 2)} "
                    f"| {fmt(cell['wrapper']['yaw_ratio_to_baseline'], 2)} "
                    f"| {fmt(cell['mekf']['tilt_ratio_worst_sea'], 2)} "
                    f"| {fmt(cell['mekf_ou3']['tilt_ratio_to_baseline'], 2)} |"
                )
        lines.append("")
    lines += [
        f"Derivative corner of the corrected arm, IMU {100 * SWEEP_DISTANCE_M:.0f} cm {AXIS_LABEL[SWEEP_AXIS]} the CG "
        "(3-D RMS / CG baseline):",
        "",
        "| Corner [Hz] | " + " | ".join(FILTER_LABELS[f] for f in FILTERS) + " |",
        "| ---: | " + " | ".join("---:" for _ in FILTERS) + " |",
    ]
    for cutoff in sorted({float(s["cutoff_hz"]) for s in sweep}):
        cells = []
        for filt in FILTERS:
            match = [s for s in sweep if s["filter"] == filt and float(s["cutoff_hz"]) == cutoff]
            cells.append(fmt(match[0]["disp_3d_ratio_to_baseline"], 3) if match else "n/a")
        lines.append(f"| {cutoff:.0f} | " + " | ".join(cells) + " |")
    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

ARM_STYLE = {
    "unmodeled": {"color": "#c1272d", "marker": "o", "linestyle": "-"},
    "exact": {"color": "#1a7f37", "marker": "^", "linestyle": ":"},
    "wrapper": {"color": "#0072b2", "marker": "s", "linestyle": "--"},
    "mekf": {"color": "#e69f00", "marker": "D", "linestyle": "-."},
    "mekf_ou3": {"color": "#7a3fb0", "marker": "v", "linestyle": (0, (1, 1))},
}


def write_ratio_grid(path: Path, summaries: list[dict[str, Any]], field: str, ylabel: str, title: str) -> None:
    plt = ou3.reproducible_pyplot()
    fig, axarr = plt.subplots(len(FILTERS), len(AXES), figsize=(7.2, 4.4), sharex=True, sharey=True)
    for row, filt in enumerate(FILTERS):
        for col, axis in enumerate(AXES):
            ax = axarr[row][col]
            for arm, style in ARM_STYLE.items():
                xs, ys = [0.0], [1.0]
                for d in DISTANCES_M:
                    match = pick(summaries, filt, arm, axis, d)
                    if match is None:
                        continue
                    xs.append(100.0 * d)
                    ys.append(float(match[field]))
                ax.plot(xs, ys, label=ARM_LABELS[arm], markersize=3, linewidth=1.0, **style)
            ax.axhline(1.0, color="#666666", linewidth=0.8)
            ax.set_xscale("symlog", linthresh=10.0)
            ax.set_xticks([0, 10, 30, 100, 200])
            ax.set_xticklabels(["0", "10", "30", "100", "200"])
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            if row == 0:
                ax.set_title(f"IMU {AXIS_WHERE[axis]} CG", fontsize=8)
            if col == 0:
                ax.set_ylabel(f"{FILTER_LABELS[filt]}\n{ylabel}", fontsize=7.5)
    fig.supxlabel("IMU offset from CG [cm]", fontsize=8, y=0.1)
    fig.suptitle(title, fontsize=9)
    handles, labels = axarr[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, fontsize=6.5, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.11, right=0.99, bottom=0.2, top=0.88, wspace=0.08, hspace=0.12)
    ou3._save(fig, path)
    plt.close(fig)


def write_sea_state_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    """Per-sea ratios at 1 m: the pooled RMS hides the calm seas."""
    plt = ou3.reproducible_pyplot()
    arms = ("unmodeled", "wrapper", "mekf")
    fig, axarr = plt.subplots(len(FILTERS), len(SEA_STATE_CASES), figsize=(7.2, 4.6), sharex=True)
    width = 0.27
    for row, filt in enumerate(FILTERS):
        frows = [r for r in rows if r["filter"] == filt]
        for col, (field, axis, title) in enumerate(SEA_STATE_CASES):
            ax = axarr[row][col]
            seas: list[tuple[float, str]] = []
            top = 1.0
            for index, arm in enumerate(arms):
                seas, ratios = ou3._per_sea_ratio(frows, arm, axis, SEA_STATE_DISTANCE_M, field)
                top = max([top] + [v for v in ratios if math.isfinite(v)])
                ax.bar([i + (index - 1) * width for i in range(len(seas))], ratios, width=width,
                       color=ARM_STYLE[arm]["color"], label=ARM_LABELS[arm])
            ax.axhline(1.0, color="#666666", linewidth=0.8)
            ax.set_ylim(0.0, top * 1.08)
            ax.set_xticks(range(len(seas)))
            ax.set_xticklabels([f"{ou3.SPECTRUM_SHORT[spectrum]} {hs:g}" for hs, spectrum in seas], rotation=90,
                               fontsize=6.5)
            ax.tick_params(labelsize=7)
            ax.grid(True, axis="y", alpha=0.3)
            if row == 0:
                ax.set_title(f"{title}, IMU {100 * SEA_STATE_DISTANCE_M:.0f} cm {AXIS_WHERE[axis]} CG", fontsize=8)
            if col == 0:
                ax.set_ylabel(f"{FILTER_LABELS[filt]}\nRMS / CG baseline", fontsize=7.5)
    fig.supxlabel("Incident spectrum and $H_s$ [m]", fontsize=8, y=0.09)
    handles, labels = axarr[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.1, right=0.99, bottom=0.22, top=0.93, wspace=0.15, hspace=0.12)
    ou3._save(fig, path)
    plt.close(fig)


def write_cutoff_plot(path: Path, sweep: list[dict[str, Any]], summaries: list[dict[str, Any]]) -> None:
    plt = ou3.reproducible_pyplot()
    fields = (
        ("disp_3d_ratio_to_baseline", "3-D displacement"),
        ("tilt_ratio_to_baseline", "max roll/pitch"),
        ("yaw_ratio_to_baseline", "yaw"),
    )
    fig, axarr = plt.subplots(1, len(fields), figsize=(7.2, 2.6), sharex=True)
    for ax, (field, title) in zip(axarr, fields):
        for filt, style in zip(FILTERS, ("-", "--")):
            pts = [s for s in sweep if s["filter"] == filt]
            ax.plot([float(s["cutoff_hz"]) for s in pts], [float(s[field]) for s in pts], linestyle=style,
                    marker="s", markersize=3, color=ARM_STYLE["wrapper"]["color"],
                    label=f"{FILTER_LABELS[filt]} wrapper")
            ref = pick(summaries, filt, "unmodeled", SWEEP_AXIS, SWEEP_DISTANCE_M)
            if ref is not None:
                ax.axhline(float(ref[field]), color=ARM_STYLE["unmodeled"]["color"], linestyle=style, linewidth=0.9,
                           label=f"{FILTER_LABELS[filt]} unmodeled")
        ax.axhline(1.0, color="#666666", linewidth=0.8)
        ax.axvline(DEFAULT_CUTOFF_HZ, color="#666666", linewidth=0.8, linestyle=":")
        ax.set_xscale("log")
        ax.set_title(title, fontsize=8)
        ax.grid(True, alpha=0.3, which="both")
        ax.tick_params(labelsize=7)
    axarr[0].set_ylabel("RMS / CG baseline", fontsize=8)
    fig.supxlabel(f"Derivative low-pass corner [Hz] (IMU {100 * SWEEP_DISTANCE_M:.0f} cm above CG; "
                  f"dotted: default {DEFAULT_CUTOFF_HZ:.0f} Hz)", fontsize=7.5, y=0.12)
    handles, labels = axarr[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=6.5, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.33, top=0.9, wspace=0.28)
    ou3._save(fig, path)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs must be >= 1")
    if not (math.isfinite(args.window_sec) and args.window_sec > 0):
        raise SystemExit("--window-sec must be positive")

    full = args.mode == "full"
    records = RECORDS if full else (RECORDS[1],)
    axes = list(AXES) if full else ["x-forward", "z-up"]
    distances = DISTANCES_M if full else (0.30, 1.00)
    cutoffs = SWEEP_CUTOFFS_HZ if full else (5.0, 15.0)
    sources = {r: ou3.find_record(args.data_dir, r) for r in records}
    dataset = input_provenance(sources.values())

    out = args.output_dir or (DEFAULT_OUTPUT_DIR if full else DEFAULT_SMOKE_OUTPUT_DIR)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    cases: list[tuple[str, ou3.Record, str, str, float, float]] = []
    sweep_cases: list[tuple[str, ou3.Record, str, str, float, float]] = []
    for filt in FILTERS:
        for record in records:
            cases.append((filt, record, "baseline", "cg", 0.0, DEFAULT_CUTOFF_HZ))
            for axis in axes:
                for distance in distances:
                    for arm in ARMS:
                        cases.append((filt, record, arm, axis, distance, DEFAULT_CUTOFF_HZ))
            for cutoff in cutoffs:
                sweep_cases.append((filt, record, "wrapper", SWEEP_AXIS, SWEEP_DISTANCE_M, cutoff))

    def launch(pool: ThreadPoolExecutor, batch) -> list[dict[str, Any]]:
        pending = {
            pool.submit(run_case, filt, rec, sources[rec], arm, axis, distance, args.window_sec, cutoff):
                (filt, rec, arm, axis, distance, cutoff)
            for filt, rec, arm, axis, distance, cutoff in batch
        }
        collected: list[dict[str, Any]] = []
        for future in as_completed(pending):
            filt, rec, arm, axis, distance, cutoff = pending[future]
            row = future.result()
            collected.append(row)
            print(f"DONE {filt} {rec.spectrum} Hs={rec.hs_m:.2f} {arm} {axis} {100 * distance:.0f}cm "
                  f"fc={cutoff:.0f} 3d={row['disp_3d_rms_m']:.6g}", flush=True)
        return collected

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows = launch(pool, cases)
        sweep_rows = launch(pool, sweep_cases)

    rows.sort(key=lambda r: (r["filter"], r["mode"], r["axis"], float(r["distance_m"]), r["spectrum"],
                             float(r["hs_m"])))
    sweep_rows.sort(key=lambda r: (r["filter"], float(r["cutoff_hz"]), r["spectrum"], float(r["hs_m"])))
    summaries = summarize(rows)
    sweep = summarize_sweep(sweep_rows, summaries)

    runs_path = out / "lever_arm_runs.csv"
    summary_path = out / "lever_arm_summary.csv"
    sweep_runs_path = out / "lever_arm_cutoff_runs.csv"
    sweep_path = out / "lever_arm_cutoff_summary.csv"
    report_path = out / "lever_arm_report.md"
    ou3.write_csv(runs_path, rows)
    ou3.write_csv(summary_path, summaries)
    ou3.write_csv(sweep_runs_path, sweep_rows)
    ou3.write_csv(sweep_path, sweep)
    report_path.write_text(markdown_report(summaries, sweep, args.window_sec, ou3.source_commit(), full))

    generated = [runs_path, summary_path, sweep_runs_path, sweep_path, report_path]
    if not args.no_plots and full:
        figures = dict(zip(FIGURES, (out / name for name in FIGURES)))
        write_ratio_grid(figures["carma_lever_arm_penalty.svg"], summaries, "disp_3d_ratio_to_baseline",
                         "3-D RMS / CG", "Pooled 3-D displacement error against IMU offset")
        write_ratio_grid(figures["carma_lever_arm_tilt.svg"], summaries, "tilt_ratio_to_baseline",
                         "tilt RMS / CG", "Pooled max roll/pitch error against IMU offset")
        write_ratio_grid(figures["carma_lever_arm_yaw.svg"], summaries, "yaw_ratio_to_baseline",
                         "yaw RMS / CG", "Pooled yaw error against IMU offset")
        write_sea_state_plot(figures["carma_lever_arm_sea_state.svg"], rows)
        write_cutoff_plot(figures["carma_lever_arm_cutoff.svg"], sweep, summaries)
        generated += list(figures.values())

    manifest = {
        "study": "CARMA IMU lever-arm installation",
        "source_commit": ou3.source_commit(),
        "simulation_data": "oceanography-waves-lib v1.2.3 vessel-rao-28ft",
        "simulation_provenance": dataset,
        "mode": args.mode,
        "scoring_window_sec": args.window_sec,
        "filters": FILTERS,
        "axes_body_ned": {a: AXES[a] for a in axes},
        "distances_m": distances,
        "arms": {a: {"input_stage_model": m, "estimator_placement": p} for a, (m, p) in ARMS.items()},
        "derivative_cutoff_hz": DEFAULT_CUTOFF_HZ,
        "cutoff_sweep_hz": cutoffs,
        "cutoff_sweep_case": {"axis": SWEEP_AXIS, "distance_m": SWEEP_DISTANCE_M},
        "model": "a_imu = a_cg + alpha x r + omega x (omega x r)",
        "files": {p.name: ou3.sha256(p) for p in generated},
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(report_path.read_text())

    if args.compare_to is not None:
        complaints = compare(args.compare_to / summary_path.name, summaries,
                             ("filter", "mode", "axis", "distance_m"), SUMMARY_COMPARE_FIELDS,
                             args.compare_tolerance)
        if sweep and (args.compare_to / sweep_path.name).exists():
            complaints += compare(args.compare_to / sweep_path.name, sweep, ("filter", "cutoff_hz"),
                                  RATIO_FIELDS, args.compare_tolerance)
        if complaints:
            print(f"\nRe-run disagrees with {args.compare_to} beyond {args.compare_tolerance:g}:", flush=True)
            for line in complaints:
                print(f"  {line}", flush=True)
            return 1
        print(f"\nRe-run agrees with {args.compare_to} within {args.compare_tolerance:g} on every published ratio.",
              flush=True)
    return 0


def compare(previous_path: Path, current: list[dict[str, Any]], keys: Sequence[str], fields: Sequence[str],
            tolerance: float) -> list[str]:
    return ou3.compare_summaries(ou3.read_csv(previous_path), current, keys, fields, tolerance)


if __name__ == "__main__":
    raise SystemExit(main())
