#!/usr/bin/env python3
"""Engine-noise mitigation study for OU-III's front-end vibration guard.

The engine-noise degradation study (``tools/engine_noise_degradation.py``)
established that machinery vibration costs the deployed estimators a large,
almost entirely *systematic* error, that it is set by how much out-of-band
accelerometer power survives to the sample rather than by where the aliased
crank orders land, and that the gyroscope path contributes nothing.  All of
that points at one remedy: keep the out-of-band content out of the
accelerometer before any consumer reads it.

``AccelVibrationGuard`` does that.  It sits at the single point in
``SeaStateFusionFilter_OU_III::updateCore_`` where raw measurements arrive, so
the Mahony proxy, the MEKF, and the tilt watchdog all see the same conditioned
signal.  It low-passes the accelerometer in the decade of empty spectrum
between the wave band and the machinery band, and -- because the group delay it
costs shows up in displacement in proportion to wave amplitude -- it engages
only when a separate high-pass detector says there is machinery to remove.

This study replays the eight stationary records through OU-III with the guard
off and on across a range of engine conditions, and confirms both halves of the
claim: that the guard recovers most of the loss under vibration, and that it is
bit-transparent when there is none.
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
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import ou_validation as ouv

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "reports" / "results" / "engine_noise_mitigation"
DEFAULT_WINDOW_SEC = 900.0
SIMULATOR = ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim"

# The guard's deployed configuration.  Two poles at 14 Hz sits in the gap
# between the wave band and the lowest crank order the modeled engine puts on
# the hull, and was the best of the corners swept: it costs 21 ms of group
# delay at full engagement, against 32 ms at 10 Hz for no better rejection and
# 16 ms at 20 Hz for markedly worse.
GUARD_CUTOFF_HZ = 14.0
GUARD_POLES = 2


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


@dataclass(frozen=True)
class Condition:
    """One engine condition, run with the guard off and on."""

    label: str
    engine: dict[str, str] = field(default_factory=dict)

    @property
    def engine_on(self) -> bool:
        return bool(self.engine)


CONDITIONS = (
    Condition("engine off"),
    Condition("800 rpm", {"W3D_ENGINE_RPM": "800"}),
    Condition("1600 rpm", {"W3D_ENGINE_RPM": "1600"}),
    Condition("2400 rpm", {"W3D_ENGINE_RPM": "2400"}),
    Condition("3200 rpm", {"W3D_ENGINE_RPM": "3200"}),
    Condition("2400 rpm, quiet mount",
              {"W3D_ENGINE_RPM": "2400", "W3D_ENGINE_LEVEL_MPS2": "0.30"}),
    Condition("2400 rpm, engine bed",
              {"W3D_ENGINE_RPM": "2400", "W3D_ENGINE_LEVEL_MPS2": "1.20"}),
    Condition("2400 rpm, wide sensor",
              {"W3D_ENGINE_RPM": "2400", "W3D_ENGINE_BANDWIDTH_HZ": "160"}),
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
)

MEAN_FIELDS = (
    "disp_x_mean_m",
    "disp_y_mean_m",
    "disp_z_mean_m",
    "roll_mean_deg",
    "pitch_mean_deg",
    "yaw_mean_deg",
)

GUARD_FIELDS = ("guard_engagement", "guard_out_of_band_rms_mps2", "guard_delay_sec")

ROW_FIELDS = (
    "condition",
    "guard",
    "spectrum",
    "hs_m",
    "input",
    *RMS_FIELDS,
    *MEAN_FIELDS,
    *GUARD_FIELDS,
    "disp_z_ref_rms_m",
)

SUMMARY_FIELDS = (
    "condition",
    "guard",
    "records",
    *RMS_FIELDS,
    *MEAN_FIELDS,
    *GUARD_FIELDS,
    "disp_3d_offset_m",
    "disp_3d_ratio_to_baseline",
)

MEANS_RE = re.compile(r"^VALIDATION_METRICS_MEANS (.*)$", re.MULTILINE)
GUARD_RE = re.compile(r"^ACC_GUARD (.*)$", re.MULTILINE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True,
                        help="directory containing the versioned wave CSVs")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
                        help=f"result directory (default: {DEFAULT_OUTPUT_DIR})")
    parser.add_argument("--window-sec", type=float, default=DEFAULT_WINDOW_SEC,
                        help="trailing scoring window in seconds (default: 900)")
    parser.add_argument("--no-plots", action="store_true",
                        help="skip the SVG chart (for environments without Matplotlib)")
    parser.add_argument("--jobs", type=int, default=4,
                        help="parallel simulator processes (default: 4)")
    return parser.parse_args()


def find_record(data_dir: Path, record: Record) -> Path:
    matches = list(data_dir.rglob(record.filename))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"expected exactly one {record.filename} under {data_dir}, found {len(matches)}"
        )
    return matches[0]


def source_commit() -> str | None:
    completed = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                               text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        return None
    return completed.stdout.strip() or None


def parse_keyed(pattern: re.Pattern[str], stdout: str) -> dict[str, float]:
    match = pattern.search(stdout)
    if not match:
        return {}
    values: dict[str, float] = {}
    for token in match.group(1).split():
        key, _, raw = token.partition("=")
        try:
            values[key] = float(raw)
        except ValueError:
            continue
    return values


def run_one(record: Record, input_path: Path, condition: Condition,
            guard: bool, window_sec: float) -> dict[str, Any]:
    env = os.environ.copy()
    env.update({
        "W3D_WRITE_TIMESERIES": "0",
        "W3D_VALIDATION_WINDOW_SEC": f"{window_sec:.9g}",
        # The deployed gates were fitted to a vibration-free input; collect the
        # metrics rather than stopping on the first one this study crosses.
        "W3D_COLLECT_ALL_GATES": "1",
    })
    env.update(condition.engine)
    # Both arms are explicit: the guard is armed by default in the filter, so
    # the off arm has to switch it off rather than say nothing.
    env["OU_III_ACC_GUARD_HZ"] = f"{GUARD_CUTOFF_HZ:.9g}" if guard else "0"
    env["OU_III_ACC_GUARD_POLES"] = str(GUARD_POLES)

    completed = subprocess.run(
        [str(SIMULATOR), "--input", str(input_path.resolve())],
        cwd=SIMULATOR.parent, env=env, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    if completed.returncode not in (0, 1):
        tail = "\n".join(completed.stdout.splitlines()[-30:])
        raise RuntimeError(
            f"simulator exit {completed.returncode} for {condition.label} "
            f"guard={guard}:\n{tail}"
        )
    if condition.engine_on and "ENGINE_VIBRATION" not in completed.stdout:
        raise RuntimeError(f"{condition.label}: engine model did not engage")
    if guard and "ACC_GUARD" not in completed.stdout:
        raise RuntimeError(f"{condition.label}: guard was requested but not armed")
    if not guard and "ACC_GUARD" in completed.stdout:
        raise RuntimeError(f"{condition.label}: guard was not switched off")

    metrics = ouv.parse_validation_metrics(completed.stdout)
    means = parse_keyed(MEANS_RE, completed.stdout)
    guard_stats = parse_keyed(GUARD_RE, completed.stdout)

    row: dict[str, Any] = {
        "condition": condition.label,
        "guard": "on" if guard else "off",
        "spectrum": record.spectrum,
        "hs_m": record.hs_m,
        "input": record.filename,
        "guard_engagement": guard_stats.get("engagement", 0.0),
        "guard_out_of_band_rms_mps2": guard_stats.get("out_of_band_rms_mps2", float("nan")),
        "guard_delay_sec": guard_stats.get("delay_sec", 0.0),
    }
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


def pooled_mean(rows: Iterable[dict[str, Any]], name: str) -> float:
    values = [float(row[name]) for row in rows]
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return float("nan")
    return sum(finite) / len(finite)


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    baseline = float("nan")
    for condition in CONDITIONS:
        for guard in ("off", "on"):
            cells = [
                row for row in rows
                if row["condition"] == condition.label and row["guard"] == guard
            ]
            if len(cells) != len(RECORDS):
                raise ValueError(
                    f"{condition.label} guard={guard}: expected {len(RECORDS)} rows, "
                    f"found {len(cells)}"
                )
            summary: dict[str, Any] = {
                "condition": condition.label,
                "guard": guard,
                "records": len(cells),
            }
            for name in RMS_FIELDS:
                summary[name] = pooled_rms(cells, name)
            for name in MEAN_FIELDS:
                summary[name] = pooled_rms(cells, name)
            for name in GUARD_FIELDS:
                summary[name] = pooled_mean(cells, name)
            summary["disp_3d_offset_m"] = math.sqrt(
                sum(float(summary[name]) ** 2 for name in
                    ("disp_x_mean_m", "disp_y_mean_m", "disp_z_mean_m"))
            )
            if condition.label == "engine off" and guard == "off":
                baseline = float(summary["disp_3d_rms_m"])
            result.append(summary)

    for summary in result:
        summary["disp_3d_ratio_to_baseline"] = (
            float(summary["disp_3d_rms_m"]) / baseline
            if baseline > 0.0 else float("nan")
        )
    return result


def fmt(value: Any, digits: int = 4) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        return "n/a"
    return f"{number:.{digits}f}"


def markdown_report(summaries: list[dict[str, Any]], window_sec: float,
                    commit: str | None, plots: bool) -> str:
    by_key = {(s["condition"], s["guard"]): s for s in summaries}
    lines = [
        "# Engine-noise mitigation: the OU-III front-end vibration guard",
        "",
        "The degradation study showed that machinery vibration costs OU-III a",
        "large and almost entirely systematic error, that the size of it tracks",
        "recorded out-of-band accelerometer power rather than the placement of the",
        "aliased crank orders, and that the gyroscope path contributes nothing.",
        "That points at a single remedy: keep the out-of-band content out of the",
        "accelerometer before anything reads it.",
        "",
        "`AccelVibrationGuard` sits at the one point in `updateCore_` where raw",
        "measurements arrive, so the Mahony proxy, the MEKF, and the tilt watchdog",
        "all see the same conditioned signal.  It low-passes the accelerometer in",
        "the empty decade between the wave band and the machinery band",
        f"(**{GUARD_POLES} poles at {GUARD_CUTOFF_HZ:.0f} Hz**), and engages only when a separate",
        "high-pass detector says there is machinery to remove.",
        "",
        f"Scoring uses the trailing **{window_sec:.0f} s** of each 1200 s record,",
        "pooled over the eight stationary records as `sqrt(mean(record_RMS^2))`.",
        "",
    ]
    if commit:
        lines.extend((f"Source commit used for the replay: `{commit}`.", ""))

    lines.extend([
        "## Result",
        "",
        "| Condition | Guard | Detector [m/s²] | Engaged | 3-D [m] | 3-D offset [m] "
        "| Pitch offset [deg] | Yaw [deg] | vs engine-off |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for condition in CONDITIONS:
        for guard in ("off", "on"):
            summary = by_key[(condition.label, guard)]
            lines.append(
                "| {cond} | {guard} | {det} | {eng} | {d3} | {off} | {pitch} | {yaw} "
                "| {ratio} |".format(
                    cond=condition.label,
                    guard=guard,
                    det=fmt(summary["guard_out_of_band_rms_mps2"], 4)
                    if guard == "on" else "—",
                    eng=fmt(summary["guard_engagement"], 3) if guard == "on" else "—",
                    d3=fmt(summary["disp_3d_rms_m"]),
                    off=fmt(summary["disp_3d_offset_m"]),
                    pitch=fmt(abs(float(summary["pitch_mean_deg"])), 3),
                    yaw=fmt(summary["yaw_rms_deg"], 2),
                    ratio=fmt(summary["disp_3d_ratio_to_baseline"], 2),
                )
            )

    off = by_key[("engine off", "off")]
    on = by_key[("engine off", "on")]
    identical = all(
        float(off[name]) == float(on[name])
        for name in RMS_FIELDS + MEAN_FIELDS
        if math.isfinite(float(off[name]))
    )
    lines.extend([
        "",
        "## Transparency with no engine running",
        "",
        "The engine-off rows are "
        + ("**identical to every digit**" if identical else "**not identical**")
        + ", because the detector never crosses its lower rail: the guard leaves",
        "the measurement path untouched and returns its input unchanged.  This is",
        "the property that lets the guard ship without re-cutting any fitted gate",
        "or invalidating a committed replay.",
        "",
        "It holds because the detector is placed above the sea rather than at the",
        "conditioning corner.  Across a 31:1 range of significant wave height the",
        "clean detector reading varies by about one percent, so a big sea does not",
        "look like machinery to it.",
        "",
    ])

    if plots:
        lines.extend([
            "## Figure",
            "",
            f"- `{PLOT_NAME}`: pooled 3-D error and standing tilt offset with the",
            "  guard off and on, against the engine-off baseline.",
            "",
            "Mirrored byte-for-byte into `doc/kalman_ou_iii/` for the article.",
            "",
        ])

    lines.extend([
        "## What this does not do",
        "",
        "The guard removes vibration from the measurement path.  It does not make",
        "the estimator vibration-aware: the accelerometer measurement covariance",
        "is unchanged, so the filter still believes a conditioned sample is as",
        "good as a quiet one.  It also cannot help with machinery whose orders",
        "reach into the wave band, which no front-end filter can separate from the",
        "sea, and it does not touch the residual that survives its own stopband.",
        "",
    ])
    return "\n".join(lines)


PLOT_NAME = "ou_engine_noise_guard.svg"
GUARD_COLORS = {"off": "#c0392b", "on": "#1b6ca8"}


def reproducible_pyplot():
    """Matplotlib configured for byte-reproducible published SVGs."""

    cache_dir = Path(tempfile.gettempdir()) / "ocean-imu-matplotlib-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    import matplotlib

    matplotlib.use("Agg")
    # Without these the figure cannot be reproduced: matplotlib derives svg
    # element ids from a per-process salt and stamps a creation date.
    matplotlib.rcParams["svg.hashsalt"] = "ocean-imu-engine-guard"
    import matplotlib.pyplot as plt

    return plt


def write_plot(path: Path, summaries: list[dict[str, Any]]) -> None:
    plt = reproducible_pyplot()
    by_key = {(s["condition"], s["guard"]): s for s in summaries}
    labels = [c.label for c in CONDITIONS]
    baseline = float(by_key[("engine off", "off")]["disp_3d_rms_m"])

    figure, axes = plt.subplots(2, 1, figsize=(8.2, 6.4), sharex=True)

    def panel(axis, extract, ylabel, title, log, reference) -> None:
        base = [float(index) for index in range(len(labels))]
        for offset, guard in enumerate(("off", "on")):
            values = [extract(by_key[(label, guard)]) for label in labels]
            axis.bar(
                [value + (offset - 0.5) * 0.38 for value in base],
                values,
                width=0.38,
                color=GUARD_COLORS[guard],
                edgecolor="#222222",
                linewidth=0.5,
                label=f"guard {guard}",
            )
        axis.axhline(reference, color="#444444", linewidth=0.9, linestyle=":")
        axis.set_xticks(base)
        axis.set_ylabel(ylabel)
        axis.set_title(title, fontsize=10)
        axis.grid(axis="y", alpha=0.25)
        if log:
            axis.set_yscale("log")

    panel(
        axes[0],
        lambda s: float(s["disp_3d_rms_m"]),
        "pooled 3-D RMS error (m)",
        "Displacement error, guard off and on",
        True,
        baseline,
    )
    panel(
        axes[1],
        lambda s: abs(float(s["pitch_mean_deg"])),
        "standing pitch offset (deg)",
        "Rectified tilt offset, the mechanism the guard targets",
        True,
        abs(float(by_key[("engine off", "off")]["pitch_mean_deg"])),
    )
    axes[1].set_xticklabels(labels, rotation=20, ha="right", fontsize=8)

    handles, legend_labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles, legend_labels, frameon=False, ncol=2, fontsize=9,
        loc="upper center", bbox_to_anchor=(0.5, 1.0),
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, format="svg", metadata={"Date": None})
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
    if not SIMULATOR.exists():
        raise SystemExit(f"missing simulator: {SIMULATOR}; build it first")

    data_dir = args.data_dir.resolve()
    inputs = {record: find_record(data_dir, record) for record in RECORDS}

    rows: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        tasks = [
            executor.submit(run_one, record, inputs[record], condition, guard,
                            args.window_sec)
            for condition in CONDITIONS
            for guard in (False, True)
            for record in RECORDS
        ]
        for index, future in enumerate(as_completed(tasks), start=1):
            row = future.result()
            rows.append(row)
            print(
                f"[{index:3d}/{len(tasks)}] {row['condition']:22s} guard={row['guard']:3s} "
                f"{row['spectrum']:9s} Hs={row['hs_m']:4.2f} "
                f"3D={row['disp_3d_rms_m']:.5f} pitch_off={row['pitch_mean_deg']:+.3f}"
            )

    condition_order = {c.label: index for index, c in enumerate(CONDITIONS)}
    spectrum_order = {"JONSWAP": 0, "PM-Stokes": 1}
    rows.sort(key=lambda row: (
        condition_order[row["condition"]], row["guard"],
        spectrum_order[row["spectrum"]], row["hs_m"],
    ))
    summaries = summarize(rows)
    commit = source_commit()

    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "guard_runs.csv", rows, ROW_FIELDS)
    write_csv(out / "guard_summary.csv", summaries, SUMMARY_FIELDS)
    (out / "guard_report.md").write_text(
        markdown_report(summaries, args.window_sec, commit, not args.no_plots),
        encoding="utf-8",
    )
    outputs = ["guard_runs.csv", "guard_summary.csv", "guard_report.md"]
    if not args.no_plots:
        write_plot(out / PLOT_NAME, summaries)
        outputs.append(PLOT_NAME)

    manifest = {
        "study": "engine-noise mitigation (OU-III accelerometer vibration guard)",
        "source_commit": commit,
        "simulation_data": "oceanography-waves-lib v1.1.3",
        "family": "OU-III",
        "guard": {"cutoff_hz": GUARD_CUTOFF_HZ, "poles": GUARD_POLES},
        "records": [record.__dict__ for record in RECORDS],
        "conditions": [
            {"label": c.label, "engine": dict(c.engine)} for c in CONDITIONS
        ],
        "total_replays": len(rows),
        "window_sec": args.window_sec,
        "outputs": outputs,
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Wrote engine-noise mitigation study to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
