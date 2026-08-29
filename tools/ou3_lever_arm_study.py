#!/usr/bin/env python3
"""OU-III IMU lever-arm installation study.

The versioned wave records describe rigid-body specific force at the vessel CG.
For a fixed body-frame IMU offset r, the accelerometer at the sensor location
contains the additional rigid-body term

    a_r = alpha x r + omega x (omega x r).

This study rewrites only acc_bx/acc_by/acc_bz, leaving the translational truth,
attitude truth, gyro truth, magnetic field, stochastic sensor model, OU-III
configuration, startup logic, adaptation, and scoring unchanged.

Two matched arms are evaluated:
  * unmodeled: the off-CG acceleration reaches OU-III unchanged;
  * exact: the same deterministic term is removed before fusion, representing
    an ideal exact lever-arm model with exact angular kinematics.

The exact arm is intentionally an upper bound.  It answers how much of the
installation penalty is deterministically recoverable when r and the rotational
kinematics are known exactly; it does not claim that differentiating a noisy
MEMS gyro reproduces this bound in deployed hardware.
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
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import ou_validation as ouv

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "reports" / "results" / "ou3_lever_arm_study"
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
MODES = ("unmodeled", "exact")

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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    p.add_argument("--window-sec", type=float, default=DEFAULT_WINDOW_SEC)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--mode", choices=("smoke", "full"), default="full")
    p.add_argument("--no-plots", action="store_true")
    return p.parse_args()


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


def _omega(row: dict[str, str]) -> tuple[float, float, float]:
    return (float(row["gyro_x"]), float(row["gyro_y"]), float(row["gyro_z"]))


def _time(row: dict[str, str]) -> float:
    return float(row["time"])


def _alpha(
    prev: dict[str, str] | None,
    cur: dict[str, str],
    nxt: dict[str, str] | None,
) -> tuple[float, float, float]:
    if prev is None and nxt is None:
        return (0.0, 0.0, 0.0)
    if prev is None:
        a, b = cur, nxt
    elif nxt is None:
        a, b = prev, cur
    else:
        a, b = prev, nxt
    assert a is not None and b is not None
    dt = _time(b) - _time(a)
    if not (dt > 0.0 and math.isfinite(dt)):
        raise ValueError("wave record time must be strictly increasing")
    wa, wb = _omega(a), _omega(b)
    return (
        (wb[0] - wa[0]) / dt,
        (wb[1] - wa[1]) / dt,
        (wb[2] - wa[2]) / dt,
    )


def _rewrite_one(
    writer: csv.DictWriter,
    row: dict[str, str],
    prev: dict[str, str] | None,
    nxt: dict[str, str] | None,
    r_body_m: Sequence[float],
    exact: bool,
) -> None:
    alpha = _alpha(prev, row, nxt)
    omega = _omega(row)
    ar = lever_acceleration(omega, alpha, r_body_m)
    acc_cg = (float(row["acc_bx"]), float(row["acc_by"]), float(row["acc_bz"]))
    acc_sensor = add(acc_cg, ar)
    acc_in = add(acc_sensor, scale(-1.0, ar)) if exact else acc_sensor
    out = dict(row)
    out["acc_bx"], out["acc_by"], out["acc_bz"] = (f"{v:.17g}" for v in acc_in)
    writer.writerow(out)


def rewrite_record(
    source: Path,
    destination: Path,
    r_body_m: Sequence[float],
    exact: bool,
) -> None:
    """Create an off-CG record, optionally followed by exact compensation.

    The transform is streaming: only the previous/current/next samples are
    retained so several 20-minute records can be prepared in parallel without
    turning the study into a memory benchmark.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open(newline="") as src, destination.open("w", newline="") as dst:
        reader = csv.DictReader(src)
        if reader.fieldnames is None:
            raise ValueError(f"{source}: missing CSV header")
        required = {
            "time",
            "acc_bx",
            "acc_by",
            "acc_bz",
            "gyro_x",
            "gyro_y",
            "gyro_z",
        }
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"{source}: missing columns {sorted(missing)}")

        writer = csv.DictWriter(dst, fieldnames=reader.fieldnames, lineterminator="\n")
        writer.writeheader()
        iterator = iter(reader)
        first = next(iterator, None)
        if first is None:
            return
        second = next(iterator, None)
        if second is None:
            _rewrite_one(writer, first, None, None, r_body_m, exact)
            return

        _rewrite_one(writer, first, None, second, r_body_m, exact)
        prev, cur = first, second
        for nxt in iterator:
            _rewrite_one(writer, cur, prev, nxt, r_body_m, exact)
            prev, cur = cur, nxt
        _rewrite_one(writer, cur, prev, None, r_body_m, exact)


def source_commit() -> str | None:
    p = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    return (p.stdout.strip() or None) if p.returncode == 0 else None


def invoke(input_path: Path, window_sec: float) -> str:
    if not BINARY.exists():
        raise FileNotFoundError(f"missing simulator binary: {BINARY}")
    env = os.environ.copy()
    env.update(
        {
            "W3D_WRITE_TIMESERIES": "0",
            "W3D_VALIDATION_WINDOW_SEC": f"{window_sec:.9g}",
            "W3D_COLLECT_ALL_GATES": "1",
        }
    )
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


def run_one(
    record: Record,
    input_path: Path,
    mode: str,
    axis: str,
    distance_m: float,
    window_sec: float,
) -> dict[str, Any]:
    stdout = invoke(input_path, window_sec)
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
    return row


def run_case(
    record: Record,
    source: Path,
    mode: str,
    axis: str,
    distance_m: float,
    unit: Sequence[float],
    temp_root: Path,
    window_sec: float,
) -> dict[str, Any]:
    if mode == "baseline":
        return run_one(record, source, mode, axis, distance_m, window_sec)
    target = temp_root / mode / axis / f"{distance_m:.2f}" / record.filename
    rewrite_record(source, target, scale(distance_m, unit), exact=(mode == "exact"))
    try:
        return run_one(record, target, mode, axis, distance_m, window_sec)
    finally:
        target.unlink(missing_ok=True)


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
    window_sec: float,
    commit: str | None,
    full: bool,
) -> str:
    lines = [
        "# OU-III IMU lever-arm installation study",
        "",
        "The wave records define the vessel center of gravity (CG).  A rigidly mounted",
        "IMU displaced by body-frame vector $r$ measures the additional acceleration",
        "`alpha × r + omega × (omega × r)`.  The **unmodeled** arm injects that term",
        "into the accelerometer while OU-III remains unchanged.  The matched **exact**",
        "arm removes the same deterministic term immediately before fusion and is an",
        "ideal exact-model upper bound.  Ordinary simulator sensor noise remains enabled",
        "in both arms.",
        "",
        "The canonical body directions are x = athwartships, y = fore-aft, and z = vertical.",
        f"Scoring uses the trailing **{window_sec:.0f} s** of each 1200 s record.",
        "",
        "The exact arm assumes exact angular kinematics.  It therefore quantifies the",
        "recoverable installation penalty; it does not claim that a noisy finite-difference",
        "gyro derivative can attain the same result on hardware.",
        "",
    ]
    if commit:
        lines += [f"Source commit: `{commit}`.", ""]
    lines += [
        "## Pooled results",
        "",
        "| Mode | Axis | Offset [cm] | 3D disp [m] | 3D / CG | Z disp [m] | Z / CG | Max roll/pitch RMS [deg] | Tilt / CG | Yaw [deg] |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for s in summaries:
        lines.append(
            f"| {s['mode']} | {s['axis']} | {100*float(s['distance_m']):.0f} | "
            f"{fmt(s['disp_3d_rms_m'])} | {fmt(s['disp_3d_ratio_to_baseline'],3)}x | "
            f"{fmt(s['disp_z_rms_m'])} | {fmt(s['disp_z_ratio_to_baseline'],3)}x | "
            f"{fmt(s['max_tilt_rms_deg'])} | {fmt(s['tilt_ratio_to_baseline'],3)}x | "
            f"{fmt(s['yaw_rms_deg'])} |"
        )
    lines += [
        "",
        "## Interpretation",
        "",
        "The comparison isolates one installation effect: rigid-body rotational acceleration",
        "at the sensor location.  No filter covariance, OU schedule, pseudo-measurement,",
        "vibration guard, startup logic, or quality threshold is retuned for the off-CG cases.",
        "The exact-model result should return to the CG baseline to numerical precision; any",
        "remaining difference is CSV rounding / replay numerics rather than an unmodeled",
        "physical term in this experiment.",
        "",
        f"Study matrix: {'8 records, 3 axes, 3 offsets' if full else 'smoke subset'}.",
    ]
    return "\n".join(lines) + "\n"


def reproducible_pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matplotlib.rcParams["svg.hashsalt"] = "ocean-imu-ou3-lever-arm"
    return plt


def write_plot(
    path: Path,
    summaries: list[dict[str, Any]],
    field: str,
    ylabel: str,
) -> None:
    plt = reproducible_pyplot()
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for axis in AXES:
        xs = [0.0]
        ys = [1.0]
        for d in DISTANCES_M:
            match = next(
                s
                for s in summaries
                if s["mode"] == "unmodeled"
                and s["axis"] == axis
                and abs(s["distance_m"] - d) < 1e-9
            )
            xs.append(100.0 * d)
            ys.append(float(match[field]))
        ax.plot(xs, ys, marker="o", label=f"unmodeled {axis}")
    exact = [s for s in summaries if s["mode"] == "exact" and s["axis"] != "cg"]
    if exact:
        by_d = []
        for d in DISTANCES_M:
            vals = [float(s[field]) for s in exact if abs(s["distance_m"] - d) < 1e-9]
            if vals:
                by_d.append((100.0 * d, sum(vals) / len(vals)))
        if by_d:
            ax.plot(
                [0.0] + [p[0] for p in by_d],
                [1.0] + [p[1] for p in by_d],
                linestyle="--",
                marker="s",
                label="exact model (axis mean)",
            )
    ax.axhline(1.0, linewidth=1.0)
    ax.set_xlabel("IMU offset from CG [cm]")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, format="svg", metadata={"Date": None})
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


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs must be >= 1")
    if not (math.isfinite(args.window_sec) and args.window_sec > 0):
        raise SystemExit("--window-sec must be positive")

    records = RECORDS if args.mode == "full" else (RECORDS[1],)
    axes = (
        AXES
        if args.mode == "full"
        else {
            "x-athwartships": AXES["x-athwartships"],
            "z-vertical": AXES["z-vertical"],
        }
    )
    distances = DISTANCES_M if args.mode == "full" else (0.10, 0.30)
    sources = {r: find_record(args.data_dir, r) for r in records}

    out = args.output_dir
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    rows: list[dict[str, Any]] = []
    cases: list[tuple[Record, str, str, float, Sequence[float]]] = []
    for record in records:
        cases.append((record, "baseline", "cg", 0.0, (0.0, 0.0, 0.0)))
        for axis, unit in axes.items():
            for distance in distances:
                for mode in MODES:
                    cases.append((record, mode, axis, distance, unit))

    with tempfile.TemporaryDirectory(prefix="ou3-lever-arm-") as temp_s:
        temp = Path(temp_s)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            pending = {
                pool.submit(
                    run_case,
                    rec,
                    sources[rec],
                    mode,
                    axis,
                    distance,
                    unit,
                    temp,
                    args.window_sec,
                ): (rec, mode, axis, distance)
                for rec, mode, axis, distance, unit in cases
            }
            for future in as_completed(pending):
                rec, mode, axis, distance = pending[future]
                row = future.result()
                rows.append(row)
                print(
                    f"DONE {rec.spectrum} Hs={rec.hs_m:.2f} {mode} {axis} "
                    f"{100*distance:.0f}cm 3d={row['disp_3d_rms_m']:.6g}",
                    flush=True,
                )

    rows.sort(
        key=lambda r: (
            r["mode"],
            r["axis"],
            float(r["distance_m"]),
            r["spectrum"],
            float(r["hs_m"]),
        )
    )
    summaries = summarize(rows)
    runs_path = out / "lever_arm_runs.csv"
    summary_path = out / "lever_arm_summary.csv"
    report_path = out / "lever_arm_report.md"
    write_csv(runs_path, rows)
    write_csv(summary_path, summaries)
    report_path.write_text(
        markdown_report(
            summaries,
            args.window_sec,
            source_commit(),
            args.mode == "full",
        )
    )

    generated = [runs_path, summary_path, report_path]
    if not args.no_plots and args.mode == "full":
        p3 = out / "ou3_lever_arm_3d.svg"
        pt = out / "ou3_lever_arm_tilt.svg"
        write_plot(
            p3,
            summaries,
            "disp_3d_ratio_to_baseline",
            "Pooled 3-D RMS / CG baseline",
        )
        write_plot(
            pt,
            summaries,
            "tilt_ratio_to_baseline",
            "Pooled max tilt RMS / CG baseline",
        )
        generated += [p3, pt]

    manifest = {
        "study": "OU-III IMU lever-arm installation",
        "source_commit": source_commit(),
        "simulation_data": "oceanography-waves-lib v1.1.3",
        "mode": args.mode,
        "scoring_window_sec": args.window_sec,
        "axes": axes,
        "distances_m": distances,
        "model": "a_imu = a_cg + alpha x r + omega x (omega x r)",
        "exact_arm": (
            "subtracts the same truth-derived rigid-body term before simulator "
            "sensor corruption/fusion; ideal upper bound"
        ),
        "files": {p.name: sha256(p) for p in generated},
    }
    manifest_path = out / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(report_path.read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
