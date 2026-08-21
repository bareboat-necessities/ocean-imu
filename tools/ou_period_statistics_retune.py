#!/usr/bin/env python3
"""Paired retuning study for the OU-III period/statistics pipeline.

Compile-time axes patch constructor defaults only in the disposable checkout;
runtime axes use existing simulator environment overrides. Every candidate is
compared with the axis baseline on identical records/seeds using paired log
ratios. The controlled 1.5 -> 4.0 m transition is also scored over blend and
recovery segments so a smoother cannot win merely by being slow.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import ou_validation as ouv  # noqa: E402

WINDOW_SEC = 900.0
DURATION_SEC = 1200.0
TRANSITION_START_SEC = 540.0
TRANSITION_END_SEC = 660.0
TRANSITION_RECOVER_END_SEC = 780.0
TRANSITION_END_HEIGHT_M = 4.0

STATIONARY = (
    ("JONSWAP", 0.27, "wave_data_jonswap_H0.270_L14.047_A30.00_P60.00.csv"),
    ("JONSWAP", 1.50, "wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv"),
    ("JONSWAP", 4.00, "wave_data_jonswap_H4.000_L112.766_A30.00_P30.00.csv"),
    ("JONSWAP", 8.50, "wave_data_jonswap_H8.500_L202.839_A-30.00_P72.00.csv"),
    ("PM-Stokes", 0.27, "wave_data_pmstokes_H0.270_L14.047_A30.00_P60.00.csv"),
    ("PM-Stokes", 1.50, "wave_data_pmstokes_H1.500_L50.710_A-30.00_P120.00.csv"),
    ("PM-Stokes", 4.00, "wave_data_pmstokes_H4.000_L112.766_A30.00_P30.00.csv"),
    ("PM-Stokes", 8.50, "wave_data_pmstokes_H8.500_L202.839_A-30.00_P72.00.csv"),
)

SEGMENTS = (
    ("start", DURATION_SEC - WINDOW_SEC, TRANSITION_START_SEC),
    ("blend", TRANSITION_START_SEC, TRANSITION_END_SEC),
    ("recover", TRANSITION_END_SEC, TRANSITION_RECOVER_END_SEC),
    ("end", TRANSITION_RECOVER_END_SEC, DURATION_SEC),
)

AXES = {
    "moment": {
        "kind": "compile",
        "values": ("3", "4", "6", "8", "10", "12"),
        "baseline": "8",
    },
    "log": {
        "kind": "compile",
        "values": ("0", "0.05", "0.10", "0.25", "0.50", """1.0""", "2.0"),
        "baseline": "0.50",
    },
    "sigma_k": {
        "kind": "runtime",
        "env": "OU_SIGMA_VAR_K_PERIODS",
        "values": ("0.5", "1", "2", "3", "4", "6", "8"),
        "baseline": "2",
    },
    "adapt_tau": {
        "kind": "runtime",
        "env": "OU_III_ADAPT_TAU_SEA_PERIODS",
        "values": ("0.10", "0.20", "0.30", "0.40", "0.60", "0.80"),
        "baseline": "0.40",
    },
    "rs_mult": {
        "kind": "runtime",
        "env": "OU_III_ADAPT_RS_MULT",
        "values": ("0.50", "0.75", "1.0", "1.5", "2.0", "3.0"),
        "baseline": "1.5",
    },
}

FIELDS = (
    "disp_z_rms_m",
    "disp_3d_rms_m",
    "roll_rms_deg",
    "pitch_rms_deg",
    "yaw_rms_deg",
    "wave_period_s",
    "period_s",
    "tau_applied_s",
    "sigma_applied_mps2",
    "accel_variance_m2ps4",
    "seg_start_disp_z_rms_m",
    "seg_blend_disp_z_rms_m",
    "seg_recover_disp_z_rms_m",
    "seg_end_disp_z_rms_m",
)


def transition_record(wave_seed: int, data_dir: Path, tmpdir: Path) -> Path:
    start_path = ouv.find_default_input(data_dir, "1.500", "50.710")
    end_path = ouv.find_default_input(data_dir, "8.500", "202.839")
    columns, start_data = ouv.read_wave_csv(start_path, DURATION_SEC)
    _, end_data = ouv.read_wave_csv(end_path, DURATION_SEC)
    generated = ouv.make_nonstationary_wave(
        columns,
        start_data,
        end_data,
        wave_seed,
        TRANSITION_END_HEIGHT_M / 8.5,
        TRANSITION_START_SEC,
        TRANSITION_END_SEC,
    )
    path = (
        tmpdir
        / f"wave{wave_seed}"
        / "wave_data_jonswap_H4.000_L202.839_A-30.00_P120.00.csv"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    ouv.write_wave_csv(path, columns, generated)
    return path


def cpp_float_literal(value: str) -> str:
    """Format a finite sweep value as a portable C++ float literal."""
    number = float(value)
    if not math.isfinite(number):
        raise SystemExit(f"non-finite sweep value: {value}")
    literal = f"{number:.8g}"
    if "." not in literal and "e" not in literal.lower():
        literal += ".0"
    return literal + "f"


def patch_default(axis: str, value: str, original_wave: str) -> None:
    path = ROOT / "src" / "tuner" / "WavePeriodEstimator.h"
    text = original_wave
    literal = cpp_float_literal(value)
    if axis == "moment":
        text, count = re.subn(
            r"float moment_horizon_periods = [0-9.]+f,",
            f"float moment_horizon_periods = {literal},",
            text,
            count=1,
        )
    elif axis == "log":
        text, count = re.subn(
            r"float log_smoothing_periods = [0-9.]+f,",
            f"float log_smoothing_periods = {literal},",
            text,
            count=1,
        )
    else:
        raise SystemExit(f"not a compile-time axis: {axis}")
    if count != 1:
        raise SystemExit(f"failed to patch {axis} default")
    path.write_text(text, encoding="utf-8")


def build_binary(out: Path | None = None) -> Path:
    work = ROOT / "tests" / "kalman_ou_iii"
    subprocess.run(
        ["make", "clean"], cwd=work, check=False, stdout=subprocess.DEVNULL
    )
    subprocess.run(["make", "build"], cwd=work, check=True)
    src = work / "kalman_ou_iii-sim"
    if not src.exists():
        raise SystemExit(f"missing {src}")
    if out is None or out.resolve() == src.resolve():
        src.chmod(0o755)
        return src
    shutil.copy2(src, out)
    out.chmod(0o755)
    return out


def run_one(
    binary: Path,
    record: Path,
    env_name: str | None,
    value: str,
    seed: int | None,
    segments: bool,
) -> dict[str, float]:
    env = dict(os.environ)
    env["W3D_WRITE_TIMESERIES"] = "0"
    env["W3D_COLLECT_ALL_GATES"] = "1"
    env["W3D_VALIDATION_WINDOW_SEC"] = str(WINDOW_SEC)
    if env_name:
        env[env_name] = value
    if segments:
        env["W3D_VALIDATION_SEGMENTS"] = ",".join(
            f"{name}:{lo:.9g}:{hi:.9g}" for name, lo, hi in SEGMENTS
        )
    if seed is not None:
        env["W3D_IMU_SEED"] = str(seed)
        env["W3D_INIT_SEED"] = str(1000 + seed)

    completed = subprocess.run(
        [str(binary), "--input", str(record.resolve())],
        cwd=ROOT / "tests" / "kalman_ou_iii",
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    try:
        parsed = ouv.parse_validation_metrics(completed.stdout)
    except ValueError:
        sys.stderr.write(completed.stdout[-4000:] + completed.stderr[-4000:])
        raise SystemExit(
            f"no VALIDATION_METRICS for {record.name} value={value} seed={seed}"
        )
    return {field: float(parsed.get(field, float("nan"))) for field in FIELDS}


def paired_logs(raw, metric, value, baseline, keys, seeds):
    logs = []
    for key in keys:
        for seed in seeds:
            num = raw[(key, value, seed)][metric]
            den = raw[(key, baseline, seed)][metric]
            if num > 0.0 and den > 0.0 and math.isfinite(num) and math.isfinite(den):
                logs.append(math.log(num / den))
    return logs


def ratio_summary(logs):
    if not logs:
        return float("nan"), float("nan"), float("nan")
    mean = sum(logs) / len(logs)
    if len(logs) < 2:
        return math.exp(mean), float("nan"), float("nan")
    var = sum((x - mean) ** 2 for x in logs) / (len(logs) - 1)
    half = 1.96 * math.sqrt(var / len(logs))
    return math.exp(mean), math.exp(mean - half), math.exp(mean + half)


def fmt_ratio(summary):
    ratio, lo, hi = summary
    if not math.isfinite(ratio):
        return "n/a"
    if math.isfinite(lo) and math.isfinite(hi):
        return (
            f"{100 * (ratio - 1):+.3f}% "
            f"[{100 * (lo - 1):+.3f},{100 * (hi - 1):+.3f}]"
        )
    return f"{100 * (ratio - 1):+.3f}%"


def finite_mean(values):
    vals = [x for x in values if math.isfinite(x)]
    return sum(vals) / len(vals) if vals else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--axis", choices=sorted(AXES), required=True)
    ap.add_argument("--values")
    ap.add_argument("--baseline")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--jobs", type=int, default=max(1, os.cpu_count() or 2))
    ap.add_argument(
        "--data-dir", type=Path, default=ROOT / "plots" / "kalman_ou_iii"
    )
    ap.add_argument("--csv", type=Path)
    ap.add_argument("--report", type=Path)
    args = ap.parse_args()

    spec = AXES[args.axis]
    values = (
        tuple(v.strip() for v in args.values.split(","))
        if args.values
        else spec["values"]
    )
    baseline = args.baseline or spec["baseline"]
    if baseline not in values:
        raise SystemExit(f"baseline {baseline} not in values {values}")
    seeds = [None] if args.seeds <= 1 else list(range(1, args.seeds + 1))

    tmp_holder = tempfile.TemporaryDirectory(prefix="ou_period_retune_")
    tmp = Path(tmp_holder.name)
    transition = transition_record(11, args.data_dir, tmp)
    scenarios = [
        (args.data_dir / filename, False) for _, _, filename in STATIONARY
    ]
    scenarios.append((transition, True))

    wave_path = ROOT / "src" / "tuner" / "WavePeriodEstimator.h"
    original_wave = wave_path.read_text(encoding="utf-8")
    binaries: dict[str, Path] = {}
    try:
        if spec["kind"] == "compile":
            for value in values:
                patch_default(args.axis, value, original_wave)
                binaries[value] = build_binary(
                    tmp / f"sim-{args.axis}-{value.replace('.', '_')}"
                )
                wave_path.write_text(original_wave, encoding="utf-8")
        else:
            wave_path.write_text(original_wave, encoding="utf-8")
            binary = build_binary()
            binaries = {value: binary for value in values}
    finally:
        wave_path.write_text(original_wave, encoding="utf-8")

    raw = {}
    pending = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for record, is_transition in scenarios:
            key = "transition" if is_transition else record.name
            for value in values:
                for seed in seeds:
                    future = pool.submit(
                        run_one,
                        binaries[value],
                        record,
                        spec.get("env"),
                        value,
                        seed,
                        is_transition,
                    )
                    pending.append((key, value, seed, future))
        for key, value, seed, future in pending:
            raw[(key, value, seed)] = future.result()

    stationary_keys = [filename for _, _, filename in STATIONARY]
    transition_keys = ["transition"]
    report_lines = [
        f"# OU-III period/statistics retune: {args.axis}",
        "",
        f"seeds={args.seeds}; baseline={baseline}; values={', '.join(values)}",
        "",
        "Paired geometric-mean ratios vs baseline; negative is better.",
        "",
        "| value | stationary Z | stationary 3D | blend Z | recover Z | mean Tz [s] | mean sigma [m/s2] |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]

    for value in values:
        st_z = ratio_summary(
            paired_logs(raw, "disp_z_rms_m", value, baseline, stationary_keys, seeds)
        )
        st_3d = ratio_summary(
            paired_logs(raw, "disp_3d_rms_m", value, baseline, stationary_keys, seeds)
        )
        tr_blend = ratio_summary(
            paired_logs(
                raw,
                "seg_blend_disp_z_rms_m",
                value,
                baseline,
                transition_keys,
                seeds,
            )
        )
        tr_recover = ratio_summary(
            paired_logs(
                raw,
                "seg_recover_disp_z_rms_m",
                value,
                baseline,
                transition_keys,
                seeds,
            )
        )
        mean_period = finite_mean(
            raw[(key, value, seed)]["wave_period_s"]
            for key in stationary_keys
            for seed in seeds
        )
        mean_sigma = finite_mean(
            raw[(key, value, seed)]["sigma_applied_mps2"]
            for key in stationary_keys
            for seed in seeds
        )
        report_lines.append(
            f"| {value} | {fmt_ratio(st_z)} | {fmt_ratio(st_3d)} | "
            f"{fmt_ratio(tr_blend)} | {fmt_ratio(tr_recover)} | "
            f"{mean_period:.4f} | {mean_sigma:.4f} |"
        )

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=["scenario", "value", "seed", *FIELDS],
                lineterminator="\n",
            )
            writer.writeheader()
            for (scenario, value, seed), metrics in sorted(
                raw.items(), key=lambda item: (item[0][0], item[0][1], str(item[0][2]))
            ):
                writer.writerow(
                    {
                        "scenario": scenario,
                        "value": value,
                        "seed": "default" if seed is None else seed,
                        **metrics,
                    }
                )

    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print("\n".join(report_lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
