#!/usr/bin/env python3
"""Paired retuning study for the OU-III sea-time/statistics pipeline.

The refactored path has three statistical horizons before the existing parameter
slew dynamics:

* WavePeriodEstimator moment horizon K_T * T_z (compile-time constructor default)
* WavePeriodEstimator canonical log(T_z) EMA K_log * T_z (compile-time default)
* SeaStateAutoTuner acceleration-moment horizon K_sigma * T_z (runtime setter)

Two downstream scheduler horizons are rechecked after the statistical defaults
are selected:

* common tau/sigma slew c_adapt * T_sea
* r_S slew m_RS * tau_target

For every candidate the harness replays identical records and seeds, reports
paired log-ratios against the axis baseline, and includes the controlled
1.5 -> 4.0 m transition with separate blend/recovery segments.  Lower error is
better.  Compile-time axes are implemented by patching constructor defaults in
the disposable Actions checkout, compiling a private simulator binary, then
restoring the source before the next candidate; production APIs remain clean.
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
        "values": ("0", "0.05", "0.10", "0.25", "0.50", "1.0", "2.0"),
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
    path = tmpdir / f"wave{wave_seed}" / "wave_data_jonswap_H4.000_L202.839_A-30.00_P120.00.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    ouv.write_wave_csv(path, columns, generated)
    return path


def patch_default(axis: str, value: str, original_wave: str, original_tuner: str) -> None:
    wave = ROOT / "src" / "tuner" / "WavePeriodEstimator.h"
    tuner = ROOT / "src" / "tuner" / "SeaStateAutoTuner.h"
    wave_text = original_wave
    tuner_text = original_tuner
    if axis == "moment":
        wave_text, n = re.subn(
            r"float moment_horizon_periods = [0-9.]+f,",
            f"float moment_horizon_periods = {float(value):.8g}f,",
            wave_text,
            count=1,
        )
        if n != 1:
            raise SystemExit("could not patch moment_horizon_periods default")
    elif axis == "log":
        wave_text, n = re.subn(
            r"float log_smoothing_periods = [0-9.]+f,",
            f"float log_smoothing_periods = {float(value):.8g}f,",
            wave_text,
            count=1,
        )
        if n != 1:
            raise SystemExit("could not patch log_smoothing_periods default")
    else:
        raise SystemExit(f"axis {axis} is not compile-time")
    wave.write_text(wave_text, encoding="utf-8")
    tuner.write_text(tuner_text, encoding="utf-8")


def restore_sources(original_wave: str, original_tuner: str) -> None:
    (ROOT / "src" / "tuner" / "WavePeriodEstimator.h").write_text(original_wave, encoding="utf-8")
    (ROOT / "src" / "tuner" / "SeaStateAutoTuner.h").write_text(original_tuner, encoding="utf-8")


def build_binary(out: Path) -> None:
    work = ROOT / "tests" / "kalman_ou_iii"
    subprocess.run(["make", "clean"], cwd=work, check=False, stdout=subprocess.DEVNULL)
    subprocess.run(["make", "build"], cwd=work, check=True)
    src = work / "kalman_ou_iii-sim"
    if not src.exists():
        raise SystemExit(f"build did not produce {src}")
    shutil.copy2(src, out)
    out.chmod(0o755)


def run_one(binary: Path, record: Path, env_name: str | None, value: str,
            seed: int | None, segments: bool) -> dict[str, float]:
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
        raise SystemExit(f"no VALIDATION_METRICS for {record.name} value={value} seed={seed}")
    return {field: float(parsed.get(field, float("nan"))) for field in FIELDS}


def ratio_summary(logs: list[float]) -> tuple[float, float, float]:
    if not logs:
        return float("nan"), float("nan"), float("nan")
    mean = sum(logs) / len(logs)
    if len(logs) < 2:
        return math.exp(mean), float("nan"), float("nan")
    var = sum((x - mean) ** 2 for x in logs) / (len(logs) - 1)
    half = 1.96 * math.sqrt(var / len(logs))
    return math.exp(mean), math.exp(mean - half), math.exp(mean + half)


def paired_logs(raw, metric, value, baseline, scenario_keys, seeds):
    out = []
    for key in scenario_keys:
        for seed in seeds:
            num = raw[(key, value, seed)][metric]
            den = raw[(key, baseline, seed)][metric]
            if num > 0.0 and den > 0.0 and math.isfinite(num) and math.isfinite(den):
                out.append(math.log(num / den))
    return out


def fmt_ratio(summary):
    r, lo, hi = summary
    if not math.isfinite(r):
        return "n/a"
    if math.isfinite(lo) and math.isfinite(hi):
        return f"{100*(r-1):+.3f}% [{100*(lo-1):+.3f},{100*(hi-1):+.3f}]"
    return f"{100*(r-1):+.3f}%"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--axis", choices=sorted(AXES), required=True)
    ap.add_argument("--values")
    ap.add_argument("--baseline")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--jobs", type=int, default=max(1, os.cpu_count() or 2))
    ap.add_argument("--data-dir", type=Path, default=ROOT / "tests" / "kalman_ou_iii")
    ap.add_argument("--csv", type=Path)
    ap.add_argument("--report", type=Path)
    args = ap.parse_args()

    spec = AXES[args.axis]
    values = tuple(v.strip() for v in args.values.split(",")) if args.values else spec["values"]
    baseline = args.baseline or spec["baseline"]
    if baseline not in values:
        raise SystemExit(f"baseline {baseline} not in values {values}")
    seeds = [None] if args.seeds <= 1 else list(range(1, args.seeds + 1))

    tmp_holder = tempfile.TemporaryDirectory(prefix="ou_period_retune_")
    tmp = Path(tmp_holder.name)
    transition = transition_record(11, args.data_dir, tmp)
    scenarios = [(record, False) for _, _, record in STATIONARY]
    scenarios.append((str(transition), True))

    wave_path = ROOT / "src" / "tuner" / "WavePeriodEstimator.h"
    tuner_path = ROOT / "src" / "tuner" / "SeaStateAutoTuner.h"
    original_wave = wave_path.read_text(encoding="utf-8")
    original_tuner = tuner_path.read_text(encoding="utf-8")

    binaries: dict[str, Path] = {}
    default_binary = ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim"
    try:
        if spec["kind"] == "compile":
            for value in values:
                patch_default(args.axis, value, original_wave, original_tuner)
                out = tmp / f"kalman_ou_iii-sim-{args.axis}-{value.replace('.', '_')}"
                build_binary(out)
                binaries[value] = out
                restore_sources(original_wave, original_tuner)
        else:
            restore_sources(original_wave, original_tuner)
            build_binary(default_binary)
            binaries = {value: default_binary for value in values}
    finally:
        restore_sources(original_wave, original_tuner)

    raw: dict[tuple[str, str, int | None], dict[str, float]] = {}
    jobs = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for scenario, segments in scenarios:
            record = Path(scenario)
            key = "transition" if segments else record.name
            for value in values:
                for seed in seeds:
                    fut = pool.submit(
                        run_one,
                        binaries[value],
                        record,
                        spec.get("env"),
                        value,
                        seed,
                        segments,
                    )
                    jobs.append((key, value, seed, fut))
        for key, value, seed, fut in jobs:
            raw[(key, value, seed)] = fut.result()

    stationary_keys = [name for _, _, name in STATIONARY]
    transition_keys = ["transition"]
    rows = []
    report_lines = [
        f"# OU-III period/statistics retune: {args.axis}",
        "",
        f"seeds={args.seeds}; baseline={baseline}; values={', '.join(values)}",
        "",
        "Ratios are paired geometric means against the baseline; negative is better.",
        "",
        "| value | stationary Z | stationary 3D | transition blend Z | transition recover Z | settled period | sigma |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]

    for value in values:
        st_z = ratio_summary(paired_logs(raw, "disp_z_rms_m", value, baseline, stationary_keys, seeds))
        st_3d = ratio_summary(paired_logs(raw, "disp_3d_rms_m", value, baseline, stationary_keys, seeds))
        tr_blend = ratio_summary(paired_logs(raw, "seg_blend_disp_z_rms_m", value, baseline, transition_keys, seeds))
        tr_recover = ratio_summary(paired_logs(raw, "seg_recover_disp_z_rms_m", value, baseline, transition_keys, seeds))

        period_vals = [raw[(key, value, seed)]["wave_period_s"] for key in stationary_keys for seed in seeds]
        sigma_vals = [raw[(key, value, seed)]["sigma_applied_mps2"] for key in stationary_keys for seed in seeds]
        mean_period = sum(x for x in period_vals if math.isfinite(x)) / max(1, sum(math.isfinite(x) for x in period_vals))
        mean_sigma = sum(x for x in sigma_vals if math.isfinite(x)) / max(1, sum(math.isfinite(x) for x in sigma_vals))
        rows.append({
            "axis": args.axis,
            "value": value,
            "baseline": baseline,
            "stationary_z_ratio": st_z[0],
            "stationary_z_lo": st_z[1],
            "stationary_z_hi": st_z[2],
            "stationary_3d_ratio": st_3d[0],
            "transition_blend_z_ratio": tr_blend[0],
            "transition_recover_z_ratio": tr_recover[0],
            "mean_wave_period_s": mean_period,
            "mean_sigma_applied_mps2": mean_sigma,
        })
        report_lines.append(
            f"| {value} | {fmt_ratio(st_z)} | {fmt_ratio(st_3d)} | {fmt_ratio(tr_blend)} | {fmt_ratio(tr_recover)} | {mean_period:.4f} | {mean_sigma:.4f} |"
        )

    # Raw per-run appendix makes the study independently auditable.
    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=["scenario", "value", "seed", *FIELDS], lineterminator="\n")
            writer.writeheader()
            for (scenario, value, seed), metrics in sorted(raw.items(), key=lambda x: (x[0][0], x[0][1], str(x[0][2]))):
                writer.writerow({"scenario": scenario, "value": value, "seed": "default" if seed is None else seed, **metrics})

    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print("\n".join(report_lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
