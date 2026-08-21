#!/usr/bin/env python3
"""Joint paired confirmation for OU-III period/statistics retuning.

Runs the production baseline and the two joint candidate tuples selected from the
five-seed one-factor fine sweeps on identical stationary records and IMU/init
seeds. Both compile-time WavePeriodEstimator horizons are changed together;
the acceleration-variance horizon is applied in the same simulator run through
its existing runtime override. Production source defaults are restored after
each disposable build and are never committed by this tool.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import ou_period_statistics_retune as base  # noqa: E402
import ou_validation as ouv  # noqa: E402

BASELINE = "baseline"
CANDIDATES = {
    BASELINE: {
        "moment": "8",
        "log": "0.50",
        "sigma_k": "2",
    },
    "joint_log_005": {
        "moment": "4",
        "log": "0.05",
        "sigma_k": "4",
    },
    "joint_log_010": {
        "moment": "4",
        "log": "0.10",
        "sigma_k": "4",
    },
}


def patched_wave_defaults(original: str, moment: str, log_horizon: str) -> str:
    """Return WavePeriodEstimator text with both compile-time horizons patched."""
    text, moment_count = re.subn(
        r"float moment_horizon_periods = [0-9.]+f,",
        f"float moment_horizon_periods = {base.cpp_float_literal(moment)},",
        original,
        count=1,
    )
    text, log_count = re.subn(
        r"float log_smoothing_periods = [0-9.]+f,",
        f"float log_smoothing_periods = {base.cpp_float_literal(log_horizon)},",
        text,
        count=1,
    )
    if moment_count != 1 or log_count != 1:
        raise SystemExit(
            "failed to patch joint WavePeriodEstimator defaults: "
            f"moment={moment_count}, log={log_count}"
        )
    return text


def run_one(
    binary: Path,
    record: Path,
    sigma_k: str,
    seed: int,
    segments: bool,
) -> dict[str, float]:
    env = dict(os.environ)
    env["W3D_WRITE_TIMESERIES"] = "0"
    env["W3D_COLLECT_ALL_GATES"] = "1"
    env["W3D_VALIDATION_WINDOW_SEC"] = str(base.WINDOW_SEC)
    env["OU_SIGMA_VAR_K_PERIODS"] = sigma_k
    if segments:
        env["W3D_VALIDATION_SEGMENTS"] = ",".join(
            f"{name}:{lo:.9g}:{hi:.9g}" for name, lo, hi in base.SEGMENTS
        )
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
            f"no VALIDATION_METRICS for {record.name} "
            f"sigma_k={sigma_k} seed={seed}"
        )
    return {
        field: float(parsed.get(field, float("nan"))) for field in base.FIELDS
    }


def finite_mean(values) -> float:
    vals = [x for x in values if math.isfinite(x)]
    return sum(vals) / len(vals) if vals else float("nan")


def candidate_label(name: str) -> str:
    cfg = CANDIDATES[name]
    return (
        f"{name} "
        f"(K_T={cfg['moment']}, K_log={cfg['log']}, K_sigma={cfg['sigma_k']})"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--jobs", type=int, default=max(1, os.cpu_count() or 2))
    ap.add_argument(
        "--data-dir", type=Path, default=ROOT / "plots" / "kalman_ou_iii"
    )
    ap.add_argument("--csv", type=Path)
    ap.add_argument("--report", type=Path)
    args = ap.parse_args()

    if args.seeds < 2:
        raise SystemExit("joint confirmation requires at least two seeds")
    seeds = list(range(1, args.seeds + 1))

    tmp_holder = tempfile.TemporaryDirectory(prefix="ou_period_joint_")
    tmp = Path(tmp_holder.name)
    transition = base.transition_record(11, args.data_dir, tmp)
    scenarios = [
        (args.data_dir / filename, False)
        for _, _, filename in base.STATIONARY
    ]
    scenarios.append((transition, True))

    wave_path = ROOT / "src" / "tuner" / "WavePeriodEstimator.h"
    original_wave = wave_path.read_text(encoding="utf-8")
    binaries: dict[str, Path] = {}
    try:
        for name, cfg in CANDIDATES.items():
            wave_path.write_text(
                patched_wave_defaults(original_wave, cfg["moment"], cfg["log"]),
                encoding="utf-8",
            )
            binaries[name] = base.build_binary(tmp / f"sim-joint-{name}")
            wave_path.write_text(original_wave, encoding="utf-8")
    finally:
        wave_path.write_text(original_wave, encoding="utf-8")

    raw: dict[tuple[str, str, int], dict[str, float]] = {}
    pending = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for record, is_transition in scenarios:
            scenario = "transition" if is_transition else record.name
            for name, cfg in CANDIDATES.items():
                for seed in seeds:
                    future = pool.submit(
                        run_one,
                        binaries[name],
                        record,
                        cfg["sigma_k"],
                        seed,
                        is_transition,
                    )
                    pending.append((scenario, name, seed, future))
        for scenario, name, seed, future in pending:
            raw[(scenario, name, seed)] = future.result()

    stationary_keys = [filename for _, _, filename in base.STATIONARY]
    transition_keys = ["transition"]

    report_lines = [
        "# OU-III joint period/statistics confirmation",
        "",
        f"seeds={args.seeds}",
        f"baseline={candidate_label(BASELINE)}",
        "candidates=" + "; ".join(
            candidate_label(name) for name in CANDIDATES if name != BASELINE
        ),
        "",
        "Paired geometric-mean ratios vs production baseline; negative is better.",
        "",
        "| candidate | stationary Z | stationary 3D | blend Z | recover Z | mean Tz [s] | mean sigma [m/s2] |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]

    for name in CANDIDATES:
        st_z = base.ratio_summary(
            base.paired_logs(
                raw, "disp_z_rms_m", name, BASELINE, stationary_keys, seeds
            )
        )
        st_3d = base.ratio_summary(
            base.paired_logs(
                raw, "disp_3d_rms_m", name, BASELINE, stationary_keys, seeds
            )
        )
        tr_blend = base.ratio_summary(
            base.paired_logs(
                raw,
                "seg_blend_disp_z_rms_m",
                name,
                BASELINE,
                transition_keys,
                seeds,
            )
        )
        tr_recover = base.ratio_summary(
            base.paired_logs(
                raw,
                "seg_recover_disp_z_rms_m",
                name,
                BASELINE,
                transition_keys,
                seeds,
            )
        )
        mean_period = finite_mean(
            raw[(key, name, seed)]["wave_period_s"]
            for key in stationary_keys
            for seed in seeds
        )
        mean_sigma = finite_mean(
            raw[(key, name, seed)]["sigma_applied_mps2"]
            for key in stationary_keys
            for seed in seeds
        )
        report_lines.append(
            f"| {candidate_label(name)} | {base.fmt_ratio(st_z)} | "
            f"{base.fmt_ratio(st_3d)} | {base.fmt_ratio(tr_blend)} | "
            f"{base.fmt_ratio(tr_recover)} | {mean_period:.4f} | "
            f"{mean_sigma:.4f} |"
        )

    # Direct paired comparison resolves the remaining 0.05 vs 0.10 choice.
    report_lines.extend(
        [
            "",
            "Direct paired ratio: joint_log_010 versus joint_log_005; negative favors 0.10.",
            "",
            "| comparison | stationary Z | stationary 3D | blend Z | recover Z |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    direct = {}
    for metric, keys in (
        ("disp_z_rms_m", stationary_keys),
        ("disp_3d_rms_m", stationary_keys),
        ("seg_blend_disp_z_rms_m", transition_keys),
        ("seg_recover_disp_z_rms_m", transition_keys),
    ):
        direct[metric] = base.ratio_summary(
            base.paired_logs(
                raw, metric, "joint_log_010", "joint_log_005", keys, seeds
            )
        )
    report_lines.append(
        "| 0.10 / 0.05 | "
        f"{base.fmt_ratio(direct['disp_z_rms_m'])} | "
        f"{base.fmt_ratio(direct['disp_3d_rms_m'])} | "
        f"{base.fmt_ratio(direct['seg_blend_disp_z_rms_m'])} | "
        f"{base.fmt_ratio(direct['seg_recover_disp_z_rms_m'])} |"
    )

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=[
                    "scenario",
                    "candidate",
                    "moment_horizon_periods",
                    "log_smoothing_periods",
                    "sigma_var_k_periods",
                    "seed",
                    *base.FIELDS,
                ],
                lineterminator="\n",
            )
            writer.writeheader()
            for (scenario, name, seed), metrics in sorted(
                raw.items(),
                key=lambda item: (item[0][0], item[0][1], item[0][2]),
            ):
                cfg = CANDIDATES[name]
                writer.writerow(
                    {
                        "scenario": scenario,
                        "candidate": name,
                        "moment_horizon_periods": cfg["moment"],
                        "log_smoothing_periods": cfg["log"],
                        "sigma_var_k_periods": cfg["sigma_k"],
                        "seed": seed,
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
