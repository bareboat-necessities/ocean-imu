#!/usr/bin/env python3
"""Full paired validation of candidate OU-III effective r_S adaptation laws.

This is an experiment runner, not a production filter change.  It patches an
inactive compile-time branch into SeaStateFusionFilter_OU_III.h inside the CI
working tree, compiles the simulator three times, and reuses the established
ou_validation.py scenario/seed machinery in shard mode so no publication bundle
is rewritten.

Arms:
  deployed        current effective law (about sigma_aw * tau^2.5)
  cubic_p3        r_S,eff = 0.1548363522 * sigma_aw * tau^3
  fitted_p2p9052  r_S,eff = 0.1667018769 * sigma_aw * tau^2.9052

Both candidate normalizations pass through the independently measured full-3D
optimum r_S,eff=1.160579 m*s at the nominal JONSWAP point
(tau=2.17904091 s, sigma_aw=0.724445343 m/s^2).
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
HEADER = REPO_ROOT / "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
MAKE_DIR = REPO_ROOT / "tests/kalman_ou_iii"
VALIDATION = REPO_ROOT / "tools/ou_validation.py"

SOURCE_FILES = (
    "wave_data_jonswap_H0.270_L14.047_A30.00_P60.00.csv",
    "wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv",
    "wave_data_jonswap_H4.000_L112.766_A30.00_P30.00.csv",
    "wave_data_jonswap_H8.500_L202.839_A-30.00_P72.00.csv",
    "wave_data_pmstokes_H0.270_L14.047_A30.00_P60.00.csv",
    "wave_data_pmstokes_H1.500_L50.710_A-30.00_P120.00.csv",
    "wave_data_pmstokes_H4.000_L112.766_A30.00_P30.00.csv",
    "wave_data_pmstokes_H8.500_L202.839_A-30.00_P72.00.csv",
)
NOMINAL_FILE = "wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv"
TRANSITION_END_SOURCE = "wave_data_jonswap_H8.500_L202.839_A-30.00_P72.00.csv"

ARMS: tuple[tuple[str, float | None, float | None], ...] = (
    ("deployed", None, None),
    ("cubic_p3", 0.1548363522, 3.0),
    ("fitted_p2p9052", 0.1667018769, 2.9052),
)

PATCH_OLD = """        if (rs_law_ == RSAdaptationLaw::Cubic) {\n            return R_S_coeff_ * sigma * tau3;\n        }\n"""
PATCH_NEW = """        if (rs_law_ == RSAdaptationLaw::Cubic) {\n#ifdef OU_III_EXPERIMENT_EFFECTIVE_RS_POWER\n            // Experiment-only branch.  The public/deployed source path remains\n            // unchanged when the macro is absent.  The wrapper stores a base\n            // r_S target and later multiplies it by sqrt(T0/TS).  Therefore a\n            // desired effective-input law K*sigma*tau^p has to be divided by\n            // that information-rate scale here.\n            const float TS_exp = pseudo_update_period_for_(tau);\n            if (std::isfinite(TS_exp) && TS_exp > 0.0f) {\n                const float info_scale_exp =\n                    std::sqrt(pseudo_update_fixed_period_s_ / TS_exp);\n                return OU_III_EXPERIMENT_EFFECTIVE_RS_COEFF * sigma *\n                    std::pow(tau, OU_III_EXPERIMENT_EFFECTIVE_RS_TAU_EXP) /\n                    info_scale_exp;\n            }\n#endif\n            return R_S_coeff_ * sigma * tau3;\n        }\n"""

PAIR_KEYS = (
    "scenario",
    "wave_phase_seed",
    "imu_noise_seed",
    "initialization_seed",
)
SEED_KEYS = (
    "wave_phase_seed",
    "imu_noise_seed",
    "initialization_seed",
)


def run(command: Sequence[str], *, env: Mapping[str, str] | None = None) -> None:
    printable = " ".join(command)
    print(f"+ {printable}", flush=True)
    subprocess.run(
        list(command),
        cwd=REPO_ROOT,
        env=None if env is None else dict(env),
        check=True,
    )


def patch_experiment_branch() -> str:
    original = HEADER.read_text(encoding="utf-8")
    count = original.count(PATCH_OLD)
    if count != 1:
        raise RuntimeError(
            f"expected exactly one Cubic-law patch point, found {count}"
        )
    HEADER.write_text(original.replace(PATCH_OLD, PATCH_NEW), encoding="utf-8")
    return original


def build_arm(coeff: float | None, exponent: float | None) -> None:
    env = os.environ.copy()
    if coeff is None:
        env.pop("CPPFLAGS", None)
    else:
        env["CPPFLAGS"] = " ".join(
            (
                "-DOU_III_EXPERIMENT_EFFECTIVE_RS_POWER=1",
                f"-DOU_III_EXPERIMENT_EFFECTIVE_RS_COEFF={coeff:.10g}f",
                f"-DOU_III_EXPERIMENT_EFFECTIVE_RS_TAU_EXP={exponent:.10g}f",
            )
        )
    run(
        ("make", "-C", str(MAKE_DIR.relative_to(REPO_ROOT)), "-B", "kalman_ou_iii-sim"),
        env=env,
    )


def validation_command(
    data_dir: Path,
    output_dir: Path,
    shard_dir: Path,
    shard_index: int,
    shard_count: int,
    jobs: int,
) -> list[str]:
    command = [
        sys.executable,
        str(VALIDATION.relative_to(REPO_ROOT)),
        "--mode", "full",
        "--data-dir", str(data_dir),
        "--output-dir", str(output_dir),
        "--families", "OU_III",
        "--adaptation-modes", "Adaptive",
        "--pmstokes-modes", "Adaptive",
        "--duration-sec", "1200",
        "--window-sec", "900",
        "--transition-start-sec", "420",
        "--transition-end-sec", "780",
        "--nonstationary-start", str(data_dir / NOMINAL_FILE),
        "--nonstationary-end", str(data_dir / TRANSITION_END_SOURCE),
        "--nonstationary-end-height", "4.0",
        "--bootstrap-resamples", "100",
        "--stats-seed", "20260317",
        "--jobs", str(jobs),
        "--skip-build",
        "--no-plots",
        "--shard-count", str(shard_count),
        "--shard-index", str(shard_index),
        "--shard-dir", str(shard_dir),
    ]
    for name in SOURCE_FILES:
        command.extend(("--stationary-input", str(data_dir / name)))
    return command


def run_arm(
    label: str,
    coeff: float | None,
    exponent: float | None,
    data_dir: Path,
    output_root: Path,
    jobs: int,
    shard_count: int,
) -> None:
    print(f"\n=== ARM {label} ===", flush=True)
    arm_dir = output_root / label
    shard_dir = arm_dir / "shards"
    shutil.rmtree(arm_dir, ignore_errors=True)
    shard_dir.mkdir(parents=True, exist_ok=True)
    build_arm(coeff, exponent)
    for shard_index in range(shard_count):
        run(
            validation_command(
                data_dir, arm_dir, shard_dir, shard_index, shard_count, jobs
            )
        )


def load_rows(shard_dir: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    paths = sorted(shard_dir.glob("ou_validation_shard_*.json"))
    if not paths:
        raise FileNotFoundError(f"no validation shards in {shard_dir}")
    for path in paths:
        with path.open(encoding="utf-8") as stream:
            rows.extend(json.load(stream))
    rows = [
        row for row in rows
        if row.get("family") == "OU_III" and row.get("mode") == "Adaptive"
    ]
    keys = [tuple(row[key] for key in PAIR_KEYS) for row in rows]
    if len(keys) != len(set(keys)):
        raise RuntimeError(f"duplicate paired rows in {shard_dir}")
    if len(rows) != 90:
        raise RuntimeError(
            f"expected 90 OU-III/Adaptive rows (9 scenarios x 10 seeds), got {len(rows)}"
        )
    return rows


def row_index(rows: Iterable[Mapping[str, Any]]) -> dict[tuple[Any, ...], Mapping[str, Any]]:
    return {tuple(row[key] for key in PAIR_KEYS): row for row in rows}


def scenarios_for(rows: Sequence[Mapping[str, Any]], spectrum: str | None = None) -> list[str]:
    names = sorted(
        {
            str(row["scenario"])
            for row in rows
            if str(row.get("scenario_kind")) == "stationary"
            and (spectrum is None or str(row.get("spectrum")) == spectrum)
        }
    )
    return names


def seed_aggregates(
    rows: Sequence[Mapping[str, Any]],
    scenarios: Sequence[str],
    metric: str,
) -> dict[tuple[Any, ...], float]:
    wanted = set(scenarios)
    grouped: dict[tuple[Any, ...], list[float]] = defaultdict(list)
    seen: dict[tuple[Any, ...], set[str]] = defaultdict(set)
    for row in rows:
        scenario = str(row["scenario"])
        if scenario not in wanted:
            continue
        value = float(row.get(metric, math.nan))
        if not math.isfinite(value):
            continue
        seed = tuple(row[key] for key in SEED_KEYS)
        grouped[seed].append(value)
        seen[seed].add(scenario)
    result: dict[tuple[Any, ...], float] = {}
    for seed, values in grouped.items():
        if seen[seed] == wanted:
            result[seed] = float(np.mean(values))
    return result


def percentile_interval(values: np.ndarray, rng: np.random.Generator, resamples: int) -> tuple[float, float]:
    if values.size == 0:
        return math.nan, math.nan
    if values.size == 1:
        return float(values[0]), float(values[0])
    indices = rng.integers(0, values.size, size=(resamples, values.size))
    means = values[indices].mean(axis=1)
    low, high = np.quantile(means, (0.025, 0.975))
    return float(low), float(high)


def paired_summary(
    baseline_rows: Sequence[Mapping[str, Any]],
    arm_rows: Sequence[Mapping[str, Any]],
    scenarios: Sequence[str],
    metric: str,
    rng: np.random.Generator,
    resamples: int,
) -> dict[str, Any]:
    base = seed_aggregates(baseline_rows, scenarios, metric)
    arm = seed_aggregates(arm_rows, scenarios, metric)
    seeds = sorted(set(base).intersection(arm))
    if not seeds:
        return {"n": 0}
    base_values = np.asarray([base[key] for key in seeds], dtype=np.float64)
    arm_values = np.asarray([arm[key] for key in seeds], dtype=np.float64)
    difference = arm_values - base_values
    relative = 100.0 * difference / base_values
    diff_low, diff_high = percentile_interval(difference, rng, resamples)
    pct_low, pct_high = percentile_interval(relative, rng, resamples)
    return {
        "n": len(seeds),
        "baseline_mean": float(base_values.mean()),
        "arm_mean": float(arm_values.mean()),
        "mean_difference": float(difference.mean()),
        "difference_ci95_low": diff_low,
        "difference_ci95_high": diff_high,
        "mean_relative_change_pct": float(relative.mean()),
        "relative_ci95_low_pct": pct_low,
        "relative_ci95_high_pct": pct_high,
        "improved_seeds": int(np.sum(difference < 0.0)),
        "worsened_seeds": int(np.sum(difference > 0.0)),
        "ties": int(np.sum(difference == 0.0)),
    }


def pairwise_max_abs_difference(
    baseline_rows: Sequence[Mapping[str, Any]],
    arm_rows: Sequence[Mapping[str, Any]],
    metric: str,
) -> float:
    baseline = row_index(baseline_rows)
    arm = row_index(arm_rows)
    values: list[float] = []
    for key in set(baseline).intersection(arm):
        left = float(baseline[key].get(metric, math.nan))
        right = float(arm[key].get(metric, math.nan))
        if math.isfinite(left) and math.isfinite(right):
            values.append(abs(right - left))
    return max(values, default=math.nan)


def gate_summary(rows: Sequence[Mapping[str, Any]]) -> dict[str, int]:
    return {
        "passes": sum(int(row.get("quality_gate_pass", 0)) for row in rows),
        "total": len(rows),
    }


def fmt(value: float, digits: int = 4) -> str:
    return "--" if not math.isfinite(value) else f"{value:.{digits}f}"


def build_report(
    all_rows: Mapping[str, Sequence[Mapping[str, Any]]],
    output_root: Path,
    resamples: int,
) -> dict[str, Any]:
    baseline = all_rows["deployed"]
    js = scenarios_for(baseline, "jonswap")
    pm = scenarios_for(baseline, "pmstokes")
    stationary = sorted(js + pm)
    transition = sorted(
        {
            str(row["scenario"])
            for row in baseline
            if str(row.get("scenario_kind")) == "nonstationary"
        }
    )
    if len(js) != 4 or len(pm) != 4 or len(transition) != 1:
        raise RuntimeError(
            f"unexpected scenario set: {len(js)} JONSWAP, {len(pm)} PM-Stokes, "
            f"{len(transition)} transition"
        )

    groups = {
        "jonswap_stationary": js,
        "pmstokes_stationary": pm,
        "all_stationary": stationary,
        "transition": transition,
    }
    metrics = (
        "disp_x_rms_m",
        "disp_y_rms_m",
        "disp_z_rms_m",
        "disp_3d_rms_m",
        "disp_z_pct_hs",
        "roll_rms_deg",
        "pitch_rms_deg",
        "yaw_rms_deg",
    )
    transition_metrics = (
        "seg_start_disp_z_rms_m",
        "seg_start_disp_3d_rms_m",
        "seg_blend_disp_z_rms_m",
        "seg_blend_disp_3d_rms_m",
        "seg_end_disp_z_rms_m",
        "seg_end_disp_3d_rms_m",
    )

    report: dict[str, Any] = {
        "protocol": {
            "arms": {
                "deployed": "current adaptive law",
                "cubic_p3": "rS_eff = 0.1548363522 sigma_aw tau^3",
                "fitted_p2p9052": "rS_eff = 0.1667018769 sigma_aw tau^2.9052",
            },
            "stationary_scenarios": stationary,
            "transition_scenario": transition[0],
            "seed_triplets": 10,
            "score_window_sec": 900,
            "transition_sec": [420, 780],
            "bootstrap_resamples": resamples,
            "pairing": "identical wave, IMU-noise, and initialization seeds",
        },
        "arms": {},
    }

    for arm_index, (label, _, _) in enumerate(ARMS):
        if label == "deployed":
            continue
        arm_rows = all_rows[label]
        rng = np.random.default_rng(20260816 + arm_index)
        arm_report: dict[str, Any] = {
            "gates": gate_summary(arm_rows),
            "groups": {},
            "scenarios": {},
            "transition_segments": {},
            "tuner_invariance": {
                "max_abs_tau_applied_s": pairwise_max_abs_difference(
                    baseline, arm_rows, "tau_applied_s"
                ),
                "max_abs_sigma_applied_mps2": pairwise_max_abs_difference(
                    baseline, arm_rows, "sigma_applied_mps2"
                ),
            },
        }
        for group_name, scenario_names in groups.items():
            arm_report["groups"][group_name] = {
                metric: paired_summary(
                    baseline, arm_rows, scenario_names, metric, rng, resamples
                )
                for metric in metrics
            }
        for scenario in stationary + transition:
            arm_report["scenarios"][scenario] = {
                metric: paired_summary(
                    baseline, arm_rows, [scenario], metric, rng, resamples
                )
                for metric in ("disp_z_pct_hs", "disp_3d_rms_m")
            }
        for metric in transition_metrics:
            arm_report["transition_segments"][metric] = paired_summary(
                baseline, arm_rows, transition, metric, rng, resamples
            )
        report["arms"][label] = arm_report

    report["baseline_gates"] = gate_summary(baseline)

    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "summary.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    lines = [
        "# Full-suite effective-r_S candidate validation",
        "",
        "Paired full OU-III Adaptive validation: four stationary JONSWAP seas, four stationary PM-Stokes seas, and the standard 1.5 m -> 4.0 m non-stationary transition; ten predeclared wave/IMU/initialization seed triplets; 1200 s records with the trailing 900 s scored.",
        "",
        "Candidate laws are applied to the effective r_S standard deviation that reaches the MEKF; the existing tau-scaled pseudo-measurement cadence is unchanged.",
        "",
        "| arm | all-8 Z [%Hs] deployed | arm | paired delta [pp] | all-8 3D [m] deployed | arm | paired change [%] | improved 3D seeds |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for label, _, _ in ARMS:
        if label == "deployed":
            continue
        group = report["arms"][label]["groups"]["all_stationary"]
        z = group["disp_z_pct_hs"]
        d3 = group["disp_3d_rms_m"]
        lines.append(
            f"| {label} | {fmt(z['baseline_mean'],3)} | {fmt(z['arm_mean'],3)} | "
            f"{fmt(z['mean_difference'],3)} [{fmt(z['difference_ci95_low'],3)}, {fmt(z['difference_ci95_high'],3)}] | "
            f"{fmt(d3['baseline_mean'],4)} | {fmt(d3['arm_mean'],4)} | "
            f"{fmt(d3['mean_relative_change_pct'],2)} [{fmt(d3['relative_ci95_low_pct'],2)}, {fmt(d3['relative_ci95_high_pct'],2)}] | "
            f"{d3['improved_seeds']}/10 |"
        )

    lines.extend(("", "## Stationary-family aggregates", ""))
    lines.extend((
        "| arm | family | delta X [%] | delta Y [%] | delta Z [%] | delta 3D [%] | delta yaw [%] |",
        "|---|---|---:|---:|---:|---:|---:|",
    ))
    for label, _, _ in ARMS:
        if label == "deployed":
            continue
        for group_name in ("jonswap_stationary", "pmstokes_stationary", "all_stationary"):
            group = report["arms"][label]["groups"][group_name]
            lines.append(
                f"| {label} | {group_name} | "
                f"{fmt(group['disp_x_rms_m']['mean_relative_change_pct'],2)} | "
                f"{fmt(group['disp_y_rms_m']['mean_relative_change_pct'],2)} | "
                f"{fmt(group['disp_z_rms_m']['mean_relative_change_pct'],2)} | "
                f"{fmt(group['disp_3d_rms_m']['mean_relative_change_pct'],2)} | "
                f"{fmt(group['yaw_rms_deg']['mean_relative_change_pct'],2)} |"
            )

    lines.extend(("", "## Per-scenario displacement", ""))
    lines.extend((
        "| arm | scenario | delta Z [%Hs pp] | delta 3D [%] | improved 3D seeds |",
        "|---|---|---:|---:|---:|",
    ))
    for label, _, _ in ARMS:
        if label == "deployed":
            continue
        for scenario in stationary + transition:
            item = report["arms"][label]["scenarios"][scenario]
            lines.append(
                f"| {label} | {scenario} | "
                f"{fmt(item['disp_z_pct_hs']['mean_difference'],3)} | "
                f"{fmt(item['disp_3d_rms_m']['mean_relative_change_pct'],2)} | "
                f"{item['disp_3d_rms_m']['improved_seeds']}/10 |"
            )

    lines.extend(("", "## Transition segments", ""))
    lines.extend((
        "| arm | segment | delta Z [%] | delta 3D [%] |",
        "|---|---|---:|---:|",
    ))
    for label, _, _ in ARMS:
        if label == "deployed":
            continue
        seg = report["arms"][label]["transition_segments"]
        for name in ("start", "blend", "end"):
            lines.append(
                f"| {label} | {name} | "
                f"{fmt(seg[f'seg_{name}_disp_z_rms_m']['mean_relative_change_pct'],2)} | "
                f"{fmt(seg[f'seg_{name}_disp_3d_rms_m']['mean_relative_change_pct'],2)} |"
            )

    lines.extend(("", "## Controls", ""))
    lines.append(
        f"Deployed quality gates: {report['baseline_gates']['passes']}/{report['baseline_gates']['total']}."
    )
    for label, _, _ in ARMS:
        if label == "deployed":
            continue
        arm = report["arms"][label]
        inv = arm["tuner_invariance"]
        lines.append(
            f"- {label}: gates {arm['gates']['passes']}/{arm['gates']['total']}; "
            f"max paired |delta tau_applied|={fmt(inv['max_abs_tau_applied_s'],9)} s; "
            f"max paired |delta sigma_applied|={fmt(inv['max_abs_sigma_applied_mps2'],9)} m/s^2."
        )

    (output_root / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n" + "\n".join(lines), flush=True)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=REPO_ROOT / "plots/kalman_ou_ii",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "reports/results/ou_rs_candidate_full_validation",
    )
    parser.add_argument("--jobs", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--shards", type=int, default=2)
    parser.add_argument("--bootstrap", type=int, default=5000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    data_dir = args.data_dir.resolve()
    output_root = args.output_dir.resolve()
    missing = [name for name in SOURCE_FILES if not (data_dir / name).is_file()]
    if missing:
        raise FileNotFoundError("missing released simulation data: " + ", ".join(missing))
    if args.shards < 2:
        raise ValueError("use at least two shards so ou_validation stops after raw rows")

    output_root.mkdir(parents=True, exist_ok=True)
    original = patch_experiment_branch()
    try:
        for label, coeff, exponent in ARMS:
            run_arm(
                label,
                coeff,
                exponent,
                data_dir,
                output_root,
                max(1, args.jobs),
                args.shards,
            )
        all_rows = {
            label: load_rows(output_root / label / "shards")
            for label, _, _ in ARMS
        }
        build_report(all_rows, output_root, args.bootstrap)
    finally:
        HEADER.write_text(original, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
