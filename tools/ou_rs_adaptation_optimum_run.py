#!/usr/bin/env python3
"""Runner for the physical-period r_S optimum experiment.

Uses the transformations/fits in ou_rs_adaptation_optimum.py, but writes every
synthetic record with a valid WaveFileNaming-compatible basename.  Under the
fixed-acceleration time stretch p_s(t)=s^2 p(t/s), the physical displacement
height and deep-water wavelength both scale as s^2 while period scales as s.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import ou_rs_adaptation_optimum as exp
import ou_validation as core


def physical_filename(stretch: float) -> str:
    s2 = stretch * stretch
    return (
        f"wave_data_jonswap_H{1.5*s2:.3f}_L{50.710*s2:.3f}_"
        "A-30.00_P120.00.csv"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stretch", type=exp.parse_floats, default=exp.parse_floats("0.85,0.925,1,1.075,1.15"))
    parser.add_argument("--rs-factors", type=exp.parse_floats, default=exp.parse_floats("0.65,0.78,0.9,1,1.12,1.28,1.5"))
    parser.add_argument("--wave-seeds", type=exp.parse_ints, default=exp.parse_ints("11,29,47"))
    parser.add_argument("--imu-seeds", type=exp.parse_ints, default=exp.parse_ints("101,211,307"))
    parser.add_argument("--init-seeds", type=exp.parse_ints, default=exp.parse_ints("1009,1103,1201"))
    parser.add_argument("--duration-sec", type=float, default=600.0)
    parser.add_argument("--window-sec", type=float, default=300.0)
    parser.add_argument("--jobs", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--output-dir", type=Path, default=exp.REPO_ROOT / "reports/results/ou_rs_adaptation_optimum")
    args = parser.parse_args()

    if not (args.duration_sec > args.window_sec > 0.0):
        parser.error("duration-sec must exceed positive window-sec")
    if not (len(args.wave_seeds) == len(args.imu_seeds) == len(args.init_seeds)):
        parser.error("wave/imu/init seed lists must have equal length")
    if min(args.stretch) <= args.duration_sec / 1200.0:
        parser.error("selected stretch/duration can run beyond the 1200 s source")
    if not exp.SOURCE.exists():
        raise FileNotFoundError(f"missing released simulation record: {exp.SOURCE}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    core.build_simulators(("OU_III",), None)
    rs_input0 = exp.RS_BASE0 * exp.information_scale(exp.TAU0)

    rows: list[dict[str, Any]] = []
    diagnostics: list[dict[str, Any]] = []

    with tempfile.TemporaryDirectory(prefix="ou-rs-period-") as td:
        temp = Path(td)
        columns, source = core.read_wave_csv(exp.SOURCE)
        generated: dict[tuple[int, float], Path] = {}

        for repetition, wave_seed in enumerate(args.wave_seeds):
            phased = core.phase_randomize_wave(columns, source, wave_seed)
            for stretch in args.stretch:
                transformed, diag = exp.time_stretch_constant_accel(
                    columns, phased, stretch, args.duration_sec
                )
                # Keep the basename parseable by WaveFileNaming.  Separate
                # directories allow different phase realizations to share it.
                path = (
                    temp
                    / f"seed_{wave_seed}"
                    / f"stretch_{stretch:.6f}"
                    / physical_filename(stretch)
                )
                core.write_wave_csv(path, columns, transformed)
                generated[(repetition, stretch)] = path
                diagnostics.append(
                    {
                        "repetition": repetition,
                        "wave_seed": wave_seed,
                        "filename": path.name,
                        "declared_height_m": 1.5 * stretch * stretch,
                        "declared_length_m": 50.710 * stretch * stretch,
                        **diag,
                    }
                )

        tasks: list[tuple[int, float, float, Path]] = [
            (repetition, stretch, factor, generated[(repetition, stretch)])
            for repetition in range(len(args.wave_seeds))
            for stretch in args.stretch
            for factor in args.rs_factors
        ]

        def run_one(task: tuple[int, float, float, Path]) -> dict[str, Any]:
            repetition, stretch, factor, path = task
            tau = exp.TAU0 * stretch
            sigma = exp.SIGMA0
            deployed_center = rs_input0 * stretch**2.5
            desired_input = deployed_center * factor
            info_scale = exp.information_scale(tau)
            base_command = desired_input / info_scale
            point = core.TuningPoint(
                tau_s=tau,
                sigma_a_mps2=sigma,
                RS_ms=base_command,
            )
            metrics, gate, code = core.run_simulator(
                "OU_III",
                path,
                args.window_sec,
                args.imu_seeds[repetition],
                args.init_seeds[repetition],
                tuning_mode="fixed",
                tuning_point=point,
                aw_cov_sync="periodic",
                write_timeseries=False,
            )
            return {
                "repetition": repetition,
                "wave_seed": args.wave_seeds[repetition],
                "imu_seed": args.imu_seeds[repetition],
                "init_seed": args.init_seeds[repetition],
                "stretch": stretch,
                "declared_height_m": 1.5 * stretch * stretch,
                "declared_length_m": 50.710 * stretch * stretch,
                "tau_filter": tau,
                "sigma_filter": sigma,
                "pseudo_period": exp.pseudo_period(tau),
                "information_scale": info_scale,
                "rs_factor": factor,
                "rs_input_deployed_center": deployed_center,
                "rs_input": desired_input,
                "rs_base_command": base_command,
                "gate_pass": int(gate),
                "return_code": code,
                "disp_z_rms_m": float(metrics["disp_z_rms_m"]),
                "disp_3d_rms_m": float(metrics["disp_3d_rms_m"]),
                "roll_rms_deg": float(metrics.get("roll_rms_deg", math.nan)),
                "pitch_rms_deg": float(metrics.get("pitch_rms_deg", math.nan)),
                "yaw_rms_deg": float(metrics.get("yaw_rms_deg", math.nan)),
                "accel_bias_3d_rms_mps2": float(metrics.get("accel_bias_3d_rms_mps2", math.nan)),
            }

        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            pending = [pool.submit(run_one, task) for task in tasks]
            for i, future in enumerate(as_completed(pending), 1):
                row = future.result()
                rows.append(row)
                print(
                    f"[{i}/{len(pending)}] s={row['stretch']:.3f} "
                    f"r={row['rs_factor']:.3f} rep={row['repetition']} "
                    f"z={row['disp_z_rms_m']:.6f} 3d={row['disp_3d_rms_m']:.6f}",
                    flush=True,
                )

    rows.sort(key=lambda r: (r["stretch"], r["rs_factor"], r["repetition"]))
    raw_path = args.output_dir / "raw.csv"
    with raw_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    z_p, z_opts = exp.exponent_from_rows(rows, "disp_z_rms_m")
    d3_p, d3_opts = exp.exponent_from_rows(rows, "disp_3d_rms_m")
    z_boot = exp.bootstrap_exponent(rows, "disp_z_rms_m", args.bootstrap, 20260815)
    d3_boot = exp.bootstrap_exponent(rows, "disp_3d_rms_m", args.bootstrap, 20260816)

    summary = {
        "definition": "effective filter-input r_S std versus physical period at fixed translational acceleration RMS",
        "source": str(exp.SOURCE.relative_to(exp.REPO_ROOT)),
        "duration_sec": args.duration_sec,
        "window_sec": args.window_sec,
        "stretches": args.stretch,
        "rs_factors_about_deployed_2p5_center": args.rs_factors,
        "seeds": [
            {"wave": w, "imu": i, "init": n}
            for w, i, n in zip(args.wave_seeds, args.imu_seeds, args.init_seeds)
        ],
        "nominal": {
            "tau": exp.TAU0,
            "sigma": exp.SIGMA0,
            "rs_base": exp.RS_BASE0,
            "rs_input": rs_input0,
        },
        "vertical": {"p_tau": z_p, **z_boot, "optima": z_opts},
        "three_d": {"p_tau": d3_p, **d3_boot, "optima": d3_opts},
        "transform_diagnostics": diagnostics,
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    def table(opts: list[dict[str, float]]) -> list[str]:
        lines = [
            "| stretch | rS optimum | deployed center | opt/deployed |",
            "|---:|---:|---:|---:|",
        ]
        for o in opts:
            lines.append(
                f"| {o['stretch']:.3f} | {o['rs_opt']:.6f} | "
                f"{o['deployed_center']:.6f} | {o['opt_over_deployed']:.4f} |"
            )
        return lines

    md = [
        "# Physical-period r_S optimum",
        "",
        "The physical JONSWAP realization is time-stretched while translational acceleration RMS is held fixed. Tau follows physical period; sigma_aw is fixed; effective filter-input r_S is independently swept.",
        "",
        f"**Vertical-MSE exponent:** p_tau = **{z_p:.4f}** (bootstrap 95% {z_boot['bootstrap_low']:.4f} to {z_boot['bootstrap_high']:.4f})",
        f"**3-D-MSE exponent:** p_tau = **{d3_p:.4f}** (bootstrap 95% {d3_boot['bootstrap_low']:.4f} to {d3_boot['bootstrap_high']:.4f})",
        "",
        "## Vertical optimum",
        "",
        *table(z_opts),
        "",
        "## 3-D optimum",
        "",
        *table(d3_opts),
        "",
        "`rS optimum` is the standard deviation actually entering the MEKF after cadence information-rate normalization.",
    ]
    (args.output_dir / "summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n" + "\n".join(md), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
