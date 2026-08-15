#!/usr/bin/env python3
"""Measure the full-MEKF r_S adaptation exponent against physical sea period.

This is deliberately different from an OU-parameter sensitivity sweep.  The
physical input record is time-stretched while translational acceleration RMS is
held fixed.  At each stretch the OU time constant is moved with the physical
period, sigma_aw is held fixed, and the *effective filter-input* r_S is swept
independently.  The optimum r_S is then fitted versus physical period.

The transformation is, axis by axis,

    p_s(t) = c s^2 p(t/s),   v_s(t) = c s v(t/s),   a_s(t) = c a(t/s),

with c chosen only to remove finite-record/interpolation RMS drift so the
acceleration RMS equals the s=1 realization exactly.  Attitude is time-stretched
at fixed angular amplitude and body IMU fields are rebuilt from the transformed
world motion and attitude.

The wrapper's Cubic law rescales the commanded/base r_S by sqrt(T_S0/T_S)
before it reaches Kalman3D_Wave_OU_III.  This experiment therefore converts a
desired filter-input r_S back to the command value before calling fixed tuning,
so the quantity being optimized is unambiguous.
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
from typing import Any, Iterable

import numpy as np

import ou_validation as core

REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE = REPO_ROOT / "plots/kalman_ou_ii/wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv"

# Current nominal noise-free operating point of that record.  This is only the
# centre of the sweep; no exponent is inferred from the old r_S law.
TAU0 = 2.17904091
SIGMA0 = 0.724445343
RS_BASE0 = 2.62343205657498

TS0 = 0.015
TAU_TS0 = 1.1
CT = TS0 / TAU_TS0
TS_MIN = 0.005
TS_MAX = 0.25


def parse_floats(text: str) -> list[float]:
    out = [float(x.strip()) for x in text.split(",") if x.strip()]
    if not out or any((not math.isfinite(x) or x <= 0.0) for x in out):
        raise argparse.ArgumentTypeError("expected positive finite comma-separated values")
    return out


def parse_ints(text: str) -> list[int]:
    out = [int(x.strip()) for x in text.split(",") if x.strip()]
    if not out:
        raise argparse.ArgumentTypeError("expected comma-separated integers")
    return out


def pseudo_period(tau: float) -> float:
    return min(max(CT * tau, TS_MIN), TS_MAX)


def information_scale(tau: float) -> float:
    return math.sqrt(TS0 / pseudo_period(tau))


def interp_columns(source_t: np.ndarray, values: np.ndarray, query_t: np.ndarray) -> np.ndarray:
    result = np.empty((query_t.size, values.shape[1]), dtype=np.float64)
    for axis in range(values.shape[1]):
        result[:, axis] = np.interp(query_t, source_t, values[:, axis])
    return result


def rms_zero_mean(x: np.ndarray) -> np.ndarray:
    centred = x - np.mean(x, axis=0, keepdims=True)
    return np.sqrt(np.mean(centred * centred, axis=0))


def time_stretch_constant_accel(
    columns: list[str],
    phase_record: np.ndarray,
    stretch: float,
    duration_sec: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Time-stretch one stationary realization at fixed acceleration amplitude."""

    n = max(2, int(round(duration_sec / core.DT_SECONDS)))
    if n > phase_record.shape[0]:
        raise ValueError("duration exceeds source record")

    time_idx = columns.index("time")
    source_t = phase_record[:, time_idx]
    target = phase_record[:n].copy()
    target_t = target[:, time_idx]
    t0 = float(source_t[0])
    query_t = t0 + (target_t - t0) / stretch
    if query_t[-1] > source_t[-1]:
        raise ValueError(
            f"stretch={stretch:g} with duration={duration_sec:g}s needs source t={query_t[-1]:g}, "
            f"but source ends at {source_t[-1]:g}"
        )

    disp_i = core._column_indices(columns, ("disp_x", "disp_y", "disp_z"))
    vel_i = core._column_indices(columns, ("vel_x", "vel_y", "vel_z"))
    acc_i = core._column_indices(columns, ("acc_x", "acc_y", "acc_z"))
    att_i = core._column_indices(columns, ("roll_deg", "pitch_deg", "yaw_deg"))

    p_src = phase_record[:, disp_i]
    v_src = phase_record[:, vel_i]
    a_src = phase_record[:, acc_i]
    th_src = phase_record[:, att_i]

    p_q = interp_columns(source_t, p_src, query_t)
    v_q = interp_columns(source_t, v_src, query_t)
    a_q = interp_columns(source_t, a_src, query_t)
    th_q = interp_columns(source_t, th_src, query_t)

    # Reference amplitude is the same phase realization at s=1 over the same
    # scored-duration precursor.  Renormalization is tiny but prevents a finite
    # record from masquerading as a sigma_aw change when only period is intended.
    ref_acc = phase_record[:n, acc_i]
    ref_acc_rms = rms_zero_mean(ref_acc)
    q_acc_rms = rms_zero_mean(a_q)
    linear_gain = np.divide(
        ref_acc_rms,
        q_acc_rms,
        out=np.ones_like(ref_acc_rms),
        where=q_acc_rms > 1e-12,
    )

    p_mean = np.mean(p_q, axis=0, keepdims=True)
    target[:, disp_i] = p_mean + (p_q - p_mean) * (stretch**2) * linear_gain
    target[:, vel_i] = v_q * stretch * linear_gain
    target[:, acc_i] = a_q * linear_gain

    # Keep angular excursion fixed while changing its period.  This preserves
    # the full MEKF's attitude/gravity-leakage mechanism instead of deleting it.
    ref_att = phase_record[:n, att_i]
    ref_att_std = rms_zero_mean(ref_att)
    q_att_mean = np.mean(th_q, axis=0, keepdims=True)
    q_att_std = rms_zero_mean(th_q)
    att_gain = np.divide(
        ref_att_std,
        q_att_std,
        out=np.ones_like(ref_att_std),
        where=q_att_std > 1e-12,
    )
    target[:, att_i] = q_att_mean + (th_q - q_att_mean) * att_gain

    target = core.rebuild_body_imu(columns, target)

    out_acc_rms = rms_zero_mean(target[:, acc_i])
    out_att_std = rms_zero_mean(target[:, att_i])
    diag = {
        "stretch": stretch,
        "acc_rms_ref": ref_acc_rms.tolist(),
        "acc_rms_out": out_acc_rms.tolist(),
        "acc_rms_ratio": (out_acc_rms / ref_acc_rms).tolist(),
        "att_std_ref": ref_att_std.tolist(),
        "att_std_out": out_att_std.tolist(),
        "att_std_ratio": np.divide(out_att_std, ref_att_std, out=np.ones_like(out_att_std), where=ref_att_std > 1e-12).tolist(),
        "linear_gain": linear_gain.tolist(),
        "att_gain": att_gain.tolist(),
    }
    return target, diag


def mean_mse(rows: Iterable[dict[str, Any]], metric: str) -> float:
    vals = np.asarray([float(r[metric]) for r in rows], dtype=np.float64)
    return float(np.mean(vals * vals))


def optimum_for_stretch(rows: list[dict[str, Any]], metric: str) -> dict[str, float]:
    factors = sorted({float(r["rs_factor"]) for r in rows})
    grouped = []
    for factor in factors:
        grp = [r for r in rows if math.isclose(float(r["rs_factor"]), factor, rel_tol=0.0, abs_tol=1e-12)]
        grouped.append((factor, float(grp[0]["rs_input"]), mean_mse(grp, metric)))

    y = np.asarray([g[2] for g in grouped])
    j = int(np.argmin(y))
    # Fit a local quadratic in log r_S. Five points when available, centred on
    # the discrete minimum. This is the Hessian approximation required by the
    # implicit-optimum derivative, not a global polynomial through the tails.
    count = min(5, len(grouped))
    lo = max(0, min(j - count // 2, len(grouped) - count))
    local = grouped[lo:lo + count]
    x = np.log(np.asarray([g[1] for g in local]))
    yy = np.asarray([g[2] for g in local])
    a, b, c = np.polyfit(x, yy, 2)
    if not (math.isfinite(a) and a > 0.0):
        xopt = math.log(grouped[j][1])
    else:
        xopt = -b / (2.0 * a)
        xopt = float(np.clip(xopt, np.min(x), np.max(x)))
    ropt = math.exp(xopt)
    return {
        "rs_opt": ropt,
        "mse_opt_fit": float(a * xopt * xopt + b * xopt + c),
        "curvature_log_rs": float(2.0 * a),
        "discrete_rs": grouped[j][1],
        "discrete_mse": grouped[j][2],
    }


def exponent_from_rows(rows: list[dict[str, Any]], metric: str) -> tuple[float, list[dict[str, float]]]:
    opts: list[dict[str, float]] = []
    for stretch in sorted({float(r["stretch"]) for r in rows}):
        grp = [r for r in rows if math.isclose(float(r["stretch"]), stretch, rel_tol=0.0, abs_tol=1e-12)]
        opt = optimum_for_stretch(grp, metric)
        deployed = float(grp[0]["rs_input_deployed_center"])
        opts.append({"stretch": stretch, "deployed_center": deployed, **opt, "opt_over_deployed": opt["rs_opt"] / deployed})
    v = np.log(np.asarray([o["stretch"] for o in opts]))
    u = np.log(np.asarray([o["rs_opt"] for o in opts]))
    slope, intercept = np.polyfit(v, u, 1)
    return float(slope), opts


def bootstrap_exponent(rows: list[dict[str, Any]], metric: str, resamples: int, seed: int) -> dict[str, float]:
    seed_ids = sorted({int(r["repetition"]) for r in rows})
    if len(seed_ids) < 2 or resamples <= 0:
        return {"bootstrap_low": math.nan, "bootstrap_high": math.nan, "bootstrap_median": math.nan}
    by_seed = {sid: [r for r in rows if int(r["repetition"]) == sid] for sid in seed_ids}
    rng = np.random.default_rng(seed)
    values: list[float] = []
    for _ in range(resamples):
        draw = rng.choice(seed_ids, size=len(seed_ids), replace=True)
        sample: list[dict[str, Any]] = []
        for new_rep, sid in enumerate(draw):
            for row in by_seed[int(sid)]:
                copy = dict(row)
                copy["repetition"] = new_rep
                sample.append(copy)
        try:
            p, _ = exponent_from_rows(sample, metric)
        except Exception:
            continue
        if math.isfinite(p):
            values.append(p)
    if not values:
        return {"bootstrap_low": math.nan, "bootstrap_high": math.nan, "bootstrap_median": math.nan}
    arr = np.asarray(values)
    low, median, high = np.quantile(arr, (0.025, 0.5, 0.975))
    return {"bootstrap_low": float(low), "bootstrap_high": float(high), "bootstrap_median": float(median)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stretch", type=parse_floats, default=parse_floats("0.85,0.925,1,1.075,1.15"))
    parser.add_argument("--rs-factors", type=parse_floats, default=parse_floats("0.65,0.78,0.9,1,1.12,1.28,1.5"))
    parser.add_argument("--wave-seeds", type=parse_ints, default=parse_ints("11,29,47"))
    parser.add_argument("--imu-seeds", type=parse_ints, default=parse_ints("101,211,307"))
    parser.add_argument("--init-seeds", type=parse_ints, default=parse_ints("1009,1103,1201"))
    parser.add_argument("--duration-sec", type=float, default=600.0)
    parser.add_argument("--window-sec", type=float, default=300.0)
    parser.add_argument("--jobs", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--output-dir", type=Path, default=REPO_ROOT / "reports/results/ou_rs_adaptation_optimum")
    args = parser.parse_args()

    if not (args.duration_sec > args.window_sec > 0.0):
        parser.error("duration-sec must exceed positive window-sec")
    if not (len(args.wave_seeds) == len(args.imu_seeds) == len(args.init_seeds)):
        parser.error("wave/imu/init seed lists must have equal length")
    if min(args.stretch) <= args.duration_sec / 1200.0:
        parser.error("selected stretch/duration can run beyond the 1200 s source")
    if not SOURCE.exists():
        raise FileNotFoundError(f"missing released simulation record: {SOURCE}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    core.build_simulators(("OU_III",), None)

    # Nominal effective r_S actually entering the MEKF, after cadence matching.
    rs_input0 = RS_BASE0 * information_scale(TAU0)

    rows: list[dict[str, Any]] = []
    diagnostics: list[dict[str, Any]] = []

    with tempfile.TemporaryDirectory(prefix="ou-rs-period-") as td:
        temp = Path(td)
        columns, source = core.read_wave_csv(SOURCE)

        generated: dict[tuple[int, float], Path] = {}
        for repetition, wave_seed in enumerate(args.wave_seeds):
            phased = core.phase_randomize_wave(columns, source, wave_seed)
            for stretch in args.stretch:
                transformed, diag = time_stretch_constant_accel(columns, phased, stretch, args.duration_sec)
                path = temp / f"period_s{stretch:.6f}_seed{wave_seed}.csv"
                core.write_wave_csv(path, columns, transformed)
                generated[(repetition, stretch)] = path
                diagnostics.append({"repetition": repetition, "wave_seed": wave_seed, **diag})

        tasks: list[tuple[int, float, float, Path]] = []
        for repetition, _ in enumerate(args.wave_seeds):
            for stretch in args.stretch:
                for factor in args.rs_factors:
                    tasks.append((repetition, stretch, factor, generated[(repetition, stretch)]))

        def run_one(task: tuple[int, float, float, Path]) -> dict[str, Any]:
            repetition, stretch, factor, path = task
            tau = TAU0 * stretch
            sigma = SIGMA0
            # Centre the independent sweep on the currently deployed *input*
            # exponent. The factor grid is independent, so the fit can reject it.
            deployed_center = rs_input0 * stretch**2.5
            desired_input = deployed_center * factor
            scale = information_scale(tau)
            base_command = desired_input / scale
            point = core.TuningPoint(tau_s=tau, sigma_a_mps2=sigma, RS_ms=base_command)
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
                "tau_filter": tau,
                "sigma_filter": sigma,
                "pseudo_period": pseudo_period(tau),
                "information_scale": scale,
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
            futures = [pool.submit(run_one, task) for task in tasks]
            for i, future in enumerate(as_completed(futures), 1):
                row = future.result()
                rows.append(row)
                print(
                    f"[{i}/{len(futures)}] s={row['stretch']:.3f} "
                    f"r={row['rs_factor']:.3f} rep={row['repetition']} "
                    f"z={row['disp_z_rms_m']:.6f} 3d={row['disp_3d_rms_m']:.6f}",
                    flush=True,
                )

    rows.sort(key=lambda r: (r["stretch"], r["rs_factor"], r["repetition"]))
    raw_path = args.output_dir / "raw.csv"
    with raw_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    z_p, z_opts = exponent_from_rows(rows, "disp_z_rms_m")
    d3_p, d3_opts = exponent_from_rows(rows, "disp_3d_rms_m")
    z_boot = bootstrap_exponent(rows, "disp_z_rms_m", args.bootstrap, 20260815)
    d3_boot = bootstrap_exponent(rows, "disp_3d_rms_m", args.bootstrap, 20260816)

    summary = {
        "definition": "effective filter-input r_S standard deviation versus physical period at fixed translational acceleration RMS",
        "source": str(SOURCE.relative_to(REPO_ROOT)),
        "duration_sec": args.duration_sec,
        "window_sec": args.window_sec,
        "stretches": args.stretch,
        "rs_factors_about_deployed_2p5_center": args.rs_factors,
        "seeds": [
            {"wave": w, "imu": i, "init": n}
            for w, i, n in zip(args.wave_seeds, args.imu_seeds, args.init_seeds)
        ],
        "nominal": {
            "tau": TAU0,
            "sigma": SIGMA0,
            "rs_base": RS_BASE0,
            "rs_input": rs_input0,
        },
        "vertical": {"p_tau": z_p, **z_boot, "optima": z_opts},
        "three_d": {"p_tau": d3_p, **d3_boot, "optima": d3_opts},
        "transform_diagnostics": diagnostics,
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    def table_lines(opts: list[dict[str, float]]) -> list[str]:
        out = ["| stretch | rS optimum | deployed center | ratio |", "|---:|---:|---:|---:|"]
        for o in opts:
            out.append(
                f"| {o['stretch']:.3f} | {o['rs_opt']:.6f} | {o['deployed_center']:.6f} | {o['opt_over_deployed']:.4f} |"
            )
        return out

    md = [
        "# Physical-period r_S optimum",
        "",
        "The physical record is time-stretched while acceleration RMS is held fixed; r_S is then independently re-optimized at every stretch.",
        "",
        f"**Vertical-MSE exponent:** p_tau = **{z_p:.4f}** (bootstrap 95% {z_boot['bootstrap_low']:.4f} to {z_boot['bootstrap_high']:.4f})",
        f"**3-D-MSE exponent:** p_tau = **{d3_p:.4f}** (bootstrap 95% {d3_boot['bootstrap_low']:.4f} to {d3_boot['bootstrap_high']:.4f})",
        "",
        "## Vertical optimum",
        "",
        *table_lines(z_opts),
        "",
        "## 3-D optimum",
        "",
        *table_lines(d3_opts),
        "",
        "`rS optimum` is the standard deviation actually entering the MEKF after the wrapper's cadence information-rate normalization.",
    ]
    (args.output_dir / "summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n" + "\n".join(md), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
