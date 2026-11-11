#!/usr/bin/env python3
"""Reproduce the comparisons and degradation studies of the CARMA method article.

Statistics use the simulator's complete 200 Hz record, never the drawing data.
Only compact, peak-preserving plot traces and per-run scores are retained.
The pinned inputs, phase randomization, seed triplets, engine model and C2
low--high--low construction are shared with the OU-III publication tools.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "plots"))
import ou_validation as ov
import ou_roundtrip_transition as rt
import engine_noise_degradation as engine
from sim_dataset import input_provenance
from vector_export import envelope_indices

OUT = ROOT / "reports/results/carma_shared_modes_validation"
CSM_SIM = ROOT / "tests/kalman_carma_shared_modes/kalman_carma_shared_modes-sim"
BINARIES = {
    "CSM": CSM_SIM,
    "Independent": CSM_SIM,
    "OU-II": ROOT / "tests/kalman_ou_ii/kalman_ou_ii-sim",
    "OU-III": ROOT / "tests/kalman_ou_iii/kalman_ou_iii-sim",
    "TFG": ROOT / "tests/kalman_tfg/kalman_tfg-sim",
    "PII": ROOT / "tests/pii_observer/pii_observer-adaptive",
    "TVG-NLO": ROOT / "tests/nlo/nlo-sim",
}
KALMAN = ("CSM", "Independent", "OU-II", "OU-III", "TFG")
TRANSITION_METHODS = ("CSM", "Independent", "OU-III")
SEEDS = tuple(zip(ov.DEFAULT_FULL_WAVE_SEEDS, ov.DEFAULT_FULL_IMU_SEEDS,
                  ov.DEFAULT_FULL_INIT_SEEDS))
FIELDS = ("disp_x_rms_m", "disp_y_rms_m", "disp_z_rms_m", "disp_3d_rms_m",
          "disp_z_pct_hs", "disp_z_ref_rms_m", "roll_rms_deg", "pitch_rms_deg",
          "yaw_rms_deg", "accel_bias_3d_rms_mps2", "gyro_bias_3d_rms_radps")
MEANS = engine.MEAN_FIELDS
SOURCES = ("tools/carma_shared_modes_validation.py", "tools/ou_validation.py",
           "tools/ou_roundtrip_transition.py", "tools/engine_noise_degradation.py",
           "tests/kalman_carma_shared_modes/kalman_carma_shared_modes-sim.cpp",
           "src/util/W3dSimCommon.cpp")


@dataclass(frozen=True)
class Cell:
    name: str
    setting: engine.Setting


def engine_cells():
    cells = [Cell("off", engine.Setting("baseline", "engine off"))]
    cells += [Cell(f"rpm{rpm}", engine.Setting("speed", str(rpm), rpm=rpm))
              for rpm in (800, 1600, 2400, 3200)]
    cells += [Cell(f"level{level:g}", engine.Setting("level", str(level),
                    rpm=2400, level_mps2=level)) for level in (0.30, 1.20)]
    cells += [Cell(f"bw{bw}", engine.Setting("bandwidth", str(bw),
                    rpm=2400, bandwidth_hz=bw)) for bw in (20, 40, 160)]
    # Same recorded vibration power as the nominal 80 Hz sensor. The probe
    # banner is independent of the sea and is evaluated before the sweep.
    cells += [Cell(f"matched{bw}", engine.Setting("matched", str(bw),
                    rpm=2400, bandwidth_hz=bw)) for bw in (20, 40, 160)]
    cells += [Cell("accel-only", engine.Setting("path", "accelerometer only",
                    rpm=2400, gyro_path=False))]
    return cells


def clean_env():
    # A user's diagnostic shell settings must not silently change the study.
    prefixes = ("W3D_", "OU_", "OU_II_", "OU_III_", "SF_", "TFG_", "NLO_", "PII_")
    return {k: v for k, v in os.environ.items() if not k.startswith(prefixes)}


def guard_env(method, guard):
    prefix = {"CSM": "W3D_CSM_GUARD", "Independent": "W3D_CSM_GUARD",
              "OU-II": "OU_II_ACC_GUARD", "OU-III": "OU_III_ACC_GUARD",
              "TFG": "TFG_ACC_GUARD"}.get(method)
    if prefix is None:
        return {}
    return {prefix + "_HZ": "0" if guard == "off" else "14",
            prefix + "_RACC_GAIN": "0.75" if guard == "deployed" else "0"}


def invoke(method, source, seeds, guard="deployed", setting=None, trace=False,
           segments=()):
    binary = BINARIES[method]
    if not binary.is_file():
        raise FileNotFoundError(f"build {binary.relative_to(ROOT)} first")
    with tempfile.TemporaryDirectory(prefix="carma-replay-") as temp:
        work = Path(temp)
        record = work / source.name
        record.symlink_to(source.resolve())
        env = clean_env()
        env.update({"W3D_IMU_SEED": str(seeds[1]), "W3D_INIT_SEED": str(seeds[2]),
                    "W3D_WRITE_TIMESERIES": "1" if trace else "0",
                    "W3D_VALIDATION_WINDOW_SEC": "900", "W3D_COLLECT_ALL_GATES": "1",
                    "W3D_AW_COV_SYNC": "periodic", "W3D_TUNING_MODE": "adaptive",
                    "W3D_VALIDATION_SEGMENTS": ",".join(f"{n}:{a}:{b}" for n,a,b in segments)})
        env.update(guard_env(method, guard))
        if method in ("CSM", "Independent"):
            env["W3D_CSM_VARIANT"] = "CSM" if method == "CSM" else "CSM_OFF"
            env["W3D_CARMA3_EXCITATION"] = "C"
        if setting is not None:
            env.update(setting.env())
        command = [str(binary)]
        if method not in ("PII", "TVG-NLO"):
            command += ["--input", str(record)]
        completed = subprocess.run(command, cwd=work, env=env, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if completed.returncode not in (0, 1):
            raise RuntimeError(f"{method}: {completed.stdout[-3000:]}")
        outputs = list(work.glob("w3d_*.csv"))
        series = None
        if trace or method in ("PII", "TVG-NLO"):
            if len(outputs) != 1:
                raise RuntimeError(f"{method}: expected one trace, got {outputs}")
            # Parse only the channels used for publication; discard raw outputs.
            wanted = lambda x: x == "time" or x.startswith(("disp_", "vel_", "acc_",
                "gyro_bias", "roll_", "pitch_", "yaw_", "tau_", "sigma_a_", "R_p0_",
                "freq_tracker", "Tp_tuner", "accel_var"))
            series = pd.read_csv(outputs[0], usecols=wanted)
        if method in ("PII", "TVG-NLO"):
            scored = series.tail(180000)
            metrics = {"samples": len(scored)}
            for axis in "xyz":
                metrics[f"disp_{axis}_rms_m"] = rms(scored[f"disp_est_{axis}"] - scored[f"disp_ref_{axis}"])
            for axis in ("roll", "pitch", "yaw"):
                delta = scored[f"{axis}_est"] - scored[f"{axis}_ref"]
                metrics[f"{axis}_rms_deg"] = rms((delta + 180) % 360 - 180)
            metrics["disp_z_ref_rms_m"] = rms(scored.disp_ref_z)
            # These vertical observers do not estimate horizontal position or
            # both sensor-bias vectors. Never score zero-filled output slots.
            for key in ("disp_x_rms_m", "disp_y_rms_m", "disp_3d_rms_m",
                        "accel_bias_3d_rms_mps2", "gyro_bias_3d_rms_radps"):
                metrics[key] = float("nan")
            if method == "TVG-NLO":
                metrics["yaw_rms_deg"] = float("nan")
        else:
            metrics = ov.parse_validation_metrics(completed.stdout)
        means = engine.parse_means(completed.stdout) if "VALIDATION_METRICS_MEANS " in completed.stdout else {}
        banner = engine.parse_banner(completed.stdout) if "ENGINE_VIBRATION " in completed.stdout else {}
        if setting is not None and setting.engine_on and not banner:
            raise RuntimeError(f"{method}: requested machinery model did not engage")
        if method in ("CSM", "Independent"):
            match = re.search(r"CSM_GUARD cutoff_hz=([\d.e+-]+) racc_gain=([\d.e+-]+)", completed.stdout)
            requested = guard_env(method, guard)
            if not match or not np.allclose([float(x) for x in match.groups()],
                    [float(requested["W3D_CSM_GUARD_HZ"]), float(requested["W3D_CSM_GUARD_RACC_GAIN"])]):
                raise RuntimeError("CARMA vibration ablation was not applied")
        return metrics, means, banner, series


def rms(values):
    return float(np.sqrt(np.mean(np.asarray(values, dtype=float) ** 2)))


def write_plot_trace(out, name, series, columns, tail=False):
    data = series.loc[series.time >= series.time.iloc[-1] - 60.0] if tail else series
    data = data[["time", *columns]]
    # Retain coordinate extrema using the same renderer helper as every vector
    # export. These files are explicitly drawing data, never scoring inputs.
    idx = envelope_indices(tuple(data[c].to_numpy() for c in data), max_points=2000)
    target = out / "traces" / f"{name}.csv.gz"
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        gz.write(data.iloc[idx].to_csv(index=False, float_format="%.9g").encode())


def store_traces(out, scenario, method, series):
    tag = f"{scenario}_{method}"
    write_plot_trace(out, tag + "_z", series, ["disp_ref_z", "disp_est_z"], tail=scenario == "nominal")
    if method != "CSM":
        return
    for kind in ("disp", "vel", "acc"):
        write_plot_trace(out, tag + "_" + kind + "_xy", series,
                         [f"{kind}_{what}_{axis}" for axis in "xy" for what in ("ref", "est")], tail=True)
        write_plot_trace(out, tag + "_" + kind + "_z", series,
                         [f"{kind}_{what}_z" for what in ("ref", "est")], tail=True)
    for axis in ("roll", "pitch", "yaw"):
        write_plot_trace(out, tag + "_" + axis, series,
                         [f"{axis}_ref", f"{axis}_est"], tail=scenario == "nominal")
    for kind in ("acc", "gyro"):
        for axis in "xyz":
            write_plot_trace(out, tag + f"_{kind}_bias_{axis}", series,
                             [f"{kind}_bias_{axis}", f"{kind}_bias_est_{axis}"])
    for col in ("tau_applied", "sigma_a_applied", "R_p0_applied", "Tp_tuner_s"):
        write_plot_trace(out, tag + "_" + col, series, [col])


def hash_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def execute(args):
    out = args.results.resolve()
    out.mkdir(parents=True, exist_ok=True)
    inputs = {r: engine.find_record(args.data_dir, r) for r in engine.RECORDS}
    protocol = {"dataset": input_provenance(inputs.values()), "seed_triplets": SEEDS[:args.seeds],
                "source_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "binary_sha256": {k: hash_file(v) for k,v in BINARIES.items()},
                "source_hashes": {p: hash_file(ROOT / p) for p in SOURCES},
                "stationary_methods": KALMAN, "observers": ("PII", "TVG-NLO"),
                "scored_seconds": 900, "sample_rate_hz": 200,
                "transition_segments": rt.roundtrip_segments(),
                "plot_geometry": "plots/vector_export.py; 2000-point drawing traces; full-rate scores"}
    manifest = out / "protocol.json"
    if manifest.exists() and json.loads(manifest.read_text()) != json.loads(json.dumps(protocol)):
        raise RuntimeError("existing results use another protocol; choose another --results directory")
    manifest.write_text(json.dumps(protocol, indent=2) + "\n")
    rows_path = out / "runs.jsonl"
    existing = [json.loads(x) for x in rows_path.read_text().splitlines()] if rows_path.exists() else []
    have = {r["key"] for r in existing}
    with tempfile.TemporaryDirectory(prefix="carma-validation-inputs-") as temp:
        cache = Path(temp)
        tasks = []
        # The vibration banner depends on the sensor model, not on vessel
        # motion. Probe just one second once, rather than replaying eight seas
        # merely to read the same nominal and bandwidth powers repeatedly.
        first_record, first_source = next(iter(inputs.items()))
        probe_columns, probe_data = ov.read_wave_csv(first_source, 1.0)
        probe = cache / "probe" / first_source.name
        probe.parent.mkdir(parents=True)
        ov.write_wave_csv(probe, probe_columns, probe_data)
        cells = engine_cells()
        nominal = next(c for c in cells if c.name == "rpm2400")
        _, _, nominal_banner, _ = invoke("OU-III", probe, SEEDS[0], "off", nominal.setting)
        target_rms = nominal_banner["recorded_rms_mps2"]
        matched = {}
        for cell in cells:
            if cell.setting.arm != "matched":
                continue
            _, _, banner, _ = invoke("OU-III", probe, SEEDS[0], "off", cell.setting)
            matched[cell.name] = engine.Setting("matched", cell.setting.label, rpm=2400,
                level_mps2=0.6 * target_rms / banner["recorded_rms_mps2"],
                bandwidth_hz=cell.setting.bandwidth_hz, matched_to_mps2=target_rms)
        for record, source in inputs.items():
            columns, data = ov.read_wave_csv(source, 1200.0)
            for seed_index, seeds in enumerate(SEEDS[:args.seeds]):
                path = cache / record.spectrum / f"{record.hs_m:g}" / str(seeds[0]) / source.name
                path.parent.mkdir(parents=True, exist_ok=True)
                ov.write_wave_csv(path, columns, ov.phase_randomize_wave(columns, data, seeds[0]))
                methods = KALMAN + (("PII", "TVG-NLO") if seed_index == 0 else ())
                for method in methods:
                    key = f"stationary/{record.filename}/{seeds[0]}/{method}"
                    trace = seed_index == 0 and record.spectrum == "JONSWAP" and record.hs_m == 1.5
                    tasks.append((key, "stationary", record, seeds, method, "deployed", None, path, trace))
                if seed_index == 0:
                    for cell in cells:
                        setting = matched.get(cell.name, cell.setting)
                        for method in TRANSITION_METHODS:
                            key = f"engine/{record.filename}/{cell.name}/{method}/off"
                            tasks.append((key, "engine", record, seeds, method, "off", setting, path, False))
                        if cell.name in ("off", "rpm800", "rpm1600", "rpm2400", "rpm3200", "level0.3", "level1.2", "bw160"):
                            for method in ("CSM", "OU-III"):
                                for guard in ("guard", "deployed"):
                                    key = f"engine/{record.filename}/{cell.name}/{method}/{guard}"
                                    tasks.append((key, "engine", record, seeds, method, guard, setting, path, False))
            del data
        low_record = next(r for r in inputs if r.spectrum == "JONSWAP" and r.hs_m == 1.5)
        high_record = next(r for r in inputs if r.spectrum == "JONSWAP" and r.hs_m == 8.5)
        columns, low = ov.read_wave_csv(inputs[low_record], 1200.0)
        _, high = ov.read_wave_csv(inputs[high_record], 1200.0)
        for seed_index, seeds in enumerate(SEEDS[:args.seeds]):
            path = cache / "roundtrip" / str(seeds[0]) / rt.ROUNDTRIP_WAVE_NAME
            path.parent.mkdir(parents=True, exist_ok=True)
            ov.write_wave_csv(path, columns, rt.make_roundtrip_wave(ov, columns, low, high, seeds[0], 4/8.5))
            for method in TRANSITION_METHODS:
                key = f"roundtrip/{seeds[0]}/{method}"
                tasks.append((key, "roundtrip", high_record, seeds, method, "deployed", None, path, seed_index == 0))
        del low, high
        pending = [t for t in tasks if t[0] not in have]
        print(f"{len(have)} completed; {len(pending)} pending", flush=True)

        def run_task(task):
            key, study, record, seeds, method, guard, setting, path, trace = task
            segments = rt.roundtrip_segments() if study == "roundtrip" else ()
            metrics, means, banner, series = invoke(method, path, seeds, guard, setting, trace, segments)
            if int(metrics["samples"]) != 180000:
                raise RuntimeError(f"{key}: expected 180000 full-rate score samples")
            if method in ("PII", "TVG-NLO"):
                metrics["disp_z_pct_hs"] = 100 * metrics["disp_z_rms_m"] / record.hs_m
            if trace:
                store_traces(out, "roundtrip" if study == "roundtrip" else "nominal", method, series)
            row = {"key": key, "study": study, "spectrum": record.spectrum, "hs_m": record.hs_m,
                   "wave_seed": seeds[0], "imu_seed": seeds[1], "init_seed": seeds[2],
                   "method": method, "guard": guard, "input": path.name,
                   "score_samples": int(metrics["samples"]),
                   "condition": key.split("/")[-3] if study == "engine" else "off",
                   **{k: metrics.get(k, float("nan")) for k in FIELDS},
                   **{k: v for k,v in metrics.items() if k.startswith("seg_")},
                   **{k: means.get(k, float("nan")) for k in MEANS},
                   **{k: banner.get(k, 0.0) for k in ("recorded_rms_mps2", "hull_rms_mps2", "firing_hz", "blade_hz", "vre_offset_mps2")}}
            if setting is not None:
                row.update(rpm=setting.rpm, level_mps2=setting.level_mps2, bandwidth_hz=setting.bandwidth_hz)
            return row

        with rows_path.open("a") as stream, ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(run_task,t): t[0] for t in pending}
            for i, future in enumerate(as_completed(futures), 1):
                row = future.result()
                stream.write(json.dumps(row, sort_keys=True) + "\n")
                stream.flush()
                if i % 10 == 0 or i == len(pending):
                    print(f"{i}/{len(pending)} {row['key']}", flush=True)
    rows = [json.loads(x) for x in rows_path.read_text().splitlines()]
    if len({r["key"] for r in rows}) != len(tasks) or len(rows) != len(tasks):
        raise RuntimeError("incomplete or duplicate run matrix")
    ov.write_csv(out / "scores.csv", sorted(rows, key=lambda r:r["key"]))
    # The append-only checkpoint is regenerable, not a second evidence copy.
    print(f"complete: {len(rows)} scored replays", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=("run", "figures", "report", "all"))
    p.add_argument("--data-dir", type=Path, default=ROOT / "tests/kalman_carma_shared_modes")
    p.add_argument("--results", type=Path, default=OUT)
    p.add_argument("--outdir", type=Path, default=ROOT / "doc/kalman_carma_shared_modes")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--seeds", type=int, default=10, choices=range(1,11))
    args = p.parse_args()
    if args.stage in ("run", "all"):
        execute(args)
    if args.stage in ("report", "figures", "all"):
        from carma_shared_modes_validation_plots import publish
        publish(args.results, args.outdir)


if __name__ == "__main__":
    main()
