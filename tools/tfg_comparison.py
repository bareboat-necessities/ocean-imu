#!/usr/bin/env python3
"""Shipping OU-III/TFG companion study; shared histories, no estimator retuning.

Full evidence is nine scenarios x ten seed triplets x two adaptive estimators.
The primary companion contrast is TFG minus OU-III in mean JONSWAP Z %Hs,
first averaging the four heights within each seed. Other endpoints are
exploratory, not multiplicity-adjusted confirmations. PM-Stokes is separate.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np

import ou_evidence_provenance as provenance
import ou_validation as ov
import sim_dataset

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports/results/tfg_comparison"
DOC = ROOT / "doc/kalman_ou_iii/w3d-tfg-comparison-generated.tex-part"
BINARY = {
    "OU_III": ROOT / "tests/kalman_ou_iii/kalman_ou_iii-sim",
    "TFG": ROOT / "tests/kalman_tfg/kalman_tfg-sim",
}
METRICS = (
    "disp_z_pct_hs", "disp_z_rms_m", "disp_3d_rms_m",
    "roll_rms_deg", "pitch_rms_deg", "yaw_rms_deg",
    "accel_3d_rms_mps2", "accel_bias_3d_rms_mps2", "gyro_bias_3d_rms_radps",
)
LABELS = (
    r"Z RMS [$\%H_s$]", "Z RMS [m]", "3-D RMS [m]",
    "Roll RMS [deg]", "Pitch RMS [deg]", "Yaw RMS [deg]",
    r"Accel. RMS [m/s$^2$]", r"Accel. bias [m/s$^2$]", "Gyro bias [rad/s]",
)
TRANSITION = "nonstationary_H1_5_to_H4_0_Tp5_7_to_11_4"
SEED_FIELDS = ("wave_phase_seed", "imu_noise_seed", "initialization_seed")
ENV_PREFIXES = ("W3D_", "TFG_", "SF_", "OU_II_", "OU_III_", "OU_SIGMA_")


def protocol(mode: str) -> dict:
    seeds = [asdict(s) for s in ov.broadcast_seed_triplets(
        ov.DEFAULT_FULL_WAVE_SEEDS, ov.DEFAULT_FULL_IMU_SEEDS, ov.DEFAULT_FULL_INIT_SEEDS
    )]
    inputs = sorted(sim_dataset.REFERENCE_CSV_SHA256)
    if mode == "smoke":
        seeds = seeds[:1]
        inputs = [n for n in inputs if "jonswap_H1.500_" in n]
    return {
        "mode": mode, "duration_sec": 1200.0 if mode == "full" else 180.0,
        "window_sec": 900.0 if mode == "full" else 60.0,
        "dt_sec": ov.DT_SECONDS, "seed_triplets": seeds, "inputs": inputs,
        "scenarios": [ov._scenario_slug(Path(n)) for n in inputs] + [TRANSITION],
        "families": list(BINARY), "tuning": "shipping adaptive; periodic a_w covariance sync",
        "bootstrap_resamples": 10000, "stats_seed": 20260317,
        "primary": "JONSWAP seed-level mean Z %Hs; difference TFG minus OU-III",
        "pooling": ov.PMSTOKES_POOLING, "phase_method": ov.WAVE_PHASE_METHOD,
        "transition_method": ov.TRANSITION_METHOD,
        "sampling_unit": "seed triplet, not record or time sample",
        "transition_normalization": "endpoint incident Hs=4 m, matching the OU validation protocol",
        "inclusion": "all complete finite replays; regression gates never filter statistical rows",
        "scope": "development RAO ensemble; exploratory companion comparison, no unseen vessel claim",
    }


def environment(family: str, seed: dict, p: dict) -> dict[str, str]:
    # Clear inherited experiments, but do not overwrite each family's tuned defaults.
    env = {k: v for k, v in os.environ.items() if not k.startswith(ENV_PREFIXES)}
    env.update(W3D_IMU_SEED=str(seed["imu_noise_seed"]),
               W3D_INIT_SEED=str(seed["initialization_seed"]),
               W3D_WRITE_TIMESERIES="0", W3D_COLLECT_ALL_GATES="1",
               W3D_VALIDATION_WINDOW_SEC=str(p["window_sec"]))
    if family == "TFG":
        env.update(TFG_TUNING="adaptive", TFG_AW_COV_SYNC="1")
    else:
        env.update(W3D_TUNING_MODE="adaptive", W3D_AW_COV_SYNC="periodic")
    return env


def run_pair(job: tuple, p: dict) -> tuple[list[dict], list[dict]]:
    scenario, seed, columns, source, endpoint, filename = job
    with tempfile.TemporaryDirectory(prefix="tfg-pair-") as temp:
        path = Path(temp) / filename
        data = (ov.phase_randomize_wave(columns, source, seed["wave_phase_seed"])
                if endpoint is None else ov.make_nonstationary_wave(
                    columns, source, endpoint, seed["wave_phase_seed"], 4.0 / 8.5,
                    0.45 * p["duration_sec"], 0.55 * p["duration_sec"]))
        ov.write_wave_csv(path, columns, data)
        digest = ov.sha256_file(path)
        rows, gates = [], []
        del data
        for family, binary in BINARY.items():
            completed = subprocess.run([str(binary), "--input", str(path)],
                cwd=binary.parent, env=environment(family, seed, p), text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=600, check=False)
            output = completed.stdout
            gate = re.findall(r"QUALITY_GATE: PASS=([01])", output)
            skipped = (p["mode"] == "smoke" and
                       "QUALITY_GATE: SKIPPED REASON=record_shorter_than_900s_window" in output)
            if completed.returncode not in (0, 1) or (not gate and not skipped):
                raise RuntimeError(f"{family}/{scenario}: exit {completed.returncode}\n{output[-6000:]}")
            metrics = ov.parse_validation_metrics(output)
            if completed.returncode != (0 if skipped or gate[-1] == "1" else 1):
                raise ValueError("exit code does not match the reported regression gate")
            if metrics["family"] != family:
                raise ValueError("wrong simulator family")
            base = {"scenario": scenario, **seed, "family": family,
                    "input_sha256": digest, "samples": int(metrics["samples"]),
                    "window_s": metrics["window_s"]}
            rows.append({**base, **{m: float(metrics[m]) for m in METRICS}})
            # These are diagnostics only: randomized histories need not satisfy
            # regression bars derived from the deterministic reference records.
            gates.append({**base, "return_code": completed.returncode,
                          "regression_gate_pass": None if skipped else gate[-1] == "1"})
        print(f"Completed {scenario}, wave seed {seed['wave_phase_seed']}", flush=True)
        return rows, gates


def validate_rows(rows: list[dict], p: dict) -> None:
    expected = {(s, *(seed[f] for f in SEED_FIELDS), family)
                for s in p["scenarios"] for seed in p["seed_triplets"] for family in BINARY}
    found = {}
    for row in rows:
        key = (row["scenario"], *(row[f] for f in SEED_FIELDS), row["family"])
        if key in found or key not in expected:
            raise ValueError(f"duplicate or unexpected replay: {key}")
        found[key] = row
        if row["samples"] != round(p["window_sec"] / p["dt_sec"]) or row["window_s"] != p["window_sec"]:
            raise ValueError(f"wrong scoring window: {key}")
        if not all(math.isfinite(row[m]) and row[m] >= 0.0 for m in METRICS):
            raise ValueError(f"invalid scored metric: {key}")
        if not re.fullmatch(r"[0-9a-f]{64}", row["input_sha256"]):
            raise ValueError("missing replay input hash")
    if set(found) != expected:
        raise ValueError(f"incomplete paired evidence: {len(found)} of {len(expected)}")
    for key, row in found.items():
        partner = found[(*key[:-1], "TFG" if key[-1] == "OU_III" else "OU_III")]
        if row["input_sha256"] != partner["input_sha256"]:
            raise ValueError("paired estimators saw different histories")


def summarize(rows: list[dict], p: dict) -> list[dict]:
    validate_rows(rows, p)
    result = []
    groups = {name: [name] for name in p["scenarios"]}
    groups.update({spectrum: [s for s in p["scenarios"] if ov.scenario_spectrum(s) == spectrum]
                   for spectrum in ("jonswap", "pmstokes")})
    for group, cases in groups.items():
        if not cases:
            continue
        for metric in METRICS:
            values = {}
            for family in BINARY:
                values[family] = np.array([
                    np.mean([row[metric] for row in rows if row["family"] == family
                             and row["scenario"] in cases
                             and all(row[f] == seed[f] for f in SEED_FIELDS)])
                    for seed in p["seed_triplets"]])
            diff = values["TFG"] - values["OU_III"]
            # Independent deterministic stream per endpoint; adding a display
            # metric cannot change an already specified bootstrap interval.
            salt = int(hashlib.sha256(f"{group}/{metric}".encode()).hexdigest()[:8], 16)
            rng = np.random.default_rng(p["stats_seed"] + salt)
            low, high = ov._bootstrap_mean_ci(diff, p["bootstrap_resamples"], rng)
            result.append({"group": group, "metric": metric,
                "ou3_mean": float(values["OU_III"].mean()),
                "tfg_mean": float(values["TFG"].mean()),
                "difference": float(diff.mean()), "ci95_low": low, "ci95_high": high,
                "seed_differences": diff.tolist(), **ov.paired_inference(diff, rng)})
    return result


def publication(summary: list[dict]) -> str:
    index = {(row["group"], row["metric"]): row for row in summary}
    primary = index[("jonswap", "disp_z_pct_hs")]
    verdict = ("TFG has lower mean normalized heave error in this ensemble."
               if primary["ci95_high"] < 0 else
               "OU--III has lower mean normalized heave error in this ensemble."
               if primary["ci95_low"] > 0 else
               "The paired interval does not resolve a normalized-heave advantage.")
    lines = ["% Generated by tools/tfg_comparison.py; full evidence only.",
             r"\begin{table}[t]\centering\footnotesize",
             r"\caption{Ten-seed JONSWAP companion comparison. Each seed averages four heights. Differences are TFG minus OU--III; negative favors TFG. Other channels are exploratory.}",
             r"\label{tab:tfg-paired}\setlength{\tabcolsep}{3pt}",
             r"\begin{tabular}{@{}lrrr@{}}\toprule",
             r"Metric & OU--III & TFG & Difference\\\midrule"]
    for metric, label in zip(METRICS, LABELS):
        row = index[("jonswap", metric)]
        lines.append(f"{label} & {row['ou3_mean']:.4g} & {row['tfg_mean']:.4g} & {row['difference']:+.3g}" + r"\\")
    lines += [r"\bottomrule\end{tabular}\end{table}",
        f"The companion JONSWAP heave difference is {primary['difference']:+.3f} "
        f"percentage points of $H_s$ (paired-bootstrap 95\\% interval "
        f"[{primary['ci95_low']:+.3f}, {primary['ci95_high']:+.3f}]; "
        f"paired Student-$t$ interval [{primary['t_ci95_low']:+.3f}, {primary['t_ci95_high']:+.3f}]). "
        f"The exact sign and sign-flip $p$-values are {primary['sign_p_value']:.4f} "
        f"and {primary['randomization_p_value']:.4f}, respectively. " + verdict,
        "These are ten paired seed differences, not forty independent seas.",
        r"\begin{table}[t]\centering\footnotesize",
        r"\caption{Companion heave contrasts by declared ensemble. The transition is a controlled crossfade, not an evolving random sea. Intervals are paired-bootstrap 95\%.}",
        r"\label{tab:tfg-ensembles}\setlength{\tabcolsep}{3pt}",
        r"\begin{tabular}{@{}lrrl@{}}\toprule",
        r"Ensemble & OU--III & TFG & Difference [CI]\\\midrule"]
    for group, label in (("jonswap", "JONSWAP"), ("pmstokes", "PM--Stokes"), (TRANSITION, "Crossfade")):
        row = index[(group, "disp_z_pct_hs")]
        lines.append(f"{label} & {row['ou3_mean']:.2f} & {row['tfg_mean']:.2f} & "
                     f"{row['difference']:+.2f} [{row['ci95_low']:+.2f}, {row['ci95_high']:+.2f}]" + r"\\")
    lines += [r"\bottomrule\end{tabular}\end{table}"]
    return "\n".join(lines) + "\n"


def source_hashes() -> dict[str, str]:
    roots = [ROOT / "src/util/W3dSimCommon.cpp"] + [b.with_suffix(".cpp") for b in BINARY.values()]
    files = set(provenance.implementation_closure(roots))
    # Simulator-only experiment hooks must not invalidate committed shipping\n    # comparison evidence. The evidence hashes the shipping estimator closure;\n    # dedicated tuning studies hash their own instrumented simulator binary.\n    files.discard(ROOT / "tests/kalman_tfg/kalman_tfg-sim.cpp")
    files.update(b.parent / "Makefile" for b in BINARY.values())
    files.update(ROOT / n for n in ("tools/tfg_comparison.py", "tools/ou_validation.py",
        "tools/ou_evidence_provenance.py", "tools/sim_dataset.py"))
    return {str(p.relative_to(ROOT)): ov.sha256_file(p) for p in sorted(files)}


def check_bundle(output: Path, require_ci: bool = False) -> None:
    manifest = json.loads((output / "manifest.json").read_text())
    if manifest["protocol"] != protocol("full"):
        raise ValueError("only the complete declared full protocol is publication evidence")
    if manifest["source_hashes"] != source_hashes():
        raise ValueError("TFG comparison source provenance is stale")
    if require_ci and not manifest.get("workflow"):
        raise ValueError("publication requires CI replay provenance")
    required = {"raw_runs.json", "raw_runs.csv", "paired_summary.json", "gate_diagnostics.json", "comparison.tex"}
    if set(manifest["result_sha256"]) != required:
        raise ValueError("incomplete evidence inventory")
    if manifest["dataset"].get("input_sha256") != sim_dataset.REFERENCE_CSV_SHA256:
        raise ValueError("reference dataset provenance differs from pinned release")
    for name, digest in manifest["result_sha256"].items():
        if ov.sha256_file(output / name) != digest:
            raise ValueError(f"changed evidence bytes: {name}")
    rows = json.loads((output / "raw_runs.json").read_text())
    expected = ov._json_safe(summarize(rows, manifest["protocol"]))
    if expected != json.loads((output / "paired_summary.json").read_text()):
        raise ValueError("paired statistics do not reproduce from retained rows")
    if publication(expected) != (output / "comparison.tex").read_text():
        raise ValueError("publication does not reproduce from paired statistics")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("full", "smoke"), default="full")
    parser.add_argument("--output-dir", type=Path, default=OUT)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "plots/kalman_ou_ii")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--require-ci", action="store_true")
    parser.add_argument("--publish", action="store_true")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if args.check:
        check_bundle(output, args.require_ci)
    else:
        if args.jobs < 1:
            parser.error("jobs must be positive")
        if args.mode == "smoke" and output == OUT:
            parser.error("smoke must use a separate output directory")
        p = protocol(args.mode)
        hashes = source_hashes()
        data_paths = [args.data_dir / n for n in sim_dataset.REFERENCE_CSV_SHA256]
        data_provenance = sim_dataset.input_provenance(data_paths)
        if not args.skip_build:
            for b in BINARY.values():
                subprocess.run(["make", "-C", str(b.parent), b.name], check=True)
        binaries = {family: ov.sha256_file(b) for family, b in BINARY.items()}
        data = {path.name: ov.read_wave_csv(path, p["duration_sec"]) for path in data_paths}
        start = next(n for n in data if "jonswap_H1.500_" in n)
        end = next(n for n in data if "jonswap_H8.500_" in n)
        jobs = [(ov._scenario_slug(Path(n)), seed, *data[n], None, n)
                for n in p["inputs"] for seed in p["seed_triplets"]]
        if data[start][0] != data[end][0]:
            raise ValueError("transition source columns differ")
        jobs += [(TRANSITION, seed, *data[start], data[end][1], "wave_data_jonswap_H4.000_L202.839_A-30.00_P120.00.csv")
                 for seed in p["seed_triplets"]]
        rows, gates = [], []
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for pair, diagnostic in pool.map(lambda job: run_pair(job, p), jobs):
                rows.extend(pair)
                gates.extend(diagnostic)
        summary = summarize(rows, p)
        if hashes != source_hashes() or binaries != {f: ov.sha256_file(b) for f, b in BINARY.items()}:
            raise ValueError("source or simulator changed during replay")
        output.mkdir(parents=True, exist_ok=True)
        for name, value in (("raw_runs.json", rows), ("paired_summary.json", summary), ("gate_diagnostics.json", gates)):
            ov.write_json(output / name, value)
        ov.write_csv(output / "raw_runs.csv", rows)
        names = ["raw_runs.json", "raw_runs.csv", "paired_summary.json", "gate_diagnostics.json"]
        if args.mode == "full":
            (output / "comparison.tex").write_text(publication(summary))
            names.append("comparison.tex")
        ov.write_json(output / "manifest.json", {
            "schema_version": 1, "protocol": p, "source_hashes": hashes,
            "binary_sha256": binaries, "dataset": data_provenance,
            "git_commit": ov.git_output("rev-parse", "HEAD"),
            "git_diff_stat": ov.git_output("diff", "--stat"),
            "workflow": provenance.workflow_metadata(),
            "build_environment": provenance.environment_metadata(),
            "generated_utc": provenance.utc_now(),
            "result_sha256": {name: ov.sha256_file(output / name) for name in names},
        })
    if args.publish:
        check_bundle(output, args.require_ci)
        DOC.write_text((output / "comparison.tex").read_text())


if __name__ == "__main__":
    main()
