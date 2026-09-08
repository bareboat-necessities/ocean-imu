#!/usr/bin/env python3
"""Replay current startup/magnetic ablations on the pinned vessel records."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from model_mismatch_ablation import FAMILIES, RECORDS, source_commit
from ou_validation import parse_validation_metrics
from sim_dataset import input_provenance

ROOT = Path(__file__).resolve().parents[1]
SEEDS = ("default", "11", "23", "202", "3033")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, required=True)
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--output-dir", type=Path,
                    default=ROOT / "reports/results/startup_ablation")
    args = ap.parse_args()
    paths = {r.filename: (args.data_dir / r.filename).resolve() for r in RECORDS}
    provenance = input_provenance(paths.values())
    tasks = []
    for family in ("OU-II", "OU-III"):
        for seed in SEEDS:
            for arm, env in (("deployed", {}), ("hard-iron-off", {"SF_MAG_CONT_HI": "0"})):
                tasks += [(family, seed, arm, env, r) for r in RECORDS]
    for arm, env in (("deployed", {}), ("hard-iron-off", {"TFG_MAG_HARD_IRON": "0"}),
                     ("refinement-off", {"TFG_MAG_REFINE": "0"}),
                     ("staged-startup", {"TFG_STARTUP_INIT": "staged"})):
        tasks += [("TFG", "default", arm, env, r) for r in RECORDS]

    def run(task):
        family, seed, arm, extra, record = task
        env = dict(os.environ, W3D_WRITE_TIMESERIES="0", W3D_COLLECT_ALL_GATES="1",
                   W3D_VALIDATION_WINDOW_SEC="900")
        if seed != "default":
            env["W3D_INIT_SEED"] = seed
        env.update(extra)
        binary = FAMILIES[family]
        result = subprocess.run([str(binary), "--input", str(paths[record.filename])],
                                cwd=binary.parent, env=env, capture_output=True, text=True)
        if result.returncode not in (0, 1) or (result.returncode == 1 and
                                                "QUALITY_GATE: PASS=0" not in result.stdout):
            raise RuntimeError(f"{family}/{arm}/{record.filename}: {result.returncode}\n"
                               + result.stderr[-2000:])
        metrics = parse_validation_metrics(result.stdout)
        return dict(metrics, family=family, seed=seed, arm=arm,
                    input=record.filename, spectrum=record.spectrum, hs_m=record.hs_m,
                    regression_exit_code=result.returncode)

    rows = []
    with ThreadPoolExecutor(args.jobs) as pool:
        for row in pool.map(run, tasks):
            rows.append(row)
            print(f"[{len(rows)}/{len(tasks)}] {row['family']} {row['arm']} "
                  f"seed={row['seed']} {row['input']}", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "startup_runs.csv"
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (args.output_dir / "manifest.json").write_text(json.dumps({
        "source_commit": source_commit(), "simulation_provenance": provenance,
        "protocol": "default sensor draw; final 900 s; paired initial-calibration seeds",
        "seeds": list(SEEDS), "replays": len(rows),
        "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "binary_sha256": {f: hashlib.sha256(p.read_bytes()).hexdigest()
                          for f, p in FAMILIES.items()},
        "files": {output.name: hashlib.sha256(output.read_bytes()).hexdigest()},
    }, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
