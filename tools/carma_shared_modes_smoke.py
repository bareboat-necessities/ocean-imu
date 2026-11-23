#!/usr/bin/env python3
"""Run and check the paired coupling-off CARMA benchmark."""

import argparse
import csv
import os
from pathlib import Path
import subprocess


def compare_coupling_off(rows):
    reference = rows["A41@C"]
    coupling_off = rows["CSM:off"]
    rms_metrics = [prefix + axis for prefix in ("rms_", "rmsv_") for axis in "xyz"]
    axis_nees_metrics = ["nees_" + axis for axis in "xyz"]
    # NEES uses an inverse position/velocity covariance, which amplifies
    # rounding between the two float realizations. The realization test
    # compares the states and covariance directly at 1e-7 in double.
    limits = ((rms_metrics + ["nees_pos3", "nees_pv6"], 1e-4),
              (axis_nees_metrics, 5e-4))
    for metrics, limit in limits:
        for metric in metrics:
            baseline = float(reference[metric])
            difference = abs(float(coupling_off[metric]) - baseline) / max(abs(baseline), 1e-12)
            if not difference < limit:
                raise ValueError(f"CSM:off / A41@C {metric}: relative difference {difference:.6g}, limit {limit:g}")


def run_smoke(bench, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    # Match the historical band-amplitude mapping on both sides of the
    # identity. The production replay exercises default band normalization.
    with output.open("w", newline="") as stream:
        subprocess.run(
            [str(bench.resolve()), "--header", "--extended", "--coupling", "--level", "full",
             "--scenario", "jonswap_fp0.29_hs0.8_g3.3+T600", "--seed", "1",
             "--filters", "A41@C,CSM,CSM:noconf,CSM:off", "--score-start", "200"],
            stdout=stream, env={**os.environ, "CSM_BAND_NORMALIZATION": "0"}, check=True,
        )
    with output.open(newline="") as stream:
        rows = {row["filter"]: row for row in csv.DictReader(stream)}
    for name in ("A41@C", "CSM", "CSM:noconf", "CSM:off"):
        if rows[name]["failed"] != "0":
            raise ValueError(f"{name}: benchmark reported a failure")
    for name in ("CSM", "CSM:noconf", "CSM:off"):
        if rows[name]["states"] != "27":
            raise ValueError(f"{name}: expected 27 states, got {rows[name]['states']}")
    for axis in "xyz":
        if not 0.0 < float(rows["CSM"]["pct_" + axis]) < 25.0:
            raise ValueError(f"CSM: RMS as a percentage of truth range is out of bounds on {axis}")
    compare_coupling_off(rows)
    print("csm smoke:", {name: round(1000 * float(row["rms_h"]), 2) for name, row in rows.items()})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("benchmark", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    run_smoke(args.benchmark, args.output_csv)
