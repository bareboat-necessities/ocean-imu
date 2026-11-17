#!/usr/bin/env python3
"""Paired OU-III Live handover component ablation on pinned, identical sensor prefixes."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from itertools import islice, product
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
import pandas as pd

from sim_dataset import ROOT, CASES, REFERENCE_CSV_SHA256_BY_CASE, digest, input_provenance

# Keep the intervention at the actual assignment site. The shipping replay and
# all its sensor, covariance, prediction, refinement and release code are reused.
ARMS = {
    "previous": 0,
    "full": 15,
    "previous_tilt": 14,
    "previous_yaw": 13,
    "previous_field": 11,
    "previous_gyro_bias": 7,
    "gravity_axis_accel_bias": 31,
}
WINDOWS = ((0, 1), (1, 10), (10, 60))
SECONDS = 110
DT = 0.005


def substitute(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f"expected one shipping source anchor: {old}")
    return text.replace(old, new)


def build(work, control=False):
    source = ROOT / "src/kalman_common"
    overlay = work / ("control" if control else "instrumented") / "kalman_common"
    overlay.mkdir(parents=True)
    learner = (source / "WaveStartup.h").read_text()
    learner = substitute(
        learner, "V gyro_bias = V::Zero();", "V gyro_bias = V::Zero();\n        V accel_bias = V::Zero();"
    )
    # With known gravity, the mean force supplies one scalar bias observation.
    # Use the same uniform force and rotation moments, carry the inferred slow
    # bias along the mean body-down direction, and retain the zero prior on
    # orthogonal directions. Skip poorly conditioned or out-of-range seeds.
    learner = substitute(
        learner,
        "out = candidate;",
        """out = candidate;
        const V axis = integrated_rotation_.transpose() * down / static_cast<float>(time_);
        if (axis.squaredNorm() > 0.25f) {
            const V bias = axis * ((gravity + down.dot(uniform)) / axis.squaredNorm());
            if (bias.allFinite() && bias.norm() <= 0.4f) out.accel_bias = bias;
        }""",
    )
    if not control:
        (overlay / "WaveStartup.h").write_text(learner)
    wrapper = (source / "ProxyStartupFusion.h").read_text()
    old = """            q_seed = learned.attitude;
            setMagWorldRef_(learned.field);
            impl_.mekf().set_initial_gyro_bias_body(learned.gyro_bias);"""
    new = """            const unsigned mask = std::stoul(std::getenv("W3D_HANDOVER_MASK"));
            if ((mask & 3u) == 3u) q_seed = learned.attitude;
            else if (mask & 1u) q_seed = boatQuatWithAbsoluteYaw(learned.attitude, yawFromBoatQuatRad(q_seed));
            else if (mask & 2u) q_seed = boatQuatWithAbsoluteYaw(q_seed, yawFromBoatQuatRad(learned.attitude));
            if (mask & 4u) setMagWorldRef_(learned.field);
            if (mask & 8u) impl_.mekf().set_initial_gyro_bias_body(learned.gyro_bias);
            if (mask & 16u) impl_.mekf().set_initial_acc_bias(learned.accel_bias);
            std::cerr << "HANDOVER_COMPONENTS " << mask << "\\n";"""
    if control:
        # Select the existing public enable flag for the two intact policies.
        # Neither the fit, the seed assignments, nor the live estimator is instrumented.
        wrapper = substitute(
            wrapper,
            "bool wave_startup_enabled = true;",
            'bool wave_startup_enabled = std::stoi(std::getenv("W3D_HANDOVER_MASK")) != 0;',
        )
    else:
        wrapper = substitute(wrapper, old, new)
    (overlay / "ProxyStartupFusion.h").write_text(wrapper)
    binary = work / ("production-control" if control else "ablation")
    eigen = ROOT / "third_party/eigen"
    if not (eigen / "Eigen/Dense").is_file():
        eigen = Path("/usr/include/eigen3")
    cmd = [
        "g++",
        "-O2",
        "-std=c++20",
        "-fno-finite-math-only",
        "-march=x86-64-v3",
        "-I" + str(overlay.parent),
        "-I" + str(ROOT / "src"),
        "-isystem",
        str(eigen),
        str(ROOT / "tests/kalman_ou_iii/kalman_ou_iii-sim.cpp"),
        str(ROOT / "src/util/W3dSimCommon.cpp"),
        "-o",
        str(binary),
    ]
    subprocess.run(cmd, check=True)
    return binary


def score(path, live):
    columns = ["time"] + [f"{name}_{side}" for name in ("roll", "pitch", "yaw") for side in ("ref", "est")]
    columns += [f"{kind}_bias{suffix}_{axis}" for kind in ("gyro", "acc") for suffix in ("", "_est") for axis in "xyz"]
    data = pd.read_csv(path, usecols=columns)
    # Match the wrapper's float clock sample by sample. Epoch zero supplies
    # its first dt, and the printed handover time has six significant digits.
    clocks = np.cumsum(np.full(len(data), DT, dtype=np.float32), dtype=np.float32)
    first = int(np.argmin(np.abs(clocks - live)))
    if abs(float(clocks[first]) - live) > 0.0001:
        raise ValueError(f"no sensor epoch matches handover {live}")
    start = float(data.time.iloc[first])
    pitch, roll = {}, {}
    for side in ("est", "ref"):
        pitch[side] = np.deg2rad(data[f"pitch_{side}"].to_numpy())
        roll[side] = np.deg2rad(data[f"roll_{side}"].to_numpy())
    dot = np.sin(pitch["est"]) * np.sin(pitch["ref"]) + np.cos(pitch["est"]) * np.cos(pitch["ref"]) * np.cos(
        roll["est"] - roll["ref"]
    )
    errors = dict(
        tilt_deg=np.rad2deg(np.arccos(np.clip(dot, -1, 1))),
        yaw_deg=(data.yaw_est.to_numpy() - data.yaw_ref.to_numpy() + 180) % 360 - 180,
    )
    for kind, unit in [("acc", "mps2"), ("gyro", "radps")]:
        errors[f"{kind}_bias_{unit}"] = np.linalg.norm(
            np.column_stack([data[f"{kind}_bias_est_{a}"] - data[f"{kind}_bias_{a}"] for a in "xyz"]), axis=1
        )
    epoch = np.rint((data.time.to_numpy() - start) / DT).astype(int)
    result = []
    for lo, hi in WINDOWS:
        keep = (epoch >= round(lo / DT)) & (epoch < round(hi / DT))
        if keep.sum() != round((hi - lo) / DT):
            raise ValueError(f"incomplete live window starting at {start}")
        row = {k: float(np.sqrt(np.mean(v[keep] ** 2))) for k, v in errors.items()}
        if not all(math.isfinite(v) for v in row.values()):
            raise ValueError("nonfinite live readout")
        result.append(dict(start_s=lo, stop_s=hi, samples=int(keep.sum()), **row))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True, help="parent of following/ and head/ datasets")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument(
        "--verify-controls",
        action="store_true",
        help="compare the two endpoint arms with the intact production policies",
    )
    parser.add_argument("--long", action="store_true", help="compare the three bias policies on the final 900 seconds")
    parser.add_argument("--seeds", type=int, nargs="+", default=[7, 41])
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    arms = {k: ARMS[k] for k in ("previous", "full", "gravity_axis_accel_bias")} if args.long else ARMS
    with tempfile.TemporaryDirectory(prefix="wave-handover-") as temp:
        work = Path(temp)
        binary = build(work)
        inputs, sources = {}, {}
        for case in CASES:
            files = [args.data.resolve() / case / n for n in REFERENCE_CSV_SHA256_BY_CASE[case]]
            sources[case] = input_provenance(files, case)
            for file in files:
                path = work / case / file.name
                path.parent.mkdir(exist_ok=True)
                if args.long:
                    path.symlink_to(file)
                else:
                    with file.open("rb") as src, path.open("wb") as dst:
                        dst.writelines(islice(src, round(SECONDS / DT) + 2))
                inputs[(case, file.name)] = path

        def run(job, executable=binary):
            (case, name), seed, (arm, mask) = job
            with tempfile.TemporaryDirectory(dir=work) as tmp:
                env = {k: v for k, v in os.environ.items() if not k.startswith(("W3D_", "SF_", "OU_", "TFG_"))}
                env.update(
                    OCEAN_IMU_SIM_CASE=case,
                    W3D_SEED=str(seed),
                    W3D_HANDOVER_MASK=str(mask),
                    W3D_WRITE_TIMESERIES="0" if args.long else "1",
                    W3D_COLLECT_ALL_GATES="1",
                )
                if args.long:
                    env["W3D_VALIDATION_WINDOW_SEC"] = "900"
                local = Path(tmp) / name
                local.symlink_to(inputs[(case, name)])
                r = subprocess.run([str(executable), "--input", name], cwd=tmp, env=env, capture_output=True, text=True)
                log = r.stdout + r.stderr
                if r.returncode not in (0, 1) or not re.search(r"^QUALITY_GATE:", log, re.M):
                    raise ValueError(log[-1600:])
                live = float(re.search(r"STARTUP live_s=([.0-9]+)", log).group(1))
                if args.long:
                    (line,) = [l for l in log.splitlines() if l.startswith("VALIDATION_METRICS ")]
                    metrics = dict(token.split("=", 1) for token in line.split()[1:])
                    if metrics["samples"] != "180000" or metrics["start_s"] != "300":
                        raise ValueError("incomplete steady-state replay")
                    keys = dict(
                        yaw_deg="yaw_rms_deg",
                        roll_deg="roll_rms_deg",
                        pitch_deg="pitch_rms_deg",
                        acc_bias_mps2="accel_bias_3d_rms_mps2",
                        gyro_bias_radps="gyro_bias_3d_rms_radps",
                        y_m="disp_y_rms_m",
                        z_m="disp_z_rms_m",
                    )
                    values = {column: float(metrics[key]) for column, key in keys.items()}
                    if not all(math.isfinite(v) for v in values.values()):
                        raise ValueError("nonfinite steady-state score")
                    rows = [dict(start_s=300, stop_s=1200, samples=180000, **values)]
                else:
                    (output,) = [p for p in Path(tmp).glob("*.csv") if p != local]
                    rows = score(output, live)
                return [
                    dict(
                        case=case,
                        input=name,
                        seed=seed,
                        arm=arm,
                        live_s=live,
                        qualified=f"HANDOVER_COMPONENTS {mask}\n" in log,
                        **row,
                    )
                    for row in rows
                ]

        rows = []
        with ThreadPoolExecutor(args.jobs) as pool:
            for part in pool.map(run, product(inputs, args.seeds, arms.items())):
                rows.extend(part)
        # A missing joint fit is kept and contributes the fallback score to every arm.
        for k, group in pd.DataFrame(rows).groupby(["case", "input", "seed", "start_s"]):
            if set(group.arm) != set(arms) or group.live_s.nunique() != 1 or group.qualified.nunique() != 1:
                raise ValueError(f"unpaired case {k}")
        if not any(row["qualified"] for row in rows):
            raise ValueError("all starts used the fallback; no handover intervention was exercised")
        control_checks = []
        control_binary_sha256 = None
        if args.verify_controls:
            control = build(work, control=True)
            control_binary_sha256 = digest(control)
            for (case, name), (arm, mask) in product(inputs, [("previous", 0), ("full", 15)]):
                if "jonswap_H4.000" not in name:
                    continue
                actual = run(((case, name), args.seeds[0], (arm, mask)), executable=control)
                expected = [
                    r for r in rows if (r["case"], r["input"], r["arm"], r["seed"]) == (case, name, arm, args.seeds[0])
                ]
                if len(expected) != len(actual):
                    raise ValueError("missing production control")
                for observed, baseline in zip(actual, expected):
                    metrics = [k for k in observed if k.endswith(("_deg", "_mps2", "_radps", "_m"))]
                    if observed["live_s"] != baseline["live_s"] or any(
                        not math.isclose(observed[k], baseline[k], abs_tol=1e-8, rel_tol=1e-6) for k in metrics
                    ):
                        raise ValueError(f"production control differs: {observed} {baseline}")
                    control_checks.append(dict(case=case, input=name, arm=arm, start_s=baseline["start_s"]))
        with (args.output / "records.csv").open("w") as out:
            w = csv.DictWriter(out, fieldnames=list(rows[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        bound_files = [
            ROOT / "src/kalman_common/WaveStartup.h",
            ROOT / "src/kalman_common/ProxyStartupFusion.h",
            ROOT / "src/kalman_ou_iii/Kalman3D_Wave_OU_III.h",
            ROOT / "tests/kalman_ou_iii/kalman_ou_iii-sim.cpp",
            ROOT / "src/util/W3dSimCommon.cpp",
        ]
        receipt = dict(
            source_sha256={str(p.relative_to(ROOT)): digest(p) for p in bound_files},
            runtime_source=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            script_sha256=digest(Path(__file__)),
            compiler=subprocess.check_output(["g++", "--version"], text=True).splitlines()[0],
            control_binary_sha256=control_binary_sha256,
            control_checks=control_checks,
            records_sha256=digest(args.output / "records.csv"),
            binary_sha256=digest(binary),
            dataset=sources,
            seeds=args.seeds,
            arms=arms,
            windows=[(300, 1200)] if args.long else WINDOWS,
            gravity_axis_accel_bias="experimental only",
        )
        (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
