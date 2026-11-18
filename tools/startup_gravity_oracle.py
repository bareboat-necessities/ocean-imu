#!/usr/bin/env python3
"""One-time exact-tilt/field interventions in the normal OU-III and TFG wave replays."""

import argparse
import csv
import gzip
import hashlib
import io
import itertools
import json
import math
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import tempfile

from sim_dataset import ROOT, REFERENCE_CSV_SHA256_BY_CASE, digest, input_provenance
from startup_source_snapshot import SourceSnapshot

POLICIES = ("previous", "deployed", "previous_exact_tilt", "exact_tilt", "exact_field", "exact_both")
MEASURED_WAVE = "measured_wave"
FILTERS = {
    "ou3": (
        "kalman_ou_iii",
        "kalman_common/ProxyStartupFusion.h",
        "q_seed",
        "use_learned",
        "        if (!q_seed.coeffs().allFinite()) return;",
    ),
    "tfg": (
        "kalman_tfg",
        "kalman_tfg/SeaStateFusionFilter_TFG.h",
        "q_bw",
        "wave_startup_used_",
        "        mekf_.initialize_from_truth(q_bw.normalized(),",
    ),
}
WINDOWS = {
    "early": (40, 60),
    "pre_refine": (60, 90),
    "refining": (90, 120),
    "post_refine": (120, 180),
    "late": (300, 1200),
}
METRICS = (
    "roll_rms_deg",
    "pitch_rms_deg",
    "yaw_rms_deg",
    "disp_x_rms_m",
    "disp_y_rms_m",
    "disp_z_rms_m",
    "disp_3d_rms_m",
    "accel_bias_3d_rms_mps2",
    "gyro_bias_3d_rms_radps",
)


def insert(source, anchor, replacement, count=1):
    if source.count(anchor) != count:
        raise ValueError(f"shipping source anchor changed: {anchor}")
    return source.replace(anchor, replacement)


def trace_updates(src, pointer, quat):
    for stage, call in (
        ("prediction", "time_update(gyro, dt);"),
        ("accelerometer", "measurement_update_acc_only(acc_in, tempC);"),
    ):
        anchor = pointer + call
        hook = f'if (startup_oracle::trace_sample) startup_oracle::report_stage("{stage}", {pointer}{quat}());'
        src = insert(src, anchor, anchor + "\n            " + hook)
    return src


def seed_measured_wave(source, family):
    if family == "ou3":
        anchor = "        wave_startup_used_ = use_learned;"
        q, gate, gyro, acc, core = "q_seed", "use_learned", "last_gyro_body_ned_", "last_acc_body_ned_", "impl_.mekf()"
        assign = "impl_.mekf().seed_oracle_wave(wave);"
    else:
        anchor = "        seedHandoffAttitudeCovariance_();"
        q, gate, gyro, acc, core = "q_bw", "wave_startup_used_", "last_gyro_body_", "last_acc_body_", "mekf_"
        assign = "mekf_.state().X.col(3) = wave;"
    block = f"""
        if ({gate} && startup_oracle::policy == "measured_wave") {{
            const Eigen::Vector3f down = {q}.conjugate() * Eigen::Vector3f::UnitZ();
            const Eigen::Vector3f turn = turn_separator_.correctionBodyMs2(down, {gyro} - learned.gyro_bias);
            const Eigen::Vector3f wave = {q} * ({acc} - {core}.get_acc_bias() - turn)
                                        + cfg_.gravity_magnitude * Eigen::Vector3f::UnitZ();
            {assign}
        }}"""
    replacement = anchor + block if family == "ou3" else block + "\n" + anchor
    return insert(source, anchor, replacement)


def build(work, families):
    fixture = ROOT / "tests/startup_learning/GravityHandoverOracle.h"
    eigen = ROOT / "third_party/eigen"
    if not (eigen / "Eigen/Dense").exists():
        eigen = Path("/usr/include/eigen3")
    cxx = os.environ.get("CXX", "g++")
    flags = [
        "-std=c++20",
        "-O3",
        "-funroll-loops",
        "-fno-finite-math-only",
        "-march=x86-64-v3",
        "-I" + str(work),
        "-I" + str(fixture.parent),
        "-I" + str(ROOT / "src"),
        "-isystem",
        str(eigen),
    ]
    sources = [fixture, Path(__file__).resolve(), ROOT / "src/util/W3dSimCommon.cpp"]
    common = sources[-1].read_text()
    anchor = "        fusion_adapter_.update(options_.dt,"
    common = insert(common, anchor, "        startup_oracle::observe(rec.time, C_wb_zu);\n" + anchor, 2)
    common_path = work / "W3dSimCommon.cpp"
    common_path.write_text('#include "GravityHandoverOracle.h"\n' + common)
    obj = work / "W3dSimCommon.o"
    subprocess.run([cxx, *flags, "-c", str(common_path), "-o", str(obj)], check=True)
    if "ou3" in families:
        impl = Path("kalman_ou_iii/SeaStateFusionFilter_OU_III.h")
        sources.append(ROOT / "src" / impl)
        core = Path("kalman_ou_iii/Kalman3D_Wave_OU_III.h")
        sources.append(ROOT / "src" / core)
        core_overlay = work / core
        core_overlay.parent.mkdir(parents=True, exist_ok=True)
        accessor = "    [[nodiscard]] Vector3 get_world_accel() const"
        core_overlay.write_text(
            insert(
                sources[-1].read_text(),
                accessor,
                "    void seed_oracle_wave(const Vector3& wave) { xext.template segment<3>(OFF_AW) = wave; }\n"
                + accessor,
            )
        )
        over = work / impl
        over.parent.mkdir(parents=True, exist_ok=True)
        over.write_text(
            '#include "GravityHandoverOracle.h"\n'
            + trace_updates((ROOT / "src" / impl).read_text(), "mekf_->", "quaternion_boat")
        )
    binaries = {}
    for family in families:
        folder, wrapper, q, learned, anchor = FILTERS[family]
        path = ROOT / "src" / wrapper
        src = seed_measured_wave(path.read_text(), family)
        src = insert(
            src, "bool wave_startup_enabled = true;", "bool wave_startup_enabled = !startup_oracle::previous();"
        )
        src = insert(
            src,
            anchor,
            f"        startup_oracle::handover({q}, mag_world_ref_uT_, {learned}, "
            "[&](const Eigen::Vector3f& ref) { setMagWorldRef_(ref); });\n" + anchor,
        )
        if family == "tfg":
            src = trace_updates(src, "mekf_.", "quaternion")
        overlay = work / wrapper
        overlay.parent.mkdir(parents=True, exist_ok=True)
        overlay.write_text('#include "GravityHandoverOracle.h"\n' + src)
        sim = ROOT / "tests" / folder / (folder + "-sim.cpp")
        sources += [path, sim]
        binaries[family] = work / family
        subprocess.run([cxx, *flags, str(sim), str(obj), "-o", str(binaries[family])], check=True)
    return binaries, {
        "sources": {str(p.relative_to(ROOT)): digest(p) for p in sources},
        "compiler": subprocess.check_output([cxx, "--version"], text=True).splitlines()[0],
        "binaries": {f: digest(p) for f, p in binaries.items()},
    }


def parse(line):
    return dict(x.split("=", 1) for x in line.split()[1:])


def replay(job):
    (family, binary), (case, path), seed, policy = job
    env = {k: v for k, v in os.environ.items() if not k.startswith(("W3D_", "OU_", "TFG_", "SF_"))}
    env.update(
        OCEAN_IMU_SIM_CASE=case,
        STARTUP_ORACLE_POLICY=policy,
        W3D_WRITE_TIMESERIES="0",
        W3D_COLLECT_ALL_GATES="1",
        W3D_VALIDATION_WINDOW_SEC="900",
        W3D_VALIDATION_SEGMENTS=",".join(f"{k}:{a}:{b}" for k, (a, b) in WINDOWS.items()),
    )
    if seed != "default":
        env["W3D_SEED"] = str(int(seed))
    with tempfile.TemporaryDirectory() as cwd:
        proc = subprocess.run([binary, "--input", path], text=True, capture_output=True, env=env, cwd=cwd)
    if proc.returncode not in (0, 1) or "QUALITY_GATE:" not in proc.stdout:
        raise RuntimeError(proc.stderr + proc.stdout[-2000:])
    lines = proc.stdout.splitlines()
    (event,) = [parse(l) for l in lines if l.startswith("GRAVITY_ORACLE ")]
    if not all(math.isfinite(float(x)) for x in event.values()) or abs(float(event["yaw_change_rad"])) > 2e-6:
        raise ValueError(f"invalid handover: {event}")
    if "exact" in policy and ("tilt" in policy or policy == "exact_both") and float(event["installed_tilt_deg"]) > 1e-4:
        raise ValueError(f"wrong tilt intervention: {event}")
    if policy in ("exact_field", "exact_both") and float(event["installed_field_error_uT"]) > 1e-5:
        raise ValueError(f"wrong field intervention: {event}")
    key = dict(family=family, case=case, record=path.name, seed=seed, policy=policy)
    rows = []
    for line in lines:
        if not line.startswith("VALIDATION_SEGMENT "):
            continue
        v = parse(line)
        a, b = WINDOWS[v["segment"]]
        if int(v["samples"]) != (b - a) * 200 or float(v["start_s"]) != a:
            raise ValueError(f"incomplete window: {line}")
        metrics = {k: float(v[k]) for k in METRICS}
        if not all(math.isfinite(x) for x in metrics.values()):
            raise ValueError(f"nonfinite score: {line}")
        rows.append(dict(**key, window=v["segment"], samples=v["samples"], **metrics))
    if [r["window"] for r in rows] != list(WINDOWS):
        raise ValueError("missing, duplicate or out-of-order segment")
    stages = [parse(line) for line in lines if line.startswith("GRAVITY_UPDATE ")]
    expected = list(itertools.product((1, 10, 200, 1000, 2000), ("prediction", "accelerometer")))
    if [(int(v["sample"]), v["stage"]) for v in stages] != expected or not all(
        math.isfinite(float(v["tilt_deg"])) for v in stages
    ):
        raise ValueError("incomplete update trace")
    return dict(**key, **event), rows, [dict(**key, **v) for v in stages]


def write_csv(path, rows):
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    path.write_bytes(gzip.compress(stream.getvalue().encode(), mtime=0))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", type=Path, required=True, help="directory with following/ and head/")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument("--seeds", nargs="+", default=["default"])
    p.add_argument("--families", nargs="+", choices=FILTERS, default=list(FILTERS))
    p.add_argument("--policies", nargs="+", choices=(*POLICIES, MEASURED_WAVE), default=list(POLICIES))
    args = p.parse_args()
    paths = [p for p in (ROOT / "src").rglob("*") if p.is_file()]
    paths.extend(p for p in (ROOT / "tests/common").iterdir() if p.suffix in {".h", ".inc"})
    paths.extend(ROOT / "tests" / folder / (folder + "-sim.cpp") for folder, *_ in FILTERS.values())
    paths.extend(
        ROOT / name
        for name in (
            "tests/startup_learning/GravityHandoverOracle.h",
            "tools/sim_dataset.py",
            "tools/sim_dataset_manifest.json",
            "tools/startup_source_snapshot.py",
        )
    )
    snapshot = SourceSnapshot(ROOT, paths + [Path(__file__).resolve()])
    files = [
        (case, args.data.resolve() / case / n) for case, names in REFERENCE_CSV_SHA256_BY_CASE.items() for n in names
    ]
    inputs = {case: input_provenance([f for c, f in files if c == case], case) for case in REFERENCE_CSV_SHA256_BY_CASE}
    args.output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="gravity-oracle-") as tmp:
        binaries, receipt = build(Path(tmp), args.families)
        snapshot.verify()
        jobs = list(itertools.product(binaries.items(), files, args.seeds, args.policies))
        events, windows, updates = [], [], []
        with ThreadPoolExecutor(args.jobs) as pool:
            for i, (event, rows, stages) in enumerate(pool.map(replay, jobs), 1):
                events.append(event)
                windows.extend(rows)
                updates.extend(stages)
                if i % len(args.policies) == 0:
                    print(
                        f"{i}/{len(jobs)} {event['family']} {event['case']} {event['record']} {event['seed']}",
                        flush=True,
                    )
        groups = itertools.groupby(events, key=lambda e: (e["family"], e["case"], e["record"], e["seed"]))
        for key, group in groups:
            if len({r["epoch_s"] for r in group}) != 1:
                raise ValueError(f"unmatched handover times: {key}")
    snapshot.verify()
    outputs = {"handover.csv.gz": events, "windows.csv.gz": windows, "updates.csv.gz": updates}
    for name, data in outputs.items():
        write_csv(args.output / name, data)
    hashes = {name: digest(args.output / name) for name in outputs}
    receipt.update(
        source_commit=snapshot.revision,
        record_limits_sha256=snapshot.hashes["tests/common/RecordReplayLimits.inc"],
        source_inventory_sha256=hashlib.sha256(json.dumps(snapshot.hashes, sort_keys=True).encode()).hexdigest(),
        inputs=inputs,
        output_sha256=hashes,
        runs=len(events),
        seeds=args.seeds,
        policies=args.policies,
        windows=WINDOWS,
    )
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
