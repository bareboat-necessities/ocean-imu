#!/usr/bin/env python3
"""Observe OU-II/III bias and literal prediction angles without changing inputs.

Builds temporary, read-only source taps. The regular simulator/test stdout is
retained separately from the JSON diagnostic. No production telemetry or API
is added. Use the same command on the baseline and candidate source trees.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import io
import tarfile
import sys
import math
import threading

REPO = Path(__file__).resolve().parents[1]
OBSERVER = r'''#pragma once
#include <algorithm>
#include <cmath>
#include <cstdio>
namespace ou_gyro_audit {
struct Stats {
    int family;
    unsigned long long predictions=0, corrections=0, nonfinite=0, projected=0;
    double max_bias=0, max_angle=0, max_dt=0, max_projection=0;
    unsigned long long over[4]={0,0,0,0};
    explicit Stats(int f): family(f) {}
    template<class V> void bias(const V& v) {
        const double n=std::hypot(double(v[0]),double(v[1]),double(v[2]));
        if (!std::isfinite(n)) { ++nonfinite; return; }
        max_bias=std::max(max_bias,n);
        const double radii[4]={0.1,0.2,0.5,1.0};
        for (int i=0;i<4;++i) if(n>radii[i]) ++over[i];
    }
    template<class V, class T> void prediction(const V& w,T dt) {
        ++predictions;
        max_dt=std::max(max_dt,double(dt));
        const double angle=double(dt)*std::hypot(double(w[0]),double(w[1]),double(w[2]));
        if (!std::isfinite(angle)) { ++nonfinite; return; }
        max_angle=std::max(max_angle,angle);
    }
    template<class V> void projection(const V& before,const V& after) {
        const double d=std::hypot(double(before[0])-double(after[0]),
            double(before[1])-double(after[1]),double(before[2])-double(after[2]));
        if (d>0) { ++projected; max_projection=std::max(max_projection,d); }
    }
    ~Stats() {
        if (!predictions && !corrections) return;
        std::fprintf(stderr,"OU_GYRO_AUDIT {\"family\":%d,\"predictions\":%llu,"
          "\"corrections\":%llu,\"nonfinite\":%llu,\"max_bias_rad_s\":%.17g,"
          "\"max_prediction_angle_rad\":%.17g,\"max_dt_s\":%.17g,"
          "\"projection_count\":%llu,\"max_projection_rad_s\":%.17g,"
          "\"candidate_exceedances\":[%llu,%llu,%llu,%llu]}\n",family,
          predictions,corrections,nonfinite,max_bias,max_angle,max_dt,projected,
          max_projection,over[0],over[1],over[2],over[3]);
    }
};
inline Stats ou2(2),ou3(3);
}
'''


def instrument(source: str, family: int) -> str:
    def once(old, new):
        nonlocal source
        if source.count(old) != 1:
            raise ValueError(f"observer anchor changed: {old}")
        source = source.replace(old, new)
    observer = f"ou_gyro_audit::ou{family}"
    once("#pragma once", '#pragma once\n#include "ou_gyro_audit_observer.h"')
    once("::applyQuaternionCorrectionFromErrorState()\n{",
         "::applyQuaternionCorrectionFromErrorState()\n{\n"
         f"    ++{observer}.corrections; {observer}.bias(gyroscope_bias());")
    once("    last_gyr_bias_corrected = gyr - gyro_bias;",
         "    last_gyr_bias_corrected = gyr - gyro_bias;\n"
         f"    {observer}.bias(gyro_bias);\n"
         f"    {observer}.prediction(last_gyr_bias_corrected, Ts);")
    projection = "        ocean_imu::kalman::ou_detail::project_gyro_bias<T>(bias);"
    if projection in source:
        once(projection, "        const Vector3 audit_before = bias;\n"+projection+
             f"\n        {observer}.projection(audit_before, Vector3(bias));")
    return source


def compare_metrics(before, after, key="run_id"):
    """Pair actual emitted metrics, retaining maxima by metric and family.

    This is descriptive evidence. Existing simulator/validation gates decide
    performance acceptance; the diagnostic never changes their thresholds.
    """
    left = {row[key]: row for row in before}
    right = {row[key]: row for row in after}
    if len(left) != len(before) or len(right) != len(after):
        raise ValueError("duplicate replay identity")
    if left.keys() != right.keys():
        raise ValueError("before/after replay identities differ")
    maxima = {}
    for identity, a in left.items():
        b = right[identity]
        for name in a.keys() & b.keys():
            try:
                x, y = float(a[name]), float(b[name])
            except (TypeError, ValueError):
                continue
            if not (math.isfinite(x) and math.isfinite(y)):
                if x == y or (math.isnan(x) and math.isnan(y)):
                    continue
                raise ValueError(f"nonfinite metric difference: {identity} {name}")
            family = a.get("family", "unknown")
            row = maxima.setdefault(family, {})
            row[name] = max(row.get(name, 0), abs(x-y))
    return {"paired_runs": len(left), "max_absolute_metric_differences": maxima}


def deterministic_metrics(path):
    return [dict(token.split("=", 1) for token in line.split()[1:] if "=" in token)
            for line in Path(path).read_text().splitlines() if line.startswith("VALIDATION_METRICS ")]


def summarize_observations(records):
    result = {}
    for record in records:
        for row in record["observations"]:
            summary = result.setdefault(str(row["family"]), {"runs": 0, "predictions": 0,
                "corrections": 0, "nonfinite": 0, "projection_count": 0,
                "max_bias_rad_s": 0, "max_prediction_angle_rad": 0,
                "max_projection_rad_s": 0, "max_dt_s": 0})
            summary["runs"] += 1
            for name in ("predictions", "corrections", "nonfinite", "projection_count"):
                summary[name] += row[name]
            for name in ("max_bias_rad_s", "max_prediction_angle_rad", "max_projection_rad_s", "max_dt_s"):
                summary[name] = max(summary[name], row[name])
    return result


def compare_audit_directories(before_dir, after_dir, part="all"):
    """Reproduce paired evidence from retained audit files, without replays."""
    directories = [Path(before_dir), Path(after_dir)]
    audits = [json.loads((directory/"audit.json").read_text()) for directory in directories]
    if audits[0]["observer_sha256"] != audits[1]["observer_sha256"]:
        raise ValueError("before/after observers differ")
    if audits[0]["candidate_radii_rad_s"] != audits[1]["candidate_radii_rad_s"]:
        raise ValueError("before/after candidate radii differ")
    result = {"candidate_radii_rad_s": audits[0]["candidate_radii_rad_s"],
              "observer_sha256": audits[0]["observer_sha256"],
              "source_sha256": [audit["source_sha256"] for audit in audits]}
    targets = [{row["target"]: row for row in audit["runs"]} if part != "replay" else {}
               for audit in audits]
    if targets[0].keys() != targets[1].keys():
        raise ValueError("before/after standalone targets differ")
    if targets[0]:
        metrics = []
        identical = {}
        for name in targets[0]:
            rows = [target[name] for target in targets]
            if any(row["exit_code"] for row in rows):
                raise ValueError(f"failed standalone target: {name}")
            if rows[0]["test_source_sha256"] != rows[1]["test_source_sha256"]:
                raise ValueError(f"standalone test inputs differ: {name}")
            outputs = [directory/(name+".stdout") for directory in directories]
            for path, row in zip(outputs, rows):
                if hashlib.sha256(path.read_bytes()).hexdigest() != row["stdout_sha256"]:
                    raise ValueError(f"standalone stdout changed: {path}")
            identical[name] = rows[0]["stdout_sha256"] == rows[1]["stdout_sha256"]
            metrics.append({"target": name, **compare_metrics(
                *[deterministic_metrics(path) for path in outputs], key="input")})
        result["standalone"] = {
            "baseline": summarize_observations(audits[0]["runs"]),
            "modified": summarize_observations(audits[1]["runs"]),
            "metrics": metrics, "stdout_identical": identical}
    if part != "standalone" and any("validation_replays" in audit for audit in audits):
        study = audits[0].get("replay_study", "validation")
        if study != audits[1].get("replay_study", "validation"):
            raise ValueError("before/after replay studies differ")
        records = [audit["validation_replays"] for audit in audits]
        if any(audit.get("validation_exit_code") != 0 for audit in audits):
            raise ValueError("incomplete replay study")
        if any(not row.get("input_immutable_verified") for rows in records for row in rows):
            raise ValueError("replay input immutability was not verified")
        def identities(rows):
            # Frozen tuning points are measured in each build and are retained
            # below. Pair the physical input, seeds and update mode exactly.
            return sorted(json.dumps([row["input_identity"], row["input_sha256"],
                row["input_rows"], {k:v for k,v in row["environment"].items()
                                    if not k.startswith("W3D_FIXED_")}], sort_keys=True)
                for row in rows)
        left, right = map(identities, records)
        if left != right:
            raise ValueError("paired physical inputs/settings differ")
        studies = [json.loads((directory/study/f"ou_{study}.json").read_text())
                   for directory in directories]
        if studies[0]["protocol"] != studies[1]["protocol"]:
            raise ValueError("before/after replay protocols differ")
        metrics = compare_metrics(*[data["raw_runs"] for data in studies])
        rows_by_id = [{row["run_id"]: row for row in data["raw_runs"]} for data in studies]
        for identity, row in rows_by_id[0].items():
            for name, value in row.items():
                if name in ("samples", "window_s", "start_s") or "disp_z_ref_rms_m" in name:
                    other = rows_by_id[1][identity][name]
                    if value != other and not (isinstance(value, float) and
                            isinstance(other, float) and math.isnan(value) and math.isnan(other)):
                        raise ValueError(f"reference motion differs: {identity} {name}")
        summaries = [summarize_observations(rows) for rows in records]
        for summary, rows in zip(summaries, records):
            for family, row in summary.items():
                row["candidate_exceedances"] = [sum(observation["candidate_exceedances"][i]
                    for record in rows for observation in record["observations"]
                    if str(observation["family"]) == family) for i in range(4)]
        result["replay"] = {"baseline": summaries[0], "modified": summaries[1],
            "metrics": metrics, "input_records_identical": True,
            "input_immutability_verified": True, "reference_motion_metrics_identical": True,
            "paired_input_records_sha256": hashlib.sha256(json.dumps(left).encode()).hexdigest(),
            "paired_replay_count": len(left), "protocol": studies[0]["protocol"],
            "frozen_tuning_points": [data.get("fixed_tuning_points", data.get("nominal_tuning_point"))
                                     for data in studies]}
    if not targets[0] and "replay" not in result:
        raise ValueError("no paired observations")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--family", choices=("ii", "iii"))
    p.add_argument("--compare", nargs=2, type=Path, metavar=("BASELINE_DIR", "CANDIDATE_DIR"),
                   help="Compare retained audits and write paired-audit.json; does not build or replay")
    p.add_argument("--compare-part", choices=("all", "standalone", "replay"), default="all")
    p.add_argument("--targets", nargs="+", default=["sim", "startup_init-test", "stationary_device-test"])
    p.add_argument("--eigen-dir", type=Path, default=REPO/"third_party/eigen")
    p.add_argument("--no-build", action="store_true")
    p.add_argument("--source-ref", help="Build src/ from this immutable git revision for a paired baseline")
    p.add_argument("--validation-mode", choices=("smoke", "full"), help="Also observe the existing paired validation replay protocol")
    p.add_argument("--validation-peer-dir", type=Path, help="Audit build directory for the other OU family; validation requires the pair")
    p.add_argument("--validation-jobs", type=int, default=2)
    p.add_argument("--validation-study", choices=("validation", "robustness"), default="validation")
    p.add_argument("--replay-only", action="store_true", help="Reuse built audit binaries and run only the requested replay study")
    args = p.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if args.compare:
        result = compare_audit_directories(*args.compare, part=args.compare_part)
        (out/"paired-audit.json").write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
        print(f"Wrote {out/'paired-audit.json'}", flush=True)
        return
    if args.family is None:
        p.error("--family is required unless --compare is used")
    src = out/"src"
    family = 2 if args.family == "ii" else 3
    folder = "kalman_ou_"+args.family
    if not args.no_build:
        if args.source_ref:
            archive = subprocess.check_output(["git", "archive", args.source_ref, "src"], cwd=REPO)
            with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
                tar.extractall(out, filter="data")
        else:
            shutil.copytree(REPO/"src", src, dirs_exist_ok=True)
        hashes = {str(path.relative_to(src)): hashlib.sha256(path.read_bytes()).hexdigest()
                  for path in sorted(src.rglob("*")) if path.is_file()}
        (out/"source-sha256.json").write_text(json.dumps(hashes, indent=2)+"\n")
        (src/"ou_gyro_audit_observer.h").write_text(OBSERVER)
        header = src/folder/("Kalman3D_Wave_OU_"+args.family.upper()+".h")
        header.write_text(instrument(header.read_text(), family))
    result = {"source_commit": subprocess.check_output(["git", "rev-parse", args.source_ref or "HEAD"], cwd=REPO, text=True).strip(),
              "source_sha256": json.loads((out/"source-sha256.json").read_text()),
              "observer_sha256": hashlib.sha256(OBSERVER.encode()).hexdigest(),
              "candidate_radii_rad_s": [0.1,0.2,0.5,1.0], "runs": []}
    for target in ([] if args.replay_only else args.targets):
        name = folder+"-sim" if target == "sim" else target
        binary = out/name
        if not args.no_build:
            cmd = ["g++", "-O3", "-std=c++20", "-funroll-loops", "-fno-finite-math-only", "-march=native",
                   "-I"+str(src), "-isystem", str(args.eigen_dir.resolve()),
                   str(REPO/"tests"/folder/(name+".cpp"))]
            if target == "sim" or (family == 3 and target in ("startup_init-test", "tuner_coupling-test", "imu_lever_arm-test")):
                cmd.append(str(src/"util/W3dSimCommon.cpp"))
            cmd += ["-o", str(binary)]
            with (out/(name+".build.log")).open("w") as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        env = dict(os.environ, W3D_WRITE_TIMESERIES="0", W3D_COLLECT_ALL_GATES="1", W3D_VALIDATION_WINDOW_SEC="900")
        with (out/(name+".stdout")).open("w") as stdout, (out/(name+".stderr")).open("w") as stderr:
            run = subprocess.run([str(binary)], cwd=REPO/"tests"/folder, env=env, stdout=stdout, stderr=stderr)
        records = [json.loads(line.removeprefix("OU_GYRO_AUDIT ")) for line in (out/(name+".stderr")).read_text().splitlines() if line.startswith("OU_GYRO_AUDIT ")]
        row = {"target": name, "exit_code": run.returncode, "observations": records,
               "test_source_sha256": hashlib.sha256((REPO/"tests"/folder/(name+".cpp")).read_bytes()).hexdigest(),
               "stdout_sha256": hashlib.sha256((out/(name+".stdout")).read_bytes()).hexdigest()}
        result["runs"].append(row)
        (out/"audit.json").write_text(json.dumps(result, indent=2)+"\n")
        print(json.dumps(row), flush=True)
        if run.returncode:
            raise SystemExit(run.returncode)
    if args.validation_mode:
        if args.validation_peer_dir is None:
            raise ValueError("paired validation requires --validation-peer-dir")
        sys.path.insert(0, str(REPO/"tools"))
        import ou_validation as validation
        key = "OU_"+args.family.upper()
        validation.FAMILY_BINARY[key] = out/(folder+"-sim")
        peer = "iii" if args.family == "ii" else "ii"
        peer_dir = args.validation_peer_dir.resolve()
        validation.FAMILY_BINARY["OU_"+peer.upper()] = peer_dir/("kalman_ou_"+peer+"-sim")
        result["validation_peer_source_sha256"] = json.loads((peer_dir/"source-sha256.json").read_text())
        binaries = {str(value) for value in validation.FAMILY_BINARY.values()}
        fingerprints = {name: hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in binaries}
        cache = out/"replay-cache"
        cache.mkdir(exist_ok=True)
        lock = threading.Lock()
        original_run = validation.subprocess.run
        original_write = validation.write_wave_csv
        inputs = {}
        released_inputs = {path.resolve() for path in (REPO/"tests"/folder).glob("wave_data_*.csv")}
        replays = []
        def checked_write(path, columns, data):
            original_write(path, columns, data)
            payload = path.read_bytes()
            rows = payload.count(b"\n")-1
            if rows != len(data):
                raise RuntimeError(f"incomplete generated input: {path}: {rows} != {len(data)} rows")
            with lock:
                inputs[str(path.resolve())] = (hashlib.sha256(payload).hexdigest(), rows)

        def observed_run(command, *a, **kw):
            is_replay = command and command[0] in binaries
            cache_file = None
            if is_replay:
                input_path = Path(command[-1])
                payload = input_path.read_bytes()
                input_hash = hashlib.sha256(payload).hexdigest()
                resolved = input_path.resolve()
                if str(resolved) not in inputs and resolved in released_inputs:
                    # Robustness calibrates on the pinned release record itself.
                    inputs[str(resolved)] = (input_hash, payload.count(b"\n")-1)
                expected_hash, input_rows = inputs[str(resolved)]
                if input_hash != expected_hash:
                    raise RuntimeError(f"generated input changed before replay: {input_path}")
                parts = input_path.parts
                temporary_index = next((i for i, part in enumerate(parts) if part.startswith("ocean-imu-ou-")), None)
                input_identity = ("dataset/"+input_path.name if temporary_index is None
                                  else "/".join(parts[temporary_index+1:]))
                environment = {k:v for k,v in kw["env"].items() if k.startswith(("W3D_", "SF_"))}
                # Only reuse complete scalar-output replays from these exact
                # binary bytes, input bytes, seeds and settings. A diagnostic
                # needing time-series files must execute again.
                if environment.get("W3D_WRITE_TIMESERIES") == "0":
                    key_data = ["immutable-input-v2", fingerprints[command[0]], input_hash,
                                input_path.name, environment, command[1:-1]]
                    key_hash = hashlib.sha256(json.dumps(key_data, sort_keys=True).encode()).hexdigest()
                    cache_file = cache/(key_hash+".json")
            if cache_file is not None and cache_file.exists():
                saved = json.loads(cache_file.read_text())
                completed = subprocess.CompletedProcess(command, saved["returncode"], saved["stdout"], saved["stderr"])
            else:
                completed = original_run(command, *a, **kw)
                if is_replay and hashlib.sha256(input_path.read_bytes()).hexdigest() != expected_hash:
                    raise RuntimeError(f"generated input changed during replay: {input_path}")
                if cache_file is not None and "OU_GYRO_AUDIT " in completed.stdout:
                    with lock:
                        temporary = cache_file.with_suffix(".tmp")
                        temporary.write_text(json.dumps({"returncode": completed.returncode,
                            "stdout": completed.stdout, "stderr": completed.stderr})+"\n")
                        temporary.replace(cache_file)
            if is_replay:
                observations = [json.loads(line.removeprefix("OU_GYRO_AUDIT "))
                    for line in completed.stdout.splitlines() if line.startswith("OU_GYRO_AUDIT ")]
                record = {"observations": observations, "exit_code": completed.returncode,
                    "input_sha256": input_hash, "input_identity": input_identity,
                    "input_rows": input_rows, "input_immutable_verified": True,
                    "binary_sha256": fingerprints[command[0]], "environment": environment}
                with lock:
                    replays.append(record)
                    with (out/"replay-progress.jsonl").open("a") as log:
                        log.write(json.dumps(record)+"\n")
            return completed
        validation.subprocess.run = observed_run
        validation.write_wave_csv = checked_write
        try:
            study = validation
            study_args = ["--families", "both"]
            if args.validation_study == "robustness":
                import ou_robustness as study
                study_args = []
            code = study.main(["--mode", args.validation_mode, *study_args,
                "--data-dir", str(REPO/"tests"/folder), "--output-dir", str(out/args.validation_study),
                "--jobs", str(args.validation_jobs), "--skip-build", "--no-plots"])
        finally:
            validation.subprocess.run = original_run
            validation.write_wave_csv = original_write
            result["validation_replays"] = replays
            (out/"audit.json").write_text(json.dumps(result, indent=2)+"\n")
        result["validation_replays"] = replays
        result["replay_study"] = args.validation_study
        result["validation_exit_code"] = code
        (out/"audit.json").write_text(json.dumps(result, indent=2)+"\n")
        if code:
            raise SystemExit(code)


if __name__ == "__main__":
    main()
