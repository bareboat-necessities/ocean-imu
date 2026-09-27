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


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--family", choices=("ii", "iii"), required=True)
    p.add_argument("--targets", nargs="+", default=["sim", "startup_init-test", "stationary_device-test"])
    p.add_argument("--eigen-dir", type=Path, default=REPO/"third_party/eigen")
    p.add_argument("--no-build", action="store_true")
    p.add_argument("--source-ref", help="Build src/ from this immutable git revision for a paired baseline")
    p.add_argument("--validation-mode", choices=("smoke", "full"), help="Also observe the existing paired validation replay protocol")
    p.add_argument("--validation-peer-dir", type=Path, help="Audit build directory for the other OU family; validation requires the pair")
    p.add_argument("--validation-jobs", type=int, default=2)
    args = p.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
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
    for target in args.targets:
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
        original_run = validation.subprocess.run
        replays = []
        def observed_run(command, *a, **kw):
            completed = original_run(command, *a, **kw)
            if command and command[0] in binaries:
                observations = [json.loads(line.removeprefix("OU_GYRO_AUDIT "))
                    for line in completed.stdout.splitlines() if line.startswith("OU_GYRO_AUDIT ")]
                replays.append({"observations": observations, "exit_code": completed.returncode,
                    "input_sha256": hashlib.sha256(Path(command[-1]).read_bytes()).hexdigest(),
                    "environment": {k:v for k,v in kw["env"].items() if k.startswith("W3D_")}})
            return completed
        validation.subprocess.run = observed_run
        try:
            code = validation.main(["--mode", args.validation_mode, "--families", "both",
                "--data-dir", str(REPO/"tests"/folder), "--output-dir", str(out/"validation"),
                "--jobs", str(args.validation_jobs), "--skip-build", "--no-plots"])
        finally:
            validation.subprocess.run = original_run
            result["validation_replays"] = replays
            (out/"audit.json").write_text(json.dumps(result, indent=2)+"\n")
        result["validation_replays"] = replays
        result["validation_exit_code"] = code
        (out/"audit.json").write_text(json.dumps(result, indent=2)+"\n")
        if code:
            raise SystemExit(code)


if __name__ == "__main__":
    main()
