"""Diagnose whether live covariance, bias release or magnetic refinement loses the startup gain."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile
from startup_handover_components import build, substitute
from sim_dataset import ROOT, REFERENCE_CSV_SHA256_BY_CASE, input_provenance, digest
from startup_source_snapshot import startup_snapshot

POLICIES = {
    "reference": {
        "deployed": {},
        "keep_seed": {"SF_MAG_REFINE": "0"},
        "keep_seed_no_cal": {"SF_MAG_REFINE": "0", "SF_MAG_CONT_HI": "0"},
        "refine_no_cal": {"SF_MAG_CONT_HI": "0"},
        "refine_60s": {"SF_MAG_REFINE_WINDOW_SEC": "60"},
    },
    "refinement": {
        "deployed": {"W3D_REFINE_JOINT": "0"},
        "joint_field": {"W3D_REFINE_JOINT": "1"},
        "joint_field_yaw": {"W3D_REFINE_JOINT": "2"},
    },
    "retention": {
        "deployed": {},
        "tilt_1deg": {"SF_PROXY_TILT_SIGMA": "0.01745329252"},
        "early_bias": {"W3D_EARLY_BIAS": "1"},
        "tilt_half_deg": {"SF_PROXY_TILT_SIGMA": "0.00872664626"},
        "tilt_half_deg_early_bias": {"SF_PROXY_TILT_SIGMA": "0.00872664626", "W3D_EARLY_BIAS": "1"},
        "tilt_1deg_early_bias": {"SF_PROXY_TILT_SIGMA": "0.01745329252", "W3D_EARLY_BIAS": "1"},
    },
    "bias_prior": {
        "deployed": {},
        "previous_prior": {"OU_ACC_BIAS_INIT_STD": "0.045295", "W3D_HANDOVER_MASK": "0"},
        "bias_prior": {"OU_ACC_BIAS_INIT_STD": "0.045295"},
        "bias_prior_early": {"OU_ACC_BIAS_INIT_STD": "0.045295", "W3D_EARLY_BIAS": "1"},
        "bias_prior_tilt": {"OU_ACC_BIAS_INIT_STD": "0.045295", "SF_PROXY_TILT_SIGMA": "0.00872664626"},
        "bias_prior_tilt_early": {
            "OU_ACC_BIAS_INIT_STD": "0.045295",
            "SF_PROXY_TILT_SIGMA": "0.00872664626",
            "W3D_EARLY_BIAS": "1",
        },
    },
    "timeline": {"previous": {"W3D_HANDOVER_MASK": "0"}, "deployed": {}},
    "frontend": {
        "deployed": {"W3D_PROXY_SEED": "0"},
        "proxy_attitude": {"W3D_PROXY_SEED": "1"},
        "proxy_gyro_bias": {"W3D_PROXY_SEED": "2"},
        "proxy_both": {"W3D_PROXY_SEED": "3"},
    },
}
SEGMENTS = "startup:40:60,before_refine:60:90,refining:90:120,after_refine:120:180,late:300:1200"


def joint_refinement(wrapper):
    old = r"""        if (stage_ != Stage::Live) {
            if (useWaveStartup_()) wave_startup_.imu(dt, gyro_body_ned, acc_body_ned, cfg_.gyro_sample_mode);"""
    new = (
        r"""        if (stage_ == Stage::Live && cfg_.mag_refine_enabled && !mag_refine_done_ &&
            t_ >= cfg_.mag_refine_start_sec && useWaveStartup_()) {
            if (!joint_refine_started_) { wave_startup_.reset(); joint_refine_started_ = true; }
            wave_startup_.imu(dt, gyro_body_ned, acc_body_ned, cfg_.gyro_sample_mode);
        }
"""
        + old
    )
    wrapper = substitute(wrapper, old, new)
    wrapper = substitute(
        wrapper,
        "if (stage_ != Stage::Live && useWaveStartup_()) wave_startup_.mag(mag_body_ned);",
        "if (useWaveStartup_() && (stage_ != Stage::Live || (joint_refine_started_ && !mag_refine_done_))) wave_startup_.mag(mag_body_ned);",
    )
    wrapper = substitute(
        wrapper,
        "bool wave_startup_used_ = false;",
        r"""bool wave_startup_used_ = false;
    bool joint_refine_started_ = false;""",
    )
    wrapper = substitute(
        wrapper,
        r"""        setMagWorldRef_(mag_world_ref_uT);

        const float mag_tilt_yaw_rad = mag_auto_tuner_.getYawGaugeCorrectionRad();""",
        r"""        WaveStartup::Seed refined;
        const bool have_joint = wave_startup_.solve(cfg_.gravity_magnitude, cfg_.wave_startup_scale_slope, refined);
        const int variant = std::stoi(std::getenv("W3D_REFINE_JOINT"));
        if (variant && have_joint) {
            mag_world_ref_uT = refined.field;
            std::cerr << "JOINT_REFERENCE_APPLIED " << t_ << "\n";
        }
        setMagWorldRef_(mag_world_ref_uT);

        const float mag_tilt_yaw_rad = variant >= 2 && have_joint
            ? -yawFromBoatQuatRad(refined.attitude)
            : mag_auto_tuner_.getYawGaugeCorrectionRad();""",
    )
    return wrapper


def early_bias(wrapper):
    wrapper = substitute(
        wrapper,
        "wave_startup_used_ = use_learned;",
        r"""wave_startup_used_ = use_learned;
        if (use_learned && std::getenv("W3D_EARLY_BIAS") && std::string(std::getenv("W3D_EARLY_BIAS")) == "1")
            impl_.setAccBiasHold(false);""",
    )
    return wrapper


def frontend_wrapper(wrapper):
    return substitute(
        wrapper,
        "wave_startup_used_ = use_learned;",
        """wave_startup_used_ = use_learned;
        if (use_learned) impl_.seedStartupProxy(learned.attitude, learned.gyro_bias,
                                               std::stoi(std::getenv("W3D_PROXY_SEED")));""",
    )


def frontend_sources(overlay):
    # Preserve observer clocks, wave statistics, raw gyro endpoints and gains.
    # The accepted sensor-only startup fit is independent of live MEKF state.
    path = Path("kalman_common/MarineWaveFrontEnd.h")
    source = (ROOT / "src" / path).read_text()
    source = substitute(
        source,
        "    void resetFrontendGyroHistory()",
        """    void seedStartupProxy(const Eigen::Quaternionf& q, const Eigen::Vector3f& b, int mask) {
        vertical_accel_comp_.seedFromWave(q, b, mask);
    }

    void resetFrontendGyroHistory()""",
    )
    (overlay / path).write_text(source)
    path = Path("tuner/VerticalAccelComplementary.h")
    source = (ROOT / "src" / path).read_text()
    source = substitute(
        source,
        "    void setGains(float two_kp, float two_ki) {",
        """    void seedFromWave(const Eigen::Quaternionf& attitude, const Eigen::Vector3f& bias, int mask) {
        if (mask & 1) {
            const Eigen::Quaternionf q = attitude.normalized();
            ahrs_.q0 = q.w(); ahrs_.q1 = q.x(); ahrs_.q2 = q.y(); ahrs_.q3 = q.z();
        }
        if (mask & 2) {
            ahrs_.integralFBx = -bias.x();
            ahrs_.integralFBy = -bias.y();
            ahrs_.integralFBz = -bias.z();
        }
    }

    void setGains(float two_kp, float two_ki) {""",
    )
    (overlay / path).parent.mkdir(parents=True, exist_ok=True)
    (overlay / path).write_text(source)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", choices=POLICIES, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--seed", type=int, help="optional independent sensor seed")
    parser.add_argument("--policies", nargs="+", help="optional subset of the selected study policies")
    args = parser.parse_args()
    policies = POLICIES[args.study]
    if args.policies:
        if not set(args.policies) <= policies.keys():
            parser.error("unknown policy for the selected study")
        policies = {key: policies[key] for key in args.policies}
    args.output.mkdir(parents=True, exist_ok=True)
    inputs = [
        (case, args.data.resolve() / case / name)
        for case, records in REFERENCE_CSV_SHA256_BY_CASE.items()
        for name in records
    ]
    provenance = {
        case: input_provenance([path for selected, path in inputs if selected == case], case)
        for case in REFERENCE_CSV_SHA256_BY_CASE
    }
    transform = {
        "refinement": joint_refinement,
        "retention": early_bias,
        "bias_prior": early_bias,
        "frontend": frontend_wrapper,
    }.get(args.study)
    sources_before_build = startup_snapshot(__file__)
    with tempfile.TemporaryDirectory(prefix="startup-live-") as directory:
        binary = build(
            Path(directory),
            control=True,
            transform_wrapper=transform,
            transform_sources=frontend_sources if args.study == "frontend" else None,
        )
        sources_before_build.verify()

        def run(job):
            (case, path), (policy, overrides) = job
            env = {
                key: value for key, value in os.environ.items() if not key.startswith(("W3D_", "SF_", "OU_", "TFG_"))
            }
            env.update(
                OCEAN_IMU_SIM_CASE=case,
                W3D_HANDOVER_MASK="15",
                W3D_VALIDATION_WINDOW_SEC="900",
                W3D_WRITE_TIMESERIES="0",
                W3D_COLLECT_ALL_GATES="1",
            )
            if args.study in ("retention", "timeline", "bias_prior", "frontend"):
                env["W3D_VALIDATION_SEGMENTS"] = SEGMENTS
            if args.seed is not None:
                env["W3D_SEED"] = str(args.seed)
            env.update(overrides)
            with tempfile.TemporaryDirectory(dir=directory) as cwd:
                result = subprocess.run(
                    [str(binary), "--input", str(path)], cwd=cwd, env=env, text=True, capture_output=True
                )
            output = result.stdout + result.stderr
            if result.returncode not in (0, 1) or "QUALITY_GATE:" not in output:
                raise ValueError(output[-2000:])
            lines = [line for line in result.stdout.splitlines() if line.startswith("VALIDATION_METRICS ")]
            if len(lines) != 1:
                raise ValueError("expected one complete validation window")
            metrics = dict(token.split("=", 1) for token in lines[0].split()[1:])
            if metrics["samples"] != "180000" or metrics["start_s"] != "300":
                raise ValueError("incomplete 900-second score")
            log = args.output / f"{case}-{policy}-{path.stem}.log"
            log.write_text(output)
            return dict(case=case, arm=policy, **metrics, joint_applied=int("JOINT_REFERENCE_APPLIED" in output))

        jobs = [(pair, policy) for pair in inputs for policy in policies.items()]
        with ThreadPoolExecutor(args.jobs) as pool:
            rows = list(pool.map(run, jobs))
        sources_before_build.verify()
        segments = []
        for (case, path), (policy, _) in jobs:
            output = (args.output / f"{case}-{policy}-{path.stem}.log").read_text()
            for line in output.splitlines():
                if line.startswith("VALIDATION_SEGMENT "):
                    metrics = dict(token.split("=", 1) for token in line.split()[1:])
                    segments.append(dict(case=case, arm=policy, **metrics))
        if segments:
            with (args.output / "segments.csv").open("w") as stream:
                writer = csv.DictWriter(stream, list(segments[0]), lineterminator="\n")
                writer.writeheader()
                writer.writerows(segments)
        with (args.output / "records.csv").open("w") as stream:
            writer = csv.DictWriter(stream, list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        receipt = dict(
            source_commit=sources_before_build.revision,
            source_sha256=sources_before_build.hashes,
            script_sha256=sources_before_build.hashes[str(Path(__file__).resolve().relative_to(ROOT))],
            binary_sha256=digest(binary),
            records_sha256=digest(args.output / "records.csv"),
            policies=policies,
            sensor_seed=args.seed,
            inputs=provenance,
            scope="Selected sensor realization, both directions, full 1200-second replay. Experimental interventions.",
        )
        (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
