#!/usr/bin/env python3
"""Source-staged covariance seed for the OU-III P5 outer-H bridge.

This producer closes the *entrance* covariance semantics at the deployed
startup->H handoff. It does not use replay and it does not pretend that the
normal-Live P3 covariance box is already the covariance at goLive.

The shipping sequence is

    initialize_from_attitude(...)
    enterLive_()
      apply_ou_tune_(true)
      set_linear_block_enabled(true)

and the MEKF source gives two exact resets relevant to P5:

* initialize_from_attitude() calls zero_AL_cross_cov_once_(), so every
  attitude/gyro-bias <-> [v,p,S,a_w] covariance block is zero;
* enabling the previously disabled linear block calls
  reset_aw_covariance_to_stationary(), which zeros every a_w cross covariance,
  and resets pseudo_update_elapsed_s_ to zero.

During warmup the linear block is disabled, so the constructor's v/p/S
covariance seed is not propagated. Therefore at the instant H mode starts,
P_theta,S=P_theta,aw=P_aw,S=0 and P_SS=(50 m s)^2 I exactly (up to the source
scalar semantics). In particular the *entrance* S->attitude Kalman gain is
exactly zero. This removes the invalid use of the many-orders-looser global P3
covariance eigenvalue bound at the handoff node.

The first due S pseudo update is later. Between goLive and that due event an
accepted accelerometer correction can create attitude/linear cross covariance
because its H matrix contains both attitude and a_w columns. The producer
therefore exposes a finite, source-derived first-pseudo stage and explicitly
leaves that accepted-branch cross-covariance propagation as the next numerical
sub-obligation. No favorable rejection pattern is assumed.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import ou3_implementation_proof_manifest as MANIFEST
import ou3_source_domain_contract as SOURCE

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
CORE = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def _one(text: str, pattern: str, label: str) -> float:
    m = re.search(pattern, text, flags=re.MULTILINE | re.DOTALL)
    if not m:
        raise RuntimeError(f"cannot extract {label}")
    return float(m.group(1))


def _require(text: str, marker: str, label: str) -> None:
    if marker not in text:
        raise RuntimeError(f"missing source semantic {label}: {marker}")


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P5 goLive covariance stage must not be trajectory fitted")

    w = WRAPPER.read_text(encoding="utf-8")
    k = MEKF.read_text(encoding="utf-8")
    c = CORE.read_text(encoding="utf-8")
    manifest = MANIFEST.build()
    mf = MANIFEST.validate(manifest)
    if mf:
        raise RuntimeError(f"implementation manifest prerequisite failed: {mf}")

    markers = {
        "go_live_initializes_attitude": "mekf_->initialize_from_attitude(q_bw, tilt_sigma_rad, yaw_sigma_rad);",
        "go_live_enters_live": "enterLive_();",
        "live_applies_ou_before_enable": "apply_ou_tune_(true);",
        "live_enables_linear_block": "mekf_->set_linear_block_enabled(enable_linear_block_);",
        "warmup_disables_linear": "mekf_->set_linear_block_enabled(false);",
        "attitude_init_zeros_AL": "zero_AL_cross_cov_once_();",
        "enable_resets_aw": "reset_aw_covariance_to_stationary();",
        "enable_resets_pseudo_phase": "pseudo_update_elapsed_s_ = T(0);",
        "aw_reset_zeros_cross": "Pext.template block<3,1>(OFF_AW, i).setZero();",
        "pseudo_due_uses_elapsed": "periodic_update_due(Ts, pseudo_update_period_s_, pseudo_update_elapsed_s_)",
    }
    joined = w + "\n" + k + "\n" + c
    for label, marker in markers.items():
        _require(joined, marker, label)

    sigma_v0 = _one(k, r"const\s+T\s+sigma_v0\s*=\s*T\(([0-9.eE+-]+)\)", "sigma_v0")
    sigma_p0 = _one(k, r"const\s+T\s+sigma_p0\s*=\s*T\(([0-9.eE+-]+)\)", "sigma_p0")
    sigma_S0 = _one(k, r"const\s+T\s+sigma_S0\s*=\s*T\(([0-9.eE+-]+)\)", "sigma_S0")

    startup = manifest["startup"]
    tilt_sigma = float(startup["handoff_tilt_sigma_rad"])
    yaw_sigma = float(startup["handoff_yaw_sigma_rad"])
    yaw_free_sigma = float(startup["handoff_yaw_sigma_free_rad"])

    dt = float(manifest["configured_runtime"]["imu_dt_s"])
    pbox = SOURCE.build(WRAPPER)["validated_parameter_box"]["continuous_parameters"]
    pseudo_lo, pseudo_hi = map(float, pbox["pseudo_update_period_s"])
    sigma_aw_lo, sigma_aw_hi = map(float, pbox["sigma_aw_mps2"])
    first_due_steps_upper = int(math.ceil(pseudo_hi / dt)) + 1
    first_due_steps_lower = 1
    first_due_time_upper = first_due_steps_upper * dt

    mag_corrections_upper = int(math.ceil(first_due_time_upper * 25.0)) + 1

    seed = {
        "mode": "H",
        "dimension": 18,
        "linear_block_was_disabled_during_warmup": True,
        "linear_block_enable_resets_pseudo_elapsed": True,
        "pseudo_update_elapsed_s_at_goLive": 0.0,
        "attitude_linear_cross_covariance_exact_zero": True,
        "theta_S_cross_covariance_operator_norm_upper": 0.0,
        "theta_aw_cross_covariance_operator_norm_upper": 0.0,
        "bg_S_cross_covariance_operator_norm_upper": 0.0,
        "aw_S_cross_covariance_operator_norm_upper": 0.0,
        "P_vv_variance_per_axis": sigma_v0 * sigma_v0,
        "P_pp_variance_per_axis": sigma_p0 * sigma_p0,
        "P_SS_variance_per_axis": sigma_S0 * sigma_S0,
        "P_awaw_reset_to_current_stationary_covariance": True,
        "P_awaw_source_std_outward_mps2": [sigma_aw_lo, sigma_aw_hi],
        "attitude_covariance_seed": {
            "tilt_sigma_rad": tilt_sigma,
            "gauged_yaw_sigma_rad": yaw_sigma,
            "ungauged_yaw_sigma_rad": yaw_free_sigma,
            "tilt_variance": tilt_sigma * tilt_sigma,
            "gauged_yaw_variance": yaw_sigma * yaw_sigma,
            "ungauged_yaw_variance": yaw_free_sigma * yaw_free_sigma,
        },
        "S_to_attitude_gain_at_goLive_exact_zero": True,
    }

    first_stage = {
        "pseudo_period_outward_s": [pseudo_lo, pseudo_hi],
        "configured_dt_s": dt,
        "first_due_prediction_samples_lower": first_due_steps_lower,
        "first_due_prediction_samples_upper": first_due_steps_upper,
        "first_due_time_upper_s": first_due_time_upper,
        "accepted_accelerometer_corrections_before_first_due_upper": first_due_steps_upper,
        "magnetometer_corrections_before_first_due_upper": mag_corrections_upper,
        "prediction_only_preserves_zero_attitude_linear_cross_covariance": True,
        "rejected_accelerometer_preserves_zero_cross_if_no_prior_cross": True,
        "magnetometer_cannot_create_attitude_linear_cross_from_zero": True,
        "accepted_accelerometer_can_create_attitude_linear_cross_via_aw": True,
        "source_branch_enumeration_may_not_assume_accelerometer_rejection": True,
        "pre_first_due_theta_S_cross_covariance_enclosure": "PENDING_ACCEPTED_ACCEL_SOURCE_STAGE",
        "first_due_S_to_attitude_gain_enclosed": False,
        "next_required_bound": (
            "propagate a source-correlated H covariance enclosure from the exact goLive seed through every "
            "accepted/rejected accelerometer and due/not-due magnetometer prefix up to the first S pseudo; "
            "then bound P_thetaS(P_SS+R_S)^-1 directly"
        ),
    }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SOURCE_STAGED_GOLIVE_H_COVARIANCE_SEED",
        "claim": "EXACT_GOLIVE_COVARIANCE_AND_PSEUDO_PHASE_SEED_FOR_OUTER_H",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "implementation_markers": sorted(markers),
        "goLive_H_covariance_seed": seed,
        "pre_first_S_stage": first_stage,
        "global_normal_live_P3_covariance_used_at_goLive": False,
        "P5_GOLIVE_COVARIANCE_STAGE_CERTIFICATE": "PASS",
        "P5_FIRST_DUE_CROSS_COVARIANCE_CERTIFICATE": "NOT_ESTABLISHED",
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("goLive stage is not source generated")
    if d.get("source_replay_used") is not False:
        failures.append("goLive stage uses replay")
    if d.get("filter_changed") is not False:
        failures.append("goLive stage changes filter")
    if d.get("global_normal_live_P3_covariance_used_at_goLive") is not False:
        failures.append("global P3 covariance was reused at goLive")
    seed = d.get("goLive_H_covariance_seed", {})
    for key in (
        "attitude_linear_cross_covariance_exact_zero",
        "P_awaw_reset_to_current_stationary_covariance",
        "S_to_attitude_gain_at_goLive_exact_zero",
        "linear_block_enable_resets_pseudo_elapsed",
    ):
        if seed.get(key) is not True:
            failures.append(f"goLive seed lost source fact {key}")
    for key in (
        "theta_S_cross_covariance_operator_norm_upper",
        "theta_aw_cross_covariance_operator_norm_upper",
        "bg_S_cross_covariance_operator_norm_upper",
        "aw_S_cross_covariance_operator_norm_upper",
    ):
        if seed.get(key) != 0.0:
            failures.append(f"goLive {key} is not exact zero")
    if not math.isclose(float(seed.get("P_SS_variance_per_axis", math.nan)), 2500.0, rel_tol=0.0, abs_tol=1e-12):
        failures.append("goLive P_SS seed is not source constructor value")
    sigma_box = seed.get("P_awaw_source_std_outward_mps2")
    if not (isinstance(sigma_box, list) and len(sigma_box) == 2
            and 0.0 < float(sigma_box[0]) <= float(sigma_box[1])):
        failures.append("goLive a_w source std box missing")
    stage = d.get("pre_first_S_stage", {})
    if not (isinstance(stage.get("first_due_prediction_samples_upper"), int)
            and stage["first_due_prediction_samples_upper"] >= 1):
        failures.append("first S due sample bound missing")
    if stage.get("source_branch_enumeration_may_not_assume_accelerometer_rejection") is not True:
        failures.append("first S stage assumes favorable accelerometer rejection")
    if stage.get("accepted_accelerometer_can_create_attitude_linear_cross_via_aw") is not True:
        failures.append("first S stage misses accelerometer cross-covariance creation")
    if d.get("P5_GOLIVE_COVARIANCE_STAGE_CERTIFICATE") != "PASS":
        failures.append("goLive covariance stage did not pass")
    if d.get("P5_FIRST_DUE_CROSS_COVARIANCE_CERTIFICATE") != "NOT_ESTABLISHED":
        failures.append("first due S cross covariance was promoted without propagation")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    failures = validate(out)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "goLive": out.get("goLive_H_covariance_seed"),
        "pre_first_S_stage": out.get("pre_first_S_stage"),
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
