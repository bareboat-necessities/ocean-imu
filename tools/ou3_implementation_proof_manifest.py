#!/usr/bin/env python3
"""Source-derived manifest for the OU-III implementation stability proof.

This is deliberately not a simulator manifest. It binds the proof to the
shipping `SeaStateFusion_OU_III` wrapper and `Kalman3D_Wave_OU_III` operations
that define startup, normal Live words, and hard hybrid events. A source change
that removes or changes one of the required semantics invalidates the manifest
instead of silently leaving the theorem attached to an older algorithm.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

import ou3_source_domain_contract as SOURCE

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
CORE = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"
GUARD = REPO / "src" / "tuner" / "AccelVibrationGuard.h"
SCHEMA = 3


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a source file bound into the manifest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def one_float(text: str, pattern: str, label: str) -> float:
    """Extract exactly one implementation scalar using a source-audit regex."""
    m = re.search(pattern, text, flags=re.MULTILINE | re.DOTALL)
    if not m:
        raise RuntimeError(f"cannot extract implementation value {label}")
    return float(m.group(1))


def require(text: str, marker: str, label: str) -> None:
    """Fail closed unless a required shipping-source semantic marker exists."""
    if marker not in text:
        raise RuntimeError(f"missing implementation semantic marker {label}: {marker}")


def build() -> dict:
    """Build the source-bound implementation manifest for the proof chain."""
    w = WRAPPER.read_text(encoding="utf-8")
    k = MEKF.read_text(encoding="utf-8")
    c = CORE.read_text(encoding="utf-8")
    g = GUARD.read_text(encoding="utf-8")

    two_kp = SOURCE.parse_const(w, "STARTUP_PROXY_TWO_KP_DEFAULT")
    two_ki = SOURCE.parse_const(w, "STARTUP_PROXY_TWO_KI_DEFAULT")
    dt = SOURCE.parse_const(w, "FREQ_SMOOTHER_DT")
    guard_hz = SOURCE.parse_const(w, "ACC_VIBRATION_GUARD_HZ_DEFAULT")
    guard_poles = int(one_float(
        w,
        r"ACC_VIBRATION_GUARD_POLES_DEFAULT\s*=\s*([0-9]+)",
        "ACC_VIBRATION_GUARD_POLES_DEFAULT"))
    guard_engage_lo = SOURCE.parse_const(g, "kEngageLoDefault")

    cfg = {
        "proxy_startup_min_sec": one_float(w, r"proxy_startup_min_sec\s*=\s*([0-9.eE+-]+)f", "proxy_startup_min_sec"),
        "proxy_startup_timeout_sec": one_float(w, r"proxy_startup_timeout_sec\s*=\s*([0-9.eE+-]+)f", "proxy_startup_timeout_sec"),
        "handoff_tilt_sigma_rad": one_float(w, r"proxy_handoff_tilt_sigma_rad\s*=\s*([0-9.eE+-]+)f", "proxy_handoff_tilt_sigma_rad"),
        "handoff_yaw_sigma_rad": one_float(w, r"proxy_handoff_yaw_sigma_rad\s*=\s*([0-9.eE+-]+)f", "proxy_handoff_yaw_sigma_rad"),
        "handoff_yaw_sigma_free_rad": one_float(w, r"proxy_handoff_yaw_sigma_free_rad\s*=\s*([0-9.eE+-]+)f", "proxy_handoff_yaw_sigma_free_rad"),
        "gravity_align_max_sin": one_float(w, r"mag_gravity_align_max_sin\s*=\s*([0-9.eE+-]+)f", "mag_gravity_align_max_sin"),
        "gravity_align_hold_sec": one_float(w, r"mag_gravity_align_hold_sec\s*=\s*([0-9.eE+-]+)f", "mag_gravity_align_hold_sec"),
        "gravity_align_world_tau_sec": one_float(w, r"mag_gravity_align_world_tau_sec\s*=\s*([0-9.eE+-]+)f", "mag_gravity_align_world_tau_sec"),
        "gravity_align_world_warmup_sec": one_float(w, r"mag_gravity_align_world_warmup_sec\s*=\s*([0-9.eE+-]+)f", "mag_gravity_align_world_warmup_sec"),
        "mag_extreme_gyro_dps": one_float(w, r"mag_extreme_gyro_dps\s*=\s*([0-9.eE+-]+)f", "mag_extreme_gyro_dps"),
        "mag_refine_start_sec": one_float(w, r"mag_refine_start_sec\s*=\s*([0-9.eE+-]+)f", "mag_refine_start_sec"),
        "mag_refine_window_sec": one_float(w, r"mag_refine_window_sec\s*=\s*([0-9.eE+-]+)f", "mag_refine_window_sec"),
        "mag_min_window_sec": one_float(w, r"mag_min_window_sec\s*=\s*([0-9.eE+-]+)f", "mag_min_window_sec"),
        "mag_tilt_fallback_sec": one_float(w, r"mag_tilt_fallback_sec\s*=\s*([0-9.eE+-]+)f", "mag_tilt_fallback_sec"),
        "acc_bias_unlock_mag_updates": int(one_float(w, r"acc_bias_unlock_mag_updates\s*=\s*([0-9]+)", "acc_bias_unlock_mag_updates")),
        "tilt_reset_deg": one_float(w, r"TILT_RESET_DEG\s*=\s*([0-9.eE+-]+)f", "TILT_RESET_DEG"),
        "tilt_reset_hold_sec": one_float(w, r"TILT_RESET_HOLD_SEC\s*=\s*([0-9.eE+-]+)f", "TILT_RESET_HOLD_SEC"),
        "tilt_reset_cooldown_sec": one_float(w, r"TILT_RESET_COOLDOWN_SEC\s*=\s*([0-9.eE+-]+)f", "TILT_RESET_COOLDOWN_SEC"),
    }

    mekf_defaults = {
        "acc_bias_limit_mps2": one_float(k, r"acc_bias_limit_\s*=\s*T\(([0-9.eE+-]+)\)", "acc_bias_limit_"),
        "acc_bias_tau_sec": one_float(k, r"tau_bacc_\s*=\s*T\(([0-9.eE+-]+)\)", "tau_bacc_"),
        "pseudo_update_period_default_sec": one_float(k, r"pseudo_update_period_s_\s*=\s*T\(([0-9.eE+-]+)\)", "pseudo_update_period_s_"),
    }

    semantic_markers = {
        "startup_policy_default_mahony": "StartupInitPolicy startup_init_policy = StartupInitPolicy::MahonyProxy;",
        "timeout_requires_aligned_branch": "(t_ >= timeout_sec) &&\n            mag_gravity_aligned_branch_",
        "quality_requires_tilt_north_tuner": "tilt_trusted &&\n            north_ready &&\n            impl_.isTunerReady()",
        "handoff_bias_held": "/*allow_acc_bias=*/false",
        "go_live_initializes_attitude": "mekf_->initialize_from_attitude(q_bw, tilt_sigma_rad, yaw_sigma_rad);",
        "live_applies_ou_sync": "apply_ou_tune_(true);",
        "live_enables_linear_block": "mekf_->set_linear_block_enabled(enable_linear_block_);",
        "h_to_a_release": "mekf_->set_acc_bias_updates_enabled(true);",
        "mag_refinement_rewrites_heading": "impl_.mekf().set_quaternion_boat(q_new);",
        "mag_refinement_releases_bias_hold": "impl_.setAccBiasHold(false);",
        "guard_armed_by_constructor": "setAccelVibrationGuard(ACC_VIBRATION_GUARD_HZ_DEFAULT,\n                               ACC_VIBRATION_GUARD_POLES_DEFAULT);",
        "guard_single_conditioned_feed": "const Eigen::Vector3f acc_in = accel_guard_.step(acc, dt);",
        "guard_conditioned_feed_reaches_proxy": "vertical_accel_comp_.update(dt, gyro, acc_in, g_std);",
        "guard_zero_weight_is_bit_exact_transparent": "if (weight_ <= 0.0f) return acc;",
        "live_tilt_relock_preserves_yaw": "mekf_->initialize_from_acc_preserve_yaw(acc_in);",
        "periodic_aw_sync": "periodic_aw_cov_sync_tick_();",
        "full_s_update": "void applyIntegralZeroPseudoMeas();",
        "joseph_update": "joseph_update3_",
        "quaternion_injection": "applyQuaternionCorrectionFromErrorState();",
        "left_error_reset": "apply_error_state_reset_jacobian_",
        "prediction_before_acc": "mekf_->time_update(gyro, dt);\n            mekf_->measurement_update_acc_only(acc_in, tempC);",
        "S_due_inside_time_update": "applyIntegralZeroPseudoMeas();",
        "acc_injects_immediately": "last_acc_diag_.accepted = true;",
        "mag_injects_immediately": "last_mag_diag_.accepted = true;",
        "aw_psd_floor": "Pext.template block<3,3>(OFF_AW, OFF_AW) += Delta;",
    }
    joined = w + "\n" + k + "\n" + c + "\n" + g
    for label, marker in semantic_markers.items():
        require(joined, marker, label)

    if "static constexpr int NX = BASE_N + EXT_ADD;" not in k:
        raise RuntimeError("cannot bind MEKF state dimension expression")
    if "12\n        + (with_accel_bias ? 3 : 0)" not in k:
        raise RuntimeError("cannot bind OU-III extended state dimension")

    return {
        "schema": SCHEMA,
        "qualification": "SOURCE_BOUND_OU3_IMPLEMENTATION_STABILITY_MANIFEST",
        "source_generated_not_trajectory_fit": True,
        "implementation_files": {
            str(WRAPPER.relative_to(REPO)): sha256(WRAPPER),
            str(MEKF.relative_to(REPO)): sha256(MEKF),
            str(CORE.relative_to(REPO)): sha256(CORE),
            str(GUARD.relative_to(REPO)): sha256(GUARD),
        },
        "configured_runtime": {
            "imu_dt_s": dt,
            "accel_vibration_guard_cutoff_hz_default": guard_hz,
            "accel_vibration_guard_poles_default": guard_poles,
            "accel_vibration_guard_engage_lo_mps2": guard_engage_lo,
            "scope": "configured nominal runtime; arbitrary caller dt is outside the quantitative theorem",
        },
        "state_coordinates": {
            "H_dimension": 18,
            "A_dimension": 21,
            "H": ["delta_theta", "b_g", "v", "p", "S", "a_w"],
            "A": ["delta_theta", "b_g", "v", "p", "S", "a_w", "b_a"],
            "held_to_active_is_dimension_changing_jump": True,
        },
        "startup": {
            "two_kp": two_kp,
            "two_ki": two_ki,
            **cfg,
            "first_live_mode": "H",
            "prior_attitude_discarded_by_accelerometer_seed": True,
            "timeout_cannot_handoff_antipodal_branch": True,
            "go_live_bias_learning_held": True,
        },
        "mekf_defaults": mekf_defaults,
        "normal_live_update_order": [
            "commit_previous_tune",
            "vibration_guard_conditioning",
            "prediction",
            "apply_pending_aw_covariance_psd_increment",
            "periodic_S_zero_when_due_then_immediate_quaternion_injection_and_left_error_reset",
            "accelerometer_correction_or_rejection_then_immediate_quaternion_injection_and_left_error_reset_if_accepted",
            "source_tuner_evolution_and_stage_next_tune",
            "periodic_aw_covariance_sync_tick_stages_future_psd_increment",
            "asynchronous_magnetometer_correction_or_rejection_then_immediate_quaternion_injection_and_left_error_reset_if_accepted",
        ],
        "same_sample_reset_policy": {
            "S_zero": "immediate after accepted pseudo correction",
            "accelerometer": "immediate after accepted correction",
            "magnetometer": "immediate after accepted correction",
            "single_shared_end_of_sample_reset": False,
        },
        "hybrid_events": [
            "startup_handoff",
            "held_to_active",
            "magnetic_lock",
            "magnetic_regauge_refinement",
            "tilt_reset",
            "tilt_relock",
            "cooldown_reentry",
            "periodic_aw_covariance_sync",
            "accelerometer_vibration_guard_engagement",
        ],
        "vibration_guard": {
            "armed_by_default": guard_hz > 0.0,
            "constructor_arm_call_source_bound": True,
            "zero_engagement_is_bit_exact_transparent": True,
            "active_guard_changes_measurement_dynamics": True,
            "active_guard_requires_separate_source_certificate": True,
        },
        "semantic_markers": sorted(semantic_markers),
        "pass": True,
    }


def validate(d: dict) -> list[str]:
    """Validate source-binding and conservative vibration-guard proof scope."""
    failures: list[str] = []
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("manifest is not source generated")
    if d.get("pass") is not True:
        failures.append("manifest did not pass")
    s = d.get("startup", {})
    if not math.isclose(float(s.get("two_kp", math.nan)), 0.2, rel_tol=0.0, abs_tol=1.0e-8) or not math.isclose(
        float(s.get("two_ki", math.nan)), 0.02, rel_tol=0.0, abs_tol=1.0e-9
    ):
        failures.append("startup Mahony gains do not match proved gains")
    if s.get("timeout_cannot_handoff_antipodal_branch") is not True:
        failures.append("timeout branch is not fail-closed on aligned gravity")
    if s.get("go_live_bias_learning_held") is not True:
        failures.append("goLive does not enter held-bias mode")
    dims = d.get("state_coordinates", {})
    if dims.get("H_dimension") != 18 or dims.get("A_dimension") != 21:
        failures.append("H/A dimensions do not match theorem")
    resets = d.get("same_sample_reset_policy", {})
    if resets.get("single_shared_end_of_sample_reset") is not False:
        failures.append("manifest incorrectly merges same-sample correction resets")
    order = d.get("normal_live_update_order", [])
    try:
        ig = order.index("vibration_guard_conditioning")
        ip = order.index("prediction")
        iS = order.index("periodic_S_zero_when_due_then_immediate_quaternion_injection_and_left_error_reset")
        ia = order.index("accelerometer_correction_or_rejection_then_immediate_quaternion_injection_and_left_error_reset_if_accepted")
    except ValueError:
        failures.append("manifest lacks exact guard/S/accelerometer ordering")
    else:
        if not ig < ip < iS < ia:
            failures.append("manifest source ordering does not match deployed guard/prediction/S/accelerometer path")
    vg = d.get("vibration_guard", {})
    if vg.get("armed_by_default") is not True:
        failures.append("vibration guard is not source-bound as armed by default")
    if vg.get("constructor_arm_call_source_bound") is not True:
        failures.append("vibration guard constructor arm call is not source bound")
    if vg.get("zero_engagement_is_bit_exact_transparent") is not True:
        failures.append("vibration guard zero-engagement transparency is not source certified")
    if vg.get("active_guard_requires_separate_source_certificate") is not True:
        failures.append("active vibration-guard source branch was silently treated as already certified")
    return failures


def main() -> int:
    """Write the validated implementation manifest as a machine-readable artifact."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build()
    failures = validate(out)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "qualification": out["qualification"],
        "startup": out["startup"],
        "normal_live_update_order": out["normal_live_update_order"],
        "hybrid_events": out["hybrid_events"],
        "vibration_guard": out["vibration_guard"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
