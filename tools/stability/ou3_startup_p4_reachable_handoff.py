#!/usr/bin/env python3
"""Shipping-reachable startup -> Live/P4 handoff contract for OU-III.

This is a theorem-facing structural lemma, not a replay and not a covariance
confidence argument.  It reads the shipping implementation and records the
state fiber that can actually be presented to a Live stability proof.

The essential coordinate choice is the estimator's natural handoff-local
translation gauge.  For a handoff time t_h define

    p^h_true(t) = p_phys(t) - p_phys(t_h),
    S^h_true(t) = integral_{t_h}^t p^h_true(s) ds.

The shipping MEKF has never propagated before t_h, and its linear state was
constructed at zero.  Hence at t_h

    v_hat = p_hat = S_hat = a_w_hat = 0,
    e_v = v_true(t_h), e_p = 0, e_S = 0, e_aw = a_w_true(t_h).

In particular the historical independent ||e_S||<=300 m*s ball is not a
shipping-reachable uncertainty at the handoff.  Large S bounds belong to the
Live working/prefix domain after integration begins, not to the entrance fiber.

The outer shipping timeout can hand over without magnetic north and without
TunerReady.  Therefore the capture theorem cannot identify goLive() with entry
into a fully-gauged/final-tuner P4 cell.  It must carry the same-history early
Live/H18 transient (or make the P4 basin wide enough to contain it) until the
actual tuner/magnetic state reaches a certified P4 section.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SHIPPING_REACHABLE_STARTUP_TO_P4_HANDOFF_V1"


def _config_float(text: str, name: str) -> float:
    m = re.search(rf"float\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f\s*;", text)
    if not m:
        raise RuntimeError(f"cannot extract shipping Config::{name}")
    return float(m.group(1))


def build() -> dict:
    wrapper = WRAPPER.read_text(encoding="utf-8")
    mekf = MEKF.read_text(encoding="utf-8")
    domain = json.loads(DOMAIN.read_text(encoding="utf-8"))

    timeout = _config_float(wrapper, "proxy_startup_timeout_sec")
    mag_settle = _config_float(wrapper, "proxy_mag_settle_sec")
    mag_window = _config_float(wrapper, "mag_min_window_sec")
    mag_fallback = _config_float(wrapper, "mag_tilt_fallback_sec")
    deployed_timeout = max(timeout, mag_settle + 2.0 * max(mag_window, 1.0) + mag_fallback)

    parity = {
        "bootstrap_uses_frontend_only": "impl_.updateFrontEnd(dt, gyro_body_ned, acc_body_ned);" in wrapper,
        "frontend_does_not_drive_mekf": "updateCore_(dt, gyro, acc, /*tempC=*/35.0f, /*drive_mekf=*/false);" in wrapper,
        "constructor_zeroes_full_state": "xext.setZero();" in mekf,
        "constructor_seeds_translation_covariance_only": all(s in mekf for s in (
            "const T sigma_v0 = T(1.0)", "const T sigma_p0 = T(20.0)", "const T sigma_S0 = T(50.0)",
            "set_initial_linear_uncertainty(sigma_v0, sigma_p0, sigma_S0);")),
        "goLive_only_reseeds_attitude_state": "mekf_->initialize_from_attitude(q_bw, tilt_sigma_rad, yaw_sigma_rad);" in wrapper,
        "goLive_enters_live": "enterLive_();" in wrapper,
        "live_entry_resets_aw_covariance_not_aw_state": "mekf_->reset_aw_covariance_to_stationary();" in wrapper,
        "aw_covariance_reset_clears_cross_covariances": "Pext.template block<3,1>(OFF_AW, i).setZero();" in mekf,
        "timeout_does_not_require_tuner_ready": "(t_ >= timeout_sec) &&\n            mag_gravity_aligned_branch_" in wrapper,
        "quality_path_requires_tuner_ready": "impl_.isTunerReady();" in wrapper,
        "quality_path_requires_north_when_mag": "const bool north_ready = !cfg_.with_mag || mag_ref_set_;" in wrapper,
        "timeout_may_have_no_yaw_gauge": "? cfg_.proxy_handoff_yaw_sigma_rad\n            : cfg_.proxy_handoff_yaw_sigma_free_rad" in wrapper,
        "scheduler_constructed_at_zero_elapsed": "T pseudo_update_elapsed_s_ = T(0);" in mekf,
    }
    failures = [k for k, v in parity.items() if not v]

    legacy = domain["startup"]["physical_handoff_coordinate_bounds"]
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "shipping_source_parity": parity,
        "shipping_source_parity_failures": failures,
        "filter_changed": False,
        "quality_gates_changed": False,
        "p3_delta_changed": False,
        "deployed_live_handoff_upper_bound_s": deployed_timeout,
        "handoff_translation_coordinate_gauge": {
            "definition": "p_true_h(t)=p_phys(t)-p_phys(t_h); S_true_h(t)=integral_[t_h,t] p_true_h(s) ds",
            "translation_invariant_dynamics_preserved": True,
            "truth_position_at_handoff_m": [0.0, 0.0, 0.0],
            "truth_integral_displacement_at_handoff_m_s": [0.0, 0.0, 0.0],
            "covariance_confidence_used_as_true_error_set": False,
        },
        "estimated_state_at_handoff": {
            "velocity_mps": [0.0, 0.0, 0.0],
            "position_m": [0.0, 0.0, 0.0],
            "integral_displacement_m_s": [0.0, 0.0, 0.0],
            "latent_acceleration_mps2": [0.0, 0.0, 0.0],
            "gyro_bias_rad_s": [0.0, 0.0, 0.0],
            "accelerometer_bias_mps2": [0.0, 0.0, 0.0],
            "linear_state_reset_at_goLive": False,
            "states_inherited_from_constructor_without_preLive_propagation": True,
        },
        "true_error_at_handoff": {
            "velocity": "e_v=v_true(t_h); requires SAME-BRMM physical velocity bound/correlation",
            "position": "e_p=0 exactly in handoff-local displacement gauge",
            "integral_displacement": "e_S=0 exactly in handoff-local integral gauge",
            "latent_acceleration": "e_aw=a_w_true(t_h); requires SAME-BRMM source bound/correlation",
            "gyro_bias": "e_bg=b_g,true(t_h) because b_g,hat=0",
            "accelerometer_bias": "e_ba=b_a,true(t_h) because b_a,hat=0 and is held in H18",
        },
        "covariance_at_handoff": {
            "translation_state_covariance_propagated_before_live": False,
            "translation_cross_covariances_zero_before_first_live_prediction": True,
            "attitude_block_reseeded_from_proxy_tilt_yaw_sigmas": True,
            "attitude_cross_covariances_dropped_by_initialize_from_attitude": True,
            "aw_marginal_reseeded_to_current_committed_stationary_covariance": True,
            "aw_cross_covariances_cleared": True,
            "scheduler_elapsed_s": 0.0,
        },
        "legacy_entry_model_audit": {
            "legacy_independent_S_radius_m_s": float(legacy["integral_displacement_error_norm_upper_m_s"]),
            "legacy_independent_position_radius_m": float(legacy["position_error_norm_upper_m"]),
            "independent_S_ball_shipping_reachable_at_handoff": False,
            "independent_position_ball_shipping_reachable_at_handoff": False,
            "S_and_position_must_move_to_working_prefix_domain": True,
            "velocity_must_not_be_zeroed_as_error": True,
            "legacy_300_m_s_may_promote_P4": False,
        },
        "live_capture_modes": {
            "quality_handoff": "proxy+tilt+north(if fitted)+TunerReady -> Live/H18",
            "timeout_handoff": "proxy+aligned-gravity-branch -> Live/H18; north/tuner may still be unready",
            "timeout_requires_full_P4_membership_at_goLive": False,
            "capture_theorem_must_cover_early_live_transient": True,
            "H18_to_A21_is_separate_hybrid_release": True,
        },
        "theorem_architecture": [
            "STARTUP_REACHABLE_TUBE",
            "SHIPPING_LIVE_HANDOFF_FIBER",
            "EARLY_LIVE_H18_CAPTURE_TUBE",
            "P4_H18_INVARIANT_BASIN",
            "H18_TO_A21_HYBRID_MAP",
            "P4_A21_INVARIANT_BASIN",
        ],
        "P4_ENTRY_MUST_BE_SHIPPING_REACHABLE": True,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_PASS": False,
        "remaining_obligations": [
            "source-uniform Mahony/proxy tilt-error reachable tube over all admitted startup BRMM/BIAS histories",
            "source-uniform timeout early-Live tuner/magnetic/scheduler reachable tube",
            "same-signal COMPLETE-BRMM coefficient/source cover and endpoint/every-prefix augmented LDLT",
            "binary32 complete Eigen/Kalman/reset additive ISS enclosure",
            "source-dependent compatible storage and first-exit retention for the widened correlated basin",
            "finite-time intersection of startup/early-Live reachable tube with that certified P4 basin",
        ],
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    f.extend(d.get("shipping_source_parity_failures", []))
    if d.get("true_error_at_handoff", {}).get("position") != "e_p=0 exactly in handoff-local displacement gauge":
        f.append("position handoff error gauge lost")
    if d.get("true_error_at_handoff", {}).get("integral_displacement") != "e_S=0 exactly in handoff-local integral gauge":
        f.append("S handoff error gauge lost")
    if d.get("legacy_entry_model_audit", {}).get("legacy_300_m_s_may_promote_P4") is not False:
        f.append("legacy independent S ball may not promote P4")
    if d.get("live_capture_modes", {}).get("capture_theorem_must_cover_early_live_transient") is not True:
        f.append("timeout early-Live capture transient omitted")
    for gate in ("P4_MOTION_PASS", "P4_PASS", "P5_PASS"):
        if d.get(gate) is not False:
            f.append(f"{gate} promoted before closure")
    return list(dict.fromkeys(f))


if __name__ == "__main__":
    d = build()
    failures = validate(d)
    print(json.dumps({**d, "validation_pass": not failures, "validation_failures": failures}, indent=2))
    raise SystemExit(1 if failures else 0)
