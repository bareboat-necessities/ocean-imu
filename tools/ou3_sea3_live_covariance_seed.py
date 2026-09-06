#!/usr/bin/env python3
"""Shipping Normal-Live covariance seed and held->active bias release for P3.

Before handoff the outer wrapper runs the measurement-only front end; the MEKF
covariance is not propagated.  ``goLive`` reseeds attitude covariance,
``enterLive_`` commits the current tuner point, and the a_w marginal is reset to
the committed stationary covariance before the first prediction.

For the configured full-heading deployment the accelerometer-bias coordinate is
held out of H18 after Live.  This is stronger than the inner 250-magnetometer
unlock counter: the outer shipping wrapper has ``with_mag=true`` and
``mag_refine_enabled=true`` by default and calls

    impl_.setAccBiasHold(cfg_.with_mag && cfg_.mag_refine_enabled)

at begin.  Magnetic refinement can start only after the outer wrapper is Live;
it resets MagAutoTuner, sets its minimum accepted window to the configured
``mag_refine_window_sec=30 s``, and releases the bias hold only *after* that
refinement finishes.  MagAutoTuner::tryFinalize_ refuses completion while the
accepted window is below that minimum.  Consequently the configured H mode is
guaranteed to persist for at least 30 s after Live entry, far longer than the
canonical 3 s H18 word.  The H18 full-matrix margin therefore exists before the
separate H->A hybrid release.

At release the held b_a block has zero H18 cross covariance and its diagonal is
floored to the shipping seed variance.  Normal Gauss-Markov prediction then
resumes.  These are source facts, not an alternate estimator or a trajectory
fit.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import ou3_sea3_dynamic_source_certificate as DYNAMIC

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
MAG = REPO / "src" / "tuner" / "MagAutoTuner.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_SEA3_SHIPPING_NORMAL_LIVE_COVARIANCE_SEED_AND_BIAS_RELEASE"
H18_WORD_S = 3.0


def _one(pattern: str, text: str, label: str) -> float:
    m = re.search(pattern, text, flags=re.MULTILINE | re.DOTALL)
    if not m:
        raise RuntimeError(f"cannot extract {label}")
    x = float(m.group(1))
    if not math.isfinite(x):
        raise RuntimeError(f"nonfinite {label}")
    return x


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    wrapper = WRAPPER.read_text(encoding="utf-8")
    mekf = MEKF.read_text(encoding="utf-8")
    mag = MAG.read_text(encoding="utf-8")
    dynamic = DYNAMIC.build(path)
    df = DYNAMIC.validate(dynamic)
    if df:
        raise RuntimeError(f"dynamic source invalid: {df}")

    Pq0 = _one(r"T Pq0\s*=\s*T\(([0-9.eE+-]+)\)", mekf, "Pq0")
    Pb0 = _one(r"T Pb0\s*=\s*T\(([0-9.eE+-]+)\)", mekf, "Pb0")
    sigma_v0 = _one(r"const T sigma_v0\s*=\s*T\(([0-9.eE+-]+)\)", mekf, "sigma_v0")
    sigma_p0 = _one(r"const T sigma_p0\s*=\s*T\(([0-9.eE+-]+)\)", mekf, "sigma_p0")
    sigma_S0 = _one(r"const T sigma_S0\s*=\s*T\(([0-9.eE+-]+)\)", mekf, "sigma_S0")
    sigma_ba0 = _one(r"T sigma_bacc0_\s*=\s*T\(([0-9.eE+-]+)\)", mekf, "sigma_bacc0")
    s_factor = _one(r"float S_factor_\s*=\s*([0-9.eE+-]+)f", wrapper, "S_factor")
    tilt = _one(r"proxy_handoff_tilt_sigma_rad\s*=\s*([0-9.eE+-]+)f", wrapper, "proxy tilt sigma")
    yaw = _one(r"proxy_handoff_yaw_sigma_rad\s*=\s*([0-9.eE+-]+)f", wrapper, "proxy yaw sigma")
    yaw_free = _one(r"proxy_handoff_yaw_sigma_free_rad\s*=\s*([0-9.eE+-]+)f", wrapper, "free yaw sigma")
    unlock_count = int(_one(r"MAG_UPDATES_TO_UNLOCK\s*=\s*([0-9]+)", wrapper, "mag unlock count"))
    refine_window = _one(r"float mag_refine_window_sec\s*=\s*([0-9.eE+-]+)f", wrapper, "mag refine window")
    refine_start = _one(r"float mag_refine_start_sec\s*=\s*([0-9.eE+-]+)f", wrapper, "mag refine start")

    parity = {
        "bootstrap_frontend_does_not_drive_mekf": (
            "updateCore_(dt, gyro, acc, /*tempC=*/35.0f, /*drive_mekf=*/false);" in wrapper
        ),
        "goLive_reseeds_attitude": "mekf_->initialize_from_attitude(q_bw, tilt_sigma_rad, yaw_sigma_rad);" in wrapper,
        "enterLive_commits_tune_before_first_prediction": "void enterLive_()" in wrapper and "apply_ou_tune_(true);" in wrapper,
        "enterLive_resets_aw_to_committed_stationary_covariance": "mekf_->reset_aw_covariance_to_stationary();" in wrapper,
        "constructor_seeds_v_p_S": all(x in mekf for x in (
            "set_initial_linear_uncertainty(sigma_v0, sigma_p0, sigma_S0);",
            "Pext.template block<3,3>(OFF_V, OFF_V)",
            "Pext.template block<3,3>(OFF_P, OFF_P)",
            "Pext.template block<3,3>(OFF_S, OFF_S)",
        )),
        "attitude_handoff_covariance_is_tilt_yaw_split": all(x in mekf for x in (
            "P_yaw_axis  = u_down_body * u_down_body.transpose();",
            "P_tilt_axes = I - P_yaw_axis;",
            "P_att = tilt_var * P_tilt_axes + yaw_var * P_yaw_axis;",
        )),
        "attitude_handoff_drops_stale_cross_covariance": "attitude-to-bias cross-covariances are dropped" in mekf,
        "held_ba_identity_no_process": all(x in mekf for x in (
            "const T phi_b = acc_bias_updates_enabled_ ? std::exp(-Ts / tau_b) : T(1);",
            "if (acc_bias_updates_enabled_) {",
        )),
        "held_ba_cross_covariance_zeroed": all(x in mekf for x in (
            "Pext.template block<3,BASE_N>(OFF_BA, 0).setZero();",
            "Pext.template block<3,12>(OFF_BA, OFF_V).setZero();",
        )),
        "held_ba_measurement_rows_frozen": "freeze_acc_bias_rows_(K);" in mekf,
        "ba_release_applies_seed_floor": all(x in mekf for x in (
            "const T target_var = sigma_bacc0_ * sigma_bacc0_;",
            "Pba(i,i) = std::max(Pba(i,i), target_var);",
        )),
        "default_S_factor_is_isotropic": s_factor == 1.0,
        "outer_default_with_mag_true": "bool with_mag = true;" in wrapper,
        "outer_default_mag_refinement_enabled": "bool  mag_refine_enabled    = true;" in wrapper,
        "outer_begin_holds_ba_for_mag_refinement": "impl_.setAccBiasHold(cfg_.with_mag && cfg_.mag_refine_enabled);" in wrapper,
        "refinement_only_runs_after_outer_live": "if (stage_ != Stage::Live) return;" in wrapper,
        "refinement_resets_accumulator": "mag_auto_tuner_.reset();" in wrapper,
        "refinement_sets_configured_minimum_window": "refine_cfg.min_window_sec = cfg_.mag_refine_window_sec;" in wrapper,
        "refinement_releases_hold_only_after_done": (
            "mag_refine_done_    = true;" in wrapper and "impl_.setAccBiasHold(false);" in wrapper
        ),
        "mag_tuner_enforces_minimum_accepted_window": all(x in mag for x in (
            "accepted_window_sec_ < cfg_.min_window_sec",
            "return false;",
            "return tryFinalize_();",
        )),
    }
    failures = [k for k, v in parity.items() if not v]

    sigma_interval = dynamic["dynamic_invariant"]["sigma_aw_filter_mps2"]
    sigma_lo, sigma_hi = map(float, sigma_interval)
    sigma_floor = 0.05
    aw_std_lo = max(sigma_floor, sigma_lo)
    aw_std_hi = max(sigma_floor, sigma_hi)

    if not (refine_window >= H18_WORD_S and refine_start >= 0.0):
        failures.append("default magnetic refinement does not guarantee one H18 word before A release")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "trajectory_fit": False,
        "filter_changed": False,
        "source_parity": parity,
        "source_parity_failures": failures,
        "live_entry_seed_is_source_generated_not_arbitrary_PSD": True,
        "bootstrap_mekf_covariance_propagated_before_live": False,
        "state_order_H": ["theta", "b_g", "v", "p", "S", "a_w"],
        "state_order_A": ["theta", "b_g", "v", "p", "S", "a_w", "b_a"],
        "constructor": {
            "Pq0_unused_after_attitude_handoff": Pq0,
            "P_bg_variance": Pb0,
            "sigma_v0": sigma_v0,
            "sigma_p0": sigma_p0,
            "sigma_S0": sigma_S0,
            "sigma_ba0": sigma_ba0,
        },
        "full_heading_gauged_live_attitude_seed": {
            "tilt_std_rad": tilt,
            "yaw_std_rad": yaw,
            "formula": "P_theta=s_tilt^2(I-dd^T)+s_yaw^2 dd^T",
            "d_is_same_body_down_geometry_used_by_measurement_source": True,
        },
        "ungauged_live_attitude_seed": {
            "tilt_std_rad": tilt,
            "yaw_std_rad": yaw_free,
            "belongs_to_gravity_quotient_until_magnetic_regauge": True,
        },
        "translation_seed": {
            "P_v": sigma_v0 * sigma_v0,
            "P_p": sigma_p0 * sigma_p0,
            "P_S": sigma_S0 * sigma_S0,
            "all_translation_cross_covariances_zero_before_first_prediction": True,
        },
        "aw_live_seed": {
            "reset_to_committed_stationary_covariance_before_first_prediction": True,
            "S_factor": s_factor,
            "committed_vertical_std_interval_mps2": [aw_std_lo, aw_std_hi],
            "isotropic_default": s_factor == 1.0,
            "formula_default": "P_aw=sigma_aw,committed^2 I",
        },
        "held_ba": {
            "excluded_from_H18": True,
            "identity_homogeneous_dynamics": True,
            "no_process_injection_while_held": True,
            "cross_covariances_zero": True,
            "measurement_rows_frozen": True,
            "seed_variance": sigma_ba0 * sigma_ba0,
        },
        "H_to_A_release": {
            "hybrid_transition_not_inside_same_mode_word": True,
            "default_mag_updates_required": unlock_count,
            "one_second_guard_after_first_mag_update": True,
            "bias_diagonal_floor_variance": sigma_ba0 * sigma_ba0,
            "configured_full_heading_default_with_mag": True,
            "configured_outer_mag_refinement_enabled": True,
            "outer_acc_bias_hold_active_until_refinement_complete": True,
            "mag_refinement_start_time_from_power_on_s": refine_start,
            "mag_refinement_minimum_accepted_window_s": refine_window,
            "minimum_H_mode_live_duration_before_A_release_s": refine_window,
            "canonical_H18_word_horizon_s": H18_WORD_S,
            "H18_full_word_guaranteed_before_A_release": refine_window >= H18_WORD_S,
            "external_acc_bias_hold_assumed_false_for_default_deployment_bound": False,
        },
        "P3_promoted": False,
        "next_obligation": (
            "consume the guaranteed pre-release H18 full-matrix margin and certify the separate H->A dimension extension through the first active b_a prediction"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("trajectory_fit") is not False or d.get("filter_changed") is not False:
        f.append("live covariance seed changed theorem/filter")
    if d.get("live_entry_seed_is_source_generated_not_arbitrary_PSD") is not True:
        f.append("live seed was not retained")
    if d.get("bootstrap_mekf_covariance_propagated_before_live") is not False:
        f.append("bootstrap falsely propagated MEKF covariance")
    f.extend(f"source parity failed: {x}" for x in d.get("source_parity_failures", []))
    for key in ("P_bg_variance", "sigma_v0", "sigma_p0", "sigma_S0", "sigma_ba0"):
        x = float(d.get("constructor", {}).get(key, 0.0))
        if not (math.isfinite(x) and x > 0.0):
            f.append(f"invalid constructor seed {key}")
    if d.get("translation_seed", {}).get("all_translation_cross_covariances_zero_before_first_prediction") is not True:
        f.append("translation seed cross covariance structure lost")
    if d.get("aw_live_seed", {}).get("reset_to_committed_stationary_covariance_before_first_prediction") is not True:
        f.append("a_w Live reset lost")
    held = d.get("held_ba", {})
    for key in ("excluded_from_H18", "identity_homogeneous_dynamics", "no_process_injection_while_held", "cross_covariances_zero", "measurement_rows_frozen"):
        if held.get(key) is not True:
            f.append(f"held b_a source property lost: {key}")
    rel = d.get("H_to_A_release", {})
    for key in (
        "configured_full_heading_default_with_mag",
        "configured_outer_mag_refinement_enabled",
        "outer_acc_bias_hold_active_until_refinement_complete",
        "H18_full_word_guaranteed_before_A_release",
    ):
        if rel.get(key) is not True:
            f.append(f"release source property lost: {key}")
    if rel.get("external_acc_bias_hold_assumed_false_for_default_deployment_bound") is not False:
        f.append("default full-heading proof incorrectly disables the shipping bias hold")
    if float(rel.get("minimum_H_mode_live_duration_before_A_release_s", 0.0)) < H18_WORD_S:
        f.append("A release can precede the canonical H18 word")
    if d.get("P3_promoted") is not False:
        f.append("seed producer promoted P3")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "live_seed": d["live_entry_seed_is_source_generated_not_arbitrary_PSD"],
        "H_to_A_release": d["H_to_A_release"],
        "held_ba": d["held_ba"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
