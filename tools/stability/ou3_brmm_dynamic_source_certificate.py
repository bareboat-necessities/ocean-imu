#!/usr/bin/env python3
"""Canonical BRMM dynamic adaptive-source certificate for OU-III.

The canonical proof no longer turns the tuner into an 800-state arbitrary
switching language.  Its source state is the shipping adaptive state plus the
effective OU scale consumed by the filter.

A crucial distinction is retained explicitly:

* ``sigma_target_raw`` / ``sigma_candidate_raw`` are the shipping tuner values.
  Once the variance EMA is ready, quiet input may drive the target down to
  ``sigma_coeff*sqrt(1e-6)``; there is no 0.05 floor on this raw state.
* ``sigma_aw_filter`` is the effective OU stationary standard deviation after
  ``apply_ou_tune_`` applies ``max(0.05, band_noise_floor, tune_.sigma_applied)``.

Conflating those coordinates would exclude the BRMM quiet continuation and is
therefore not admissible for P4.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from ou3_interval import Interval
import ou3_source_domain_contract as SOURCE
import ou3_brmm_physical_admissibility as PHYSICAL
import ou3_brmm_wave_period_frontend as FRONTEND
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[2]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
ESTIMATOR = REPO / "src" / "tuner" / "WavePeriodEstimator.h"
TUNER = REPO / "src" / "tuner" / "SeaStateAutoTuner.h"
LIMITS = REPO / "src" / "tuner" / "SeaStateAdaptationLimits.h"
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 3
QUALIFICATION = "OU3_BRMM_DYNAMIC_ADAPTIVE_SOURCE_CERTIFICATE"


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _point(x: float) -> Interval:
    return Interval.point(float(x))


def _member_float(text: str, name: str) -> float:
    m = re.search(rf"\b{name}\s*=\s*([0-9.eE+-]+)f\b", text)
    if not m:
        raise RuntimeError(f"cannot extract deployed member {name}")
    return float(m.group(1))


def _clamp(lo: float, hi: float, a: float, b: float) -> tuple[float, float]:
    if not (a <= b and lo <= hi):
        raise ValueError("invalid clamp interval")
    return max(a, min(b, lo)), max(a, min(b, hi))


def _alpha_upper(dt: float, horizon_lower: float) -> float:
    x = _point(dt) / _point(horizon_lower)
    e = VT.exp_interval(-x)
    return up(1.0 - e.lo)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("BRMM dynamic source may not be trajectory fitted")

    wrapper = WRAPPER.read_text(encoding="utf-8")
    estimator = ESTIMATOR.read_text(encoding="utf-8")
    tuner = TUNER.read_text(encoding="utf-8")
    limits = LIMITS.read_text(encoding="utf-8")

    physical = PHYSICAL.build(path)
    pf = PHYSICAL.validate(physical)
    frontend = FRONTEND.build(REPO)
    ff = FRONTEND.validate(frontend)
    if pf or ff:
        raise RuntimeError(f"BRMM prerequisites failed: physical={pf}, frontend={ff}")

    names = (
        "MIN_TUNE_FREQ_HZ", "MAX_TUNE_FREQ_HZ", "MIN_TAU_S", "MAX_TAU_S",
        "MAX_SIGMA_A", "MIN_R_S", "MAX_R_S", "ADAPT_TAU_SEA_PERIODS",
        "ADAPT_RS_MULT", "ADAPT_EVERY_SECS", "PSEUDO_UPDATE_TAU_RATIO_DEFAULT",
        "PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT", "PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT",
        "FREQ_SMOOTHER_DT",
    )
    c = {name: SOURCE.parse_const(wrapper, name) for name in names}
    hmin = SOURCE.parse_const(limits, "kDynamicEmaHorizonMinSec")
    hmax = SOURCE.parse_const(limits, "kDynamicEmaHorizonMaxSec")
    smin = SOURCE.parse_const(limits, "kDynamicEmaTimeScaleMinSec")
    smax = SOURCE.parse_const(limits, "kDynamicEmaTimeScaleMaxSec")

    tau_coeff = _member_float(wrapper, "tau_coeff_")
    sigma_coeff = _member_float(wrapper, "sigma_coeff_")
    tau_initial = _member_float(wrapper, "tau_applied")
    sigma_initial = _member_float(wrapper, "sigma_applied")
    rs_initial = _member_float(wrapper, "RS_applied")

    f_lo, f_hi = c["MIN_TUNE_FREQ_HZ"], c["MAX_TUNE_FREQ_HZ"]
    tau_target_lo, tau_target_hi = _clamp(
        tau_coeff * 0.5 / f_hi, tau_coeff * 0.5 / f_lo,
        c["MIN_TAU_S"], c["MAX_TAU_S"],
    )
    tau_applied_lo = min(tau_initial, tau_target_lo)
    tau_applied_hi = max(tau_initial, tau_target_hi)

    # After variance readiness shipping executes
    #   var_wave=max(var_total-var_noise,0); var_wave=max(var_wave,1e-6)
    #   sigma_target=min(sigma_coeff*sqrt(var_wave),MAX_SIGMA_A).
    # Therefore quiet continuation permits 0.001*sigma_coeff.  The separate
    # 0.05 floor appears only later in apply_ou_tune_ when the filter consumes
    # the candidate.
    sigma_target_lo = min(c["MAX_SIGMA_A"], max(0.0, sigma_coeff) * 1.0e-3)
    sigma_target_hi = c["MAX_SIGMA_A"]
    sigma_candidate_lo = min(sigma_initial, sigma_target_lo)
    sigma_candidate_hi = max(sigma_initial, sigma_target_hi)
    sigma_filter_lo = 0.05
    sigma_filter_hi = max(0.05, sigma_candidate_hi)

    rs_lo = min(rs_initial, c["MIN_R_S"])
    rs_hi = max(rs_initial, c["MAX_R_S"])

    sea_time_lo, sea_time_hi = _clamp(0.5 / f_hi, 0.5 / f_lo, smin, smax)
    common_h_lo, common_h_hi = _clamp(
        c["ADAPT_TAU_SEA_PERIODS"] * sea_time_lo,
        c["ADAPT_TAU_SEA_PERIODS"] * sea_time_hi,
        hmin, hmax,
    )
    tau_scale_lo, tau_scale_hi = _clamp(tau_target_lo, tau_target_hi, smin, smax)
    rs_h_lo, rs_h_hi = _clamp(
        c["ADAPT_RS_MULT"] * tau_scale_lo,
        c["ADAPT_RS_MULT"] * tau_scale_hi,
        hmin, hmax,
    )

    dt = c["FREQ_SMOOTHER_DT"]
    alpha_common = _alpha_upper(dt, common_h_lo)
    alpha_rs = _alpha_upper(dt, rs_h_lo)
    per_sample_tau = up(alpha_common * (tau_target_hi - tau_target_lo))
    per_sample_sigma = up(alpha_common * (sigma_target_hi - sigma_target_lo))
    per_sample_rs = up(alpha_rs * (rs_hi - rs_lo))

    commit_gap_samples = int(math.ceil(c["ADAPT_EVERY_SECS"] / dt)) + 2
    commit_gap_s_upper = up(commit_gap_samples * dt)
    tau_commit_jump = up(commit_gap_samples * per_sample_tau)
    sigma_commit_jump = up(commit_gap_samples * per_sample_sigma)
    rs_commit_jump = up(commit_gap_samples * per_sample_rs)
    pseudo_commit_jump = up(c["PSEUDO_UPDATE_TAU_RATIO_DEFAULT"] * tau_commit_jump)

    parity = {
        "single_existing_wave_period_estimator_used": bool(
            re.search(r"\bWavePeriodEstimator\s+wave_period_\s*;", wrapper)
        ),
        "normal_live_measured_period_selector_present": "wave_period_.hasUsablePeriod()" in wrapper,
        "tau_sigma_candidates_smoothed_each_valid_sample": (
            "tune_.tau_applied   += alpha" in wrapper and
            "tune_.sigma_applied += alpha" in wrapper
        ),
        "RS_candidate_smoothed_each_valid_sample": "tune_.RS_applied    += alpha_RS" in wrapper,
        "active_schedule_commit_is_next_sample_predictable": (
            "online_tune_apply_pending_ = true" in wrapper and
            "void apply_pending_online_tune_()" in wrapper
        ),
        "pseudo_cadence_is_same_tau_lipschitz_image": (
            "const float requested = pseudo_update_tau_ratio_ * tau" in wrapper
        ),
        "dynamic_horizon_ceiling_35s_present": "kDynamicEmaHorizonMaxSec = 35.0f" in limits,
        "one_way_usable_period_latch_present": "if (usable_period_) return;" in estimator,
        "quiet_variance_may_reach_zero_after_ready": (
            "return std::max(0.0f, A_sq.get() - mu * mu);" in tuner
        ),
        "raw_sigma_target_uses_1e_minus6_variance_floor": (
            "var_wave = std::max(var_wave, 1e-6f);" in wrapper
            and "sigma_wave = std::sqrt(var_wave);" in wrapper
        ),
        "effective_filter_sigma_has_separate_0p05_floor": (
            "const float sigma_floor = std::max(0.05f, band_noise_floor_sigma_());" in wrapper
            and "const float sZ = std::max(sigma_floor, tune_.sigma_applied);" in wrapper
        ),
    }
    parity_failures = [name for name, ok in parity.items() if not ok]

    live = domain["normal_live"]
    status = "PASS" if not parity_failures else "FAIL"
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "theorem_is_conditional_on_admitted_BRMM_event": True,
        "deterministic_finite_window_BRMM_realization_closed_here": False,
        "old_P2_800_state_graph_consumed": False,
        "old_P2_history_word_enumeration_consumed": False,
        "adaptive_state": [
            "tau_candidate_active", "sigma_candidate_raw", "R_S_candidate_active",
            "sigma_aw_filter_effective", "pseudo_update_period", "pending_commit_progress",
        ],
        "raw_and_effective_sigma_are_distinct_coordinates": True,
        "source_parity": parity,
        "source_parity_failures": parity_failures,
        "reference_model_metadata": {
            "required_for_shipping_clamp_invariant": False,
            "sea_modes_max": physical["sea_modes_max"],
            "total_Hs_upper_m": physical["repository_total_Hs_upper_m"],
            "height_period_coupling_retained": physical["three_partition_contract"][
                "independent_H_r_and_T_p_rectangular_extrema_forbidden"
            ],
            "frontend_source_parity_pass": all(frontend["source_parity"].values()),
            "frontend_full_SEA0_promoted": frontend["SEA0_full_certificate_promoted"],
        },
        "normal_live_contract": {
            "accelerometer_update_required_each_valid_sample": live[
                "accelerometer_update_required_each_valid_imu_sample_after_live_entry"
            ],
            "accelerometer_rejection_in_scope": live[
                "accelerometer_rejection_in_normal_live_scope"
            ],
            "vector_PE_recurrence_window_s": live["vector_pe_recurrence_window_s"],
            "quiet_case_admitted": True,
        },
        "dynamic_invariant": {
            "tuning_frequency_hz": [down(f_lo), up(f_hi)],
            "tau_target_s": [down(tau_target_lo), up(tau_target_hi)],
            "tau_applied_s": [down(tau_applied_lo), up(tau_applied_hi)],
            "sigma_target_raw_mps2": [down(sigma_target_lo), up(sigma_target_hi)],
            "sigma_candidate_raw_mps2": [down(sigma_candidate_lo), up(sigma_candidate_hi)],
            "sigma_aw_filter_mps2": [down(sigma_filter_lo), up(sigma_filter_hi)],
            "R_S_applied": [down(rs_lo), up(rs_hi)],
            "pseudo_update_period_s": [
                down(c["PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT"]),
                up(c["PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT"]),
            ],
            "common_tau_sigma_horizon_s": [down(common_h_lo), up(common_h_hi)],
            "R_S_horizon_s": [down(rs_h_lo), up(rs_h_hi)],
        },
        "validated_rate_and_jump_bounds": {
            "dt_s": dt,
            "tau_sigma_alpha_per_sample_upper": alpha_common,
            "R_S_alpha_per_sample_upper": alpha_rs,
            "tau_candidate_abs_delta_per_sample_upper_s": per_sample_tau,
            "sigma_candidate_abs_delta_per_sample_upper_mps2": per_sample_sigma,
            "R_S_candidate_abs_delta_per_sample_upper": per_sample_rs,
            "active_commit_gap_samples_upper": commit_gap_samples,
            "active_commit_gap_s_upper": commit_gap_s_upper,
            "tau_active_abs_jump_per_commit_upper_s": tau_commit_jump,
            "sigma_active_abs_jump_per_commit_upper_mps2": sigma_commit_jump,
            "R_S_active_abs_jump_per_commit_upper": rs_commit_jump,
            "pseudo_period_active_abs_jump_per_commit_upper_s": pseudo_commit_jump,
            "target_jump_may_span_full_box": True,
            "proof_relies_on_implemented_smoothing_not_unproved_sea_parameter_derivatives": True,
        },
        "P2_DYNAMIC_SOURCE_CERTIFICATE": status,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "next_obligation": (
            "propagate raw tuner candidates and the separately floored effective OU sigma through the same-history P4 source cover; do not replace them by one sigma coordinate"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit", "theorem_is_conditional_on_admitted_BRMM_event",
        "raw_and_effective_sigma_are_distinct_coordinates",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "deterministic_finite_window_BRMM_realization_closed_here",
        "old_P2_800_state_graph_consumed", "old_P2_history_word_enumeration_consumed",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    for name in d.get("source_parity_failures", []):
        f.append(f"shipping source parity failed: {name}")
    if d.get("P2_DYNAMIC_SOURCE_CERTIFICATE") != "PASS":
        f.append("dynamic source certificate did not pass")
    live = d.get("normal_live_contract", {})
    if live.get("accelerometer_update_required_each_valid_sample") is not True:
        f.append("normal Live accelerometer recurrence lost")
    if live.get("accelerometer_rejection_in_scope") is not False:
        f.append("rejected accelerometer branch re-entered theorem")
    if live.get("quiet_case_admitted") is not True:
        f.append("quiet BRMM continuation was lost")
    inv=d.get("dynamic_invariant",{})
    for name, bounds in inv.items():
        if isinstance(bounds, list) and len(bounds) == 2:
            lo, hi = map(float, bounds)
            if not (math.isfinite(lo) and math.isfinite(hi) and 0.0 < lo <= hi):
                f.append(f"invalid dynamic interval {name}")
    raw=inv.get("sigma_candidate_raw_mps2",[math.nan,math.nan])
    eff=inv.get("sigma_aw_filter_mps2",[math.nan,math.nan])
    if not (float(raw[0]) < 0.05 <= float(eff[0])):
        f.append("raw/effective sigma lower-bound distinction missing")
    if d.get("P3_PROMOTED") is not False or d.get("P4_PROMOTED") is not False:
        f.append("dynamic source stage promoted downstream theorem")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P2_DYNAMIC_SOURCE_CERTIFICATE"],
        "source_parity_failures": d["source_parity_failures"],
        "invariant": d["dynamic_invariant"],
        "rate_bounds": d["validated_rate_and_jump_bounds"],
        "failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
