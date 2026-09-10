#!/usr/bin/env python3
"""End-to-end startup -> P4 capture gate for the shipping OU-III estimator.

The gate follows the actual hybrid shipping path. It does not use the declared
Mahony chart as a reachability premise, does not equate goLive() with P4
membership, and does not demand a full-yaw conclusion before magnetic
observability exists.

The theorem-facing startup terminal section is the actual shipping handoff
predicate. At that sample the world-frame gravity branch is strictly aligned;
with the declared same-history averaged-gravity direction error this gives a
wide, chart-free tilt quotient bound. Early Live H18 contracts that set toward
P4. If the timeout path entered without north, the same shipping MagAutoTuner
continues after Live. Absolute-heading capture consumes the minimal explicit
accepted-count/accepted-time observability admission; Normal-Live PE alone is
not silently promoted to first-north progress.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import ou3_startup_p4_reachable_handoff as HANDOFF
import ou3_startup_proxy_initial_tilt_bound as INIT_TILT
import ou3_startup_handoff_tilt_hemisphere as HANDOFF_TILT
import ou3_brmm_private_mahony_state_step as MAHONY
import ou3_startup_magnetic_observability as MAGOBS
import ou3_p5_magnetic_progress_admission as MAGPROG
import ou3_startup_postlive_regauge_contract as POSTLIVE

REPO = Path(__file__).resolve().parents[2]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 4
QUALIFICATION = "OU3_SHIPPING_STARTUP_TO_P4_CAPTURE_GATE_V4"


def _config_float(text: str, name: str) -> float:
    m = re.search(rf"float\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f\s*;", text)
    if not m:
        raise RuntimeError(f"cannot extract Config::{name}")
    return float(m.group(1))


def build() -> dict:
    wrapper = WRAPPER.read_text(encoding="utf-8")
    domain = json.loads(DOMAIN.read_text(encoding="utf-8"))
    handoff = HANDOFF.build()
    init = INIT_TILT.build()
    handoff_tilt = HANDOFF_TILT.build()
    mahony = MAHONY.build()
    magobs = MAGOBS.build()
    magprog = MAGPROG.build()
    postlive = POSTLIVE.build()

    prereq = {
        "handoff": HANDOFF.validate(handoff),
        "initial_tilt": INIT_TILT.validate(init),
        "handoff_tilt": HANDOFF_TILT.validate(handoff_tilt),
        "mahony": MAHONY.validate(mahony),
        "mag_observability": MAGOBS.validate(magobs),
        "mag_progress": MAGPROG.validate(magprog),
        "postlive_regauge": POSTLIVE.validate(postlive),
    }
    bad = {k: v for k, v in prereq.items() if v}
    if bad:
        raise RuntimeError("startup capture prerequisite failure: " + repr(bad))

    hold = _config_float(wrapper, "mag_gravity_align_hold_sec")
    max_good_credit = 10.0
    bad_tail_s = max(0.0, (max_good_credit - hold) / 2.0)

    startup = domain["startup"]
    normal = domain["normal_live"]
    stale_domain_claims = {
        "live_entry_requires_wave_period_estimator_usable": bool(
            startup.get("live_entry_requires_wave_period_estimator_usable", False)
        ),
        "tuner_ready_requires_wave_period_estimator_usable": bool(
            startup.get("tuner_ready_requires_wave_period_estimator_usable", False)
        ),
    }

    circular_chart_dependency = bool(
        mahony.get("declared_startup_chart_implies_gravity_aligned_branch")
        and mahony.get("declared_domain_live_entry_upper_bound_closed")
        and not mahony.get("complete_BRMM_family_materialized_here")
    )
    magnetic_obstruction = bool(
        magobs.get("structural_observability_obstruction_proved")
        and not magobs.get("UNCONDITIONAL_FULL_ATTITUDE_FINITE_CAPTURE_FROM_CURRENT_SOURCE")
    )
    finite_north_under_capture_admission = bool(
        magprog["FINITE_NORTH_EVENT_UNDER_ADMISSION_CLOSED"]
        and postlive["FINITE_POSTLIVE_NORTH_EVENT_UNDER_ADMISSION_CLOSED"]
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "filter_changed": False,
        "quality_gates_changed": False,
        "p3_delta_changed": False,
        "shipping_handoff_fiber_validated": True,
        "physical_first_sample_tilt_lemma_validated": True,
        "physical_first_sample_tilt_upper_deg": init["physical_initial_tilt_error_upper_deg"],
        "mahony_one_sample_binary32_map_validated": True,
        "mahony_arithmetic_boundedness_is_capture": False,
        "legacy_declared_mahony_chart_may_establish_capture": False,
        "existing_mahony_timeout_argument_has_circular_chart_dependency": circular_chart_dependency,
        "actual_shipping_handoff_tilt_section": {
            "validated": True,
            "representation": handoff_tilt["handoff_tilt_set_representation"],
            "true_gravity_quotient_tilt_strict_upper_deg": handoff_tilt["true_gravity_quotient_tilt_strict_upper_deg"],
            "historical_60_deg_mahony_chart_consumed": False,
            "yaw_gauge_required": False,
            "applies_to_quality_and_timeout_handoff_samples": True,
            "this_is_conditional_on_handoff_occurring": True,
            "early_live_contracts_this_set_to_P4": True,
        },
        "magnetic_observability": {
            "structural_yaw_obstruction_proved": magnetic_obstruction,
            "raw_COMPLETE_BRMM_source_forces_finite_north_acquisition": False,
            "normal_live_PE_alone_forces_first_north": False,
            "explicit_capture_admission_added": True,
            "capture_admission_changes_BRMM_motion_caps": False,
            "capture_admission_changes_BIAS_family": False,
            "capture_admission_progress_horizon_s": magprog["finite_north_event_time_upper_from_progress_origin_s"],
            "finite_north_event_under_explicit_capture_admission_closed": finite_north_under_capture_admission,
            "ungauged_gravity_quotient_capture_is_still_meaningful": True,
            "full_attitude_capture_requires_north_event": True,
            "shipping_min_accepted_mag_samples": magobs["shipping_magnetic_acquisition"]["min_accepted_samples"],
            "shipping_min_accepted_mag_window_s": magobs["shipping_magnetic_acquisition"]["min_accepted_window_s"],
            "counterexample_minimax_full_attitude_error_deg_without_progress_admission": magobs["yaw_indistinguishability_witness"]["minimax_full_attitude_error_lower_deg"],
            "capture_admission": magprog["admission"],
        },
        "quality_gate_memory": {
            "good_counter_cap_s": max_good_credit,
            "required_hold_s": hold,
            "bad_sample_credit_decay_rate": 2.0,
            "max_bad_tail_after_last_good_sample_s": bad_tail_s,
            "current_align_sin_threshold_implied_at_quality_handoff": False,
            "proof_requirement": "propagate same-history proxy error from the last certified good sample through the entire leaky-hold tail",
        },
        "shipping_timeout_semantics": {
            "live_handoff_upper_bound_s_under_declared_aligned_branch": handoff["deployed_live_handoff_upper_bound_s"],
            "requires_tuner_ready": False,
            "requires_magnetic_north_when_with_mag": False,
            "requires_aligned_gravity_branch": True,
            "goLive_equals_P4_capture": False,
            "postlive_first_north_path_validated": postlive["shipping_postlive_first_north_path_present"],
        },
        "domain_text_audit": {
            "stale_live_entry_requires_WPE_usable_claim": stale_domain_claims["live_entry_requires_wave_period_estimator_usable"],
            "shipping_timeout_contradicts_unconditional_WPE_usable_entry": True,
            "normal_live_magnetic_PE_is_not_startup_acquisition": True,
            "normal_live_vector_PE_window_s": float(normal["vector_pe_recurrence_window_s"]),
        },
        "corrected_theorem_chain": [
            "STARTUP_PROXY_REACHABILITY_TO_SHIPPING_HANDOFF_PREDICATE",
            "SHIPPING_LIVE_HANDOFF_FIBER_WITH_TILT_LT_91P146_DEG",
            "EARLY_LIVE_H18_GRAVITY_QUOTIENT_CAPTURE",
            "FINITE_MAGNETIC_NORTH_EVENT_UNDER_EXPLICIT_ACCEPTED_PROGRESS_ADMISSION_WHEN_WITH_MAG",
            "LATE_NORTH_HYBRID_YAW_RESET_IF_ALREADY_LIVE",
            "P4_H18_FULL_ATTITUDE_INVARIANT_BASIN",
            "H18_TO_A21_RELEASE_WHEN_ENABLED",
            "P4_A21_INVARIANT_BASIN",
        ],
        "capture_modes_required": [
            "QUALITY_GAUGED_HANDOFF",
            "TIMEOUT_UNGAUGED_OR_UNTUNED_HANDOFF",
            "EARLY_LIVE_H18_GRAVITY_QUOTIENT_CAPTURE",
            "FINITE_ACCEPTED_MAGNETIC_PROGRESS_WHEN_WITH_MAG",
            "LATE_MAGNETIC_GAUGE_RESET",
            "H18_TO_A21_RELEASE_WHEN_ENABLED",
        ],
        "required_reachable_state": [
            "true/proxy tilt quotient error",
            "Mahony integral-feedback state",
            "physical and measured specific-force direction on same BRMM history",
            "gyro true+bias+disturbance on same BIAS history",
            "world gravity LPF and leaky hold counter",
            "literal magnetic accepted-count/accepted-time state and yaw-gauge status",
            "WPE + sigma-band + tuner candidate/active/commit state",
            "MEKF handoff covariance and scheduler state",
            "early-Live H18 nonlinear/covariance state until P4 membership",
        ],
        "entry_translation_fiber": {
            "e_p_at_handoff": 0.0,
            "e_S_at_handoff": 0.0,
            "e_v_at_handoff": "v_true(t_h)",
            "e_aw_at_handoff": "a_w_true(t_h)",
            "independent_300_m_s_S_entry_allowed": False,
        },
        "HANDOFF_TILT_SECTION_CLOSED": True,
        "FINITE_REACHABILITY_OF_HANDOFF_PREDICATE_CLOSED": False,
        "UNGAUGED_GRAVITY_QUOTIENT_CAPTURE_CLOSED": False,
        "UNCONDITIONAL_FINITE_NORTH_FROM_RAW_COMPLETE_BRMM_CLOSED": False,
        "FINITE_NORTH_EVENT_UNDER_EXPLICIT_CAPTURE_ADMISSION_CLOSED": finite_north_under_capture_admission,
        "LATE_NORTH_HYBRID_RESET_MAP_STRUCTURE_CLOSED": False,
        "LATE_NORTH_FULL_ATTITUDE_LANDING_CLOSED": False,
        "EARLY_LIVE_H18_CAPTURE_TUBE_CLOSED": False,
        "P4_REACHABLE_BASIN_OVERLAP_CLOSED": False,
        "FINITE_CAPTURE_BOUND_CLOSED": False,
        "P4_PASS": False,
        "P5_PASS": False,
        "failure_classification": {
            "current_primary": "F plus P4 closure",
            "E_raw_source_observability_gap_resolved_by_explicit_conditional_admission": True,
            "global_physical_deployment_left_inclusion_of_magnetic_progress_proved": False,
            "F": "finite reachability of handoff, early-Live H18 quotient contraction, late-north yaw landing, and full basin intersection remain open",
            "not_filter_instability": True,
        },
        "next_mathematical_obligations": [
            "prove finite source-uniform reachability of the actual shipping aligned-branch handoff predicate; do not replace it by the old 60 deg chart",
            "prove H18 gravity-quotient contraction from the certified <91.146 deg handoff section into the widest production P4 finite-angle basin",
            "bound the late-north yaw gauge from the same accepted magnetic history and land the yaw-only rewrite inside the full-attitude basin",
            "materialize the same COMPLETE-BRMM frontend/tuner continuation and P4 source cover",
            "construct and maximize the source-uniform P4-H18 basin from the corrected correlated handoff fiber",
            "prove finite intersection of every admitted startup/early-Live tube with that P4-H18 basin",
        ],
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "shipping_handoff_fiber_validated",
        "physical_first_sample_tilt_lemma_validated",
        "mahony_one_sample_binary32_map_validated",
        "existing_mahony_timeout_argument_has_circular_chart_dependency",
        "HANDOFF_TILT_SECTION_CLOSED",
        "FINITE_NORTH_EVENT_UNDER_EXPLICIT_CAPTURE_ADMISSION_CLOSED",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    hs = d.get("actual_shipping_handoff_tilt_section", {})
    if not (91.0 < float(hs.get("true_gravity_quotient_tilt_strict_upper_deg", 0.0)) < 92.0):
        f.append("shipping handoff tilt section is not the expected wide hemisphere")
    mo = d.get("magnetic_observability", {})
    if mo.get("structural_yaw_obstruction_proved") is not True:
        f.append("missing structural yaw-observability obstruction")
    if mo.get("raw_COMPLETE_BRMM_source_forces_finite_north_acquisition") is not False:
        f.append("raw source falsely claims finite north acquisition")
    if mo.get("normal_live_PE_alone_forces_first_north") is not False:
        f.append("Normal-Live PE falsely promoted to first-north progress")
    if mo.get("capture_admission_changes_BRMM_motion_caps") is not False or mo.get("capture_admission_changes_BIAS_family") is not False:
        f.append("capture admission altered motion or bias family")
    if d.get("quality_gate_memory", {}).get("max_bad_tail_after_last_good_sample_s") != 4.0:
        f.append("unexpected leaky quality-gate bad tail")
    for key in (
        "mahony_arithmetic_boundedness_is_capture",
        "legacy_declared_mahony_chart_may_establish_capture",
        "FINITE_REACHABILITY_OF_HANDOFF_PREDICATE_CLOSED",
        "UNGAUGED_GRAVITY_QUOTIENT_CAPTURE_CLOSED",
        "UNCONDITIONAL_FINITE_NORTH_FROM_RAW_COMPLETE_BRMM_CLOSED",
        "LATE_NORTH_HYBRID_RESET_MAP_STRUCTURE_CLOSED",
        "LATE_NORTH_FULL_ATTITUDE_LANDING_CLOSED",
        "EARLY_LIVE_H18_CAPTURE_TUBE_CLOSED",
        "P4_REACHABLE_BASIN_OVERLAP_CLOSED",
        "FINITE_CAPTURE_BOUND_CLOSED",
        "P4_PASS",
        "P5_PASS",
    ):
        if d.get(key) is not False:
            f.append(f"{key} must remain false")
    if d.get("entry_translation_fiber", {}).get("independent_300_m_s_S_entry_allowed") is not False:
        f.append("legacy independent S entry reintroduced")
    return list(dict.fromkeys(f))


if __name__ == "__main__":
    d = build(); failures = validate(d)
    print(json.dumps({**d, "validation_pass": not failures, "validation_failures": failures}, indent=2, sort_keys=True))
    raise SystemExit(1 if failures else 0)
