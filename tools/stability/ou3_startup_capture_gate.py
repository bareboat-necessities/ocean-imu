#!/usr/bin/env python3
"""End-to-end startup -> P4 capture gate for the shipping OU-III estimator.

The gate follows the actual hybrid shipping path. It does not use the declared
Mahony chart as a reachability premise, does not equate goLive() with P4
membership, and does not demand a full-yaw conclusion before magnetic
observability exists.

The theorem-facing startup terminal section is now the *actual shipping handoff
predicate*: at the sample that enters Live, the world-frame gravity branch is
strictly aligned. Together with the declared same-history averaged-gravity
direction error this gives a wide, chart-free tilt quotient bound. Early Live
H18, not the private Mahony proxy, is responsible for contracting that set into
the eventual P4 basin.
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

REPO = Path(__file__).resolve().parents[2]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 3
QUALIFICATION = "OU3_SHIPPING_STARTUP_TO_P4_CAPTURE_GATE_V3"


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

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "filter_changed": False,
        "quality_gates_changed": False,
        "p3_delta_changed": False,
        "shipping_handoff_fiber_validated": not HANDOFF.validate(handoff),
        "physical_first_sample_tilt_lemma_validated": not INIT_TILT.validate(init),
        "physical_first_sample_tilt_upper_deg": init["physical_initial_tilt_error_upper_deg"],
        "mahony_one_sample_binary32_map_validated": not MAHONY.validate(mahony),
        "mahony_arithmetic_boundedness_is_capture": False,
        "legacy_declared_mahony_chart_may_establish_capture": False,
        "existing_mahony_timeout_argument_has_circular_chart_dependency": circular_chart_dependency,
        "actual_shipping_handoff_tilt_section": {
            "validated": not HANDOFF_TILT.validate(handoff_tilt),
            "representation": handoff_tilt["handoff_tilt_set_representation"],
            "true_gravity_quotient_tilt_strict_upper_deg": handoff_tilt["true_gravity_quotient_tilt_strict_upper_deg"],
            "historical_60_deg_mahony_chart_consumed": handoff_tilt["historical_60_deg_mahony_chart_consumed"],
            "yaw_gauge_required": handoff_tilt["yaw_gauge_required"],
            "applies_to_quality_and_timeout_handoff_samples": handoff_tilt["applies_to_quality_and_timeout_handoff_samples"],
            "this_is_conditional_on_handoff_occurring": True,
            "early_live_contracts_this_set_to_P4": True,
        },
        "magnetic_observability": {
            "structural_yaw_obstruction_proved": magnetic_obstruction,
            "current_source_forces_finite_north_acquisition": False,
            "ungauged_gravity_quotient_capture_is_still_meaningful": True,
            "full_attitude_capture_requires_north_event": True,
            "shipping_min_accepted_mag_samples": magobs["shipping_magnetic_acquisition"]["min_accepted_samples"],
            "shipping_min_accepted_mag_window_s": magobs["shipping_magnetic_acquisition"]["min_accepted_window_s"],
            "counterexample_minimax_full_attitude_error_deg": magobs["yaw_indistinguishability_witness"]["minimax_full_attitude_error_lower_deg"],
            "weakest_source_extension": magobs["weakest_required_source_extension"],
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
            "MAGNETIC_NORTH_ACQUISITION_EVENT_WHEN_WITH_MAG",
            "LATE_NORTH_HYBRID_YAW_RESET_IF_ALREADY_LIVE",
            "P4_H18_FULL_ATTITUDE_INVARIANT_BASIN",
            "H18_TO_A21_RELEASE_WHEN_ENABLED",
            "P4_A21_INVARIANT_BASIN",
        ],
        "capture_modes_required": [
            "QUALITY_GAUGED_HANDOFF",
            "TIMEOUT_UNGAUGED_OR_UNTUNED_HANDOFF",
            "EARLY_LIVE_H18_GRAVITY_QUOTIENT_CAPTURE",
            "LATE_MAGNETIC_GAUGE_RESET",
            "H18_TO_A21_RELEASE_WHEN_ENABLED",
        ],
        "required_reachable_state": [
            "true/proxy tilt quotient error",
            "Mahony integral-feedback state",
            "physical and measured specific-force direction on same BRMM history",
            "gyro true+bias+disturbance on same BIAS history",
            "world gravity LPF and leaky hold counter",
            "magnetic acquisition state and yaw-gauge status",
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
        "HANDOFF_TILT_SECTION_CLOSED": bool(handoff_tilt["HANDOFF_TILT_HEMISPHERE_BOUND_CLOSED"]),
        "FINITE_REACHABILITY_OF_HANDOFF_PREDICATE_CLOSED": False,
        "UNGAUGED_GRAVITY_QUOTIENT_CAPTURE_CLOSED": False,
        "FINITE_NORTH_ACQUISITION_FROM_CURRENT_SOURCE_CLOSED": False,
        "LATE_NORTH_HYBRID_RESET_CLOSED": False,
        "EARLY_LIVE_H18_CAPTURE_TUBE_CLOSED": False,
        "P4_REACHABLE_BASIN_OVERLAP_CLOSED": False,
        "FINITE_CAPTURE_BOUND_CLOSED": False,
        "P4_PASS": False,
        "P5_PASS": False,
        "failure_classification": {
            "current_primary": "E/F",
            "E": "current COMPLETE-BRMM declaration does not force finite startup magnetic acquisition and its global source-family materialization remains open",
            "F": "finite reachability of the shipping handoff predicate plus early-Live quotient capture and late-north reset into the eventual P4 basin remain to be established",
            "not_filter_instability": True,
        },
        "next_mathematical_obligations": [
            "prove finite source-uniform reachability of the actual shipping aligned-branch handoff predicate; do not replace it by the old 60 deg chart",
            "prove H18 gravity-quotient contraction from the certified <91.146 deg handoff section into the widest production P4 finite-angle basin",
            "add/derive the weakest finite accepted-magnetic acquisition progress condition for with_mag=true",
            "certify the exact late-north shipping yaw-reset map into the full-attitude basin",
            "materialize the same COMPLETE-BRMM frontend/tuner continuation",
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
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    hs = d.get("actual_shipping_handoff_tilt_section", {})
    if hs.get("validated") is not True:
        f.append("shipping handoff tilt section is not validated")
    if not (91.0 < float(hs.get("true_gravity_quotient_tilt_strict_upper_deg", 0.0)) < 92.0):
        f.append("shipping handoff tilt section is not the expected wide hemisphere")
    if hs.get("historical_60_deg_mahony_chart_consumed") is not False or hs.get("yaw_gauge_required") is not False:
        f.append("handoff tilt section reintroduced old chart/yaw")
    mo = d.get("magnetic_observability", {})
    if mo.get("structural_yaw_obstruction_proved") is not True:
        f.append("missing structural yaw-observability obstruction")
    if mo.get("current_source_forces_finite_north_acquisition") is not False:
        f.append("current source falsely claims finite north acquisition")
    if mo.get("counterexample_minimax_full_attitude_error_deg") != 45.0:
        f.append("unexpected yaw counterexample lower bound")
    if d.get("quality_gate_memory", {}).get("max_bad_tail_after_last_good_sample_s") != 4.0:
        f.append("unexpected leaky quality-gate bad tail")
    for key in (
        "mahony_arithmetic_boundedness_is_capture",
        "legacy_declared_mahony_chart_may_establish_capture",
        "FINITE_REACHABILITY_OF_HANDOFF_PREDICATE_CLOSED",
        "UNGAUGED_GRAVITY_QUOTIENT_CAPTURE_CLOSED",
        "FINITE_NORTH_ACQUISITION_FROM_CURRENT_SOURCE_CLOSED",
        "LATE_NORTH_HYBRID_RESET_CLOSED",
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
    return f


if __name__ == "__main__":
    d = build(); failures = validate(d)
    print(json.dumps({**d, "validation_pass": not failures, "validation_failures": failures}, indent=2))
    raise SystemExit(1 if failures else 0)
