#!/usr/bin/env python3
"""End-to-end startup -> P4 capture gate for the shipping OU-III estimator.

This gate intentionally rejects two historical shortcuts:

1. treating the declared Mahony chart as if startup dynamics had proved entry
   into / retention of that chart; and
2. treating ``goLive()`` as synonymous with P4 capture.

The outer wrapper can take its 150 s timeout while north and TunerReady are
still absent.  Thus the theorem is hybrid and the capture interval may extend
into early Live/H18.  P5 can promote only after a SAME-COMPLETE-BRMM reachable
startup tube intersects a source-uniform P4-H18 invariant basin and every
allowed handoff path is covered.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import ou3_startup_p4_reachable_handoff as HANDOFF
import ou3_startup_proxy_initial_tilt_bound as INIT_TILT
import ou3_brmm_private_mahony_state_step as MAHONY

REPO = Path(__file__).resolve().parents[2]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SHIPPING_STARTUP_TO_P4_CAPTURE_GATE_V1"


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
    mahony = MAHONY.build()

    hold = _config_float(wrapper, "mag_gravity_align_hold_sec")
    # good_sec is capped at 10 s and decays at 2 seconds of credit per second
    # of a bad sample.  Therefore a saturated hold can remain >= hold for this
    # long after the last good sample.  This is shipping semantics, not a proof
    # pessimism invented here.
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

    # The current private Mahony module closes arithmetic boundedness but its
    # 150 s branch proof consumes the *declared* 60 deg chart.  It therefore
    # cannot be used as the missing reachable-chart/capture implication.
    circular_chart_dependency = bool(
        mahony.get("declared_startup_chart_implies_gravity_aligned_branch")
        and mahony.get("declared_domain_live_entry_upper_bound_closed")
        and not mahony.get("complete_BRMM_family_materialized_here")
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
        "quality_gate_memory": {
            "good_counter_cap_s": max_good_credit,
            "required_hold_s": hold,
            "bad_sample_credit_decay_rate": 2.0,
            "max_bad_tail_after_last_good_sample_s": bad_tail_s,
            "current_align_sin_threshold_implied_at_quality_handoff": False,
            "proof_requirement": "propagate same-history proxy error from the last certified good sample through the entire leaky-hold tail",
        },
        "shipping_timeout_semantics": {
            "live_handoff_upper_bound_s_under_aligned_branch": handoff["deployed_live_handoff_upper_bound_s"],
            "requires_tuner_ready": False,
            "requires_magnetic_north_when_with_mag": False,
            "requires_aligned_gravity_branch": True,
            "goLive_equals_P4_capture": False,
        },
        "domain_text_audit": {
            "stale_live_entry_requires_WPE_usable_claim": stale_domain_claims["live_entry_requires_wave_period_estimator_usable"],
            "shipping_timeout_contradicts_unconditional_WPE_usable_entry": True,
            "normal_live_magnetic_PE_is_not_yet_a_startup_PE_proof": True,
            "normal_live_vector_PE_window_s": float(normal["vector_pe_recurrence_window_s"]),
        },
        "capture_modes_required": [
            "QUALITY_GAUGED_HANDOFF",
            "TIMEOUT_UNGAUGED_OR_UNTUNED_HANDOFF",
            "EARLY_LIVE_H18_TO_P4_H18",
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
        "STARTUP_REACHABLE_TUBE_CLOSED": False,
        "EARLY_LIVE_H18_CAPTURE_TUBE_CLOSED": False,
        "P4_REACHABLE_BASIN_OVERLAP_CLOSED": False,
        "FINITE_CAPTURE_BOUND_CLOSED": False,
        "P4_PASS": False,
        "P5_PASS": False,
        "failure_classification": {
            "current_primary": "E/F",
            "E": "COMPLETE-BRMM startup source/observability qualification and materialized same-history source tube remain incomplete",
            "F": "shipping Mahony/proxy finite-time capture into the eventual P4 basin has not yet been established",
            "not_filter_instability": True,
        },
        "next_mathematical_obligations": [
            "derive a source-uniform discrete quotient storage/ISS inequality for the exact binary32 Mahony PI recurrence without assuming the 60 deg chart",
            "cover the leaky gravity-hold tail and timeout aligned-branch path",
            "materialize the same COMPLETE-BRMM startup frontend/tuner/magnetic continuation",
            "construct the widened source-uniform P4-H18 basin from the corrected correlated handoff fiber",
            "prove finite intersection of every startup/early-Live reachable tube with that P4-H18 basin",
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
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    if d.get("quality_gate_memory", {}).get("max_bad_tail_after_last_good_sample_s") != 4.0:
        f.append("unexpected leaky quality-gate bad tail")
    for key in (
        "mahony_arithmetic_boundedness_is_capture",
        "legacy_declared_mahony_chart_may_establish_capture",
        "STARTUP_REACHABLE_TUBE_CLOSED",
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
