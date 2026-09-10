#!/usr/bin/env python3
"""Startup magnetic observability theorem/counterexample for OU-III P5.

This is a structural proof, not a replay.  It answers whether the current
COMPLETE-BRMM declaration is sufficient to force finite-time *full-attitude*
capture after the shipping timeout handoff.

If no accepted magnetic acquisition is required, two physical histories that
differ only by a constant world-yaw gauge produce the same gyro/accelerometer
history.  The private Mahony observer is gravity-only and the timeout handoff
does not require north.  Until a magnetic packet is accepted, every gravity
quotient update remains yaw invariant.  Consequently no estimator driven by
those observations can guarantee membership of a strict full-attitude P4 basin
in finite time for both histories.

The theorem does not say the filter is unstable.  It identifies the weakest
missing observability premise: for with_mag=true the same COMPLETE-BRMM history
must supply enough accepted, nondegenerate magnetic packets to close the
shipping MagAutoTuner window in finite time, followed by the already-declared
Normal-Live accepted-vector recurrence.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SOURCE = REPO / "tools" / "stability" / "ou3_brmm_complete_source.py"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MAG_TUNER = REPO / "src" / "tuner" / "MagAutoTuner.h"
SCHEMA = 1
QUALIFICATION = "OU3_STARTUP_MAGNETIC_OBSERVABILITY_V1"


def _cfg_float(text: str, name: str) -> float:
    m = re.search(rf"(?:float|int)\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f?\s*;", text)
    if not m:
        raise RuntimeError(f"cannot extract Config::{name}")
    return float(m.group(1))


def build() -> dict:
    domain = json.loads(DOMAIN.read_text(encoding="utf-8"))
    source_text = SOURCE.read_text(encoding="utf-8")
    wrapper = WRAPPER.read_text(encoding="utf-8")
    mag_tuner = MAG_TUNER.read_text(encoding="utf-8")

    startup = domain["startup"]
    normal = domain["normal_live"]
    p4_candidates = [float(x) for x in domain["certificate_search"]["p4_complete_word_full_attitude_candidate_deg"]]
    widest_p4_deg = max(p4_candidates)

    min_samples = int(_cfg_float(wrapper, "mag_min_samples"))
    min_window_s = _cfg_float(wrapper, "mag_min_window_sec")
    proxy_timeout_s = _cfg_float(wrapper, "proxy_startup_timeout_sec")
    pe_window_s = float(normal["vector_pe_recurrence_window_s"])

    timeout_requires_north = (
        "const bool ready_by_timeout" in wrapper
        and "(t_ >= timeout_sec)" in wrapper
        and "mag_gravity_aligned_branch_;" in wrapper
        and "north_ready" not in wrapper.split("const bool ready_by_timeout", 1)[1].split("if (!ready_by_quality", 1)[0]
    ) is False
    # Spell the shipping fact positively as well, to make accidental parser
    # weakening visible in evidence.
    timeout_allows_ungauged = (
        "const bool ready_by_timeout" in wrapper
        and "mag_gravity_aligned_branch_;" in wrapper
    )

    startup_has_packet_recurrence = any(
        key in startup for key in (
            "accepted_magnetic_packet_recurrence_window_s",
            "magnetic_acquisition_recurrence_window_s",
            "magnetic_north_acquisition_time_upper_s",
        )
    )
    normal_has_vector_recurrence = pe_window_s > 0.0
    source_admission_open = '"BRMM_SOURCE_ADMISSION_PASS": False' in source_text
    finite_source_family_open = '"finite_window_family_materialized": False' in source_text

    # Use a 90 degree gauge pair.  Identical gyro/accelerometer observations
    # cannot distinguish the two histories before magnetic information enters.
    # For any single estimated yaw, triangle inequality on S1/SO(3) says at
    # least one history has geodesic attitude error >= separation/2.
    yaw_pair_separation_deg = 90.0
    minimax_error_lower_deg = yaw_pair_separation_deg / 2.0
    defeats_current_p4_candidates = minimax_error_lower_deg > widest_p4_deg

    mag_tuner_parity = {
        "minimum_accepted_samples_enforced": "accepted_count_ < std::max(1, cfg_.min_samples)" in mag_tuner,
        "minimum_accepted_window_enforced": "accepted_window_sec_ < cfg_.min_window_sec" in mag_tuner,
        "horizontal_field_nondegeneracy_enforced": "horiz_frac < cfg_.min_horizontal_fraction" in mag_tuner,
        "ready_only_after_finalize": "ready_ =" in mag_tuner and "return tryFinalize_();" in mag_tuner,
    }

    obstruction = (
        timeout_allows_ungauged
        and not startup_has_packet_recurrence
        and defeats_current_p4_candidates
        and all(mag_tuner_parity.values())
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "filter_changed": False,
        "quality_gates_changed": False,
        "current_source_contract": {
            "BRMM_source_admission_is_open": source_admission_open,
            "finite_source_family_materialization_is_open": finite_source_family_open,
            "startup_accepted_magnetic_packet_recurrence_declared": startup_has_packet_recurrence,
            "normal_live_vector_recurrence_declared": normal_has_vector_recurrence,
            "normal_live_vector_recurrence_window_s": pe_window_s,
        },
        "shipping_magnetic_acquisition": {
            "min_accepted_samples": min_samples,
            "min_accepted_window_s": min_window_s,
            "proxy_startup_timeout_s": proxy_timeout_s,
            "timeout_requires_north": timeout_requires_north,
            "timeout_allows_ungauged_live": timeout_allows_ungauged,
            "mag_auto_tuner_parity": mag_tuner_parity,
        },
        "yaw_indistinguishability_witness": {
            "history_A": "any admitted same-history BRMM/BIAS realization with no accepted magnetic packet before capture",
            "history_B": "history_A transformed by a constant 90 deg world-yaw gauge rotation",
            "gyro_and_accelerometer_observations_identical_before_mag": True,
            "private_mahony_gravity_quotient_observations_identical": True,
            "timeout_branch_behavior_identical": True,
            "yaw_pair_separation_deg": yaw_pair_separation_deg,
            "minimax_full_attitude_error_lower_deg": minimax_error_lower_deg,
            "widest_current_P4_candidate_deg": widest_p4_deg,
            "defeats_every_current_strict_full_attitude_P4_candidate": defeats_current_p4_candidates,
        },
        "UNCONDITIONAL_FULL_ATTITUDE_FINITE_CAPTURE_FROM_CURRENT_SOURCE": False,
        "structural_observability_obstruction_proved": obstruction,
        "failure_class": "E_SOURCE_OBSERVABILITY",
        "not_filter_instability": True,
        "weakest_required_source_extension": {
            "scope": "only when with_mag=true and theorem conclusion requires absolute/full-heading P4",
            "same_history_required": True,
            "finite_startup_acquisition": (
                "there exists finite T_north such that the shipping magnetic acquisition receives "
                f"at least {min_samples} accepted nondegenerate packets spanning at least {min_window_s:g} s, "
                "with finite/norm-consistent samples and nondegenerate horizontal mean, so mag_ref_set_ becomes true"
            ),
            "post_acquisition_recurrence": f"accepted magnetic/vector PE recurrence no slower than {pe_window_s:g} s (existing Normal-Live premise)",
            "does_not_require_heading_maneuver_when_hard_iron_estimation_is_disabled": True,
            "does_not_shrink_BRMM_motion_caps": True,
        },
        "theorem_split_required": {
            "gravity_quotient_capture_before_north": True,
            "full_attitude_P4_capture_only_after_north_gauge_event": True,
            "late_north_event_is_hybrid_yaw_reset": True,
        },
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("structural_observability_obstruction_proved") is not True:
        f.append("yaw observability obstruction not proved")
    if d.get("UNCONDITIONAL_FULL_ATTITUDE_FINITE_CAPTURE_FROM_CURRENT_SOURCE") is not False:
        f.append("false unconditional full-attitude capture promotion")
    s = d.get("shipping_magnetic_acquisition", {})
    if s.get("min_accepted_samples") != 128:
        f.append("shipping magnetic min-sample parity changed")
    if s.get("min_accepted_window_s") != 15.0:
        f.append("shipping magnetic window parity changed")
    if s.get("proxy_startup_timeout_s") != 150.0:
        f.append("shipping startup timeout parity changed")
    if s.get("timeout_requires_north") is not False or s.get("timeout_allows_ungauged_live") is not True:
        f.append("shipping timeout north semantics changed")
    w = d.get("yaw_indistinguishability_witness", {})
    if w.get("minimax_full_attitude_error_lower_deg") != 45.0:
        f.append("unexpected yaw witness minimax bound")
    if w.get("defeats_every_current_strict_full_attitude_P4_candidate") is not True:
        f.append("yaw witness no longer defeats current P4 candidate family")
    return f


if __name__ == "__main__":
    d = build()
    failures = validate(d)
    print(json.dumps({**d, "validation_pass": not failures, "validation_failures": failures}, indent=2))
    raise SystemExit(1 if failures else 0)
