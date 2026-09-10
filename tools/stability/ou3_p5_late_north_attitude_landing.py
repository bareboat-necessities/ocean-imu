#!/usr/bin/env python3
"""Late-north yaw rewrite -> widest P4 attitude-sector landing lemma.

This module composes two already-declared/certified facts:

1. the shipping late-north rewrite is a world-Z rotation and therefore leaves
   gravity-quotient tilt exactly invariant; and
2. the deployment theorem already declares an internal heading-gauge error
   bound of 10 degrees.

It does *not* shrink the P4 basin.  P4 keeps its widest current 30 degree full
SO(3) candidate.  P5 merely targets a 19 degree tilt-only capture core before
(or while waiting for) north.  By the SO(3) triangle inequality, after the
shipping yaw-only rewrite

    d_SO3(R_new,R_true) <= theta_tilt + |delta_yaw| < 19 deg + 10 deg = 29 deg,

leaving a full degree of strict margin inside the 30 degree P4 attitude sector.
The 19 degree core is a capture target, not a P4 entrance restriction.

This closes only the attitude coordinate of the late-north landing.  The same
history must still satisfy the P4 translation/bias/tuner/covariance entrance
fiber, and early-Live finite-time contraction to the 19 degree core remains a
separate P5 obligation.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import ou3_p5_late_north_yaw_reset as RESET

REPO = Path(__file__).resolve().parents[2]
DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P5_LATE_NORTH_ATTITUDE_LANDING_V1"
CAPTURE_TILT_CORE_DEG = 19.0


def build() -> dict:
    domain = json.loads(DOMAIN.read_text(encoding="utf-8"))
    reset = RESET.build()
    if RESET.validate(reset):
        raise RuntimeError("late-north reset structure prerequisite failed")

    startup = domain["startup"]
    candidates = [float(x) for x in domain["certificate_search"]["p4_complete_word_full_attitude_candidate_deg"]]
    widest_p4_deg = max(candidates)
    yaw_bound_rad = float(startup["internal_heading_gauge_error_upper_rad"])
    yaw_bound_deg = math.degrees(yaw_bound_rad)
    landing_upper_deg = CAPTURE_TILT_CORE_DEG + yaw_bound_deg
    strict_margin_deg = widest_p4_deg - landing_upper_deg

    structure = bool(
        reset["LATE_NORTH_HYBRID_RESET_MAP_STRUCTURE_CLOSED"]
        and reset["gravity_quotient_tilt_error_exactly_invariant"]
        and reset["linear_navigation_state_exactly_unchanged"]
        and reset["gyro_bias_state_exactly_unchanged"]
        and reset["accelerometer_bias_state_exactly_unchanged"]
        and reset["covariance_exactly_unchanged_by_setter"]
    )
    attitude_landing = bool(
        structure
        and yaw_bound_deg > 0.0
        and CAPTURE_TILT_CORE_DEG > 0.0
        and landing_upper_deg < widest_p4_deg
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "filter_changed": False,
        "P4_attitude_basin_shrunk": False,
        "P4_widest_full_attitude_candidate_deg": widest_p4_deg,
        "P5_tilt_capture_core_deg": CAPTURE_TILT_CORE_DEG,
        "P5_tilt_core_is_not_P4_basin_radius": True,
        "declared_internal_heading_gauge_error_upper_rad": yaw_bound_rad,
        "declared_internal_heading_gauge_error_upper_deg": yaw_bound_deg,
        "heading_bound_is_declared_theorem_premise_not_derived_here": True,
        "SO3_triangle_inequality_consumed": True,
        "late_north_reset_structure_closed": structure,
        "post_reset_full_attitude_strict_upper_deg": landing_upper_deg,
        "strict_margin_to_widest_P4_attitude_sector_deg": strict_margin_deg,
        "LATE_NORTH_ATTITUDE_LANDING_INTO_WIDEST_P4_SECTOR_CLOSED": attitude_landing,
        "other_P4_coordinates_landing_closed_here": False,
        "finite_time_reach_of_19deg_tilt_core_closed_here": False,
        "full_P4_membership_closed_here": False,
        "P4_PASS": False,
        "P5_PASS": False,
        "next_obligation": (
            "prove finite early-Live H18 gravity-quotient contraction from the <91.146 degree shipping handoff section into the 19 degree P5 tilt core while retaining the correlated translation/bias/tuner/covariance entrance fiber"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for k in (
        "P5_tilt_core_is_not_P4_basin_radius",
        "heading_bound_is_declared_theorem_premise_not_derived_here",
        "SO3_triangle_inequality_consumed",
        "late_north_reset_structure_closed",
        "LATE_NORTH_ATTITUDE_LANDING_INTO_WIDEST_P4_SECTOR_CLOSED",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "filter_changed", "P4_attitude_basin_shrunk",
        "other_P4_coordinates_landing_closed_here",
        "finite_time_reach_of_19deg_tilt_core_closed_here",
        "full_P4_membership_closed_here", "P4_PASS", "P5_PASS",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    if abs(float(d.get("declared_internal_heading_gauge_error_upper_deg", 0.0)) - 10.0) > 1e-12:
        f.append("declared heading-gauge error is not 10 degrees")
    if float(d.get("P4_widest_full_attitude_candidate_deg", 0.0)) != 30.0:
        f.append("widest P4 candidate changed")
    if not (float(d.get("post_reset_full_attitude_strict_upper_deg", 99.0)) < 30.0):
        f.append("late-north attitude landing does not fit widest P4 sector")
    if not (float(d.get("strict_margin_to_widest_P4_attitude_sector_deg", 0.0)) > 0.0):
        f.append("late-north attitude landing lacks strict margin")
    return f


if __name__ == "__main__":
    d = build()
    failures = validate(d)
    print(json.dumps({**d, "validation_pass": not failures, "validation_failures": failures}, indent=2, sort_keys=True))
    raise SystemExit(1 if failures else 0)
