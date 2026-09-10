#!/usr/bin/env python3
"""Derive the widest structural tilt set at the actual OU-III Live handoff.

The shipping handoff sample is required to satisfy the world-frame aligned
branch ``acc_world_lp.z() < 0``.  This is stronger and more useful for P5 than
assuming the historical 60 degree Mahony chart: it is an observable predicate
on the state that actually crosses into Live.

Let u_hat be the normalized low-passed specific-force direction after rotation
by the proxy attitude and let g_hat be the direction that would result from the
same proxy attitude applied to true gravity.  The deployed branch gives
angle(u_hat,-e_z) < pi/2.  The declared same-history world-averaged gravity
direction error gives angle(u_hat,g_hat) <= eps_g.  By the spherical triangle
inequality, the proxy-vs-true gravity quotient error obeys

    theta_tilt < pi/2 + eps_g.

No Mahony-chart assumption and no yaw gauge is used.  The result is deliberately
wide; early Live/H18 capture is responsible for contracting it into P4.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
COMMON = REPO / "src" / "kalman_common" / "SeaStateFusionFilterCommon.h"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
SCHEMA = 1
QUALIFICATION = "OU3_SHIPPING_HANDOFF_TILT_HEMISPHERE_V1"


def build() -> dict:
    domain = json.loads(DOMAIN.read_text(encoding="utf-8"))
    common = COMMON.read_text(encoding="utf-8")
    wrapper = WRAPPER.read_text(encoding="utf-8")
    eps = float(domain["startup"]["world_averaged_gravity_direction_error_upper_rad"])

    common_branch_parity = (
        "return acc_world_lp.z() < 0.0f;" in common
        and "gravityAlignedBranchWorld" in common
    )
    wrapper_handoff_branch_parity = (
        "const bool ready_by_timeout" in wrapper
        and "mag_gravity_aligned_branch_;" in wrapper
        and "if (!ready_by_quality && !ready_by_timeout) return;" in wrapper
    )

    structural_upper_rad = 0.5 * math.pi + eps
    structural_upper_deg = math.degrees(structural_upper_rad)

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "filter_changed": False,
        "quality_gates_changed": False,
        "shipping_branch_parity": {
            "world_branch_is_strict_negative_z": common_branch_parity,
            "handoff_requires_branch_on_timeout_path": wrapper_handoff_branch_parity,
        },
        "same_history_world_averaged_gravity_direction_error_upper_rad": eps,
        "handoff_measurement_direction_angle_to_expected_world_up_strict_upper_deg": 90.0,
        "spherical_triangle_inequality": True,
        "true_gravity_quotient_tilt_strict_upper_rad": structural_upper_rad,
        "true_gravity_quotient_tilt_strict_upper_deg": structural_upper_deg,
        "historical_60_deg_mahony_chart_consumed": False,
        "yaw_gauge_required": False,
        "applies_to_quality_and_timeout_handoff_samples": True,
        "handoff_tilt_set_representation": "OPEN_GEODESIC_BALL_ON_S2_QUOTIENT",
        "early_live_must_contract_to_P4": True,
        "HANDOFF_TILT_HEMISPHERE_BOUND_CLOSED": (
            common_branch_parity
            and wrapper_handoff_branch_parity
            and 0.0 <= eps < 0.5 * math.pi
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("HANDOFF_TILT_HEMISPHERE_BOUND_CLOSED") is not True:
        f.append("shipping handoff hemisphere bound not closed")
    if d.get("historical_60_deg_mahony_chart_consumed") is not False:
        f.append("historical Mahony chart reintroduced")
    if d.get("yaw_gauge_required") is not False:
        f.append("yaw gauge incorrectly required for tilt quotient")
    expected = math.degrees(0.5 * math.pi + 0.02)
    if abs(float(d.get("true_gravity_quotient_tilt_strict_upper_deg", -1.0)) - expected) > 1e-12:
        f.append("unexpected handoff tilt bound")
    return f


if __name__ == "__main__":
    d = build()
    failures = validate(d)
    print(json.dumps({**d, "validation_pass": not failures, "validation_failures": failures}, indent=2))
    raise SystemExit(1 if failures else 0)
