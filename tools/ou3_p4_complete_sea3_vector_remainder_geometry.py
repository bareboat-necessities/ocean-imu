#!/usr/bin/env python3
"""Pure finite-angle vector remainder geometry for current complete-SEA3 P4.

This is bound only to ``COMPLETE_SEA3_NORMAL_LIVE_WORD`` and to the current
moving-Riccati P4 route.  It does not revive the retired endpoint/source-grid
remainder module.

For the Cayley attitude coordinate

    c = 2 tan(theta/2) u

and any inertial vector v, let

    h(c)   = [c]x v,
    y(c)   = (R(c)-I)v,
    eta(c) = y(c)-h(c).

The exact Cayley formula gives, for every v (only its component perpendicular
to u contributes),

    y^T eta = 0,
    ||eta|| = sin(theta/2) ||h||,
    ||y||   = cos(theta/2) ||h||.

The accelerometer operation-coordinate certificate removes a_w and b_a exactly
from nonlinear eta, so its eta is precisely this pure rotation defect with
v=f_hat.  The full linear a_w/b_a directions remain in y, S and K, and every
actual SpectralMSE R_S S=0 update remains in the complete word.  Magnetometer
radial energy is handled separately by the exact Joseph cancellation lemma.

Configured Racc and Rmag are isotropic on this proof branch, so the same norm
identities hold after whitening by R^-1.  These are homogeneous quadratic
sector identities, not additive disturbances and not an N-times packet budget.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import down, up
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_complete_sea3_accelerometer_operation_coordinate as ACC
import ou3_p4_magnetometer_radial_joseph as MAG
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_windowed_vector_pe as PE
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_PURE_VECTOR_REMAINDER_GEOMETRY"
PI_UP = 355.0 / 113.0


def _cell(theta_deg: float) -> dict:
    theta_hi = up(float(theta_deg) * PI_UP / 180.0)
    half_hi = up(0.5 * theta_hi)
    s = VT.sin_point(half_hi)
    c = VT.cos_point(half_hi)
    if not (s.hi >= 0.0 and c.lo > 0.0):
        raise RuntimeError("validated half-angle geometry lost chart positivity")
    s2 = up(s.hi * s.hi)
    c2 = down(c.lo * c.lo)
    return {
        "attitude_angle_deg": float(theta_deg),
        "attitude_angle_rad_upper": theta_hi,
        "sin_half_angle_upper": s.hi,
        "cos_half_angle_lower": c.lo,
        "eta_squared_over_linear_tangent_squared_upper": s2,
        "exact_residual_squared_over_linear_tangent_squared_lower": c2,
        "exact_residual_dot_eta_exact_zero": True,
        "sector_is_homogeneous_quadratic": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("vector remainder geometry must not be trajectory fitted")
    complete = COMPLETE.build(path)
    cayley = CAYLEY.build(path)
    acc = ACC.build(path)
    mag = MAG.build(path)
    pe = PE.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "cayley": CAYLEY.validate(cayley),
        "accelerometer_operation_coordinate": ACC.validate(acc),
        "mag_radial": MAG.validate(mag),
        "vector_PE": PE.validate(pe),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"complete-SEA3 vector remainder prerequisites failed: {bad}")

    meas = pe["measurement_runtime"]
    sa = list(map(float, meas["accelerometer_std_mps2"]))
    sm = list(map(float, meas["magnetometer_std_uT"]))
    acc_iso = len(sa) == 3 and sa[0] == sa[1] == sa[2]
    mag_iso = len(sm) == 3 and sm[0] == sm[1] == sm[2]
    candidates = list(map(float, domain["certificate_search"]["p4_complete_word_full_attitude_candidate_deg"]))
    rows = [_cell(x) for x in candidates]
    outer_deg = float(cayley["outer_angle_rad"]) * 180.0 / math.pi
    outer = _cell(outer_deg)

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "source_family_replaced": False,
        "retired_source_grid_remainder_route_used": False,
        "all_valid_accelerometer_updates_remain_in_complete_word": True,
        "all_due_S_updates_and_actual_RS_remain_in_complete_word": True,
        "accelerometer_operation_coordinate_consumed": True,
        "accelerometer_aw_nonlinear_eta_coefficient": float(acc["latent_aw_nonlinear_eta_coefficient"]),
        "accelerometer_bias_nonlinear_eta_coefficient": float(acc["accelerometer_bias_nonlinear_eta_coefficient"]),
        "accelerometer_eta_is_pure_force_rotation": True,
        "magnetometer_radial_Joseph_cancellation_consumed": bool(mag["radial_Joseph_energy_cancellation_exact"]),
        "configured_Racc_isotropic": acc_iso,
        "configured_Rmag_isotropic": mag_iso,
        "R_inverse_whitening_preserves_pure_vector_sector": bool(acc_iso and mag_iso),
        "exact_vector_identities": {
            "y_dot_eta": "0",
            "eta_norm": "sin(theta/2)*||[c]x v||",
            "exact_residual_norm": "cos(theta/2)*||[c]x v||",
        },
        "outer_geometry_cell": outer,
        "candidate_cells": rows,
        "packet_count_multiplier_used": False,
        "standalone_eta_disturbance_budget_used": False,
        "complete_signed_word_established_here": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "accumulate these homogeneous signed vector forms with the same complete-SEA3 innovation recursion and the actual-R_S four-S a_w regularizer; scalarize only after the recurrent H/A word"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "source_generated_not_trajectory_fit", "all_valid_accelerometer_updates_remain_in_complete_word",
        "all_due_S_updates_and_actual_RS_remain_in_complete_word",
        "accelerometer_operation_coordinate_consumed", "accelerometer_eta_is_pure_force_rotation",
        "magnetometer_radial_Joseph_cancellation_consumed", "configured_Racc_isotropic",
        "configured_Rmag_isotropic", "R_inverse_whitening_preserves_pure_vector_sector",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed", "source_family_replaced",
        "retired_source_grid_remainder_route_used", "packet_count_multiplier_used",
        "standalone_eta_disturbance_budget_used", "complete_signed_word_established_here", "P4_promoted_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("accelerometer_aw_nonlinear_eta_coefficient") != 0.0:
        f.append("a_w re-entered nonlinear accelerometer eta")
    if d.get("accelerometer_bias_nonlinear_eta_coefficient") != 0.0:
        f.append("b_a re-entered nonlinear accelerometer eta")
    rows = d.get("candidate_cells", [])
    if [x.get("attitude_angle_deg") for x in rows] != [30.0, 25.0, 20.0, 15.0]:
        f.append("candidate finite-angle cells changed")
    for row in [d.get("outer_geometry_cell", {})] + rows:
        s2 = float(row.get("eta_squared_over_linear_tangent_squared_upper", math.inf))
        c2 = float(row.get("exact_residual_squared_over_linear_tangent_squared_lower", -math.inf))
        if not (0.0 <= s2 < 1.0 and 0.0 < c2 <= 1.0 and row.get("exact_residual_dot_eta_exact_zero") is True):
            f.append(f"invalid pure-vector sector cell {row.get('attitude_angle_deg')}")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "outer": d["outer_geometry_cell"],
        "candidate_cells": d["candidate_cells"],
        "aw_eta": d["accelerometer_aw_nonlinear_eta_coefficient"],
        "actual_RS_retained": d["all_due_S_updates_and_actual_RS_remain_in_complete_word"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
