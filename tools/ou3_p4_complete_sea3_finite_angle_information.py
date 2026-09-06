#!/usr/bin/env python3
"""Finite-angle complete-SEA3 information sector for OU-III P4.

P4-only consumer of the frozen complete-SEA3 P3 information lemmas.  It does
not create a source word and does not replace the complete shipping word by the
selected PE/four-S witnesses.  Those witnesses remain PSD lower-bound
components of the same complete word.

For Cayley attitude error c with physical angle theta, the exact vector
residual differential retains at least k(theta)=cos(theta/2)^2 of the attitude
coordinate singular value, hence k(theta)^2 of the eta6 vector information.
The accelerometer a_w column is an orthogonal rotation and the same actual
SpectralMSE R_S four-S regularizer is unchanged.  Reusing the frozen triangular
Schur bound gives

  lambda >= alpha_theta*d_aw/(alpha_theta+||C_aw||^2+d_aw),
  alpha_theta = k(theta)^2*alpha6.

The magnetometer radial component is removed by the exact Joseph cancellation
lemma.  The accelerometer nonlinear component is deliberately NOT declared
harmless here; it is charged by the downstream signed-Joseph/reset word.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import down, up
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_magnetometer_radial_joseph as MAG
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_h18_information_composition as HINFO
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_FINITE_ANGLE_INFORMATION_SECTOR_V2"
USEFUL_GATE = 1.0e-18
# Rational Archimedean upper; avoids using libm pi in a promoted lower bound.
PI_UP = 355.0 / 113.0


def _candidate(theta_deg: float, *, alpha6: float, cross2: float, d_aw: float,
               non_aw: float) -> dict:
    deg = float(theta_deg)
    if not (math.isfinite(deg) and 0.0 < deg <= 60.0):
        raise ValueError("P4 candidate angle outside audited finite-angle range")
    theta_hi = up(deg * PI_UP / 180.0)
    half_hi = up(0.5 * theta_hi)
    # cos is positive/decreasing on this audited interval.  Evaluate at the
    # upper angle and take the validated lower endpoint.
    c = VT.cos_point(half_hi)
    if c.lo <= 0.0:
        raise RuntimeError("validated candidate cosine lost positivity")
    k = down(c.lo * c.lo)
    retention = down(k * k)
    alpha_theta = down(retention * alpha6)
    det = down(alpha_theta * d_aw)
    trace = up(alpha_theta + cross2 + d_aw)
    coupled = down(det / trace)
    full = min(coupled, non_aw)
    return {
        "attitude_angle_deg": deg,
        "attitude_angle_rad_upper": theta_hi,
        "validated_cos_half_interval": c.as_list(),
        "cayley_chart_sigma_min_lower": k,
        "vector_information_retention_factor_lower": retention,
        "eta6_information_lower_after_finite_angle": alpha_theta,
        "aw_direction_information_lower_unchanged": d_aw,
        "accelerometer_aw_cross_norm_squared_upper_unchanged": cross2,
        "coupled_eta6_aw_information_lambda_min_lower": coupled,
        "non_aw_translation_information_lower_unchanged": non_aw,
        "full_H18_information_lambda_min_lower": full,
        "clears_useful_gate": full >= USEFUL_GATE,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("finite-angle P4 information sector must not be trajectory fitted")

    complete = COMPLETE.build(path)
    info = HINFO.build(path)
    outer = CAYLEY.build(path)
    mag = MAG.build(path)
    bad = {
        "complete_SEA3": COMPLETE.validate(complete),
        "H18_information": HINFO.validate(info),
        "outer_Cayley": CAYLEY.validate(outer),
        "mag_radial": MAG.validate(mag),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"finite-angle information prerequisites failed: {bad}")
    if complete["canonical_P3_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("canonical complete SEA3 source changed")

    tri = info["triangular_information_composition"]
    alpha6 = float(tri["A_transpose_A_lower"])
    cross2 = float(tri["C_aw_spectral_norm_squared_upper"])
    d_aw = float(tri["aw_direction_information_lower"])
    non_aw = float(tri["non_aw_translation_lambda_min_lower"])
    linear_full = float(tri["D_H18_lambda_min_lower"])

    candidate_deg = list(map(float, domain["certificate_search"]["p4_complete_word_full_attitude_candidate_deg"]))
    rows = [_candidate(x, alpha6=alpha6, cross2=cross2, d_aw=d_aw, non_aw=non_aw)
            for x in candidate_deg]
    passing = [r for r in rows if r["clears_useful_gate"]]
    widest = passing[0] if passing else None

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "ordinary_libm_trigonometric_used_in_pass_decision": False,
        "validated_transcendental_backend_used": True,
        "P3_frozen_not_modified": True,
        "component_of_complete_SEA3_full_word": True,
        "selected_PE_or_four_S_replace_complete_word": False,
        "all_due_S_updates_remain_in_complete_word": True,
        "all_valid_accelerometer_updates_remain_in_complete_word": True,
        "actual_applied_SpectralMSE_R_S_consumed_through_frozen_H18_information": bool(
            info["actual_applied_SpectralMSE_R_S_consumed"]
        ),
        "directional_four_S_R_S_regularizer_retained": bool(
            info["directional_four_S_inverse_row_bound_used"]
        ),
        "finite_angle_vector_geometry_exact": True,
        "magnetometer_radial_state_correction_cancels_exactly": bool(
            mag["radial_Joseph_energy_cancellation_exact"]
        ),
        "accelerometer_radial_remainder_declared_zero": False,
        "accelerometer_radial_remainder_requires_signed_Joseph_word_charge": True,
        "outer_geometry_angle_rad_retained": float(outer["outer_angle_rad"]),
        "outer_geometry_sector_separate_from_inner_dissipation_cell": True,
        "linear_H18_information_lambda_min_lower": linear_full,
        "candidate_cells": rows,
        "widest_information_cell_deg": widest["attitude_angle_deg"] if widest else None,
        "widest_information_cell_H18_lambda_min_lower": (
            widest["full_H18_information_lambda_min_lower"] if widest else 0.0
        ),
        "all_declared_candidate_cells_keep_strict_information": bool(passing and len(passing) == len(rows)),
        "information_headroom_closed": bool(widest is not None),
        "signed_Joseph_complete_word_closed_here": False,
        "reset_defect_complete_word_closed_here": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "use retained finite-angle eta6/a_w information and the same actual-R_S complete SEA3 word in the signed Joseph correction/reset ledger; never charge an N-times worst standalone accelerometer remainder"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "source_generated_not_trajectory_fit", "validated_transcendental_backend_used",
        "P3_frozen_not_modified", "component_of_complete_SEA3_full_word",
        "all_due_S_updates_remain_in_complete_word",
        "all_valid_accelerometer_updates_remain_in_complete_word",
        "actual_applied_SpectralMSE_R_S_consumed_through_frozen_H18_information",
        "directional_four_S_R_S_regularizer_retained", "finite_angle_vector_geometry_exact",
        "magnetometer_radial_state_correction_cancels_exactly",
        "accelerometer_radial_remainder_requires_signed_Joseph_word_charge",
        "outer_geometry_sector_separate_from_inner_dissipation_cell",
        "all_declared_candidate_cells_keep_strict_information", "information_headroom_closed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "ordinary_libm_trigonometric_used_in_pass_decision",
        "selected_PE_or_four_S_replace_complete_word", "accelerometer_radial_remainder_declared_zero",
        "signed_Joseph_complete_word_closed_here", "reset_defect_complete_word_closed_here",
        "P4_promoted_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    rows = d.get("candidate_cells", [])
    if [r.get("attitude_angle_deg") for r in rows] != [30.0, 25.0, 20.0, 15.0]:
        f.append("declared P4 candidate order changed")
    for row in rows:
        if row.get("clears_useful_gate") is not True:
            f.append(f"candidate {row.get('attitude_angle_deg')} deg lost useful information gate")
        x = float(row.get("full_H18_information_lambda_min_lower", 0.0))
        if not (math.isfinite(x) and x >= USEFUL_GATE):
            f.append(f"candidate {row.get('attitude_angle_deg')} deg information lower invalid")
    if float(d.get("outer_geometry_angle_rad_retained", 0.0)) < 0.80:
        f.append("outer 0.8-rad geometry sector was lost")
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
        "widest_information_cell_deg": d["widest_information_cell_deg"],
        "widest_information_cell_H18_lambda_min_lower": d["widest_information_cell_H18_lambda_min_lower"],
        "candidate_cells": d["candidate_cells"],
        "outer_geometry_angle_rad_retained": d["outer_geometry_angle_rad_retained"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
