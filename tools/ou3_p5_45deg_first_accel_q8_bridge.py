#!/usr/bin/env python3
"""Sign-complete first-Live accelerometer bridge from the 45 deg P5 entrance.

PR #445 gives P5 a deployment-scale 45 deg SO(3) entrance while the eventual
P4 complete-word sector may be narrower (30/25/20/15 deg).  The exact-source
first-accelerometer range certificate proves the first shipping correction is
inside the already validated six-radian deployed-quaternion helper, but a
correction-range bound alone does not prove that the attitude remains inside
the narrow P4 sector.

This producer proves the correct intermediate statement without assuming that
the accelerometer correction has a favorable sign.  For a pre-update Cayley
vector c with ||c||<=q and a shipping correction d with ||d||<=dmax, the
homogeneous scalar of the exact left quaternion product is

    W = 2 cos(||d||/2) - v_d^T c,
    ||v_d|| = sin(||d||/2).

Because every first-source correction has dmax < pi and the half-angle lies in
the monotone first quadrant,

    W >= 2 cos(dmax/2) - q sin(dmax/2).

The normalized product quaternion has scalar W/sqrt(4+q^2).  A Cayley norm
strictly below 8 is therefore equivalent to the algebraic condition

    17 W^2 > 4 + q^2,

provided W>0.  This avoids atan/tan in the promoted calculation.  sin/cos are
outward enclosed by the validated exact-rational Taylor backend.

The same correction norm also bounds the shipping reset congruence
G=I+0.5[d]_x by ||G||_2^2 <= 1+dmax^2/4.  Joseph gives P+<=P- in Loewner order,
so both the accepted/reset branch and the identity/rejected branch have a
finite source-complete covariance multiplier.  This stage does not claim that
the first packet returns to 30 deg, and it does not promote P4/P5; it establishes
a finite q<8 transient chart from which source-correlated recapture can be
continued.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_candidate_first_accel_exact_source as FIRST
import ou3_p4_p5_entrance_search_domain as ENTRANCE
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_startup_stability_certificate as P1
import ou3_validated_transcendentals as VT

DEFAULT_DOMAIN = FIRST.DEFAULT_DOMAIN
SCHEMA = 1
Q_TRANSIENT = 8.0
Q8_SQUARED_QUAT_SCALAR_DENOM = 17.0  # q<8 iff scalar^2 > 1/17


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _q8_product_bound(q: float, dmax: float) -> dict:
    q = float(q)
    dmax = float(dmax)
    if not (math.isfinite(q) and q >= 0.0):
        raise ValueError("finite nonnegative pre-update Cayley radius required")
    if not (math.isfinite(dmax) and 0.0 <= dmax < math.pi):
        raise ValueError("first correction must lie below pi for monotone half-angle bound")

    half_hi = up(0.5 * dmax)
    s = VT.sin_point(half_hi)
    c = VT.cos_point(half_hi)
    # On [0,half_hi] subset [0,pi/2], cos decreases and sin increases.
    qsin = up(q * s.hi)
    Wlo = down(down(2.0 * c.lo) - qsin)
    q2 = up(q * q)
    rhs = up(4.0 + q2)
    W2lo = down(Wlo * Wlo) if Wlo > 0.0 else 0.0
    lhs = down(Q8_SQUARED_QUAT_SCALAR_DENOM * W2lo)
    q8_safe = Wlo > 0.0 and lhs > rhs

    qplus_upper = math.inf
    if Wlo > 0.0:
        # scalar_product^2 >= Wlo^2/(4+q^2), hence
        # q_plus^2 = 4/scalar^2-4 <= 4(4+q^2)/Wlo^2-4.
        ratio = up(up(4.0 * rhs) / W2lo) if W2lo > 0.0 else math.inf
        qplus2 = up(max(0.0, ratio - 4.0)) if math.isfinite(ratio) else math.inf
        if math.isfinite(qplus2):
            qplus_upper = P1.sqrt_interval_point(qplus2).hi

    return {
        "pre_update_q_upper": q,
        "correction_norm_upper_rad": dmax,
        "half_correction_upper_rad": half_hi,
        "validated_sin_half_upper": s.hi,
        "validated_cos_half_lower": c.lo,
        "homogeneous_product_scalar_lower": Wlo,
        "q8_test_lhs_lower_17W2": lhs,
        "q8_test_rhs_upper_4plusq2": rhs,
        "post_update_q_upper_from_scalar": qplus_upper,
        "inside_q8_for_every_correction_direction": q8_safe,
    }


def _source_max_correction_for_q(
    qpred: float, *, gravity: float, ba_H: float, source_rows: list[dict]
) -> tuple[float, float, dict | None]:
    max_d = 0.0
    max_resid = 0.0
    worst = None
    rotational = up(gravity * qpred)
    for row in source_rows:
        aw = float(row["predicted_aw_error_norm_upper_mps2"])
        k = float(row["Ktheta_norm_upper"])
        residual = up(rotational + up(aw + ba_H))
        d = up(k * residual)
        if d > max_d:
            max_d = d
            max_resid = residual
            worst = {
                "source_phase_cell": row["source_phase_cell"],
                "pseudo_phase": row["pseudo_phase"],
                "tau_s": row["tau_s"],
                "sigma_aw_mps2": row["sigma_aw_mps2"],
                "R_S_filter_std": row["R_S_filter_std"],
                "predicted_aw_error_norm_upper_mps2": aw,
                "Ktheta_norm_upper": k,
                "rotational_residual_norm_upper_mps2": rotational,
                "combined_residual_norm_upper_mps2": residual,
                "correction_norm_upper_rad": d,
            }
    return max_d, max_resid, worst


def _reset_bounds(dmax: float) -> dict:
    d2 = up(dmax * dmax)
    multiplier = up(1.0 + up(0.25 * d2))
    op = P1.sqrt_interval_point(multiplier).hi
    return {
        "reset_operator_norm_upper": op,
        "reset_covariance_spectral_multiplier_upper": multiplier,
        "accepted_Joseph_posterior_Loewner_below_prior": True,
        "rejected_branch_is_identity": True,
        "accepted_or_rejected_covariance_multiplier_upper": max(1.0, multiplier),
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2) -> dict:
    path = Path(domain_path).resolve()
    dom = json.loads(path.read_text(encoding="utf-8"))
    if dom.get("trajectory_fit") is not False:
        raise RuntimeError("45deg q8 bridge must not be trajectory fitted")
    if dom.get("configured_runtime", {}).get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("45deg q8 bridge requires lever arm disabled")

    first = FIRST.build(path, source_pieces=source_pieces)
    entrance = ENTRANCE.build(path)
    failures = [f"first-accel: {x}" for x in FIRST.validate(first)]
    failures += [f"entrance: {x}" for x in ENTRANCE.validate(entrance)]
    if first.get("P4_CANDIDATE_FIRST_ACCEL_EXACT_SOURCE_CERTIFICATE") != "PASS":
        failures.append("exact-source first-accelerometer prerequisite did not pass")

    h = float(FULL._source_cell()["dt_s"])
    gravity = float(dom["startup"]["gravity_mps2"])
    ba_H = float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])

    # The first-source gain/covariance is independent of the physical candidate
    # attitude radius.  Reuse its source-correlated rows and replace only the
    # exact rotational residual bound g*q by the requested entrance q.
    candidate_rows = first.get("candidate_rows", [])
    if not candidate_rows:
        failures.append("candidate first-accelerometer rows are missing")
        source_rows = []
    else:
        source_rows = candidate_rows[0].get("source_rows", [])
    if not source_rows:
        failures.append("source-phase first-accelerometer rows are missing")

    q45 = float(entrance["P5_entrance"]["attitude_geometry"]["cayley_norm_upper"])
    q45pred = RG._q_after_first_prediction(q45, dom, h)
    d45, residual45, worst45 = _source_max_correction_for_q(
        q45pred, gravity=gravity, ba_H=ba_H, source_rows=source_rows)
    p5_product = _q8_product_bound(q45pred, d45)
    p5_reset = _reset_bounds(d45)

    rows = []
    for row in candidate_rows:
        qpred = float(row["post_prediction_q_upper"])
        dmax = float(row["max_first_accelerometer_correction_norm_upper_rad"])
        prod = _q8_product_bound(qpred, dmax)
        rows.append({
            "angle_deg": float(row["angle_deg"]),
            "pre_update_q_upper": qpred,
            "correction_norm_upper_rad": dmax,
            "product": prod,
            "reset": _reset_bounds(dmax),
        })

    p5_safe = p5_product["inside_q8_for_every_correction_direction"] is True
    candidates_safe = bool(rows) and all(r["product"]["inside_q8_for_every_correction_direction"] for r in rows)
    if not p5_safe:
        failures.append("45deg entrance first accepted correction is not certified inside q8")
    if not candidates_safe:
        failures.append("a P4 candidate first accepted correction is not certified inside q8")

    passed = not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_45DEG_FIRST_ACCEL_SIGN_COMPLETE_Q8_AND_RESET_BRIDGE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "shipping_accelerometer_has_no_assumed_NIS_rejection_gate": True,
        "accepted_correction_sign_assumed_favorable": False,
        "all_correction_directions_covered_by_dot_product_extremum": True,
        "validated_sin_cos_used": True,
        "atan_tan_used_in_promoted_q8_test": False,
        "deployed_correction_limit_rad": float(first["deployed_correction_limit_rad"]),
        "deployed_correction_limit_increased": False,
        "transient_cayley_chart_q_upper": Q_TRANSIENT,
        "P5_45deg_entrance_first_accel": {
            "entrance_q_upper": q45,
            "post_prediction_q_upper": q45pred,
            "max_combined_residual_norm_upper_mps2": residual45,
            "max_first_accelerometer_correction_norm_upper_rad": d45,
            "worst_source_phase_child": worst45,
            "product": p5_product,
            "reset": p5_reset,
            "accepted_or_rejected_branch_family_inside_q8": p5_safe,
        },
        "P4_candidate_first_accel_rows": rows,
        "all_P4_candidate_first_accel_branches_inside_q8": candidates_safe,
        "H_dimension": 18,
        "A_dimension": 21,
        "H_first_accel_range_bound_covers_A_attitude_gain": first.get("A_first_prefix_attitude_gain_bounded_by_H_gain") is True,
        "Joseph_reset_covariance_spectral_bound_applies_to_H_and_A": True,
        "detailed_post_reset_cross_covariance_propagated_here": False,
        "returned_to_30deg_P4_sector_here": False,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P5_FINITE_CAPTURE_TO_P4_ESTABLISHED_HERE": False,
        "P5_45DEG_FIRST_ACCEL_Q8_BRIDGE_CERTIFICATE": "PASS" if passed else "NOT_ESTABLISHED",
        "next_obligation": (
            "propagate the source-correlated H/A Joseph/reset covariance and state children from this finite q<8 first-packet bridge to sample 1, then prove finite recapture into the 30deg P4 candidate without replacing the signed/source correlations by a global cube"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "shipping_accelerometer_has_no_assumed_NIS_rejection_gate",
        "all_correction_directions_covered_by_dot_product_extremum",
        "validated_sin_cos_used",
        "Joseph_reset_covariance_spectral_bound_applies_to_H_and_A",
        "H_first_accel_range_bound_covers_A_attitude_gain",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "accepted_correction_sign_assumed_favorable",
        "atan_tan_used_in_promoted_q8_test", "deployed_correction_limit_increased",
        "detailed_post_reset_cross_covariance_propagated_here", "returned_to_30deg_P4_sector_here",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE", "P5_FINITE_CAPTURE_TO_P4_ESTABLISHED_HERE",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction range changed")
    if float(d.get("transient_cayley_chart_q_upper", 0.0)) != 8.0:
        f.append("transient q chart changed")
    p5 = d.get("P5_45deg_entrance_first_accel", {})
    if p5.get("accepted_or_rejected_branch_family_inside_q8") is not True:
        f.append("45deg first-packet branch family is not inside q8")
    prod = p5.get("product", {})
    if not float(prod.get("homogeneous_product_scalar_lower", -math.inf)) > 0.0:
        f.append("45deg product scalar can cross Cayley antipode")
    if not float(prod.get("q8_test_lhs_lower_17W2", -math.inf)) > float(prod.get("q8_test_rhs_upper_4plusq2", math.inf)):
        f.append("45deg algebraic q8 margin is not positive")
    if d.get("all_P4_candidate_first_accel_branches_inside_q8") is not True:
        f.append("candidate first-packet branches are not all inside q8")
    if [float(r.get("angle_deg", -1.0)) for r in d.get("P4_candidate_first_accel_rows", [])] != [30.0,25.0,20.0,15.0]:
        f.append("candidate ladder changed")
    if d.get("H_dimension") != 18 or d.get("A_dimension") != 21:
        f.append("H/A dimensions changed")
    st = d.get("P5_45DEG_FIRST_ACCEL_Q8_BRIDGE_CERTIFICATE")
    if st not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid q8 bridge status")
    if st == "PASS" and f:
        f.append("PASS carries validation failures")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(x.domain.resolve(), source_pieces=x.source_pieces)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    p = d["P5_45deg_entrance_first_accel"]
    print(json.dumps({
        "status": d["P5_45DEG_FIRST_ACCEL_Q8_BRIDGE_CERTIFICATE"],
        "q45_pred": p["post_prediction_q_upper"],
        "d45_max": p["max_first_accelerometer_correction_norm_upper_rad"],
        "W_lower": p["product"]["homogeneous_product_scalar_lower"],
        "qplus_upper": p["product"]["post_update_q_upper_from_scalar"],
        "q8_lhs": p["product"]["q8_test_lhs_lower_17W2"],
        "q8_rhs": p["product"]["q8_test_rhs_upper_4plusq2"],
        "reset_cov_multiplier": p["reset"]["reset_covariance_spectral_multiplier_upper"],
        "candidate_qplus": [
            [r["angle_deg"], r["product"]["post_update_q_upper_from_scalar"]]
            for r in d["P4_candidate_first_accel_rows"]
        ],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
