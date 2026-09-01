#!/usr/bin/env python3
"""Exact-source first-Live accelerometer range for finite-angle P4 candidates.

The generic candidate range deliberately allowed any normal-Live predicted
specific-force magnitude and any relative yaw-covariance/force alignment.  That
is not source reachable at the first Live packet.

The source-audited P5 first-accelerometer certificate proves that at this packet

* the MEKF a_w, v, p and S means are still zero;
* enabling the linear block resets a_w covariance, not its mean;
* the first homogeneous prediction keeps those means zero;
* an immediately-due S=0 pseudo update has zero mean residual;
* hence the accelerometer Jacobian force is exactly R_wb(-g), with norm g;
* the goLive yaw covariance rank-one axis is the same body image of world down.

Thus the structured first-accelerometer gain has alignment x=0 and force m=g.
For a candidate Cayley radius q after the first 5 ms transport, the gravity
rotational residual is bounded directly by g*q_tangent <= g*q.  The latent
linear and finite-rotation-cross terms are combined before taking a norm,
||e_aw+(R^T-I)e_aw||=||e_aw||.  The H-mode physical bias-error allowance also
contains the declared A-mode entrance ball, while the A covariance adds only an
isotropic PSD innovation term at this first packet.

This certifies only the correction *range* of the first Live accelerometer
packet.  It does not yet propagate the signed correction through Joseph/reset
or promote complete-word P4 dissipation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_candidate_first_accel_range_v2 as ASTRUCT
import ou3_p4_candidate_first_accel_range_v3 as TANGENT
import ou3_p4_candidate_full_word as CAND
import ou3_p4_p5_entrance_search_domain as ENTRANCE
import ou3_p5_first_accel_exact_source_v2 as EXACT
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_first_accel_rotation_gauge_v3 as RG3
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = CAND.DEFAULT_DOMAIN
SCHEMA = 4
LIMIT = 6.0


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2) -> dict:
    path = Path(domain_path).resolve()
    dom = json.loads(path.read_text(encoding="utf-8"))
    if dom.get("trajectory_fit") is not False:
        raise RuntimeError("candidate exact-source first accel must not be trajectory fitted")
    if dom.get("configured_runtime", {}).get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("candidate exact-source first accel requires lever arm disabled")

    semantics, missing = EXACT._source_semantics()
    entrance = ENTRANCE.build(path)
    vector = VECTOR.build()
    failures = [f"source semantic missing: {x}" for x in missing]
    failures += [f"entrance: {x}" for x in ENTRANCE.validate(entrance)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]

    RG3._install_backend(path, source_pieces)
    FULL3._install_backend()
    src_phases = RG._source_phase_children(source_pieces)
    if not src_phases:
        failures.append("no first-prefix source phase cells")

    h = float(FULL._source_cell()["dt_s"])
    gravity = float(dom["startup"]["gravity_mps2"])
    m = Interval.outward_bounds(gravity, gravity)
    x_aligned = Interval.outward_bounds(0.0, 0.0)
    tilt, yaw, eps = RG._attitude_covariance_epsilon(path, h)
    vc = vector["configured_measurement_bounds"]
    Racc = FULL._R_diag(float(vc["acc_measurement_std_mps2"]))
    racc_var = Racc[0][0]
    ba_H = float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    ba_A = float(dom["normal_live"]["active_accelerometer_bias_state_norm_upper_mps2"])
    if not (0.0 <= ba_A <= ba_H):
        failures.append("H bias-error allowance no longer contains A entrance ball")

    astruct = ASTRUCT._structural_a_check(path, dom, src_phases[0][0]) if src_phases else {}
    if astruct.get("A_bias_innovation_addition_isotropic_PSD") is not True:
        failures.append("A first-prefix isotropic innovation addition not established")
    if astruct.get("first_prefix_theta_aw_S_to_ba_cross_exact_zero") is not True:
        failures.append("A first-prefix forbidden ba cross covariance")

    rows = []
    widest = None
    for crow in entrance["P4_complete_word_search"]["candidate_rows"]:
        angle = float(crow["angle_deg"])
        q0 = float(crow["cayley_norm_upper"])
        qpred = RG._q_after_first_prediction(q0, dom, h)
        # No separate tilt assumption is needed: the gravity-tangent Cayley
        # component is bounded by the full candidate Cayley norm.
        c_tangent = qpred
        rotational = FULL.up(gravity * c_tangent)

        max_k = 0.0
        max_resid = 0.0
        max_d = 0.0
        min_margin = math.inf
        over = 0
        first_over = None
        source_rows = []
        for si, (src, phase) in enumerate(src_phases):
            P0 = FULL._initial_covariance(src, path)
            F, Q, _ = FULL._transition_and_Q(src, dom)
            Pp = FULL._psd_tighten(FULL.matrix_add(FULL.matrix_mul(FULL.matrix_mul(F, P0), FULL.matrix_transpose(F)), Q))
            _pss, _psa, paw_pred = RG._scalar_axis_structure(Pp)
            aw_pred, eS_pred = RG._prediction_norms(src, dom)
            if phase == "due":
                paw, aw_norm = RG._due_paw_and_error_norm(Pp, src, aw_pred, eS_pred)
            else:
                paw, aw_norm = paw_pred, aw_pred

            k, _kh, detail = TANGENT._tangent_structured_gain_bounds(
                tilt=tilt, yaw=yaw, eps=eps, x=x_aligned, m=m,
                paw=paw, racc_var=racc_var)
            residual = FULL.up(rotational + FULL.up(aw_norm + ba_H))
            d = FULL.up(k * residual)
            margin = LIMIT - d
            max_k = max(max_k, k)
            max_resid = max(max_resid, residual)
            max_d = max(max_d, d)
            min_margin = min(min_margin, margin)
            r = {
                "source_phase_cell": si,
                "pseudo_phase": phase,
                "tau_s": src["tau_s"].as_list(),
                "sigma_aw_mps2": src["sigma_aw_mps2"].as_list(),
                "R_S_filter_std": src["R_S_filter_std"].as_list(),
                "P_aw_variance_interval": paw.as_list(),
                "predicted_aw_error_norm_upper_mps2": aw_norm,
                "candidate_cayley_tangent_norm_upper": c_tangent,
                "rotational_residual_norm_upper_mps2": rotational,
                "combined_residual_norm_upper_mps2": residual,
                "Ktheta_norm_upper": k,
                "correction_norm_upper_rad": d,
                "gain_detail": detail,
            }
            source_rows.append(r)
            if not math.isfinite(d) or d > LIMIT:
                over += 1
                if first_over is None:
                    first_over = r

        safe = bool(source_rows) and over == 0
        if safe and widest is None:
            widest = angle
        rows.append({
            "angle_deg": angle,
            "candidate_q_upper": q0,
            "post_prediction_q_upper": qpred,
            "candidate_tangent_q_upper": c_tangent,
            "first_force_magnitude_mps2": gravity,
            "yaw_force_alignment_x": 0.0,
            "evaluated_source_phase_cells": len(source_rows),
            "children_above_6rad": over,
            "max_Ktheta_norm_upper": max_k,
            "max_combined_residual_norm_upper_mps2": max_resid,
            "max_first_accelerometer_correction_norm_upper_rad": max_d,
            "minimum_correction_range_margin_rad": min_margin,
            "first_accelerometer_range_safe": safe,
            "first_unclosed_child": first_over,
            "source_rows": source_rows,
        })

    all_safe = bool(rows) and all(r["first_accelerometer_range_safe"] for r in rows) and not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_CANDIDATE_FIRST_ACCEL_EXACT_SOURCE_GEOMETRY",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "source_semantics": semantics,
        "first_accel_aw_mean_exact_zero_before_measurement": True,
        "first_due_S_mean_correction_exact_zero": True,
        "first_accel_specific_force_magnitude_exact_gravity": True,
        "first_accel_specific_force_magnitude_mps2": gravity,
        "first_accel_yaw_covariance_axis_aligned_with_force_axis": True,
        "yaw_alignment_x_equals_zero": True,
        "candidate_full_q_used_as_conservative_gravity_tangent_q": True,
        "latent_linear_plus_rotation_cross_combined_before_norm": True,
        "independent_latent_rotation_cross_norm_added": False,
        "V12D_tangent_PSD_resolvent_used": True,
        "PSD_remainder_axial_noise_floor_inverse_used": False,
        "H_bias_error_bound_contains_A": ba_H >= ba_A,
        "A_first_prefix_attitude_gain_bounded_by_H_gain": astruct.get("A_bias_innovation_addition_isotropic_PSD") is True,
        "A_structure_proved_before_generic_PSD_boxing": True,
        "A_mode_structure": astruct,
        "deployed_correction_limit_rad": LIMIT,
        "deployed_correction_limit_increased": False,
        "candidate_rows": rows,
        "widest_candidate_first_accel_range_safe_deg": widest,
        "all_candidate_first_accelerometer_ranges_safe": all_safe,
        "P4_CANDIDATE_FIRST_ACCEL_EXACT_SOURCE_CERTIFICATE": "PASS" if all_safe else "NOT_ESTABLISHED",
        "signed_correction_Joseph_reset_propagated_here": False,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "PROPAGATE_EXACT_SOURCE_FIRST_ACCEL_CORRECTION_THROUGH_SIGNED_JOSEPH_RESET"
            if all_safe else
            "REFINE_EXACT_SOURCE_FIRST_ACCEL_RESIDUAL_DIRECTION"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "first_accel_aw_mean_exact_zero_before_measurement",
        "first_due_S_mean_correction_exact_zero",
        "first_accel_specific_force_magnitude_exact_gravity",
        "first_accel_yaw_covariance_axis_aligned_with_force_axis",
        "yaw_alignment_x_equals_zero",
        "candidate_full_q_used_as_conservative_gravity_tangent_q",
        "latent_linear_plus_rotation_cross_combined_before_norm",
        "V12D_tangent_PSD_resolvent_used",
        "H_bias_error_bound_contains_A",
        "A_first_prefix_attitude_gain_bounded_by_H_gain",
        "A_structure_proved_before_generic_PSD_boxing",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "independent_latent_rotation_cross_norm_added",
        "PSD_remainder_axial_noise_floor_inverse_used", "deployed_correction_limit_increased",
        "signed_correction_Joseph_reset_propagated_here", "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
        "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction range changed")
    if [float(r.get("angle_deg", -1.0)) for r in d.get("candidate_rows", [])] != [30.0,25.0,20.0,15.0]:
        f.append("candidate ladder changed")
    st = d.get("P4_CANDIDATE_FIRST_ACCEL_EXACT_SOURCE_CERTIFICATE")
    if st == "PASS":
        if d.get("all_candidate_first_accelerometer_ranges_safe") is not True:
            f.append("PASS without all candidate ranges safe")
        if any(float(r["max_first_accelerometer_correction_norm_upper_rad"]) > 6.0 for r in d["candidate_rows"]):
            f.append("PASS above deployed correction range")
    elif st == "NOT_ESTABLISHED":
        if not any(r.get("first_unclosed_child") is not None for r in d.get("candidate_rows", [])) and not f:
            f.append("nonclosure lacks source witness")
    else:
        f.append("invalid candidate exact-source status")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--output", type=Path, required=True)
    x=ap.parse_args()
    d=build(x.domain.resolve(), source_pieces=x.source_pieces)
    vf=validate(d); d["validation_pass"]=not vf; d["validation_failures"]=vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps({
        "status":d["P4_CANDIDATE_FIRST_ACCEL_EXACT_SOURCE_CERTIFICATE"],
        "widest_safe_deg":d["widest_candidate_first_accel_range_safe_deg"],
        "force_mps2":d["first_accel_specific_force_magnitude_mps2"],
        "rows":[{
            "angle_deg":r["angle_deg"],"qpred":r["post_prediction_q_upper"],
            "K":r["max_Ktheta_norm_upper"],"residual":r["max_combined_residual_norm_upper_mps2"],
            "d":r["max_first_accelerometer_correction_norm_upper_rad"],
            "margin":r["minimum_correction_range_margin_rad"],"safe":r["first_accelerometer_range_safe"]
        } for r in d["candidate_rows"]],
        "validation_failures":vf},indent=2,sort_keys=True))
    return 0 if not vf else 2


if __name__=="__main__": raise SystemExit(main())
