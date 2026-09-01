#!/usr/bin/env python3
"""Analytic first-accelerometer correction range for finite-angle P4 candidates.

The generic interval prefix backend deliberately encloses the accelerometer
orientation matrix entrywise.  At the first Live vector update that loses two
source facts which are available without replay:

* J_aw=R_wb is orthogonal, so the rotation gauge sends it exactly to I;
* the finite attitude residual is one rotated-vector difference, not a linear
  term plus an independently charged nonlinear remainder.

For ||c||<=q after the first 5 ms transport step,

    ||(R(c)-I) f|| <= 2 q/sqrt(4+q^2) ||f||.

The latent term is combined before taking a norm,

    e_aw + (R(c)^T-I)e_aw = R(c)^T e_aw,

so its norm is exactly ||e_aw||.  The remaining accelerometer-bias error is
additive.  The first-prefix attitude gain is bounded by the existing analytic
rank-two rotation-gauge certificate, retaining the source-varying yaw axis,
force magnitude, source/tuner cells, optional first S pseudo phase, and the
small PSD attitude-covariance remainder.

The H-mode bound also covers the first A-mode attitude correction: the A seed
adds an isotropic positive accelerometer-bias covariance to the innovation,
while the first-prefix theta/ba cross covariance is zero.  Thus the numerator
of K_theta is unchanged and the canonical innovation denominator only grows.
The deterministic H-mode 0.5 m/s^2 bias-error norm also contains the declared
A-mode 0.45 m/s^2 state-error ball.

This is a range certificate only.  It does not replace the signed correction
cell needed for Joseph/reset propagation and does not promote P4 dissipation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_candidate_full_word as CAND
import ou3_p4_p5_entrance_search_domain as ENTRANCE
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_first_accel_rotation_gauge_v3 as RG3
import ou3_p5_first_accel_structured_gain as SG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = CAND.DEFAULT_DOMAIN
SCHEMA = 1
DEPLOYED_CORRECTION_LIMIT_RAD = 6.0


def _rotation_residual_gain_upper(q: float) -> float:
    q = float(q)
    if not (math.isfinite(q) and q >= 0.0):
        raise ValueError("finite nonnegative Cayley radius required")
    den = FULL.down(math.sqrt(FULL.down(4.0 + FULL.down(q * q))))
    if den <= 0.0:
        raise RuntimeError("finite-angle residual denominator lost positivity")
    return FULL.up(FULL.up(2.0 * q) / den)


def _a_mode_first_prefix_isotropic_bias_check(domain_path: Path, domain: dict, src: dict) -> dict:
    CAND._configure_mode("A")
    F, Q, _Rstep, _ba_process = CAND._transition_and_Q("A", src, domain)
    P0 = CAND._initial_covariance("A", src, domain_path)
    Pp = FULL._psd_tighten(FULL.matrix_add(FULL.matrix_mul(FULL.matrix_mul(F, P0), FULL.matrix_transpose(F)), Q))
    ba_diag = [Pp[i][i] for i in CAND.H.BA]
    isotropic = all(x.lo == ba_diag[0].lo and x.hi == ba_diag[0].hi for x in ba_diag[1:])
    cross_zero = True
    for i in CAND.H.BA:
        for j in list(CAND.H.TH) + list(CAND.H.AW) + list(CAND.H.SS):
            z = Pp[i][j]
            if z.lo != 0.0 or z.hi != 0.0:
                cross_zero = False
                break
    CAND._configure_mode("H")
    return {
        "A_bias_predicted_variance_interval": ba_diag[0].as_list(),
        "A_bias_innovation_addition_isotropic_PSD": isotropic and ba_diag[0].lo >= 0.0,
        "first_prefix_theta_aw_S_to_ba_cross_exact_zero": cross_zero,
    }


def build(
    domain_path: Path = DEFAULT_DOMAIN,
    *,
    source_pieces: int = 2,
    alignment_pieces: int = 16,
    force_magnitude_pieces: int = 4,
) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("candidate first-accel range domain must not be trajectory fitted")
    if domain.get("configured_runtime", {}).get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("candidate first-accel range requires lever arm disabled")

    RG3._install_backend(domain_path, source_pieces)
    FULL3._install_backend()
    entrance = ENTRANCE.build(domain_path)
    vector = VECTOR.build()
    failures = [f"entrance: {x}" for x in ENTRANCE.validate(entrance)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]

    src_phases = RG._source_phase_children(source_pieces)
    if not src_phases:
        failures.append("no first-prefix source phase cells")
    xcells = SG._linear_cells(alignment_pieces)
    live = domain["normal_live"]
    force_cells = RG._geom_ranges(
        float(live["specific_force_norm_lower_mps2"]),
        float(live["specific_force_norm_upper_mps2"]),
        force_magnitude_pieces,
    )
    h = float(FULL._source_cell()["dt_s"])
    tilt, yaw, eps = RG._attitude_covariance_epsilon(domain_path, h)
    vc = vector["configured_measurement_bounds"]
    Racc = FULL._R_diag(float(vc["acc_measurement_std_mps2"]))
    racc_var = Racc[0][0]

    # H-mode physical bias error is the larger of the two declared first-prefix
    # bias-error sets, so one residual bound covers H and A.
    ba_H = float(domain["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    ba_A = float(domain["normal_live"]["active_accelerometer_bias_state_norm_upper_mps2"])
    if not (0.0 <= ba_A <= ba_H):
        failures.append("H-mode accelerometer-bias error no longer contains A-mode entrance ball")

    a_check = _a_mode_first_prefix_isotropic_bias_check(domain_path, domain, src_phases[0][0]) if src_phases else {}
    if a_check.get("A_bias_innovation_addition_isotropic_PSD") is not True:
        failures.append("A first-prefix bias covariance is not an isotropic PSD innovation addition")
    if a_check.get("first_prefix_theta_aw_S_to_ba_cross_exact_zero") is not True:
        failures.append("A first-prefix bias gained forbidden cross covariance")

    rows = []
    widest_safe = None
    for crow in entrance["P4_complete_word_search"]["candidate_rows"]:
        angle = float(crow["angle_deg"])
        q0 = float(crow["cayley_norm_upper"])
        qpred = RG._q_after_first_prediction(q0, domain, h)
        rot_gain = _rotation_residual_gain_upper(qpred)
        total = 0
        over = 0
        max_d = 0.0
        min_margin = math.inf
        max_k = 0.0
        max_residual = 0.0
        first_over = None

        for si, (src, phase) in enumerate(src_phases):
            P0 = FULL._initial_covariance(src, domain_path)
            F, Q, _ = FULL._transition_and_Q(src, domain)
            Pp = FULL._psd_tighten(FULL.matrix_add(FULL.matrix_mul(FULL.matrix_mul(F, P0), FULL.matrix_transpose(F)), Q))
            _pss, _psa, paw_pred = RG._scalar_axis_structure(Pp)
            aw_pred, eS_pred = RG._prediction_norms(src, domain)
            if phase == "due":
                paw, aw_norm = RG._due_paw_and_error_norm(Pp, src, aw_pred, eS_pred)
            else:
                paw, aw_norm = paw_pred, aw_pred

            for xi, x in enumerate(xcells):
                for mi, m in enumerate(force_cells):
                    k, _kh, detail = SG._structured_gain_bounds(
                        tilt=tilt,
                        yaw=yaw,
                        eps=eps,
                        x=x,
                        m=m,
                        paw=paw,
                        racc_var=racc_var,
                    )
                    rotational = FULL.up(rot_gain * m.hi)
                    residual = FULL.up(rotational + FULL.up(aw_norm + ba_H))
                    d = FULL.up(k * residual)
                    total += 1
                    max_k = max(max_k, k)
                    max_residual = max(max_residual, residual)
                    max_d = max(max_d, d)
                    margin = DEPLOYED_CORRECTION_LIMIT_RAD - d
                    min_margin = min(min_margin, margin)
                    if not math.isfinite(d) or d > DEPLOYED_CORRECTION_LIMIT_RAD:
                        over += 1
                        if first_over is None:
                            first_over = {
                                "source_phase_cell": si,
                                "pseudo_phase": phase,
                                "alignment_cell": xi,
                                "alignment_x_tangent_yaw_fraction": x.as_list(),
                                "force_magnitude_cell": mi,
                                "force_magnitude_mps2": m.as_list(),
                                "tau_s": src["tau_s"].as_list(),
                                "sigma_aw_mps2": src["sigma_aw_mps2"].as_list(),
                                "R_S_filter_std": src["R_S_filter_std"].as_list(),
                                "predicted_aw_error_norm_upper_mps2": aw_norm,
                                "finite_rotation_residual_norm_upper_mps2": rotational,
                                "combined_residual_norm_upper_mps2": residual,
                                "Ktheta_norm_upper": k,
                                "correction_norm_upper_rad": d,
                                "gain_detail": detail,
                            }

        safe = total > 0 and over == 0
        if safe and widest_safe is None:
            widest_safe = angle
        rows.append({
            "angle_deg": angle,
            "candidate_q_upper": q0,
            "post_prediction_q_upper": qpred,
            "finite_rotation_residual_gain_upper": rot_gain,
            "evaluated_children": total,
            "children_above_validated_correction_limit": over,
            "max_Ktheta_norm_upper": max_k,
            "max_combined_residual_norm_upper_mps2": max_residual,
            "max_first_accelerometer_correction_norm_upper_rad": max_d,
            "minimum_correction_range_margin_rad": min_margin,
            "first_accelerometer_range_safe": safe,
            "first_unclosed_child": first_over,
        })

    all_safe = bool(rows) and all(r["first_accelerometer_range_safe"] for r in rows) and not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_CANDIDATE_FIRST_ACCEL_ANALYTIC_RANGE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "deployed_correction_limit_rad": DEPLOYED_CORRECTION_LIMIT_RAD,
        "deployed_correction_limit_increased": False,
        "candidate_ball_norm_used_without_cartesian_cover_inflation": True,
        "rotation_gauge_sets_J_aw_to_identity": True,
        "analytic_rank_two_gain_used": True,
        "finite_rotation_residual_used_directly": True,
        "finite_rotation_residual_bound": "||(R-I)f||<=2q/sqrt(4+q^2)||f||",
        "latent_linear_plus_rotation_cross_combined_before_norm": True,
        "latent_combined_norm_identity": "||e+(R^T-I)e||=||e||",
        "independent_accelerometer_eta_penalty_used": False,
        "H_bias_error_norm_upper_mps2": ba_H,
        "A_bias_error_norm_upper_mps2": ba_A,
        "H_bias_error_bound_contains_A": ba_H >= ba_A,
        "A_first_prefix_attitude_gain_bounded_by_H_gain": not failures and a_check.get("A_bias_innovation_addition_isotropic_PSD") is True,
        "A_mode_structure": a_check,
        "candidate_rows": rows,
        "widest_candidate_first_accel_range_safe_deg": widest_safe,
        "all_candidate_first_accelerometer_ranges_safe": all_safe,
        "P4_CANDIDATE_FIRST_ACCEL_RANGE_CERTIFICATE": "PASS" if all_safe else "NOT_ESTABLISHED",
        "signed_correction_Joseph_reset_propagated_here": False,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "INTERSECT_FIRST_ACCEL_SIGNED_CORRECTION_CELL_WITH_CERTIFIED_NORM_RANGE_AND_PROPAGATE_JOSEPH_RESET"
            if all_safe else
            "REFINE_CANDIDATE_FIRST_ACCEL_FORCE_ALIGNMENT_OR_SOURCE_COVARIANCE"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "candidate_ball_norm_used_without_cartesian_cover_inflation",
        "rotation_gauge_sets_J_aw_to_identity",
        "analytic_rank_two_gain_used",
        "finite_rotation_residual_used_directly",
        "latent_linear_plus_rotation_cross_combined_before_norm",
        "H_bias_error_bound_contains_A",
        "A_first_prefix_attitude_gain_bounded_by_H_gain",
    ):
        if d.get(k) is not True:
            failures.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "independent_accelerometer_eta_penalty_used", "signed_correction_Joseph_reset_propagated_here",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE", "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(k) is not False:
            failures.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        failures.append("deployed correction range changed")
    rows = d.get("candidate_rows", [])
    if [float(r.get("angle_deg", -1.0)) for r in rows] != [30.0, 25.0, 20.0, 15.0]:
        failures.append("candidate ladder changed")
    if d.get("P4_CANDIDATE_FIRST_ACCEL_RANGE_CERTIFICATE") == "PASS":
        if d.get("all_candidate_first_accelerometer_ranges_safe") is not True:
            failures.append("PASS without all candidate first-accel ranges safe")
        if any(float(r["max_first_accelerometer_correction_norm_upper_rad"]) > 6.0 for r in rows):
            failures.append("PASS with correction bound above deployed range")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--alignment-pieces", type=int, default=16)
    ap.add_argument("--force-magnitude-pieces", type=int, default=4)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(
        args.domain.resolve(),
        source_pieces=args.source_pieces,
        alignment_pieces=args.alignment_pieces,
        force_magnitude_pieces=args.force_magnitude_pieces,
    )
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P4_CANDIDATE_FIRST_ACCEL_RANGE_CERTIFICATE"],
        "widest_safe_deg": d["widest_candidate_first_accel_range_safe_deg"],
        "rows": [{
            "angle_deg": r["angle_deg"],
            "q_pred": r["post_prediction_q_upper"],
            "max_Ktheta": r["max_Ktheta_norm_upper"],
            "max_residual": r["max_combined_residual_norm_upper_mps2"],
            "max_d": r["max_first_accelerometer_correction_norm_upper_rad"],
            "margin": r["minimum_correction_range_margin_rad"],
            "safe": r["first_accelerometer_range_safe"],
        } for r in d["candidate_rows"]],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
