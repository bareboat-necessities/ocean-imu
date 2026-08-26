#!/usr/bin/env python3
"""One-plus-two block resolvent refinement for OU-III P5 sample 1.

V12B correctly restored the unconditional shipping innovation floor
``S' >= R_a I`` but still bounded the perturbation of the *full* six-state
Kalman gain and multiplied that numerator by ``1/R_a``.  The 24^3 run showed
that this is still far too coarse: the physical V10 attitude correction is
small while a tiny covariance/Jacobian perturbation is amplified by the raw
measurement-noise inverse.

The nominal sample-1 innovation in V7's un-Rx gauge is exactly a scalar block
plus a 2x2 block.  Both have cancellation-free positive identities.  For the
scalar block V7 gives ``S_x>0``.  For the 2x2 SPD block,

    lambda_min(S_yz) >= det(S_yz) / trace(S_yz),

with V7's positive determinant identity.  Hence a source-cell nominal inverse
bound is available without interval LDLT.  If ``e=||Delta S||`` and
``||S0^-1|| e < 1``, the exact resolvent gives

    ||S'^-1|| <= ||S0^-1|| / (1-||S0^-1|| e).

The unconditional ``1/R_a`` bound remains a fallback, so no positivity
assumption is introduced.

Only the attitude rows of the gain are needed for the q<8 capture gate.  With
``C_theta=P_theta H_theta^T+P_theta,aw`` and perturbation bounds
``||Delta P_block||<=dP``, ``||Delta H_theta||<=dH``,

    ||Delta C_theta||
      <= dP ||H_theta|| + ||P_theta|| dH + dP dH + dP.

The last term is the theta/aw cross-block perturbation.  Therefore

    ||Delta K_theta||
      <= (||Delta C_theta|| + ||K_theta|| ||Delta S||) ||S'^-1||.

All V11 PSD/S residual and covariance perturbations are retained.  No source
bound, filter parameter, deployed six-radian limit, q-chart target, or theorem
promotion is changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_structured_full_gain_v11 as V11

DEFAULT_DOMAIN = V11.DEFAULT_DOMAIN
SCHEMA = 1202
RANGE = V11.RANGE
FULL = V11.FULL


def _sum_up(*xs: float) -> float:
    y = 0.0
    for x in xs:
        y = FULL.up(y + float(x))
    return y


def _nominal_block_inverse_upper(*, a: Interval, Y: Interval, c0: Interval,
                                 alpha: Interval, qaw: Interval,
                                 b: Interval, bz: Interval,
                                 det_first: Interval, d: Interval,
                                 fy: Interval, fz: Interval,
                                 r: Interval) -> dict:
    """Cancellation-free inverse bound for V7's scalar plus 2x2 innovation."""
    h = FULL.I(0.5) * d
    Bt = alpha.square() * b + qaw
    Bz = alpha.square() * bz + qaw
    delta = alpha.square() * det_first + a * qaw
    if not (a.lo > 0.0 and Y.lo > 0.0 and Bz.lo > 0.0 and delta.lo > 0.0 and r.lo > 0.0):
        raise RuntimeError("positive V12C canonical covariance/noise floors required")

    U = fz - h * fy
    V = -(h * fz) - fy
    Nu = a * U + alpha * c0
    Sx = Nu.square() / a + delta / a + Y * V.square() + r
    if not Sx.lo > 0.0:
        raise RuntimeError("V12C scalar innovation lost positive floor")

    cx = -(alpha * c0)
    q = cx - a * fz
    A = delta + a * r
    det = fy.square() * A + (q.square() + A) * (Bz + r) / a
    trace = (q.square() + A) / a + a * fy.square() + Bz + r
    if not (A.lo > 0.0 and det.lo > 0.0 and trace.lo > 0.0):
        raise RuntimeError("V12C two-by-two innovation lost positivity")
    lam2 = FULL.down(det.lo / trace.hi)
    lam0 = min(Sx.lo, lam2)
    if not lam0 > 0.0:
        raise RuntimeError("V12C nominal block eigenvalue floor nonpositive")
    inv0 = FULL.up(1.0 / FULL.down(lam0))
    return {
        "scalar_innovation_lower": Sx.lo,
        "two_by_two_determinant_lower": det.lo,
        "two_by_two_trace_upper": trace.hi,
        "two_by_two_lambda_min_lower": lam2,
        "nominal_block_lambda_min_lower": lam0,
        "nominal_block_inverse_operator_upper": inv0,
    }


def _attitude_gain_perturbation(*, dP: float, dH: float,
                                ptheta_norm: float, htheta_norm: float,
                                ktheta_norm: float, dS: float,
                                nominal_inverse_upper: float,
                                r_floor: float) -> dict:
    vals = (dP, dH, ptheta_norm, htheta_norm, ktheta_norm,
            dS, nominal_inverse_upper, r_floor)
    if not all(math.isfinite(float(x)) for x in vals):
        raise RuntimeError("finite V12C gain inputs required")
    if min(dP, dH, ptheta_norm, htheta_norm, ktheta_norm, dS,
           nominal_inverse_upper) < 0.0 or not r_floor > 0.0:
        raise RuntimeError("invalid V12C gain inputs")

    dCtheta = _sum_up(
        FULL.up(dP * htheta_norm),
        FULL.up(ptheta_norm * dH),
        FULL.up(dP * dH),
        dP,
    )
    beta = FULL.up(nominal_inverse_upper * dS)
    noise_inv = FULL.up(1.0 / FULL.down(r_floor))
    if beta < 1.0:
        denom = FULL.down(1.0 - beta)
        structured_inv = FULL.up(nominal_inverse_upper / denom)
        actual_inv = min(noise_inv, structured_inv)
        backend = "one_plus_two_neumann_intersect_noise_floor"
    else:
        structured_inv = math.inf
        actual_inv = noise_inv
        backend = "noise_floor_fallback"
    numerator = _sum_up(dCtheta, FULL.up(ktheta_norm * dS))
    dk = FULL.up(numerator * actual_inv)
    return {
        "sample1_attitude_cross_covariance_perturbation_upper": dCtheta,
        "nominal_inverse_times_innovation_perturbation_upper": beta,
        "structured_actual_innovation_inverse_operator_upper": structured_inv,
        "noise_floor_actual_innovation_inverse_operator_upper": noise_inv,
        "actual_innovation_inverse_operator_upper": actual_inv,
        "actual_inverse_backend": backend,
        "sample1_attitude_gain_operator_perturbation_upper": dk,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24) -> dict:
    V11.FULL3._install_backend()
    path = Path(domain_path).resolve()
    dom = json.loads(path.read_text(encoding="utf-8"))
    core = V11.V10.build(path, source_pieces=source_pieces,
                         source_cell_index=source_cell_index,
                         p_pieces=p_pieces, tangent_pieces=tangent_pieces,
                         axial_pieces=axial_pieces)
    first = V11.FIRST.build(path, source_pieces=source_pieces)
    vec = V11.VECTOR.build()
    failures = [f"V10: {x}" for x in V11.V10.validate(core)]
    failures += [f"first: {x}" for x in V11.FIRST.validate(first)]
    failures += [f"vector: {x}" for x in V11.VECTOR.validate(vec)]
    if core.get("P5_SAMPLE1_COMBINED_PERPENDICULAR_RESIDUAL_V10") != "PASS":
        failures.append("V10 prerequisite did not pass")

    src, phase = V11.RG._source_phase_children(source_pieces)[source_cell_index]
    if phase != "due":
        failures.append("V12C focused perturbation requires first due source cell")
    if source_cell_index != 0:
        failures.append("V12C S perturbation helper currently certified for source cell 0")

    fr = first["source_cells"][source_cell_index]
    p_all = Interval.outward_bounds(*map(float, fr["P_aw_variance_interval"]))
    rho0 = float(fr["combined_useful_residual_norm_upper_mps2"])
    aw_pre = float(fr["predicted_aw_error_norm_upper_mps2"])
    hstep = float(src["dt_s"])
    g = float(dom["startup"]["gravity_mps2"])
    tilt, yaw, eps = V11.RG._attitude_covariance_epsilon(path, hstep)
    t = Interval.outward_bounds(tilt, FULL.up(tilt + eps))
    Y = Interval.outward_bounds(yaw, FULL.up(yaw + eps))
    r = FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))[0][0]
    F, Q, _ = FULL._transition_and_Q(src, dom)
    alpha = F[15][15]
    alpha_hi = alpha.hi
    qaw = Q[15][15]
    pcells = V11.SUB.parts(p_all.lo, p_all.hi, p_pieces)
    sb = V11._sample1_s_bounds(path, src, dom, first, vec)

    rows = []
    bad = None
    worst = None
    unclosed = 0
    fallback = 0
    max_d = max_dP = max_dH = max_dS = max_drho = max_dk = 0.0
    max_inv = max_beta = 0.0
    min_nominal_lambda = math.inf

    for base in core["rows"]:
        pi = int(base["p_cell"])
        p = pcells[pi]
        rt = Interval.outward_bounds(*map(float, base["first_tangent_residual_magnitude_mps2"]))
        rz = Interval.outward_bounds(*map(float, base["first_axial_residual_mps2"]))
        d = Interval.outward_bounds(*map(float, base["first_attitude_correction_rad"]))
        D = FULL.I(g * g) * t + p + r
        a = t * (p + r) / D
        c0 = -(FULL.I(g) * t * p / D)
        b = p * (FULL.I(g * g) * t + r) / D
        bz = p * r / (p + r)
        det_first = t * p * r / D
        fy = -(alpha * (p / D) * rt)
        fz = FULL.I(g) + alpha * (p / (p + r)) * rz

        psd = V11._first_psd_perturbation(
            t=t, Y=Y, p=p, r=r, g=g, eps=eps, rho0=rho0,
            dhi=max(0.0, d.hi), rt=rt, rz=rz, alpha_hi=alpha_hi,
            aw_pre=aw_pre)
        vecnorm = FULL.up(
            g + FULL.up(alpha_hi * (
                aw_pre
                + psd["first_nominal_aw_correction_norm_upper_mps2"]
                + psd["first_aw_x_correction_upper_mps2"])))
        drS = FULL.up(
            sb["sample1_S_aw_correction_upper_mps2"]
            + FULL.up(sb["sample1_S_attitude_correction_upper_rad"] * vecnorm))
        drho = FULL.up(psd["PSD_induced_sample1_residual_perturbation_upper_mps2"] + drS)
        dP = FULL.up(
            psd["sample1_reduced_covariance_PSD_perturbation_upper"]
            + sb["sample1_S_total_reduced_covariance_perturbation_upper"])
        dH = drho

        Pn, Hn, _Sn = V11._nominal_sample1_matrices(
            t=t, Y=Y, p=p, r=r, g=g, alpha=alpha, qaw=qaw,
            d=d, fy=fy, fz=fz)
        Ptheta = [row[:3] for row in Pn[:3]]
        ptheta_norm = V11._op(Ptheta)
        htheta_norm = V11.V5._norm2_upper(fy.abs_upper(), fz.abs_upper())
        hnorm = V11._op(Hn)
        pnorm = V11._op(Pn)

        # Same perturbation of S as V12B, but its inverse is bounded from the
        # exact nominal 1+2 blocks before intersecting with the unconditional
        # shipping noise-floor inverse.
        dS = _sum_up(
            FULL.up(FULL.up(hnorm * hnorm) * dP),
            FULL.up(FULL.up(2.0 * hnorm * pnorm) * dH),
            FULL.up(pnorm * FULL.up(dH * dH)),
            FULL.up(FULL.up(2.0 * hnorm * dP) * dH),
            FULL.up(dP * FULL.up(dH * dH)),
        )
        inv = _nominal_block_inverse_upper(
            a=a, Y=Y, c0=c0, alpha=alpha, qaw=qaw, b=b, bz=bz,
            det_first=det_first, d=d, fy=fy, fz=fz, r=r)
        k0 = max(float(base["Ktheta_perpendicular_block_upper"]),
                 float(base["Ktheta_parallel_block_upper"]))
        gp = _attitude_gain_perturbation(
            dP=dP, dH=dH, ptheta_norm=ptheta_norm,
            htheta_norm=htheta_norm, ktheta_norm=k0, dS=dS,
            nominal_inverse_upper=float(inv["nominal_block_inverse_operator_upper"]),
            r_floor=r.lo)

        rho = float(base["sample1_full_residual_norm_upper_mps2"])
        d10 = float(base["combined_directional_correction_norm_upper_rad"])
        dk = float(gp["sample1_attitude_gain_operator_perturbation_upper"])
        d12c = FULL.up(
            d10
            + FULL.up(FULL.up(k0 * drho)
                      + FULL.up(dk * FULL.up(rho + drho))))
        closed = math.isfinite(d12c) and d12c < RANGE
        row = {
            "p_cell": pi,
            "tangent_residual_cell": base["tangent_residual_cell"],
            "axial_residual_cell": base["axial_residual_cell"],
            "V10_directional_correction_upper_rad": d10,
            "PSD_residual_perturbation_upper_mps2": psd["PSD_induced_sample1_residual_perturbation_upper_mps2"],
            "S_residual_perturbation_upper_mps2": drS,
            "total_residual_perturbation_upper_mps2": drho,
            "PSD_reduced_covariance_perturbation_upper": psd["sample1_reduced_covariance_PSD_perturbation_upper"],
            "S_reduced_covariance_perturbation_upper": sb["sample1_S_total_reduced_covariance_perturbation_upper"],
            "total_reduced_covariance_perturbation_upper": dP,
            "nominal_reduced_covariance_operator_upper": pnorm,
            "nominal_attitude_covariance_operator_upper": ptheta_norm,
            "sample1_H_operator_upper": hnorm,
            "sample1_Htheta_operator_upper": htheta_norm,
            "sample1_H_perturbation_upper": dH,
            "sample1_innovation_perturbation_upper": dS,
            "actual_innovation_noise_floor_lower": r.lo,
            **inv,
            **gp,
            "V12C_correction_norm_upper_rad": d12c,
            "inside_9rad_range": closed,
            **psd,
        }
        rows.append(row)
        if gp["actual_inverse_backend"] == "noise_floor_fallback":
            fallback += 1
        max_d = max(max_d, d12c)
        max_dP = max(max_dP, dP)
        max_dH = max(max_dH, dH)
        max_dS = max(max_dS, dS)
        max_drho = max(max_drho, drho)
        max_dk = max(max_dk, dk)
        max_inv = max(max_inv, float(gp["actual_innovation_inverse_operator_upper"]))
        max_beta = max(max_beta, float(gp["nominal_inverse_times_innovation_perturbation_upper"]))
        min_nominal_lambda = min(min_nominal_lambda, float(inv["nominal_block_lambda_min_lower"]))
        if worst is None or d12c > worst["V12C_correction_norm_upper_rad"]:
            worst = row
        if not closed:
            unclosed += 1
            if bad is None:
                bad = row

    ok = bool(rows) and unclosed == 0 and bad is None and not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_ONE_PLUS_TWO_ATTITUDE_RESOLVENT_V12C",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V10_canonical_core_retained": True,
        "V11_PSD_and_S_perturbation_magnitudes_retained": True,
        "nominal_one_plus_two_innovation_structure_used": True,
        "scalar_positive_innovation_identity_used": True,
        "two_by_two_positive_determinant_identity_used": True,
        "two_by_two_lambda_floor_det_over_trace_used": True,
        "attitude_gain_rows_only_bounded_for_q_gate": True,
        "actual_shipping_covariance_PSD_used": True,
        "actual_innovation_noise_floor_fallback_retained": True,
        "sample1_force_Jacobian_perturbation_included": True,
        "first_attitude_PSD_cross_axis_remainder_included": True,
        "second_prediction_attitude_process_remainder_included": True,
        "sample1_S_covariance_update_included": True,
        "sample1_S_attitude_injection_included": True,
        "sample1_S_aw_mean_correction_included": True,
        "sample1_S_solver_identity_branch_contained_as_zero_perturbation": True,
        "generic_interval_innovation_LDLT_floor_used": False,
        "broad_sample1_3x3_interval_inverse_reintroduced": False,
        "temporal_force_slew_assumed": False,
        "complete_sample1_branch_closed_here": False,
        "signed_cayley_q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "validated_deployed_quaternion_range_rad": RANGE,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "attitude_covariance_remainder_spectral_upper": eps,
        "sample1_S_perturbation_bounds": sb,
        "evaluated_joint_cells": len(rows),
        "unclosed_joint_cells": unclosed,
        "noise_floor_inverse_fallback_cells": fallback,
        "minimum_nominal_block_lambda_lower": min_nominal_lambda,
        "max_nominal_inverse_times_innovation_perturbation_upper": max_beta,
        "max_actual_innovation_inverse_operator_upper": max_inv,
        "max_total_residual_perturbation_upper_mps2": max_drho,
        "max_total_reduced_covariance_perturbation_upper": max_dP,
        "max_sample1_H_perturbation_upper": max_dH,
        "max_sample1_innovation_perturbation_upper": max_dS,
        "max_sample1_attitude_gain_operator_perturbation_upper": max_dk,
        "max_V12C_correction_norm_upper_rad": max_d,
        "first_unclosed_joint_cell": bad,
        "worst_joint_cell": worst,
        "P5_SAMPLE1_ONE_PLUS_TWO_ATTITUDE_RESOLVENT_V12C": "PASS" if ok else "NOT_ESTABLISHED",
        "next_obligation": (
            "SIGNED_RADIAL_SUBDIVIDE_AND_CAYLEY_COMPOSE_SAMPLE1_INSIDE_Q8"
            if ok else "REFINE_FIRST_PSD_OFFAXIS_DIRECTION_OR_BLOCKWISE_INNOVATION_PERTURBATION"
        ),
        "failures": failures,
        "rows": rows,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    for k in (
        "source_generated_not_trajectory_fit", "V10_canonical_core_retained",
        "V11_PSD_and_S_perturbation_magnitudes_retained",
        "nominal_one_plus_two_innovation_structure_used",
        "scalar_positive_innovation_identity_used",
        "two_by_two_positive_determinant_identity_used",
        "two_by_two_lambda_floor_det_over_trace_used",
        "attitude_gain_rows_only_bounded_for_q_gate",
        "actual_shipping_covariance_PSD_used",
        "actual_innovation_noise_floor_fallback_retained",
        "sample1_force_Jacobian_perturbation_included",
        "first_attitude_PSD_cross_axis_remainder_included",
        "second_prediction_attitude_process_remainder_included",
        "sample1_S_covariance_update_included",
        "sample1_S_attitude_injection_included",
        "sample1_S_aw_mean_correction_included",
        "sample1_S_solver_identity_branch_contained_as_zero_perturbation",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "generic_interval_innovation_LDLT_floor_used",
        "broad_sample1_3x3_interval_inverse_reintroduced", "temporal_force_slew_assumed",
        "complete_sample1_branch_closed_here", "signed_cayley_q8_composed_here",
        "q8_word_promoted_here", "whole_word_promoted_here", "N_H_words_set_here",
        "deployed_correction_limit_increased",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if int(d.get("evaluated_joint_cells", 0)) <= 0:
        f.append("no V12C cells")
    for k in (
        "minimum_nominal_block_lambda_lower",
        "max_actual_innovation_inverse_operator_upper",
        "max_sample1_attitude_gain_operator_perturbation_upper",
        "max_V12C_correction_norm_upper_rad",
    ):
        if not math.isfinite(float(d.get(k, math.nan))):
            f.append(f"nonfinite {k}")
    if float(d.get("minimum_nominal_block_lambda_lower", 0.0)) <= 0.0:
        f.append("nominal one-plus-two block lost positive eigenvalue floor")
    st = d.get("P5_SAMPLE1_ONE_PLUS_TWO_ATTITUDE_RESOLVENT_V12C")
    w = d.get("first_unclosed_joint_cell")
    if st == "PASS" and (w is not None or int(d.get("unclosed_joint_cells", -1)) != 0):
        f.append("V12C PASS retains unclosed cell")
    if st == "NOT_ESTABLISHED" and w is None and not f:
        f.append("missing V12C witness")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=24)
    ap.add_argument("--tangent-pieces", type=int, default=24)
    ap.add_argument("--axial-pieces", type=int, default=24)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(x.domain, source_pieces=x.source_pieces,
              source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
              tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_ONE_PLUS_TWO_ATTITUDE_RESOLVENT_V12C"],
        "cells": d["evaluated_joint_cells"],
        "unclosed": d["unclosed_joint_cells"],
        "noise_floor_fallback_cells": d["noise_floor_inverse_fallback_cells"],
        "min_nominal_lambda": d["minimum_nominal_block_lambda_lower"],
        "max_beta": d["max_nominal_inverse_times_innovation_perturbation_upper"],
        "max_actual_inverse": d["max_actual_innovation_inverse_operator_upper"],
        "max_drho": d["max_total_residual_perturbation_upper_mps2"],
        "max_dP": d["max_total_reduced_covariance_perturbation_upper"],
        "max_dH": d["max_sample1_H_perturbation_upper"],
        "max_dS": d["max_sample1_innovation_perturbation_upper"],
        "max_dKtheta": d["max_sample1_attitude_gain_operator_perturbation_upper"],
        "max_d": d["max_V12C_correction_norm_upper_rad"],
        "first_unclosed": d["first_unclosed_joint_cell"],
        "worst": d["worst_joint_cell"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
