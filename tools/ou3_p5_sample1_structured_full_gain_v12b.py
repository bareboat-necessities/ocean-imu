#!/usr/bin/env python3
"""PSD-preserving refinement of the OU-III P5 sample-1 perturbation closure.

V12 attempted to improve V11 by certifying a lower eigenvalue of an interval
3x3 innovation enclosure and then subtracting ||Delta S|| from that bound.  The
24^3 certificate showed that this is the wrong abstraction: the interval
innovation can lose its positive pivot even though every *real* shipping
innovation is

    S' = H' P' H'^T + R_a >= R_a I,

because Joseph updates, left-error reset congruences, prediction, and the S
measurement all leave the actual covariance P' positive semidefinite.

This producer retains the complete V11 PSD/S perturbation magnitudes and V10's
source-correlated 1+2 nominal correction.  It changes only the final gain
perturbation calculus.  For C=P H^T, S=H P H^T+R and the actual perturbed
C',S', the exact resolvent identity is

    K' - K = (Delta C - K Delta S) S'^(-1),

so ||S'^(-1)|| <= 1/R_a is available *without* a small-denominator condition.
The small PSD/S mean perturbations also change the accelerometer attitude
Jacobian through the predicted force; that Jacobian perturbation is charged
explicitly here rather than silently holding H fixed.

For ||Delta P||<=e, ||Delta H||<=h_e, ||P||<=p and ||H||<=h,

    ||Delta C|| <= e h + p h_e + e h_e,
    ||Delta S|| <= h^2 e + 2 h p h_e + p h_e^2
                   + 2 h e h_e + e h_e^2.

The actual innovation noise floor then gives a finite, source-safe gain
perturbation for every child.  No source domain, filter parameter, correction
limit, q-chart gate, or theorem promotion is changed here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_structured_full_gain_v11 as V11

DEFAULT_DOMAIN = V11.DEFAULT_DOMAIN
SCHEMA = 1201
RANGE = V11.RANGE
FULL = V11.FULL


def _sum_up(*xs: float) -> float:
    y = 0.0
    for x in xs:
        y = FULL.up(y + float(x))
    return y


def _gain_perturbation_psd_floor(*, dP: float, dH: float, pnorm: float,
                                 hnorm: float, ktheta_norm: float,
                                 r_floor: float) -> dict:
    vals = (dP, dH, pnorm, hnorm, ktheta_norm, r_floor)
    if not all(math.isfinite(float(x)) for x in vals):
        raise RuntimeError("finite PSD-resolvent inputs required")
    if min(dP, dH, pnorm, hnorm, ktheta_norm) < 0.0 or not r_floor > 0.0:
        raise RuntimeError("invalid PSD-resolvent signs/floor")

    dC = _sum_up(
        FULL.up(dP * hnorm),
        FULL.up(pnorm * dH),
        FULL.up(dP * dH),
    )
    dS = _sum_up(
        FULL.up(FULL.up(hnorm * hnorm) * dP),
        FULL.up(FULL.up(2.0 * hnorm * pnorm) * dH),
        FULL.up(pnorm * FULL.up(dH * dH)),
        FULL.up(FULL.up(2.0 * hnorm * dP) * dH),
        FULL.up(dP * FULL.up(dH * dH)),
    )
    numerator = _sum_up(dC, FULL.up(ktheta_norm * dS))
    dk = FULL.up(numerator / FULL.down(r_floor))
    return {
        "sample1_H_perturbation_upper": dH,
        "sample1_cross_covariance_perturbation_upper": dC,
        "sample1_innovation_perturbation_upper": dS,
        "actual_innovation_inverse_operator_upper": FULL.up(1.0 / FULL.down(r_floor)),
        "sample1_gain_operator_perturbation_upper": dk,
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
        failures.append("V12B focused perturbation requires first due source cell")
    if source_cell_index != 0:
        failures.append("V12B S perturbation helper currently certified for source cell 0")

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
    max_d = max_dP = max_dH = max_dS = max_dC = max_drho = max_dk = 0.0

    for base in core["rows"]:
        pi = int(base["p_cell"])
        p = pcells[pi]
        rt = Interval.outward_bounds(*map(float, base["first_tangent_residual_magnitude_mps2"]))
        rz = Interval.outward_bounds(*map(float, base["first_axial_residual_mps2"]))
        d = Interval.outward_bounds(*map(float, base["first_attitude_correction_rad"]))
        D = FULL.I(g * g) * t + p + r
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

        # H_theta=-[f_hat]x, hence ||Delta H_theta||_2=||Delta f_hat||_2.
        # The same source mean/gauge perturbations used in drho bound Delta f.
        dH = drho
        Pn, Hn, _Sn = V11._nominal_sample1_matrices(
            t=t, Y=Y, p=p, r=r, g=g, alpha=alpha, qaw=qaw,
            d=d, fy=fy, fz=fz)
        hnorm = V11._op(Hn)
        pnorm = V11._op(Pn)
        k0 = max(float(base["Ktheta_perpendicular_block_upper"]),
                 float(base["Ktheta_parallel_block_upper"]))
        gp = _gain_perturbation_psd_floor(
            dP=dP, dH=dH, pnorm=pnorm, hnorm=hnorm,
            ktheta_norm=k0, r_floor=r.lo)

        rho = float(base["sample1_full_residual_norm_upper_mps2"])
        d10 = float(base["combined_directional_correction_norm_upper_rad"])
        dk = float(gp["sample1_gain_operator_perturbation_upper"])
        d12b = FULL.up(
            d10
            + FULL.up(FULL.up(k0 * drho)
                      + FULL.up(dk * FULL.up(rho + drho))))
        closed = math.isfinite(d12b) and d12b < RANGE
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
            "sample1_H_operator_upper": hnorm,
            "actual_innovation_noise_floor_lower": r.lo,
            **gp,
            "V12B_correction_norm_upper_rad": d12b,
            "inside_9rad_range": closed,
            **psd,
        }
        rows.append(row)
        max_d = max(max_d, d12b)
        max_dP = max(max_dP, dP)
        max_dH = max(max_dH, dH)
        max_dS = max(max_dS, float(gp["sample1_innovation_perturbation_upper"]))
        max_dC = max(max_dC, float(gp["sample1_cross_covariance_perturbation_upper"]))
        max_drho = max(max_drho, drho)
        max_dk = max(max_dk, dk)
        if worst is None or d12b > worst["V12B_correction_norm_upper_rad"]:
            worst = row
        if not closed:
            unclosed += 1
            if bad is None:
                bad = row

    ok = bool(rows) and unclosed == 0 and bad is None and not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_PSD_S_PSD_PRESERVING_RESOLVENT_V12B",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V10_canonical_core_retained": True,
        "V11_PSD_and_S_perturbation_magnitudes_retained": True,
        "actual_shipping_covariance_PSD_used": True,
        "actual_innovation_noise_floor_Ra_used_without_subtraction": True,
        "exact_gain_resolvent_identity_used": True,
        "sample1_force_Jacobian_perturbation_included": True,
        "first_attitude_PSD_cross_axis_remainder_included": True,
        "second_prediction_attitude_process_remainder_included": True,
        "sample1_S_covariance_update_included": True,
        "sample1_S_attitude_injection_included": True,
        "sample1_S_aw_mean_correction_included": True,
        "sample1_S_solver_identity_branch_contained_as_zero_perturbation": True,
        "generic_interval_innovation_LDLT_floor_used": False,
        "perturbation_subtracted_from_measurement_noise_floor": False,
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
        "max_total_residual_perturbation_upper_mps2": max_drho,
        "max_total_reduced_covariance_perturbation_upper": max_dP,
        "max_sample1_H_perturbation_upper": max_dH,
        "max_sample1_cross_covariance_perturbation_upper": max_dC,
        "max_sample1_innovation_perturbation_upper": max_dS,
        "max_sample1_gain_operator_perturbation_upper": max_dk,
        "max_V12B_correction_norm_upper_rad": max_d,
        "first_unclosed_joint_cell": bad,
        "worst_joint_cell": worst,
        "P5_SAMPLE1_PSD_S_PSD_RESOLVENT_V12B": "PASS" if ok else "NOT_ESTABLISHED",
        "next_obligation": (
            "SIGNED_RADIAL_SUBDIVIDE_AND_CAYLEY_COMPOSE_SAMPLE1_INSIDE_Q8"
            if ok else "REFINE_PSD_RESOLVENT_NUMERATOR_BY_ONE_PLUS_TWO_BLOCK_AT_FIRST_WITNESS"
        ),
        "failures": failures,
        "rows": rows,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    for k in (
        "source_generated_not_trajectory_fit", "V10_canonical_core_retained",
        "V11_PSD_and_S_perturbation_magnitudes_retained",
        "actual_shipping_covariance_PSD_used",
        "actual_innovation_noise_floor_Ra_used_without_subtraction",
        "exact_gain_resolvent_identity_used", "sample1_force_Jacobian_perturbation_included",
        "first_attitude_PSD_cross_axis_remainder_included",
        "second_prediction_attitude_process_remainder_included",
        "sample1_S_covariance_update_included", "sample1_S_attitude_injection_included",
        "sample1_S_aw_mean_correction_included",
        "sample1_S_solver_identity_branch_contained_as_zero_perturbation",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "generic_interval_innovation_LDLT_floor_used",
        "perturbation_subtracted_from_measurement_noise_floor",
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
        f.append("no V12B cells")
    for k in (
        "max_total_residual_perturbation_upper_mps2",
        "max_total_reduced_covariance_perturbation_upper",
        "max_sample1_H_perturbation_upper",
        "max_sample1_gain_operator_perturbation_upper",
        "max_V12B_correction_norm_upper_rad",
    ):
        if not math.isfinite(float(d.get(k, math.nan))):
            f.append(f"nonfinite {k}")
    st = d.get("P5_SAMPLE1_PSD_S_PSD_RESOLVENT_V12B")
    w = d.get("first_unclosed_joint_cell")
    if st == "PASS" and (w is not None or int(d.get("unclosed_joint_cells", -1)) != 0):
        f.append("V12B PASS retains unclosed cell")
    if st == "NOT_ESTABLISHED" and w is None and not f:
        f.append("missing V12B witness")
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
        "status": d["P5_SAMPLE1_PSD_S_PSD_RESOLVENT_V12B"],
        "cells": d["evaluated_joint_cells"],
        "unclosed": d["unclosed_joint_cells"],
        "max_drho": d["max_total_residual_perturbation_upper_mps2"],
        "max_dP": d["max_total_reduced_covariance_perturbation_upper"],
        "max_dH": d["max_sample1_H_perturbation_upper"],
        "max_dS": d["max_sample1_innovation_perturbation_upper"],
        "max_dK": d["max_sample1_gain_operator_perturbation_upper"],
        "max_d": d["max_V12B_correction_norm_upper_rad"],
        "first_unclosed": d["first_unclosed_joint_cell"],
        "worst": d["worst_joint_cell"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
