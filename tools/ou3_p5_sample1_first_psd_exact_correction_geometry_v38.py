#!/usr/bin/env python3
"""V38: exact canonical first-correction geometry for the V36 PSD cone.

V36 proves that after diagonal absorption the omitted first-attitude PSD
remainder is a symmetric zero-diagonal matrix

    O = [[0,a,b],[a,0,c],[b,c,0]],  |a|,|b|,|c| <= eps/2.

V36 still converts O to one operator ball before bounding the first Kalman
correction.  In the canonical gravity gauge the first accelerometer residual is
[0,-r_t,r_z], and the O-induced innovation perturbation has exact tangent form

    Delta S_t = [[0,-g^2 a],[-g^2 a,0]],

with zero axial row/column.  Let D_x,D_y be the two nominal tangent innovation
diagonals, each lying in D=g^2 t+p+r, u=g^2 a, and
Delta=D_x D_y-u^2.  The exact second gain column gives

  |delta d_x| <= r_t g^3 a^2 (p+r)/(D_x Delta),
  |delta d_y| <= r_t g |a| (p+r)/Delta,
  |delta d_z| <= r_t g (|b| D_y + |c u|)/Delta.

The a_w tangent correction difference is similarly

  |delta a_x| <= r_t p |u|/Delta,
  |delta a_y| <= r_t p u^2/(D_x Delta).

V38 uses these source-directional bounds only for the state/residual and reset
direction perturbations.  V36's full gain-operator bounds remain in force for
the Joseph covariance-update perturbation, so no covariance direction is
silently discarded.  The resulting helper is audited at the authoritative
V12D witness before any q8 lift.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_first_psd_offdiagonal_cone_v36 as V36

DEFAULT_DOMAIN = V36.DEFAULT_DOMAIN
SCHEMA = 3800
V12D = V36.V12D
FULL = V36.FULL
V11 = V36.V11
V30 = V36.V30


def _mul_up(*xs: float) -> float:
    y = 1.0
    for x in xs:
        y = FULL.up(y * float(x))
    return y


def _sum_up(*xs: float) -> float:
    y = 0.0
    for x in xs:
        y = FULL.up(y + float(x))
    return y


def _first_psd_perturbation_exact_correction(*, t, Y, p, r, g: float,
                                              eps: float, rho0: float,
                                              dhi: float, rt, rz,
                                              alpha_hi: float,
                                              aw_pre: float) -> dict:
    base = V36._first_psd_perturbation_psd_cone(
        t=t, Y=Y, p=p, r=r, g=g, eps=eps, rho0=rho0, dhi=dhi,
        rt=rt, rz=rz, alpha_hi=alpha_hi, aw_pre=aw_pre)

    e = FULL.up(0.5 * float(eps))
    D = FULL.I(g * g) * t + p + r
    if not D.lo > 0.0:
        raise RuntimeError("positive first tangent innovation floor required")
    u = _mul_up(g * g, e)
    u2 = FULL.up(u * u)
    dprod_lo = FULL.down(D.lo * D.lo)
    det_lo = FULL.down(dprod_lo - u2)
    if not det_lo > 0.0:
        raise RuntimeError("PSD-cone tangent determinant floor lost")

    rt_hi = rt.abs_upper()
    pr_hi = (p + r).hi
    p_hi = max(0.0, p.hi)
    gabs = abs(float(g))

    # Difference from the nominal x-directed attitude correction.  The x term
    # is second order in a; y is linear in a; z is linear in b plus the a*c
    # coupling.  D_x,D_y are independently enclosed by the same interval D.
    dx_num = _mul_up(rt_hi, gabs ** 3, FULL.up(e * e), pr_hi)
    dx_den = FULL.down(D.lo * det_lo)
    if not dx_den > 0.0:
        raise RuntimeError("PSD-cone x-correction denominator floor lost")
    dx = FULL.up(dx_num / dx_den)
    dy = FULL.up(_mul_up(rt_hi, gabs, e, pr_hi) / det_lo)
    dz = FULL.up(
        _mul_up(rt_hi, gabs, e, _sum_up(D.hi, u)) / det_lo)
    dd = FULL.up(math.sqrt(FULL.up(
        FULL.up(dx * dx) + FULL.up(FULL.up(dy * dy) + FULL.up(dz * dz)))))

    # Exact tangent a_w correction difference induced only through S_t^{-1}.
    dax = FULL.up(_mul_up(rt_hi, p_hi, u) / det_lo)
    day_num = _mul_up(rt_hi, p_hi, u2)
    day = FULL.up(day_num / dx_den)
    daw = FULL.up(math.sqrt(FULL.up(FULL.up(dax * dax) + FULL.up(day * day))))

    # Keep V36's full gain-operator/Joseph covariance perturbation.  Only the
    # source-directional mean correction and reset-direction terms are sharper.
    dxaw = float(base["first_nominal_aw_correction_norm_upper_mps2"])
    vec = FULL.up(gabs + FULL.up(alpha_hi * FULL.up(aw_pre + FULL.up(dxaw + daw))))
    drho = FULL.up(FULL.up(alpha_hi * daw) + FULL.up(dd * vec))

    dPplus = float(base["first_posterior_covariance_perturbation_upper"])
    Tnom = FULL.up(math.sqrt(FULL.up(1.0 + FULL.up(0.25 * dhi * dhi))))
    Gactual = FULL.up(math.sqrt(FULL.up(
        1.0 + FULL.up(0.25 * FULL.up((dhi + dd) * (dhi + dd))))))
    dT = FULL.up(FULL.up(dd * Gactual) + FULL.up(0.5 * dd))
    Pprior = max(Y.hi, t.hi, p.hi)
    dirterm = FULL.up(
        FULL.up(FULL.up(2.0 * Tnom * dT) + FULL.up(dT * dT)) * Pprior)
    after_reset = FULL.up(FULL.up(Tnom * Tnom * dPplus) + dirterm)
    dP1 = FULL.up(after_reset + eps)

    out = dict(base)
    out.update({
        "first_PSD_exact_canonical_residual_direction_used": True,
        "first_PSD_exact_tangent_2x2_inverse_geometry_used": True,
        "first_PSD_axial_residual_has_zero_attitude_effect": True,
        "first_PSD_offdiagonal_entry_abs_upper": e,
        "first_PSD_tangent_offdiagonal_innovation_abs_upper": u,
        "first_PSD_tangent_determinant_lower": det_lo,
        "first_PSD_exact_delta_attitude_x_abs_upper_rad": dx,
        "first_PSD_exact_delta_attitude_y_abs_upper_rad": dy,
        "first_PSD_exact_delta_attitude_z_abs_upper_rad": dz,
        "first_PSD_exact_delta_aw_x_abs_upper_mps2": dax,
        "first_PSD_exact_delta_aw_y_abs_upper_mps2": day,
        "first_offaxis_attitude_correction_upper_rad": dd,
        "first_aw_x_correction_upper_mps2": daw,
        "PSD_induced_sample1_residual_perturbation_upper_mps2": drho,
        "reset_gauge_transform_perturbation_upper": dirterm,
        "sample1_reduced_covariance_PSD_perturbation_upper": dP1,
        "V36_full_gain_operator_retained_for_Joseph_covariance": True,
    })
    return out


def _witness(core: dict) -> dict:
    return V30._witness_row(core)


def _metrics(row: dict) -> dict:
    keys = (
        "first_offaxis_attitude_correction_upper_rad",
        "first_aw_x_correction_upper_mps2",
        "PSD_induced_sample1_residual_perturbation_upper_mps2",
        "first_posterior_covariance_perturbation_upper",
        "reset_gauge_transform_perturbation_upper",
        "sample1_reduced_covariance_PSD_perturbation_upper",
        "total_residual_perturbation_upper_mps2",
        "total_reduced_covariance_perturbation_upper",
        "sample1_innovation_perturbation_upper",
        "sample1_attitude_cross_covariance_perturbation_upper",
        "sample1_attitude_gain_operator_perturbation_upper",
        "V12C_correction_norm_upper_rad",
    )
    return {k: float(row[k]) for k in keys}


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24) -> dict:
    path = Path(domain_path).resolve()
    original = V12D._first_psd_perturbation_tangent

    V12D._first_psd_perturbation_tangent = V36._first_psd_perturbation_psd_cone
    try:
        baseline = V12D.build(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    finally:
        V12D._first_psd_perturbation_tangent = original
    failures = [f"V36 baseline: {x}" for x in V12D.validate(baseline)]
    b = _witness(baseline)

    V12D._first_psd_perturbation_tangent = _first_psd_perturbation_exact_correction
    try:
        refined = V12D.build(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    finally:
        V12D._first_psd_perturbation_tangent = original
    failures += [f"V38 refined: {x}" for x in V12D.validate(refined)]
    r = _witness(refined)

    bm = _metrics(b)
    rm = _metrics(r)
    for key in bm:
        if rm[key] > FULL.up(bm[key]):
            failures.append(f"refined {key} exceeded V36 parent")
    for k in (
        "first_PSD_exact_canonical_residual_direction_used",
        "first_PSD_exact_tangent_2x2_inverse_geometry_used",
        "first_PSD_axial_residual_has_zero_attitude_effect",
        "V36_full_gain_operator_retained_for_Joseph_covariance",
    ):
        if r.get(k) is not True:
            failures.append(f"refined row lost {k}")

    strict_keys = (
        "first_offaxis_attitude_correction_upper_rad",
        "PSD_induced_sample1_residual_perturbation_upper_mps2",
        "reset_gauge_transform_perturbation_upper",
        "sample1_reduced_covariance_PSD_perturbation_upper",
        "sample1_attitude_gain_operator_perturbation_upper",
    )
    strict = all(rm[k] < bm[k] for k in strict_keys)
    ratios = {k: (0.0 if bm[k] == 0.0 else rm[k] / bm[k]) for k in strict_keys}

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_FIRST_PSD_EXACT_CORRECTION_GEOMETRY_V38",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V36_PSD_cone_parent_revalidated": True,
        "canonical_first_residual_direction_used": True,
        "exact_2x2_tangent_inverse_geometry_used": True,
        "axial_residual_zero_PSD_attitude_effect_used": True,
        "V36_Joseph_gain_operator_parent_retained": True,
        "baseline_witness": bm,
        "refined_witness": rm,
        "refinement_ratios": ratios,
        "strict_witness_refinement": strict,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_FIRST_PSD_EXACT_CORRECTION_GEOMETRY_V38": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_V38_EXACT_FIRST_PSD_CORRECTION_THROUGH_V37_CURRENT_SUBBOX_COVER"
            if strict and not failures else
            "DERIVE_EXACT_FIRST_PSD_RESET_COMPONENT_MATRIX_AT_Q8_WITNESS"),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_FIRST_PSD_EXACT_CORRECTION_GEOMETRY_V38":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V36_PSD_cone_parent_revalidated",
        "canonical_first_residual_direction_used",
        "exact_2x2_tangent_inverse_geometry_used",
        "axial_residual_zero_PSD_attitude_effect_used",
        "V36_Joseph_gain_operator_parent_retained",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "q8_composed_here",
        "q8_word_promoted_here", "whole_word_promoted_here",
        "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if d.get("P5_SAMPLE1_FIRST_PSD_EXACT_CORRECTION_GEOMETRY_V38") not in (
            "PASS", "NOT_ESTABLISHED"):
        f.append("invalid V38 status")
    return list(dict.fromkeys(f))


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
    d = build(
        x.domain, source_pieces=x.source_pieces,
        source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
        tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_FIRST_PSD_EXACT_CORRECTION_GEOMETRY_V38"],
        "baseline": d.get("baseline_witness"),
        "refined": d.get("refined_witness"),
        "ratios": d.get("refinement_ratios"),
        "strict": d.get("strict_witness_refinement"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
