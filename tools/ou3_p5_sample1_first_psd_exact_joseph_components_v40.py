#!/usr/bin/env python3
"""V40: exact attitude-supported Joseph transport for the first PSD remainder.

V38 makes the first-PSD mean correction source-directional and exact in the
canonical gravity gauge, but deliberately retains V36's covariance transport

    ||(I-K'H) E (I-K'H)^T|| <= ||I-K'H||^2 ||E||,

with a full 6x6 operator norm.  The omitted PSD remainder E is not a generic
6x6 perturbation: after the diagonal part is absorbed into V12D's t/Y
intervals it is supported only on the attitude block and has the form

    O = [[0,a,b],[a,0,c],[b,c,0]],  |a|,|b|,|c| <= eps/2.

For the nominal first accelerometer update, only the first three columns B0 of
A0=I-KH act on O.  In the canonical gravity gauge those columns are mutually
orthogonal.  The two tangent columns have common norm

    s = sqrt(((p+r)/D)^2 + (g p/D)^2),  D=g^2 t+p+r,

and the yaw column has norm one.  In the normalized column basis the nominal
transport is therefore the symmetric zero-diagonal 3x3 matrix with entry
bounds eps/2*s^2, eps/2*s, eps/2*s.  Its spectral norm is bounded by its maximum
absolute row sum,

    eps/2 * max(s^2+s, 2s).

The gain perturbation changes only B through -Delta K H_theta, hence

    ||Delta B|| <= g ||Delta K||.

V40 adds the two cross terms and quadratic term explicitly against the PSD
operator bound eps.  This replaces only V36/V38's first posterior covariance
transport; V38's exact residual/correction geometry and every downstream
reset/process/source guard remain unchanged.

This is a witness-level refinement precondition for a later V39 current-subbox
lift.  It does not compose q<8 or promote sample 1/P5.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_first_psd_exact_correction_geometry_v38 as V38

DEFAULT_DOMAIN = V38.DEFAULT_DOMAIN
SCHEMA = 4000
V12D = V38.V12D
FULL = V38.FULL
V30 = V38.V30


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


def _first_psd_perturbation_exact_joseph_components(*, t, Y, p, r, g: float,
                                                     eps: float, rho0: float,
                                                     dhi: float, rt, rz,
                                                     alpha_hi: float,
                                                     aw_pre: float) -> dict:
    base = V38._first_psd_perturbation_exact_correction(
        t=t, Y=Y, p=p, r=r, g=g, eps=eps, rho0=rho0, dhi=dhi,
        rt=rt, rz=rz, alpha_hi=alpha_hi, aw_pre=aw_pre)

    D = FULL.I(g * g) * t + p + r
    if not D.lo > 0.0:
        raise RuntimeError("positive first tangent innovation floor required")

    # Exact nominal attitude-supported columns of A0=I-KH.  The tangent
    # columns have disjoint support and equal norm; the yaw column is unit.
    beta = (p + r) / D
    gamma = FULL.I(abs(float(g))) * p / D
    beta_hi = beta.abs_upper()
    gamma_hi = gamma.abs_upper()
    s2 = FULL.up(FULL.up(beta_hi * beta_hi) + FULL.up(gamma_hi * gamma_hi))
    s = FULL.up(math.sqrt(s2))
    b0 = max(1.0, s)

    e = FULL.up(0.5 * float(eps))
    eop = FULL.up(float(eps))
    nominal_row12 = _mul_up(e, _sum_up(FULL.up(s * s), s))
    nominal_row3 = _mul_up(2.0, e, s)
    nominal_transport = max(nominal_row12, nominal_row3)

    dkth = float(base["first_gain_theta_perturbation_upper"])
    dkaw = float(base["first_gain_aw_perturbation_upper"])
    dk = FULL.up(math.sqrt(FULL.up(
        FULL.up(dkth * dkth) + FULL.up(dkaw * dkaw))))
    delta_B = _mul_up(abs(float(g)), dk)
    cross = _mul_up(2.0, delta_B, eop, b0)
    quadratic = _mul_up(delta_B, delta_B, eop)
    dPplus = _sum_up(nominal_transport, cross, quadratic)

    parent_dPplus = float(base["first_posterior_covariance_perturbation_upper"])
    if dPplus > FULL.up(parent_dPplus):
        raise RuntimeError("exact Joseph component transport exceeded V38 parent")

    dd = float(base["first_offaxis_attitude_correction_upper_rad"])
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
        "first_PSD_Joseph_attitude_support_used": True,
        "first_PSD_Joseph_nominal_columns_orthogonal": True,
        "first_PSD_Joseph_zero_diagonal_component_matrix_used": True,
        "first_PSD_Joseph_tangent_column_norm_upper": s,
        "first_PSD_Joseph_nominal_B_operator_upper": b0,
        "first_PSD_Joseph_nominal_component_transport_upper": nominal_transport,
        "first_PSD_Joseph_deltaK_operator_upper": dk,
        "first_PSD_Joseph_deltaB_operator_upper": delta_B,
        "first_PSD_Joseph_cross_transport_upper": cross,
        "first_PSD_Joseph_quadratic_transport_upper": quadratic,
        "first_PSD_Joseph_parent_covariance_transport_upper": parent_dPplus,
        "first_posterior_covariance_perturbation_upper": dPplus,
        "reset_gauge_transform_perturbation_upper": dirterm,
        "sample1_reduced_covariance_PSD_perturbation_upper": dP1,
        "V38_exact_mean_correction_geometry_retained": True,
        "V36_full_gain_operator_retired_only_for_PSD_covariance_transport": True,
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

    V12D._first_psd_perturbation_tangent = V38._first_psd_perturbation_exact_correction
    try:
        baseline = V12D.build(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    finally:
        V12D._first_psd_perturbation_tangent = original
    failures = [f"V38 baseline: {x}" for x in V12D.validate(baseline)]
    b = _witness(baseline)

    V12D._first_psd_perturbation_tangent = _first_psd_perturbation_exact_joseph_components
    try:
        refined = V12D.build(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    finally:
        V12D._first_psd_perturbation_tangent = original
    failures += [f"V40 refined: {x}" for x in V12D.validate(refined)]
    r = _witness(refined)

    bm = _metrics(b)
    rm = _metrics(r)
    for key in bm:
        if rm[key] > FULL.up(bm[key]):
            failures.append(f"refined {key} exceeded V38 parent")

    for k in (
        "first_PSD_Joseph_attitude_support_used",
        "first_PSD_Joseph_nominal_columns_orthogonal",
        "first_PSD_Joseph_zero_diagonal_component_matrix_used",
        "V38_exact_mean_correction_geometry_retained",
        "V36_full_gain_operator_retired_only_for_PSD_covariance_transport",
    ):
        if r.get(k) is not True:
            failures.append(f"refined row lost {k}")

    strict_keys = (
        "first_posterior_covariance_perturbation_upper",
        "sample1_reduced_covariance_PSD_perturbation_upper",
        "total_reduced_covariance_perturbation_upper",
        "sample1_innovation_perturbation_upper",
        "sample1_attitude_gain_operator_perturbation_upper",
    )
    strict = all(rm[k] < bm[k] for k in strict_keys)
    ratios = {k: (0.0 if bm[k] == 0.0 else rm[k] / bm[k]) for k in strict_keys}

    detail = {
        "tangent_column_norm_upper": float(r["first_PSD_Joseph_tangent_column_norm_upper"]),
        "nominal_component_transport_upper": float(r["first_PSD_Joseph_nominal_component_transport_upper"]),
        "deltaB_operator_upper": float(r["first_PSD_Joseph_deltaB_operator_upper"]),
        "cross_transport_upper": float(r["first_PSD_Joseph_cross_transport_upper"]),
        "quadratic_transport_upper": float(r["first_PSD_Joseph_quadratic_transport_upper"]),
        "parent_covariance_transport_upper": float(r["first_PSD_Joseph_parent_covariance_transport_upper"]),
        "refined_covariance_transport_upper": float(r["first_posterior_covariance_perturbation_upper"]),
    }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_FIRST_PSD_EXACT_JOSEPH_COMPONENTS_V40",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V38_exact_correction_parent_revalidated": True,
        "attitude_supported_Joseph_columns_used": True,
        "zero_diagonal_PSD_component_matrix_used": True,
        "deltaK_cross_terms_retained": True,
        "baseline_witness": bm,
        "refined_witness": rm,
        "refinement_ratios": ratios,
        "joseph_component_detail": detail,
        "strict_witness_refinement": strict,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_FIRST_PSD_EXACT_JOSEPH_COMPONENTS_V40": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_V40_EXACT_JOSEPH_COMPONENTS_THROUGH_V39_CURRENT_SUBBOX_COVER"
            if strict and not failures else
            "REFINE_FIRST_PSD_RESET_TRANSPORT_COMPONENT_MATRIX_AT_V39_FIRST_OPEN_SUBBOX"),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_FIRST_PSD_EXACT_JOSEPH_COMPONENTS_V40":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V38_exact_correction_parent_revalidated",
        "attitude_supported_Joseph_columns_used",
        "zero_diagonal_PSD_component_matrix_used",
        "deltaK_cross_terms_retained",
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
    if d.get("P5_SAMPLE1_FIRST_PSD_EXACT_JOSEPH_COMPONENTS_V40") not in (
            "PASS", "NOT_ESTABLISHED"):
        f.append("invalid V40 status")
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
        "status": d["P5_SAMPLE1_FIRST_PSD_EXACT_JOSEPH_COMPONENTS_V40"],
        "detail": d.get("joseph_component_detail"),
        "ratios": d.get("refinement_ratios"),
        "strict": d.get("strict_witness_refinement"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
