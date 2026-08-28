#!/usr/bin/env python3
"""V36: sharpen the first-PSD off-diagonal remainder by matrix order.

The canonical first-attitude source covariance is

    P_theta = diag(tilt, tilt, yaw) + E,    0 <= E <= eps I.

V12D already absorbs the diagonal entries of E into
`t=[tilt,tilt+eps]` and `Y=[yaw,yaw+eps]`.  Only the symmetric zero-diagonal
remainder O remains.  For every i!=j, the 2x2 principal minors of E and
`eps I-E` imply

    |E_ij| <= min(sqrt(E_ii E_jj),
                  sqrt((eps-E_ii)(eps-E_jj))) <= eps/2.

Each row of O therefore has absolute sum at most eps, so
`||O||_2 <= sqrt(||O||_1 ||O||_inf) <= eps`.

V12D currently uses `2 eps` for this operator bound.  V36 changes only that
source lemma to `eps`, then re-runs the unchanged V12D tangent-resolvent
construction at the authoritative first q8 witness.  This is a diagnostic
precondition for a later q8 lift; it does not compose q8 or promote sample 1.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_v12d_witness_perturbation_decomposition_v35 as V35

DEFAULT_DOMAIN = V35.DEFAULT_DOMAIN
SCHEMA = 3600
V12D = V35.V12D
FULL = V12D.FULL
V11 = V12D.V11
V30 = V35.V30


def _psd_offdiag_operator_upper(eps: float) -> float:
    if not (math.isfinite(float(eps)) and float(eps) >= 0.0):
        raise ValueError("finite nonnegative PSD remainder required")
    return FULL.up(float(eps))


def _first_psd_perturbation_psd_cone(*, t, Y, p, r, g: float,
                                     eps: float, rho0: float, dhi: float,
                                     rt, rz, alpha_hi: float,
                                     aw_pre: float) -> dict:
    eoff = _psd_offdiag_operator_upper(eps)
    D = FULL.I(g * g) * t + p + r
    dS = FULL.up(g * g * eoff)
    tangent_floor = FULL.down(D.lo - dS)
    if not tangent_floor > 0.0:
        raise RuntimeError("first tangent innovation floor lost")
    inv_tangent = FULL.up(1.0 / tangent_floor)

    ktheta_tangent = (FULL.I(g) * t / D).abs_upper()
    kaw_tangent = (p / D).abs_upper()
    dCtheta = FULL.up(g * eoff)
    dkth = FULL.up(
        FULL.up(dCtheta + FULL.up(ktheta_tangent * dS)) * inv_tangent)
    dkaw = FULL.up(FULL.up(kaw_tangent * dS) * inv_tangent)
    dkfull = FULL.up(math.sqrt(FULL.up(dkth * dkth + FULL.up(dkaw * dkaw))))
    dd = FULL.up(dkth * rho0)
    daw = FULL.up(dkaw * rho0)

    D0 = FULL.I(g * g) * t + p + r
    kawt = p / D0
    kz = p / (p + r)
    awt = kawt * rt
    az = kz * rz
    dxaw = V11.V5._norm2_upper(awt.abs_upper(), az.abs_upper())
    vec = FULL.up(g + FULL.up(alpha_hi * FULL.up(aw_pre + FULL.up(dxaw + daw))))
    drho = FULL.up(FULL.up(alpha_hi * daw) + FULL.up(dd * vec))

    h0 = FULL.up(math.sqrt(FULL.up(g * g + 1.0)))
    A0 = V11._first_A_norm(t, p, r, g)
    Amax = FULL.up(A0 + FULL.up(dkfull * h0))
    dPplus = FULL.up(FULL.up(Amax * Amax) * eoff)

    Tnom = FULL.up(math.sqrt(FULL.up(1.0 + FULL.up(0.25 * dhi * dhi))))
    Gactual = FULL.up(
        math.sqrt(FULL.up(1.0 + FULL.up(0.25 * FULL.up((dhi + dd) * (dhi + dd))))))
    dT = FULL.up(FULL.up(dd * Gactual) + FULL.up(0.5 * dd))
    Pprior = max(Y.hi, t.hi, p.hi)
    dirterm = FULL.up(
        FULL.up(FULL.up(2.0 * Tnom * dT) + FULL.up(dT * dT)) * Pprior)
    after_reset = FULL.up(FULL.up(Tnom * Tnom * dPplus) + dirterm)
    dP1 = FULL.up(after_reset + eps)
    return {
        "first_attitude_offdiagonal_operator_upper": eoff,
        "first_PSD_offdiagonal_operator_from_matrix_order": True,
        "first_PSD_diagonal_absorbed_in_t_Y_intervals": True,
        "first_PSD_offdiagonal_entry_abs_upper_eps_over_2": FULL.up(0.5 * eps),
        "first_PSD_innovation_perturbation_tangent_only": True,
        "first_PSD_innovation_axial_row_column_exact_zero": True,
        "first_nominal_tangent_innovation_lower": D.lo,
        "first_tangent_innovation_perturbation_upper": dS,
        "first_perturbed_tangent_innovation_lower": tangent_floor,
        "first_perturbed_tangent_inverse_operator_upper": inv_tangent,
        "first_nominal_Ktheta_tangent_operator_upper": ktheta_tangent,
        "first_nominal_Kaw_tangent_operator_upper": kaw_tangent,
        "first_gain_theta_perturbation_upper": dkth,
        "first_gain_aw_perturbation_upper": dkaw,
        "first_offaxis_attitude_correction_upper_rad": dd,
        "first_aw_x_correction_upper_mps2": daw,
        "first_nominal_aw_correction_norm_upper_mps2": dxaw,
        "PSD_induced_sample1_residual_perturbation_upper_mps2": drho,
        "first_covariance_update_A_norm_upper": A0,
        "first_covariance_update_A_perturbed_norm_upper": Amax,
        "first_posterior_covariance_perturbation_upper": dPplus,
        "reset_gauge_transform_perturbation_upper": dirterm,
        "sample1_reduced_covariance_PSD_perturbation_upper": dP1,
    }


def _witness(core: dict) -> dict:
    return V30._witness_row(core)


def _metrics(row: dict) -> dict:
    keys = (
        "first_attitude_offdiagonal_operator_upper",
        "first_gain_theta_perturbation_upper",
        "first_gain_aw_perturbation_upper",
        "first_offaxis_attitude_correction_upper_rad",
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
    baseline = V12D.build(
        path, source_pieces=source_pieces, source_cell_index=source_cell_index,
        p_pieces=p_pieces, tangent_pieces=tangent_pieces,
        axial_pieces=axial_pieces)
    failures = [f"baseline V12D: {x}" for x in V12D.validate(baseline)]
    b = _witness(baseline)

    original = V12D._first_psd_perturbation_tangent
    V12D._first_psd_perturbation_tangent = _first_psd_perturbation_psd_cone
    try:
        refined = V12D.build(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    finally:
        V12D._first_psd_perturbation_tangent = original
    failures += [f"refined V12D: {x}" for x in V12D.validate(refined)]
    r = _witness(refined)

    bm = _metrics(b); rm = _metrics(r)
    for key in bm:
        if rm[key] > FULL.up(bm[key]):
            failures.append(f"refined {key} exceeded V12D parent")
    if r.get("first_PSD_offdiagonal_operator_from_matrix_order") is not True:
        failures.append("refined row lost PSD matrix-order flag")
    if r.get("first_PSD_diagonal_absorbed_in_t_Y_intervals") is not True:
        failures.append("refined row lost diagonal-absorption flag")

    strict_keys = [
        "first_attitude_offdiagonal_operator_upper",
        "first_offaxis_attitude_correction_upper_rad",
        "PSD_induced_sample1_residual_perturbation_upper_mps2",
        "sample1_reduced_covariance_PSD_perturbation_upper",
        "sample1_attitude_gain_operator_perturbation_upper",
    ]
    strict = all(rm[k] < bm[k] for k in strict_keys)
    ratios = {
        k: (0.0 if bm[k] == 0.0 else rm[k] / bm[k])
        for k in strict_keys
    }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_FIRST_PSD_OFFDIAGONAL_CONE_V36",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V12D_baseline_revalidated": True,
        "PSD_remainder_order_interval_0_E_epsI_used": True,
        "PSD_diagonal_absorbed_in_existing_t_Y_intervals": True,
        "PSD_principal_minor_offdiag_abs_le_eps_over_2": True,
        "PSD_zero_diagonal_remainder_operator_le_eps": True,
        "baseline_witness": bm,
        "refined_witness": rm,
        "refinement_ratios": ratios,
        "strict_witness_refinement": strict,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_FIRST_PSD_OFFDIAGONAL_CONE_V36": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_V36_PSD_OFFDIAGONAL_CONE_THROUGH_V34_CURRENT_SUBBOX_COVER"
            if strict and not failures else
            "DERIVE_EXACT_FIRST_PSD_RESET_COMPONENT_MATRIX_AT_Q8_WITNESS"
        ),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_FIRST_PSD_OFFDIAGONAL_CONE_V36":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit", "V12D_baseline_revalidated",
        "PSD_remainder_order_interval_0_E_epsI_used",
        "PSD_diagonal_absorbed_in_existing_t_Y_intervals",
        "PSD_principal_minor_offdiag_abs_le_eps_over_2",
        "PSD_zero_diagonal_remainder_operator_le_eps",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed", "q8_composed_here",
              "q8_word_promoted_here", "whole_word_promoted_here",
              "N_H_words_set_here"):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if d.get("P5_SAMPLE1_FIRST_PSD_OFFDIAGONAL_CONE_V36") not in (
            "PASS", "NOT_ESTABLISHED"):
        f.append("invalid V36 status")
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
    d = build(x.domain, source_pieces=x.source_pieces,
              source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
              tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_FIRST_PSD_OFFDIAGONAL_CONE_V36"],
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
