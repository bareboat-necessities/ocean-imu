#!/usr/bin/env python3
"""Tangent-channel refinement of the OU-III P5 sample-1 V12C closure.

The authoritative V12C 24^3 run identified the remaining dominant proof loss:
V11 bounded the first accelerometer covariance-remainder gain perturbation using
``p_aw + R_a`` as a scalar innovation floor.  That is the *axial* innovation
channel, but the omitted first attitude covariance remainder enters the
accelerometer innovation only through the gravity-tangent two-channel block.

At the first normal-Live accelerometer packet, in the exact gravity gauge,

    H_theta = -[g e3]_x,
    P_theta = diag(t,t,Y) + E,
    P_aw = p I,
    P_theta,aw = 0.

The nominal innovation is therefore

    S0 = diag(D,D,p+R_a),  D = g^2 t + p + R_a.

For the omitted symmetric attitude covariance remainder ``E`` with
``||E||<=e``, the innovation perturbation is

    Delta S = H_theta E H_theta^T,

which has an exact zero third row/column and
``||Delta S|| <= g^2 e``.  Consequently only the tangent inverse is relevant,
with

    ||S_tan'^-1|| <= 1/(D_lower-g^2 e).

The first attitude-gain perturbation follows from the exact resolvent

    Delta K_theta = (Delta C_theta-K_theta Delta S) S_tan'^-1,
    ||Delta C_theta|| <= g e.

Likewise the a_w-gain perturbation is tangent only; its nominal tangent gain is
``p/D`` and its numerator has no independent covariance term.  This removes the
incorrect amplification by the axial ``(p+R_a)^-1`` while retaining every V11
PSD/S remainder, every V12C sample-1 perturbation term, and the unconditional
shipping source family.

The producer reuses V12C after temporarily replacing only V11's first-PSD
remainder helper with this algebraically sharper source bound.  No filter,
source domain, deployed correction limit, Cayley target, or theorem promotion
is changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_structured_full_gain_v12c as V12C

DEFAULT_DOMAIN = V12C.DEFAULT_DOMAIN
SCHEMA = 1203
RANGE = V12C.RANGE
FULL = V12C.FULL
V11 = V12C.V11


def _first_psd_perturbation_tangent(*, t, Y, p, r, g: float,
                                    eps: float, rho0: float, dhi: float,
                                    rt, rz, alpha_hi: float,
                                    aw_pre: float) -> dict:
    """Refine V11 first-PSD perturbation using the exact tangent innovation."""
    eoff = FULL.up(2.0 * eps)
    D = FULL.I(g * g) * t + p + r
    dS = FULL.up(g * g * eoff)
    tangent_floor = FULL.down(D.lo - dS)
    if not tangent_floor > 0.0:
        raise RuntimeError("first tangent innovation floor lost")
    inv_tangent = FULL.up(1.0 / tangent_floor)

    # Nominal first attitude and a_w tangent gains.  The omitted covariance
    # remainder is attitude-only, hence the axial innovation/gain is unchanged.
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


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24) -> dict:
    original = V11._first_psd_perturbation
    V11._first_psd_perturbation = _first_psd_perturbation_tangent
    try:
        core = V12C.build(
            domain_path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    finally:
        V11._first_psd_perturbation = original

    failures = list(V12C.validate(core))
    rows = core.get("rows", [])
    if rows:
        for row in rows:
            if row.get("first_PSD_innovation_perturbation_tangent_only") is not True:
                failures.append("row lost tangent-only first PSD innovation structure")
                break
            if row.get("first_PSD_innovation_axial_row_column_exact_zero") is not True:
                failures.append("row lost exact axial zero in first PSD innovation perturbation")
                break

    status = core.get("P5_SAMPLE1_ONE_PLUS_TWO_ATTITUDE_RESOLVENT_V12C")
    out = dict(core)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_FIRST_PSD_TANGENT_CHANNEL_REFINEMENT_V12D",
        "V12C_one_plus_two_resolvent_retained": True,
        "V11_first_PSD_generic_axial_noise_floor_retired": True,
        "first_PSD_tangent_innovation_exact_structure_used": True,
        "first_PSD_axial_innovation_perturbation_exact_zero": True,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "P5_SAMPLE1_FIRST_PSD_TANGENT_REFINEMENT_V12D": (
            "PASS" if status == "PASS" and not failures else "NOT_ESTABLISHED"
        ),
        "next_obligation": (
            "SIGNED_RADIAL_SUBDIVIDE_AND_CAYLEY_COMPOSE_SAMPLE1_INSIDE_Q8"
            if status == "PASS" and not failures
            else "REFINE_REMAINING_SAMPLE1_BLOCK_PERTURBATION_OR_FIRST_RESET_DIRECTION"
        ),
        "failures": failures,
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V12C_one_plus_two_resolvent_retained",
        "V11_first_PSD_generic_axial_noise_floor_retired",
        "first_PSD_tangent_innovation_exact_structure_used",
        "first_PSD_axial_innovation_perturbation_exact_zero",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed",
        "deployed_correction_limit_increased", "complete_sample1_branch_closed_here",
        "signed_cayley_q8_composed_here", "q8_word_promoted_here",
        "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    st = d.get("P5_SAMPLE1_FIRST_PSD_TANGENT_REFINEMENT_V12D")
    if st not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid V12D status")
    if st == "PASS" and int(d.get("unclosed_joint_cells", -1)) != 0:
        f.append("V12D PASS retains unclosed cells")
    if st == "NOT_ESTABLISHED" and d.get("first_unclosed_joint_cell") is None and not f:
        f.append("missing V12D witness")
    # Deduplicate failures inherited through the wrapped V12C validator.
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
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_FIRST_PSD_TANGENT_REFINEMENT_V12D"],
        "cells": d["evaluated_joint_cells"],
        "unclosed": d["unclosed_joint_cells"],
        "fallbacks": d["noise_floor_inverse_fallback_cells"],
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
