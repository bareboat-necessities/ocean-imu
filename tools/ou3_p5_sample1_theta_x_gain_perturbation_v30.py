#!/usr/bin/env python3
"""V30: row-specific theta-x gain-perturbation resolvent at the first q8 witness.

V29 proves that almost the entire V12D correction remainder comes from the
source-uniform gain-perturbation ball ``Delta K_theta (r+Delta r)``.  The V12D
resolvent used the full nominal attitude-gain norm in

    ||Delta K_theta|| <= (dC_theta + ||K_theta|| dS) ||S'^{-1||.

For a single attitude row the same derivation is valid rowwise.  Since
``||e_i^T Delta C_theta||_2 <= ||Delta C_theta||_2`` and
``||e_i^T K_theta||_2`` is bounded by the corresponding exact one-plus-two
nominal gain block, the theta-x row satisfies

    ||Delta K_x||_2 <= (dC_theta + k_parallel dS) ||S'^{-1||,

whereas the y/z rows retain the existing full V12D bound.  V30 injects only
this sharper x-row gain-perturbation ball into V29's directional remainder.
Everything else, including V28 signed residuals, V23 current box, V12D y/z and
radial parents, V16/V15/V18 q8 checks, source domain, estimator, 6-rad shipping
limit, and promotion state, is unchanged.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_directional_v12d_remainder_v29 as V29

DEFAULT_DOMAIN = V29.DEFAULT_DOMAIN
SCHEMA = 3000
FULL = V29.FULL
Q_TARGET = V29.Q_TARGET
WITNESS = V29.WITNESS


def _witness_row(core: dict) -> dict:
    for row in core.get("rows", []):
        ids = (int(row["p_cell"]), int(row["tangent_residual_cell"]),
               int(row["axial_residual_cell"]))
        if ids == WITNESS:
            return row
    raise RuntimeError("V30 first-q8 witness row not found")


def _theta_x_gain_perturbation_upper(vr: dict, base: dict) -> dict:
    dC = float(vr["sample1_attitude_cross_covariance_perturbation_upper"])
    dS = float(vr["sample1_innovation_perturbation_upper"])
    inv = float(vr["actual_innovation_inverse_operator_upper"])
    kpar = float(base["Ktheta_parallel_block_upper"])
    dk_parent = float(vr["sample1_attitude_gain_operator_perturbation_upper"])
    vals = (dC, dS, inv, kpar, dk_parent)
    if not all(math.isfinite(x) and x >= 0.0 for x in vals):
        raise RuntimeError("invalid V30 row-resolvent inputs")
    numerator = FULL.up(dC + FULL.up(kpar * dS))
    dk_x = FULL.up(numerator * inv)
    if dk_x > FULL.up(dk_parent):
        raise RuntimeError("theta-x row gain perturbation exceeded V12D parent")
    return {
        "sample1_attitude_cross_covariance_perturbation_upper": dC,
        "sample1_innovation_perturbation_upper": dS,
        "actual_innovation_inverse_operator_upper": inv,
        "nominal_theta_x_gain_row_norm_upper": kpar,
        "theta_x_gain_perturbation_operator_upper": dk_x,
        "V12D_full_attitude_gain_perturbation_operator_upper": dk_parent,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6,
          current_component_pieces: int = 4) -> dict:
    path = Path(domain_path).resolve()
    V12D = V29.V28.V27.V23.V22.V21B.V21.V12D
    V10 = V12D.V11.V10
    v12 = V12D.build(path, source_pieces=source_pieces,
                     source_cell_index=source_cell_index, p_pieces=p_pieces,
                     tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    core = V10.build(path, source_pieces=source_pieces,
                     source_cell_index=source_cell_index, p_pieces=p_pieces,
                     tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    failures = [f"V12D: {x}" for x in V12D.validate(v12)]
    failures += [f"V10: {x}" for x in V10.validate(core)]
    try:
        vr = _witness_row(v12); base = _witness_row(core)
        row_detail = _theta_x_gain_perturbation_upper(vr, base)
        dkx = float(row_detail["theta_x_gain_perturbation_operator_upper"])
        rho_plus = FULL.up(float(base["sample1_full_residual_norm_upper_mps2"])
                           + float(vr["total_residual_perturbation_upper_mps2"]))
        gain_ball_x = FULL.up(dkx * rho_plus)
    except Exception as exc:
        failures.append(f"theta-x gain perturbation: {exc}")
        row_detail = None; gain_ball_x = math.inf

    original_caps = V29._directional_perturbation_caps
    def refined_caps(*, k_perp: float, k_parallel: float,
                     drho: float, dk: float, rho: float) -> dict:
        parent = original_caps(k_perp=k_perp, k_parallel=k_parallel,
                               drho=drho, dk=dk, rho=rho)
        if not math.isfinite(gain_ball_x):
            return parent
        ex = FULL.up(FULL.up(k_parallel * drho) + gain_ball_x)
        if ex > FULL.up(float(parent["x_correction_perturbation_abs_upper_rad"])):
            raise RuntimeError("V30 theta-x perturbation exceeded V29 parent")
        out = dict(parent)
        out.update({
            "theta_x_gain_perturbation_ball_upper_rad": gain_ball_x,
            "V29_full_gain_perturbation_ball_upper_rad": parent["gain_perturbation_ball_upper_rad"],
            "x_correction_perturbation_abs_upper_rad": ex,
        })
        return out

    V29._directional_perturbation_caps = refined_caps
    try:
        parent = V29.build(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces,
            current_component_pieces=current_component_pieces)
    finally:
        V29._directional_perturbation_caps = original_caps

    failures += [f"V29: {x}" for x in V29.validate(parent)]
    if parent.get("P5_SAMPLE1_DIRECTIONAL_V12D_REMAINDER_V29") != "PASS":
        failures.append("V29 directional remainder prerequisite did not pass")
    caps = parent.get("directional_perturbation_detail", {})
    used_x = float(caps.get("x_correction_perturbation_abs_upper_rad", math.inf))
    old_x = float(caps.get("V29_full_gain_perturbation_ball_upper_rad",
                          parent.get("previous_isotropic_V12D_correction_perturbation_upper_rad", math.inf)))
    x_strict = math.isfinite(used_x) and used_x < float(
        parent.get("previous_isotropic_V12D_correction_perturbation_upper_rad", math.inf))

    closed = bool(parent.get("first_open_subbox_closed_inside_q8") and not failures)
    out = dict(parent)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_THETA_X_GAIN_PERTURBATION_V30",
        "V29_directional_V12D_parent_retained": True,
        "V12D_resolvent_reused_rowwise": True,
        "theta_x_DeltaC_row_bounded_by_full_DeltaC_operator_norm": True,
        "theta_x_nominal_K_row_uses_parallel_block_bound": True,
        "theta_yz_gain_perturbation_parent_unchanged": True,
        "theta_x_gain_perturbation_detail": row_detail,
        "theta_x_correction_perturbation_strictly_refined": x_strict,
        "focused_first_open_subbox_closed_by_theta_x_gain_refinement": closed,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_THETA_X_GAIN_PERTURBATION_V30": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_V30_THETA_X_GAIN_REFINEMENT_OVER_ALL_V23_CURRENT_SUBBOXES"
            if closed else
            "DERIVE_DIRECTIONAL_DELTA_C_THETA_OR_DELTA_S_COUPLING_AT_FIRST_OPEN_SUBBOX"),
        "failures": list(dict.fromkeys(failures)),
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA: f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_THETA_X_GAIN_PERTURBATION_V30":
        f.append("qualification mismatch")
    for k in ("source_generated_not_trajectory_fit", "V29_directional_V12D_parent_retained",
              "V12D_resolvent_reused_rowwise",
              "theta_x_DeltaC_row_bounded_by_full_DeltaC_operator_norm",
              "theta_x_nominal_K_row_uses_parallel_block_bound",
              "theta_yz_gain_perturbation_parent_unchanged"):
        if d.get(k) is not True: f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed", "deployed_correction_limit_increased",
              "q8_composed_here", "q8_word_promoted_here", "whole_word_promoted_here",
              "N_H_words_set_here"):
        if d.get(k) is not False: f.append(f"{k} is not false")
    rd = d.get("theta_x_gain_perturbation_detail") or {}
    dx = float(rd.get("theta_x_gain_perturbation_operator_upper", math.inf))
    dp = float(rd.get("V12D_full_attitude_gain_perturbation_operator_upper", -math.inf))
    if not (math.isfinite(dx) and dx >= 0.0 and math.isfinite(dp) and dp >= 0.0 and dx <= FULL.up(dp)):
        f.append("invalid theta-x gain perturbation refinement")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    if d.get("P5_SAMPLE1_THETA_X_GAIN_PERTURBATION_V30") not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid V30 status")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=24)
    ap.add_argument("--tangent-pieces", type=int, default=24)
    ap.add_argument("--axial-pieces", type=int, default=24)
    ap.add_argument("--residual-x-pieces", type=int, default=6)
    ap.add_argument("--parallel-pieces", type=int, default=6)
    ap.add_argument("--current-component-pieces", type=int, default=4)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(x.domain, source_pieces=x.source_pieces,
              source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
              tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces,
              residual_x_pieces=x.residual_x_pieces,
              parallel_pieces=x.parallel_pieces,
              current_component_pieces=x.current_component_pieces)
    vf = validate(d); d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_THETA_X_GAIN_PERTURBATION_V30"],
        "theta_x_gain": d.get("theta_x_gain_perturbation_detail"),
        "perturbation": d.get("directional_perturbation_detail"),
        "q_current": d.get("current_q_upper"),
        "radial_lower": d.get("directional_radial_lower_rad"),
        "radial_upper": d.get("directional_radial_upper_rad"),
        "geodesic_q": d.get("geodesic_q_upper"),
        "product_W": d.get("product_abs_W_lower"),
        "product_q": d.get("product_q_upper"),
        "closed_q8": d.get("first_open_subbox_closed_inside_q8"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
