#!/usr/bin/env python3
"""V35 diagnostic: decompose the V12D perturbation at the first q8 witness.

This producer does not tighten or compose any q8 cell.  It records which part
of the certified V12D perturbation budget dominates the authoritative first
source-cell-0 witness after V34: first-PSD transport versus accepted S=0
pseudo-update effects, and scalar-x versus two-by-two nominal innovation
floors.  The output is intended to choose the next rigorous refinement rather
than to promote sample 1 or P5.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_theta_x_gain_perturbation_v30 as V30

DEFAULT_DOMAIN = V30.DEFAULT_DOMAIN
SCHEMA = 3500
V12D = V30.V29.V28.V27.V23.V22.V21B.V21.V12D


def _frac(part: float, total: float) -> float:
    p = float(part); t = float(total)
    if not (math.isfinite(p) and p >= 0.0 and math.isfinite(t) and t >= 0.0):
        return math.nan
    return 0.0 if t == 0.0 else p / t


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24) -> dict:
    path = Path(domain_path).resolve()
    core = V12D.build(
        path, source_pieces=source_pieces,
        source_cell_index=source_cell_index, p_pieces=p_pieces,
        tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    failures = [f"V12D: {x}" for x in V12D.validate(core)]
    if core.get("P5_SAMPLE1_FIRST_PSD_TANGENT_REFINEMENT_V12D") != "PASS":
        failures.append("V12D prerequisite did not pass")

    try:
        row = V30._witness_row(core)
    except Exception as exc:
        failures.append(f"witness row: {exc}")
        row = {}

    def num(key: str, default=math.nan) -> float:
        try:
            return float(row.get(key, default))
        except Exception:
            return math.nan

    dP_psd = num("PSD_reduced_covariance_perturbation_upper")
    dP_s = num("S_reduced_covariance_perturbation_upper")
    dP = num("total_reduced_covariance_perturbation_upper")
    dr_psd = num("PSD_residual_perturbation_upper_mps2")
    dr_s = num("S_residual_perturbation_upper_mps2")
    dr = num("total_residual_perturbation_upper_mps2")
    scalar_floor = num("scalar_innovation_lower")
    yz_floor = num("two_by_two_lambda_min_lower")
    nominal_floor = num("nominal_block_lambda_min_lower")

    detail = {
        "p_cell": row.get("p_cell"),
        "tangent_residual_cell": row.get("tangent_residual_cell"),
        "axial_residual_cell": row.get("axial_residual_cell"),
        "PSD_residual_perturbation_upper_mps2": dr_psd,
        "S_residual_perturbation_upper_mps2": dr_s,
        "total_residual_perturbation_upper_mps2": dr,
        "PSD_residual_fraction": _frac(dr_psd, dr),
        "S_residual_fraction": _frac(dr_s, dr),
        "PSD_reduced_covariance_perturbation_upper": dP_psd,
        "S_reduced_covariance_perturbation_upper": dP_s,
        "total_reduced_covariance_perturbation_upper": dP,
        "PSD_covariance_fraction": _frac(dP_psd, dP),
        "S_covariance_fraction": _frac(dP_s, dP),
        "sample1_H_perturbation_upper": num("sample1_H_perturbation_upper"),
        "sample1_innovation_perturbation_upper": num("sample1_innovation_perturbation_upper"),
        "sample1_attitude_cross_covariance_perturbation_upper": num(
            "sample1_attitude_cross_covariance_perturbation_upper"),
        "sample1_attitude_gain_operator_perturbation_upper": num(
            "sample1_attitude_gain_operator_perturbation_upper"),
        "scalar_innovation_lower": scalar_floor,
        "two_by_two_lambda_min_lower": yz_floor,
        "nominal_block_lambda_min_lower": nominal_floor,
        "scalar_inverse_upper": (math.inf if not scalar_floor > 0.0 else 1.0 / scalar_floor),
        "two_by_two_inverse_upper": (math.inf if not yz_floor > 0.0 else 1.0 / yz_floor),
        "nominal_block_inverse_operator_upper": num("nominal_block_inverse_operator_upper"),
        "actual_innovation_inverse_operator_upper": num("actual_innovation_inverse_operator_upper"),
        "actual_inverse_backend": row.get("actual_inverse_backend"),
        "nominal_inverse_times_innovation_perturbation_upper": num(
            "nominal_inverse_times_innovation_perturbation_upper"),
        "first_PSD_innovation_perturbation_tangent_only": row.get(
            "first_PSD_innovation_perturbation_tangent_only"),
        "first_PSD_innovation_axial_row_column_exact_zero": row.get(
            "first_PSD_innovation_axial_row_column_exact_zero"),
        "first_posterior_covariance_perturbation_upper": num(
            "first_posterior_covariance_perturbation_upper"),
        "reset_gauge_transform_perturbation_upper": num(
            "reset_gauge_transform_perturbation_upper"),
        "sample1_reduced_covariance_PSD_perturbation_upper": num(
            "sample1_reduced_covariance_PSD_perturbation_upper"),
    }

    finite_keys = (
        "total_residual_perturbation_upper_mps2",
        "total_reduced_covariance_perturbation_upper",
        "sample1_innovation_perturbation_upper",
        "sample1_attitude_gain_operator_perturbation_upper",
        "scalar_innovation_lower",
        "two_by_two_lambda_min_lower",
        "actual_innovation_inverse_operator_upper",
    )
    for key in finite_keys:
        x = detail.get(key)
        if not isinstance(x, (int, float)) or not math.isfinite(float(x)):
            failures.append(f"invalid diagnostic field {key}")
    if isinstance(detail.get("scalar_innovation_lower"), (int, float)) and not detail["scalar_innovation_lower"] > 0.0:
        failures.append("scalar innovation floor nonpositive")
    if isinstance(detail.get("two_by_two_lambda_min_lower"), (int, float)) and not detail["two_by_two_lambda_min_lower"] > 0.0:
        failures.append("two-by-two innovation floor nonpositive")

    # Choose only a diagnostic direction; no theorem state changes here.
    cov_dom = "PSD" if dP_psd >= dP_s else "S"
    res_dom = "PSD" if dr_psd >= dr_s else "S"
    inv_dom = "SCALAR_X" if scalar_floor <= yz_floor else "TWO_BY_TWO_YZ"
    next_obligation = (
        "REFINE_V12D_PSD_COMPONENT_MATRIX_AT_FIRST_Q8_WITNESS"
        if cov_dom == "PSD" and res_dom == "PSD" else
        "REFINE_V12D_ACCEPTED_S_PSEUDO_UPDATE_COMPONENT_MATRIX_AT_FIRST_Q8_WITNESS"
        if cov_dom == "S" and res_dom == "S" else
        "SPLIT_V12D_PSD_AND_S_COMPONENT_MATRICES_AT_FIRST_Q8_WITNESS"
    )

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_V12D_WITNESS_PERTURBATION_DECOMPOSITION_V35",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V12D_parent_revalidated": True,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "witness_perturbation_detail": detail,
        "dominant_covariance_component": cov_dom,
        "dominant_residual_component": res_dom,
        "nominal_inverse_limiting_block": inv_dom,
        "P5_SAMPLE1_V12D_WITNESS_PERTURBATION_DECOMPOSITION_V35": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": next_obligation,
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_V12D_WITNESS_PERTURBATION_DECOMPOSITION_V35":
        f.append("qualification mismatch")
    for k in ("source_generated_not_trajectory_fit", "V12D_parent_revalidated"):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed", "q8_composed_here",
              "q8_word_promoted_here", "whole_word_promoted_here",
              "N_H_words_set_here"):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if d.get("P5_SAMPLE1_V12D_WITNESS_PERTURBATION_DECOMPOSITION_V35") not in (
            "PASS", "NOT_ESTABLISHED"):
        f.append("invalid V35 status")
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
        "status": d["P5_SAMPLE1_V12D_WITNESS_PERTURBATION_DECOMPOSITION_V35"],
        "detail": d.get("witness_perturbation_detail"),
        "dominant_covariance": d.get("dominant_covariance_component"),
        "dominant_residual": d.get("dominant_residual_component"),
        "limiting_inverse_block": d.get("nominal_inverse_limiting_block"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
