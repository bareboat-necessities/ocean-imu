#!/usr/bin/env python3
"""V25: use the certified V10 structured block gain in V24's radial ball.

V24 correctly preserves the physical sample-1 residual nuisance and the V12D
PSD/S correction remainder as Euclidean balls, but reconstructed ||K_theta||
from interval component formulas.  At the first q8 witness that reconstruction
is 2.51034 and is far looser than the structured block norm already certified
by the V10 row used by V21 itself.

V10 stores independent rigorous bounds for the two orthogonal attitude-gain
blocks.  For the exact one-plus-two structure,

    ||K_theta||_2 <= max(Ktheta_perpendicular_block_upper,
                         Ktheta_parallel_block_upper).

V25 rebuilds the deterministic V10 first-witness row, validates it, and runs the
unchanged V24 radial-ball/subdivision calculation with exactly this certified
scalar in place of V24's later component-interval reconstruction.  No current
box, residual algebra, correction component enclosure, source domain, estimator,
shipping limit, q target, or promotion state is changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_radial_nuisance_ball_subdivision_v24 as V24

DEFAULT_DOMAIN = V24.DEFAULT_DOMAIN
SCHEMA = 2500
Q_TARGET = V24.Q_TARGET
WITNESS_P_CELL = 0
WITNESS_TANGENT_CELL = 0
WITNESS_AXIAL_CELL = 19


def _certified_v10_gain_norm(core: dict) -> tuple[float, dict]:
    rows = core.get("rows", [])
    for row in rows:
        ids = (int(row["p_cell"]), int(row["tangent_residual_cell"]),
               int(row["axial_residual_cell"]))
        if ids != (WITNESS_P_CELL, WITNESS_TANGENT_CELL, WITNESS_AXIAL_CELL):
            continue
        kp = float(row["Ktheta_perpendicular_block_upper"])
        ka = float(row["Ktheta_parallel_block_upper"])
        if not (math.isfinite(kp) and kp >= 0.0 and math.isfinite(ka) and ka >= 0.0):
            raise RuntimeError("invalid V10 structured gain block bound")
        return max(kp, ka), {
            "Ktheta_perpendicular_block_upper": kp,
            "Ktheta_parallel_block_upper": ka,
            "p_cell": WITNESS_P_CELL,
            "tangent_residual_cell": WITNESS_TANGENT_CELL,
            "axial_residual_cell": WITNESS_AXIAL_CELL,
        }
    raise RuntimeError("V10 first-q8-witness row not found")


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6,
          current_component_pieces: int = 4) -> dict:
    path = Path(domain_path).resolve()
    V10 = V24.V23.V22.V21B.V21.V12D.V11.V10
    core = V10.build(
        path, source_pieces=source_pieces, source_cell_index=source_cell_index,
        p_pieces=p_pieces, tangent_pieces=tangent_pieces,
        axial_pieces=axial_pieces)
    failures = [f"V10: {x}" for x in V10.validate(core)]
    try:
        k0, detail = _certified_v10_gain_norm(core)
    except Exception as exc:
        failures.append(f"V10 gain row: {exc}")
        k0 = math.inf
        detail = None

    original_gain_norm = V24._gain_operator_norm
    if math.isfinite(k0):
        V24._gain_operator_norm = lambda _parent: k0
    try:
        parent = V24.build(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces,
            current_component_pieces=current_component_pieces)
    finally:
        V24._gain_operator_norm = original_gain_norm

    failures += [f"V24: {x}" for x in V24.validate(parent)]
    if parent.get("P5_SAMPLE1_RADIAL_NUISANCE_BALL_SUBDIVISION_V24") != "PASS":
        failures.append("V24 radial nuisance-ball prerequisite did not pass")
    used = float(parent.get("gain_operator_norm_upper", math.inf))
    if math.isfinite(k0) and used != k0:
        failures.append("V25 did not use certified V10 structured gain norm")

    focused_closed = bool(
        parent.get("focused_first_witness_signed_subcell_closed_by_radial_ball")
        and not failures)
    out = dict(parent)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_V10_GAIN_RADIAL_BALL_SUBDIVISION_V25",
        "V24_radial_nuisance_ball_parent_retained": True,
        "V10_structured_gain_row_revalidated": True,
        "V10_orthogonal_block_operator_norm_used": True,
        "V24_component_gain_norm_reconstruction_superseded": True,
        "V10_structured_gain_detail": detail,
        "V10_structured_gain_operator_norm_upper": k0,
        "focused_first_witness_signed_subcell_closed_by_V10_gain_ball": focused_closed,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_V10_GAIN_RADIAL_BALL_SUBDIVISION_V25": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_V25_V10_GAIN_RADIAL_BALL_ROUTE_OVER_FIRST_BASE_ROW"
            if focused_closed else
            "DERIVE_SOURCE_CORRELATED_SAMPLE1_AW_ERROR_COMPONENTS_AT_FIRST_OPEN_SUBBOX"),
        "failures": list(dict.fromkeys(failures)),
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_V10_GAIN_RADIAL_BALL_SUBDIVISION_V25":
        f.append("qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "V24_radial_nuisance_ball_parent_retained",
        "V10_structured_gain_row_revalidated",
        "V10_orthogonal_block_operator_norm_used",
        "V24_component_gain_norm_reconstruction_superseded",
        "all_candidate_current_subboxes_accounted",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "q8_composed_here", "q8_word_promoted_here", "whole_word_promoted_here",
        "N_H_words_set_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    k = d.get("V10_structured_gain_operator_norm_upper")
    used = d.get("gain_operator_norm_upper")
    if not (isinstance(k, (int, float)) and math.isfinite(float(k)) and float(k) >= 0.0):
        f.append("invalid V10 structured gain norm")
    elif not (isinstance(used, (int, float)) and float(used) == float(k)):
        f.append("V24 radial route did not use V10 structured gain norm")
    if d.get("focused_first_witness_signed_subcell_closed_by_V10_gain_ball") is True \
            and int(d.get("open_current_subboxes", -1)) != 0:
        f.append("V25 claims focused closure with open current subboxes")
    if d.get("P5_SAMPLE1_V10_GAIN_RADIAL_BALL_SUBDIVISION_V25") not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid V25 status")
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
    d = build(
        x.domain, source_pieces=x.source_pieces,
        source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
        tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces,
        residual_x_pieces=x.residual_x_pieces,
        parallel_pieces=x.parallel_pieces,
        current_component_pieces=x.current_component_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_V10_GAIN_RADIAL_BALL_SUBDIVISION_V25"],
        "V10_gain_detail": d.get("V10_structured_gain_detail"),
        "V10_gain_norm": d.get("V10_structured_gain_operator_norm_upper"),
        "previous_component_gain_norm": 2.510342917347482,
        "nuisance_correction_ball": d.get("physical_nuisance_correction_ball_upper_rad"),
        "V12D_correction_ball": d.get("V12D_correction_perturbation_ball_upper_rad"),
        "radial_refined": d.get("radial_ball_strictly_refined_current_subboxes"),
        "candidate": d.get("candidate_current_subboxes"),
        "compatible": d.get("compatible_current_subboxes"),
        "closed": d.get("closed_current_subboxes"),
        "open": d.get("open_current_subboxes"),
        "min_radial_upper": d.get("minimum_compatible_radial_upper_rad"),
        "max_radial_upper": d.get("maximum_compatible_radial_upper_rad"),
        "min_best_q": d.get("minimum_best_q_upper"),
        "max_best_q": d.get("maximum_best_q_upper"),
        "focused_closed": d.get("focused_first_witness_signed_subcell_closed_by_V10_gain_ball"),
        "first_open": d.get("first_open_current_subbox"),
        "worst_open": d.get("worst_open_current_subbox"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
