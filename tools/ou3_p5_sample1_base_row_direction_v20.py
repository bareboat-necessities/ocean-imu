#!/usr/bin/env python3
"""V20: source-correlated base-row direction audit for the sample-1 P5 q<8 route.

V19 showed that subdividing the *final* current/correction y-z product is too
late: 1.19 million sub-pair evaluations closed no additional q<8 cells.  The
remaining dependency loss occurs earlier.  V10 already fixes a canonical SO(2)
gravity gauge for each first-accelerometer residual row, but V14/V18 rebuild the
pre-correction Cayley state from the symmetric box

    c_x,c_y in [-c_t,c_t],  c_z in [-q,q],

without asking whether that attitude box can generate the same V10 residual
row.

For the canonical gauge, positive first x-axis correction corresponds to the
first residual tangent component on the negative y axis.  Write the source row
as

    r0 = [0,-r_t,r_z] = y_R(c) + n,
    ||n|| <= ||e_aw^-|| + ||b_a||,

where the exact gravity rotational residual for the Cayley vector c is

    y_R/g = [ (2 c_x c_z + 4 c_y)/D,
              (2 c_y c_z - 4 c_x)/D,
             -2(c_x^2+c_y^2)/D ],
    D = 4 + ||c||^2.

This identity has the same small-angle Jacobian -[g e3]_x used by the shipping
accelerometer update.  V20 partitions the existing source Cayley box, evaluates
that exact residual interval on every subbox, and discards a subbox only when
its minimum possible distance from the V10 residual row is strictly larger than
the independently certified nuisance norm.  The P1 full-Cayley and gravity-
tangent balls are retained as additional compatibility filters.

The stage is intentionally diagnostic: it proves a source-compatible directional
cover for each V10 (tangent, axial) base row but does not yet compose those
subboxes through the first reset or claim q<8.  It changes no estimator, source
domain, correction limit, or theorem-promotion state.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, hull
import ou3_p5_sample1_rotation_gauge_refinement_v2 as SUB
import ou3_p5_sample1_structured_full_gain_v10 as V10

DEFAULT_DOMAIN = V10.DEFAULT_DOMAIN
SCHEMA = 2000
FULL = V10.FULL
FIRST_Q8_WITNESS_TANGENT_CELL = 0
FIRST_Q8_WITNESS_AXIAL_CELL = 19


def _minimum_abs(x: Interval) -> float:
    if x.lo <= 0.0 <= x.hi:
        return 0.0
    return min(abs(x.lo), abs(x.hi))


def _box_norm_lower(v) -> float:
    s = 0.0
    for x in v:
        a = _minimum_abs(x)
        s = FULL.down(s + FULL.down(a * a))
    if s <= 0.0:
        return 0.0
    return FULL.down(math.sqrt(s))


def _interval_distance_lower(a: Interval, b: Interval) -> float:
    """Lower bound on distance between two real intervals."""
    if a.hi >= b.lo and b.hi >= a.lo:
        return 0.0
    gap = b.lo - a.hi if a.hi < b.lo else a.lo - b.hi
    return max(0.0, FULL.down(gap))


def _vector_box_distance_lower(a, b) -> float:
    if len(a) != len(b):
        raise ValueError("equal vector dimensions required")
    s = 0.0
    for x, y in zip(a, b):
        d = _interval_distance_lower(x, y)
        s = FULL.down(s + FULL.down(d * d))
    if s <= 0.0:
        return 0.0
    return FULL.down(math.sqrt(s))


def _gravity_residual_from_cayley(c, gravity: float):
    """Exact interval y_R=(R(c)-I) g e3 in the canonical gravity gauge."""
    if len(c) != 3:
        raise ValueError("three-component Cayley vector required")
    cx, cy, cz = c
    two = FULL.I(2.0)
    four = FULL.I(4.0)
    den = four + cx.square() + cy.square() + cz.square()
    if not den.lo > 0.0:
        raise RuntimeError("Cayley gravity residual lost positive denominator")
    g = FULL.I(float(gravity))
    yx = g * (two * cx * cz + four * cy) / den
    yy = g * (two * cy * cz - four * cx) / den
    yz = -(g * two * (cx.square() + cy.square()) / den)
    return [yx, yy, yz]


def _candidate_cayley_boxes(qpre: float, ctan: float, *,
                            tangent_direction_pieces: int,
                            yaw_direction_pieces: int):
    if tangent_direction_pieces < 2 or yaw_direction_pieces < 2:
        raise ValueError("V20 direction subdivision counts must be >=2")
    cx_parts = SUB.parts(-ctan, ctan, tangent_direction_pieces)
    cy_parts = SUB.parts(-ctan, ctan, tangent_direction_pieces)
    cz_parts = SUB.parts(-qpre, qpre, yaw_direction_pieces)
    out = []
    rejected_tangent = rejected_full = 0
    for cx in cx_parts:
        for cy in cy_parts:
            if _box_norm_lower((cx, cy)) > ctan:
                rejected_tangent += len(cz_parts)
                continue
            for cz in cz_parts:
                box = (cx, cy, cz)
                if _box_norm_lower(box) > qpre:
                    rejected_full += 1
                    continue
                out.append(box)
    return out, rejected_tangent, rejected_full, len(cx_parts) * len(cy_parts) * len(cz_parts)


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          tangent_direction_pieces: int = 4, yaw_direction_pieces: int = 4) -> dict:
    path = Path(domain_path).resolve()
    dom = json.loads(path.read_text(encoding="utf-8"))
    core = V10.build(
        path,
        source_pieces=source_pieces,
        source_cell_index=source_cell_index,
        p_pieces=p_pieces,
        tangent_pieces=tangent_pieces,
        axial_pieces=axial_pieces,
    )
    failures = [f"V10: {x}" for x in V10.validate(core)]
    if core.get("P5_SAMPLE1_COMBINED_PERPENDICULAR_RESIDUAL_V10") != "PASS":
        failures.append("V10 base-row prerequisite did not pass")

    first = V10.FIRST.build(path, source_pieces=source_pieces)
    failures += [f"first: {x}" for x in V10.FIRST.validate(first)]
    source_children = V10.RG._source_phase_children(source_pieces)
    if source_cell_index < 0 or source_cell_index >= len(source_children):
        raise ValueError("invalid V20 source cell index")
    _src, phase = source_children[source_cell_index]
    if phase != "due":
        failures.append("V20 focused family requires a first-due source cell")

    fr = first["source_cells"][source_cell_index]
    qpre = float(first["post_prediction_full_cayley_norm_upper"])
    ctan = float(first["post_prediction_cayley_tangent_norm_upper"])
    aw_pre = float(fr["predicted_aw_error_norm_upper_mps2"])
    ba = float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    nuisance = FULL.up(aw_pre + ba)
    gravity = float(dom["startup"]["gravity_mps2"])
    if not all(math.isfinite(x) and x >= 0.0 for x in (qpre, ctan, aw_pre, ba, nuisance, gravity)):
        failures.append("invalid finite nonnegative V20 source bounds")

    candidates, ball_reject_tan, ball_reject_full, candidate_grid = _candidate_cayley_boxes(
        qpre, ctan,
        tangent_direction_pieces=tangent_direction_pieces,
        yaw_direction_pieces=yaw_direction_pieces,
    )
    residual_cache = [(_gravity_residual_from_cayley(box, gravity), box) for box in candidates]
    if not candidates:
        failures.append("V20 Cayley direction subdivision has no ball-compatible boxes")

    # V10 repeats each (tangent,axial) source row over p cells.  Hull any tiny
    # outward differences across p so the directional audit is source complete.
    base_rows = {}
    for row in core.get("rows", []):
        key = (int(row["tangent_residual_cell"]), int(row["axial_residual_cell"]))
        rt = Interval.outward_bounds(*map(float, row["first_tangent_residual_magnitude_mps2"]))
        rz = Interval.outward_bounds(*map(float, row["first_axial_residual_mps2"]))
        if key in base_rows:
            old_rt, old_rz = base_rows[key]
            base_rows[key] = (hull(old_rt, rt), hull(old_rz, rz))
        else:
            base_rows[key] = (rt, rz)

    parent_cx = Interval.outward_bounds(-ctan, ctan)
    parent_cy = Interval.outward_bounds(-ctan, ctan)
    parent_cz = Interval.outward_bounds(-qpre, qpre)
    rows = []
    total_compatibility_checks = total_rejected = 0
    refined_rows = empty_rows = 0
    min_survival_fraction = 1.0
    max_survival_fraction = 0.0
    first_refined = first_empty = None
    witness_row = None

    for key in sorted(base_rows):
        rt, rz = base_rows[key]
        # The V7/V10 canonical gain has K_theta,x,y<0 at gravity alignment, so
        # positive x-axis correction d corresponds to residual tangent -r_t e_y.
        target = [FULL.I(0.0), Interval.outward_bounds(-rt.hi, -rt.lo), rz]
        survivors = []
        max_rejected_distance = 0.0
        for residual, box in residual_cache:
            total_compatibility_checks += 1
            dmin = _vector_box_distance_lower(residual, target)
            if dmin > nuisance:
                total_rejected += 1
                max_rejected_distance = max(max_rejected_distance, dmin)
                continue
            survivors.append(box)

        fraction = len(survivors) / len(candidates) if candidates else 0.0
        min_survival_fraction = min(min_survival_fraction, fraction)
        max_survival_fraction = max(max_survival_fraction, fraction)
        row = {
            "tangent_residual_cell": key[0],
            "axial_residual_cell": key[1],
            "first_tangent_residual_magnitude_mps2": rt.as_list(),
            "first_axial_residual_mps2": rz.as_list(),
            "candidate_cayley_boxes": len(candidates),
            "surviving_cayley_boxes": len(survivors),
            "rejected_cayley_boxes": len(candidates) - len(survivors),
            "survival_fraction": fraction,
            "max_proved_incompatible_distance_lower_mps2": max_rejected_distance,
            "source_row_incompatible": not survivors,
        }
        if not survivors:
            empty_rows += 1
            row["cayley_component_hull"] = None
            if first_empty is None:
                first_empty = dict(row)
        else:
            hx = hull(*(b[0] for b in survivors))
            hy = hull(*(b[1] for b in survivors))
            hz = hull(*(b[2] for b in survivors))
            # Count a refinement only when the surviving hull moves an endpoint
            # inward.  Outward-rounded subcell overlap may make an unchanged
            # full cover a few ulps wider than the parent and must not count as
            # proof tightening.
            refined = (
                hx.lo > parent_cx.lo or hx.hi < parent_cx.hi
                or hy.lo > parent_cy.lo or hy.hi < parent_cy.hi
                or hz.lo > parent_cz.lo or hz.hi < parent_cz.hi
            )
            row["cayley_component_hull"] = [hx.as_list(), hy.as_list(), hz.as_list()]
            row["directionally_refined"] = refined
            if refined:
                refined_rows += 1
                if first_refined is None:
                    first_refined = dict(row)
        if key == (FIRST_Q8_WITNESS_TANGENT_CELL, FIRST_Q8_WITNESS_AXIAL_CELL):
            witness_row = dict(row)
        rows.append(row)

    if witness_row is None:
        failures.append("V20 did not find the V19 first-q8-witness base row")
    audited = bool(rows) and not failures
    witness_useful = bool(witness_row) and (
        witness_row.get("source_row_incompatible") is True
        or witness_row.get("directionally_refined") is True
    )
    next_obligation = (
        "PROPAGATE_V20_BASE_ROW_DIRECTION_CELLS_THROUGH_V18B_Q8_COMPOSITION"
        if witness_useful else
        "DERIVE_TIGHTER_JOINT_ATTITUDE_RESIDUAL_MAP_AT_FIRST_Q8_WITNESS"
    )
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_SOURCE_CORRELATED_BASE_ROW_DIRECTION_V20",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V10_direct_first_residual_coordinate_prerequisite_retained": True,
        "V10_canonical_SO2_gravity_gauge_retained": True,
        "positive_x_correction_uses_negative_y_tangent_residual_branch": True,
        "exact_cayley_gravity_residual_used": True,
        "first_aw_and_bias_nuisance_combined_before_direction_rejection": True,
        "P1_full_cayley_ball_retained": True,
        "P1_gravity_tangent_cayley_ball_retained": True,
        "base_row_direction_subdivision_is_source_complete": True,
        "source_cell_index": int(source_cell_index),
        "post_prediction_full_cayley_norm_upper": qpre,
        "post_prediction_cayley_tangent_norm_upper": ctan,
        "pre_first_aw_error_norm_upper_mps2": aw_pre,
        "accelerometer_bias_error_norm_upper_mps2": ba,
        "combined_first_residual_nuisance_norm_upper_mps2": nuisance,
        "gravity_mps2": gravity,
        "tangent_direction_pieces": int(tangent_direction_pieces),
        "yaw_direction_pieces": int(yaw_direction_pieces),
        "candidate_cayley_grid_boxes": int(candidate_grid),
        "ball_compatible_candidate_cayley_boxes": len(candidates),
        "tangent_ball_rejected_grid_boxes": int(ball_reject_tan),
        "full_ball_rejected_grid_boxes": int(ball_reject_full),
        "evaluated_base_direction_rows": len(rows),
        "source_incompatible_base_direction_rows": int(empty_rows),
        "directionally_refined_base_direction_rows": int(refined_rows),
        "total_direction_compatibility_checks": int(total_compatibility_checks),
        "total_direction_incompatible_subboxes": int(total_rejected),
        "minimum_survival_fraction": float(min_survival_fraction),
        "maximum_survival_fraction": float(max_survival_fraction),
        "first_directionally_refined_base_row": first_refined,
        "first_source_incompatible_base_row": first_empty,
        "first_v19_q8_witness_base_direction_row": witness_row,
        "base_direction_rows": rows,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_SOURCE_CORRELATED_BASE_ROW_DIRECTION_V20": "PASS" if audited else "NOT_ESTABLISHED",
        "next_obligation": next_obligation,
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_SOURCE_CORRELATED_BASE_ROW_DIRECTION_V20":
        failures.append("qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "V10_direct_first_residual_coordinate_prerequisite_retained",
        "V10_canonical_SO2_gravity_gauge_retained",
        "positive_x_correction_uses_negative_y_tangent_residual_branch",
        "exact_cayley_gravity_residual_used",
        "first_aw_and_bias_nuisance_combined_before_direction_rejection",
        "P1_full_cayley_ball_retained",
        "P1_gravity_tangent_cayley_ball_retained",
        "base_row_direction_subdivision_is_source_complete",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "q8_composed_here", "q8_word_promoted_here", "whole_word_promoted_here",
        "N_H_words_set_here",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        failures.append("deployed correction limit changed")
    rows = int(d.get("evaluated_base_direction_rows", -1))
    refined = int(d.get("directionally_refined_base_direction_rows", -1))
    empty = int(d.get("source_incompatible_base_direction_rows", -1))
    checks = int(d.get("total_direction_compatibility_checks", -1))
    rejected = int(d.get("total_direction_incompatible_subboxes", -1))
    if not (rows > 0 and rows >= refined >= 0 and rows >= empty >= 0):
        failures.append("invalid V20 base-row accounting")
    if not (checks >= rejected >= 0):
        failures.append("invalid V20 compatibility accounting")
    lo = float(d.get("minimum_survival_fraction", math.nan))
    hi = float(d.get("maximum_survival_fraction", math.nan))
    if not (0.0 <= lo <= hi <= 1.0):
        failures.append("invalid V20 survival fractions")
    st = d.get("P5_SAMPLE1_SOURCE_CORRELATED_BASE_ROW_DIRECTION_V20")
    if st == "PASS":
        if d.get("first_v19_q8_witness_base_direction_row") is None:
            failures.append("V20 PASS lacks first-q8-witness base row")
    elif st == "NOT_ESTABLISHED":
        if not failures:
            failures.append("V20 NOT_ESTABLISHED lacks validation failure")
    else:
        failures.append("invalid V20 status")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=24)
    ap.add_argument("--tangent-pieces", type=int, default=24)
    ap.add_argument("--axial-pieces", type=int, default=24)
    ap.add_argument("--tangent-direction-pieces", type=int, default=4)
    ap.add_argument("--yaw-direction-pieces", type=int, default=4)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(
        args.domain,
        source_pieces=args.source_pieces,
        source_cell_index=args.source_cell_index,
        p_pieces=args.p_pieces,
        tangent_pieces=args.tangent_pieces,
        axial_pieces=args.axial_pieces,
        tangent_direction_pieces=args.tangent_direction_pieces,
        yaw_direction_pieces=args.yaw_direction_pieces,
    )
    vf = validate(out)
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_SAMPLE1_SOURCE_CORRELATED_BASE_ROW_DIRECTION_V20"],
        "base_rows": out["evaluated_base_direction_rows"],
        "candidate_boxes": out["ball_compatible_candidate_cayley_boxes"],
        "refined_rows": out["directionally_refined_base_direction_rows"],
        "empty_rows": out["source_incompatible_base_direction_rows"],
        "checks": out["total_direction_compatibility_checks"],
        "rejected_subboxes": out["total_direction_incompatible_subboxes"],
        "min_survival_fraction": out["minimum_survival_fraction"],
        "max_survival_fraction": out["maximum_survival_fraction"],
        "first_refined": out["first_directionally_refined_base_row"],
        "first_empty": out["first_source_incompatible_base_row"],
        "first_q8_witness_row": out["first_v19_q8_witness_base_direction_row"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
