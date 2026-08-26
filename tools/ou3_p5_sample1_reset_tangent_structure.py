#!/usr/bin/env python3
"""Low-dimensional P5 sample-1 reset/yaw covariance structure diagnostic.

At the first accepted accelerometer prefix the exact source gauge has force
+g e3 and the two observable tangent channels are identical scalar Kalman
updates.  Ignoring only the already separately bounded tiny PSD attitude
remainder, one tangent/aw channel has prior covariance diag(t,p), H=[g,1].
The exact first posterior tangent variance is

    a = t (p+r) / (g^2 t+p+r).

The gravity-parallel yaw variance Y is not observed by that packet.  Rotational
symmetry lets the first attitude correction be placed on +x with magnitude d.
The shipping left-error covariance reset then has

    G_yz = [[1,-d/2],[d/2,1]],

and the proof/body gauge rotates that same yz plane by R_x(d).  Consequently
all reset-induced tangent/yaw covariance is contained in the exact 2x2 map

    L_yz = R(d) G_yz,
    P_yz+ = L_yz diag(a,Y) L_yz^T.

This producer encloses that map over source-derived p and d cells and reports
how much the second tangent attitude variance can grow solely from the reset
and yaw seed.  It is a diagnostic prerequisite: process noise, tangent-force
misalignment, sample-1 S due/not-due, and the second measurement are not yet
included, so it cannot promote P5 or N_H_words.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_mul, matrix_transpose
import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_sample1_rotation_gauge_refinement_v2 as SUB
import ou3_validated_transcendentals as VT
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = FIRST.DEFAULT_DOMAIN
SCHEMA = 1


def I(x: float) -> Interval:
    return FULL.I(float(x))


def _rx_yz(d: Interval):
    # d lies in [0,1.45], inside the validated monotone sin/cos point range.
    s = Interval(VT.sin_point(d.lo).lo, VT.sin_point(d.hi).hi)
    c = Interval(VT.cos_point(d.hi).lo, VT.cos_point(d.lo).hi)
    return [[c, -s], [s, c]]


def _g_yz(d: Interval):
    h = I(0.5)
    return [[I(1.0), -h*d], [h*d, I(1.0)]]


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 32, d_pieces: int = 32) -> dict:
    path = Path(domain_path).resolve()
    first = FIRST.build(path, source_pieces=source_pieces)
    vector = VECTOR.build()
    failures = [f"first: {x}" for x in FIRST.validate(first)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]
    src, phase = RG._source_phase_children(source_pieces)[source_cell_index]
    if phase != "due":
        failures.append("expected source-cell-0 due witness")

    fr = first["source_cells"][source_cell_index]
    p_all = Interval.outward_bounds(*map(float, fr["P_aw_variance_interval"]))
    dmax = float(fr["correction_norm_upper_rad"])
    h = float(src["dt_s"])
    g = float(json.loads(path.read_text())["startup"]["gravity_mps2"])
    tilt, yaw, eps = RG._attitude_covariance_epsilon(path, h)
    racc = FULL._R_diag(float(vector["configured_measurement_bounds"]["acc_measurement_std_mps2"]))[0][0]

    # The separately proved 0<=E<=eps I remainder is retained here by widening
    # both source attitude variances before the exact reset/gauge map.
    t = Interval.outward_bounds(tilt, FULL.up(tilt + eps))
    Y = Interval.outward_bounds(yaw, FULL.up(yaw + eps))
    pcells = SUB.parts(p_all.lo, p_all.hi, p_pieces)
    dcells = SUB.parts(0.0, dmax, d_pieces)

    rows = []
    max_y = 0.0
    max_z = 0.0
    max_yz = 0.0
    max_ratio = 0.0
    first_bad = None
    for pi, p in enumerate(pcells):
        D = I(g*g)*t + p + racc
        if D.lo <= 0.0:
            failures.append(f"p cell {pi}: first tangent innovation lost positivity")
            continue
        a = t * (p + racc) / D
        for di, d in enumerate(dcells):
            L = matrix_mul(_rx_yz(d), _g_yz(d))
            P0 = [[a, I(0.0)], [I(0.0), Y]]
            P = matrix_mul(matrix_mul(L, P0), matrix_transpose(L))
            vy = P[0][0]
            vz = P[1][1]
            cyz = P[0][1]
            ratio = FULL.up(vy.hi / max(a.lo, math.nextafter(0.0, math.inf)))
            finite = all(math.isfinite(x) for x in (vy.hi, vz.hi, cyz.abs_upper(), ratio))
            row = {
                "p_cell": pi,
                "d_cell": di,
                "P_aw_variance": p.as_list(),
                "first_correction_rad": d.as_list(),
                "first_posterior_tangent_variance": a.as_list(),
                "yaw_variance_before_reset": Y.as_list(),
                "post_reset_gauge_tangent_variance": vy.as_list(),
                "post_reset_gauge_yaw_variance": vz.as_list(),
                "post_reset_gauge_tangent_yaw_covariance": cyz.as_list(),
                "tangent_variance_multiplier_upper": ratio,
                "finite": finite,
            }
            rows.append(row)
            max_y = max(max_y, vy.hi)
            max_z = max(max_z, vz.hi)
            max_yz = max(max_yz, cyz.abs_upper())
            max_ratio = max(max_ratio, ratio)
            if not finite and first_bad is None:
                first_bad = row

    ok = bool(rows) and first_bad is None and not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_EXACT_RESET_TANGENT_YAW_STRUCTURE_DIAGNOSTIC",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "first_tangent_scalar_posterior_used": True,
        "first_correction_axis_fixed_by_rotational_symmetry": True,
        "shipping_left_error_reset_exact_yz_block_used": True,
        "proof_gauge_exact_yz_rotation_used": True,
        "tiny_attitude_PSD_remainder_retained_as_variance_widening": True,
        "process_noise_included_here": False,
        "sample1_tangent_force_misalignment_included_here": False,
        "sample1_S_due_not_due_included_here": False,
        "second_accelerometer_gain_computed_here": False,
        "complete_sample1_branch_closed_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "evaluated_cells": len(rows),
        "first_accel_correction_norm_upper_rad": dmax,
        "attitude_tilt_variance_interval": t.as_list(),
        "attitude_yaw_variance_interval": Y.as_list(),
        "attitude_PSD_remainder_upper": eps,
        "max_post_reset_gauge_tangent_variance": max_y,
        "max_post_reset_gauge_yaw_variance": max_z,
        "max_post_reset_gauge_tangent_yaw_covariance_abs": max_yz,
        "max_tangent_variance_multiplier_upper": max_ratio,
        "first_nonfinite_cell": first_bad,
        "P5_SAMPLE1_RESET_TANGENT_STRUCTURE_DIAGNOSTIC": "PASS" if ok else "NOT_ESTABLISHED",
        "next_obligation": "ADD_ONE_STEP_PROCESS_AND_SAMPLE1_S_BRANCH_TO_STRUCTURED_TANGENT_YAW_BLOCK",
        "failures": failures,
        "rows": rows,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit", "first_tangent_scalar_posterior_used",
        "first_correction_axis_fixed_by_rotational_symmetry",
        "shipping_left_error_reset_exact_yz_block_used", "proof_gauge_exact_yz_rotation_used",
        "tiny_attitude_PSD_remainder_retained_as_variance_widening",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "process_noise_included_here",
        "sample1_tangent_force_misalignment_included_here", "sample1_S_due_not_due_included_here",
        "second_accelerometer_gain_computed_here", "complete_sample1_branch_closed_here",
        "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if int(d.get("evaluated_cells", 0)) <= 0:
        f.append("no reset/yaw cells")
    for k in (
        "max_post_reset_gauge_tangent_variance", "max_post_reset_gauge_yaw_variance",
        "max_post_reset_gauge_tangent_yaw_covariance_abs", "max_tangent_variance_multiplier_upper",
    ):
        x = float(d.get(k, math.nan))
        if not (math.isfinite(x) and x >= 0.0):
            f.append(f"invalid {k}")
    if not f and d.get("P5_SAMPLE1_RESET_TANGENT_STRUCTURE_DIAGNOSTIC") != "PASS":
        f.append("reset tangent structure diagnostic did not pass")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=32)
    ap.add_argument("--d-pieces", type=int, default=32)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain, source_pieces=a.source_pieces, source_cell_index=a.source_cell_index,
              p_pieces=a.p_pieces, d_pieces=a.d_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_RESET_TANGENT_STRUCTURE_DIAGNOSTIC"],
        "cells": d["evaluated_cells"],
        "dmax": d["first_accel_correction_norm_upper_rad"],
        "max_tangent_var": d["max_post_reset_gauge_tangent_variance"],
        "max_yaw_var": d["max_post_reset_gauge_yaw_variance"],
        "max_cross": d["max_post_reset_gauge_tangent_yaw_covariance_abs"],
        "max_multiplier": d["max_tangent_variance_multiplier_upper"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
