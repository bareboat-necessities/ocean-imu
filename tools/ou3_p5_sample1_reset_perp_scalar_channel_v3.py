#!/usr/bin/env python3
"""Source-state-correlated V3 of the reset-sensitive sample-1 scalar channel.

V2 fixed the covariance dependency but still multiplied every cell by the raw
normal-Live row bound 30.5 m/s^2.  That Cartesian product is not the source
state reached after the first accepted accelerometer update.

In the ideal gravity-gauged first tangent block, SO(2) symmetry places the
first attitude correction on +e1.  The perpendicular attitude correction is
then exactly zero, hence the corresponding first accelerometer residual row is
zero.  For that row the exact finite measurement identity gives

    0 = y_R,0 + u_aw,0 + b_a,0,

where u_aw is the physical latent-acceleration error expressed in the first
body gauge.  Thus |u_aw,0| is bounded by the first gravity chord residual plus
the bias bound; it is not an independent 30 m/s^2 measurement.

The deployed reset/proof gauge is an x-axis rotation, so it leaves this aw_x
component unchanged.  One next OU prediction multiplies it by alpha; the
source-enclosed 5-ms body rotation only mixes in the other aw components by its
small off-diagonal bound.  The sample-1 aligned-force rotational residual is
bounded by the gravity-direction chord after: first handoff, the actual first
correction cell, and one bounded prediction.  Therefore every V2 covariance
cell receives a source-correlated residual bound rather than a global raw
measurement bound.

The small first attitude-PSD cross-axis remainder can make the canonical
perpendicular first correction nonzero; that term is intentionally not hidden
here and remains the next perturbation obligation.  Sample-1 tangent-force,
body-rate details beyond the already included chord transport, and the proved
1.83e-12-rad S attitude injection remain explicit later perturbations.  This is
not a complete sample-1/P5 promotion.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_p5_sample1_reset_perp_scalar_channel_v2 as V2
import ou3_validated_transcendentals as VT

DEFAULT_DOMAIN = FIRST.DEFAULT_DOMAIN
SCHEMA = 3
RANGE = 9.0


def _chord_from_cos_lower(c: float) -> float:
    if not (-1.0 <= c <= 1.0):
        raise ValueError("cosine lower outside [-1,1]")
    return min(2.0, FULL.up(math.sqrt(max(0.0, FULL.up(2.0 * FULL.up(1.0 - c))))))


def _correction_chord_upper(d_hi: float) -> float:
    if not (0.0 <= d_hi < math.pi):
        raise ValueError("first correction outside monotone half-angle range")
    s = VT.sin_point(FULL.up(0.5 * d_hi))
    return min(2.0, FULL.up(2.0 * s.hi))


def _interval_abs_upper(a) -> float:
    return max(abs(float(a[0])), abs(float(a[1])))


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 32,
          d_pieces: int = 32, axial_pieces: int = 32) -> dict:
    FULL3._install_backend()
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    first = FIRST.build(path, source_pieces=source_pieces)
    core = V2.build(path, source_pieces=source_pieces,
                    source_cell_index=source_cell_index,
                    p_pieces=p_pieces, d_pieces=d_pieces,
                    axial_pieces=axial_pieces)
    failures = [f"first: {x}" for x in FIRST.validate(first)]
    failures += [f"V2: {x}" for x in V2.validate(core)]

    src, phase = RG._source_phase_children(source_pieces)[source_cell_index]
    if phase != "due":
        failures.append("expected due source witness")
    fr = first["source_cells"][source_cell_index]
    g = float(domain["startup"]["gravity_mps2"])
    ba = float(domain["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    aw_pred = float(fr["predicted_aw_error_norm_upper_mps2"])
    h = float(src["dt_s"])
    tilt, _yaw, eps = RG._attitude_covariance_epsilon(path, h)
    t_lo = float(tilt)
    if not t_lo > 0.0:
        failures.append("nonpositive first tangent covariance floor")

    F, _Q, Rstep = FULL._transition_and_Q(src, domain)
    alpha_hi = float(F[15][15].hi)
    off = max(Rstep[0][1].abs_upper(), Rstep[0][2].abs_upper())

    chord0 = _chord_from_cos_lower(float(first["post_prediction_true_gravity_cosine_lower"]))
    pred_angle = float(first["first_prediction_transport_angle_upper_rad"])
    pred_chord = _correction_chord_upper(pred_angle)
    first_perp_aw_component = FULL.up(FULL.up(g * chord0) + ba)

    rows = []
    bad = None
    max_rho = 0.0
    max_corr = 0.0
    max_chord = 0.0
    max_eaw_x = 0.0
    min_headroom = math.inf

    for base in core["rows"]:
        p = base["P_aw_variance"]
        d = base["first_correction_rad"]
        az = base["first_axial_aw_correction_mps2"]
        m = base["sample1_aligned_signed_force_mps2"]
        k = float(base["Ktheta_abs_upper"])

        d_hi = float(d[1])
        chord = min(2.0, FULL.up(chord0 + FULL.up(_correction_chord_upper(d_hi) + pred_chord)))
        m_abs = _interval_abs_upper(m)
        rotational = FULL.up(m_abs * chord)

        # In the ideal first tangent block K_aw/K_theta=p/(g*t).  Canonical
        # d=e1 therefore changes only aw_y in the tangent plane; aw_x stays
        # constrained by the zero perpendicular first residual above.
        beta_hi = FULL.up(float(p[1]) / FULL.down(g * t_lo))
        tangent_aw_update = FULL.up(beta_hi * d_hi)
        axial_aw_update = _interval_abs_upper(az)
        update_other = FULL.up(math.sqrt(FULL.up(tangent_aw_update*tangent_aw_update + axial_aw_update*axial_aw_update)))
        other_error_norm = FULL.up(aw_pred + update_other)
        mixed = FULL.up(FULL.up(math.sqrt(2.0) * off) * other_error_norm)
        eaw_x = FULL.up(alpha_hi * FULL.up(first_perp_aw_component + mixed))

        rho = FULL.up(rotational + FULL.up(eaw_x + ba))
        corr = FULL.up(k * rho)
        headroom = FULL.down(RANGE - corr)
        closed = math.isfinite(corr) and corr < RANGE
        max_rho = max(max_rho, rho)
        max_corr = max(max_corr, corr)
        max_chord = max(max_chord, chord)
        max_eaw_x = max(max_eaw_x, eaw_x)
        min_headroom = min(min_headroom, headroom)
        row = {
            "p_cell": base["p_cell"],
            "d_cell": base["d_cell"],
            "axial_cell": base["axial_cell"],
            "P_aw_variance": p,
            "first_correction_rad": d,
            "first_axial_aw_correction_mps2": az,
            "sample1_aligned_signed_force_mps2": m,
            "Ktheta_abs_upper": k,
            "pre_first_gravity_chord_upper": chord0,
            "post_first_plus_prediction_gravity_chord_upper": chord,
            "first_perpendicular_aw_component_abs_upper_mps2": first_perp_aw_component,
            "body_rotation_offdiag_abs_upper": off,
            "post_prediction_perpendicular_aw_error_abs_upper_mps2": eaw_x,
            "sample1_rotational_residual_abs_upper_mps2": rotational,
            "sample1_source_correlated_perp_residual_abs_upper_mps2": rho,
            "correction_norm_upper_rad": corr,
            "range_headroom_rad_lower": headroom,
            "inside_9rad_range": closed,
        }
        rows.append(row)
        if not closed and bad is None:
            bad = row

    ok = bool(rows) and bad is None and not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_RESET_PERP_SCALAR_CHANNEL_SOURCE_STATE_RESIDUAL_V3",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "positive_determinant_covariance_gain_from_v2_used": True,
        "canonical_first_perpendicular_residual_zero_used": True,
        "first_exact_gravity_chord_residual_used": True,
        "first_perpendicular_aw_error_inferred_from_zero_residual": True,
        "reset_gauge_preserves_perpendicular_aw_x_component": True,
        "one_step_OU_decay_included": True,
        "one_step_body_rotation_aw_mixing_included": True,
        "sample1_force_cell_used_in_rotational_residual": True,
        "global_raw_30p5_residual_multiplier_used": False,
        "temporal_force_slew_assumed": False,
        "first_attitude_PSD_cross_axis_remainder_included": False,
        "sample1_tangent_force_perturbation_included": False,
        "sample1_S_attitude_correction_included": False,
        "complete_sample1_branch_closed_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "validated_deployed_quaternion_range_rad": RANGE,
        "evaluated_joint_cells": len(rows),
        "pre_first_gravity_chord_upper": chord0,
        "one_step_prediction_gravity_chord_increment_upper": pred_chord,
        "body_rotation_offdiag_abs_upper": off,
        "max_post_first_plus_prediction_gravity_chord_upper": max_chord,
        "max_post_prediction_perpendicular_aw_error_abs_upper_mps2": max_eaw_x,
        "max_source_correlated_perp_residual_abs_upper_mps2": max_rho,
        "max_correction_norm_upper_rad": max_corr,
        "minimum_range_headroom_rad_lower": min_headroom,
        "first_unclosed_joint_cell": bad,
        "P5_SAMPLE1_RESET_PERP_SCALAR_CHANNEL_V3": "PASS" if ok else "NOT_ESTABLISHED",
        "next_obligation": (
            "ADD_FIRST_PSD_CROSS_AXIS_TANGENT_FORCE_AND_TINY_S_PERTURBATIONS"
            if ok else
            "SUBDIVIDE_GRAVITY_CAYLEY_DIRECTION_AND_FIRST_RESIDUAL_VECTOR"
        ),
        "failures": failures,
        "rows": rows,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    for k in (
        "source_generated_not_trajectory_fit",
        "positive_determinant_covariance_gain_from_v2_used",
        "canonical_first_perpendicular_residual_zero_used",
        "first_exact_gravity_chord_residual_used",
        "first_perpendicular_aw_error_inferred_from_zero_residual",
        "reset_gauge_preserves_perpendicular_aw_x_component",
        "one_step_OU_decay_included",
        "one_step_body_rotation_aw_mixing_included",
        "sample1_force_cell_used_in_rotational_residual",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "global_raw_30p5_residual_multiplier_used",
        "temporal_force_slew_assumed", "first_attitude_PSD_cross_axis_remainder_included",
        "sample1_tangent_force_perturbation_included", "sample1_S_attitude_correction_included",
        "complete_sample1_branch_closed_here", "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if int(d.get("evaluated_joint_cells", 0)) <= 0:
        f.append("no cells")
    for k in ("max_source_correlated_perp_residual_abs_upper_mps2", "max_correction_norm_upper_rad"):
        if not math.isfinite(float(d.get(k, math.nan))):
            f.append(f"nonfinite {k}")
    st = d.get("P5_SAMPLE1_RESET_PERP_SCALAR_CHANNEL_V3")
    w = d.get("first_unclosed_joint_cell")
    if st == "PASS" and w is not None:
        f.append("PASS retains witness")
    if st == "NOT_ESTABLISHED" and w is None:
        f.append("missing witness")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=32)
    ap.add_argument("--d-pieces", type=int, default=32)
    ap.add_argument("--axial-pieces", type=int, default=32)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain, source_pieces=a.source_pieces, source_cell_index=a.source_cell_index,
              p_pieces=a.p_pieces, d_pieces=a.d_pieces, axial_pieces=a.axial_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_RESET_PERP_SCALAR_CHANNEL_V3"],
        "cells": d["evaluated_joint_cells"],
        "chord0": d["pre_first_gravity_chord_upper"],
        "chord_max": d["max_post_first_plus_prediction_gravity_chord_upper"],
        "max_eaw_x": d["max_post_prediction_perpendicular_aw_error_abs_upper_mps2"],
        "max_rho": d["max_source_correlated_perp_residual_abs_upper_mps2"],
        "max_d": d["max_correction_norm_upper_rad"],
        "headroom": d["minimum_range_headroom_rad_lower"],
        "first_unclosed": d["first_unclosed_joint_cell"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
