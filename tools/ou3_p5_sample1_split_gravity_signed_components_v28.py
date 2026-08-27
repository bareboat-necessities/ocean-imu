#!/usr/bin/env python3
"""V28: split tangent/axial gravity-decay remainders in V27's signed residual.

V27 is the first focused route that materially tightens the signed product at
the remaining sample-1 q8 witness: its product q upper falls from 9.6827 to
8.6464.  It intentionally used V10's perpendicular gravity-decay bound as one
common remainder for the tangent and axial residual components.  That is
rigorous on the certified startup chart but discards exact gravity geometry.

The source-audited first-accelerometer certificate supplies

    c_g <= cos(theta)

for the true post-prediction gravity misalignment.  Therefore independently

    |y_R,t| <= g sqrt(1-c_g^2),
    |y_R,z| <= g (1-c_g),

and under one homogeneous OU step the corresponding uncancelled pieces are
multiplied by at most 1-alpha_lower.  V28 replaces only V27's common gravity
remainder by these two source-certified component bounds.  Bias-difference,
transport/series mismatch, the exact first residual cell, V21 signed gains,
V12D PSD/S correction ball, V23 correction parent, and V16/V15/V18 q8 checks
are unchanged.

This remains a focused source-cell-0 diagnostic.  No estimator, source domain,
source branch, 6-rad shipping limit, q target, theorem promotion, or N_H_words
state is changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_signed_post_first_aw_components_v27 as V27

DEFAULT_DOMAIN = V27.DEFAULT_DOMAIN
SCHEMA = 2800
FULL = V27.FULL
Q_TARGET = V27.Q_TARGET


def _gravity_component_decay_bounds(*, cosine_lower: float,
                                    alpha_lower: float,
                                    gravity: float) -> dict:
    c = float(cosine_lower); a = float(alpha_lower); g = float(gravity)
    if not (math.isfinite(c) and -1.0 <= c <= 1.0):
        raise ValueError("gravity cosine lower must lie in [-1,1]")
    if not (math.isfinite(a) and 0.0 < a <= 1.0):
        raise ValueError("OU decay lower must lie in (0,1]")
    if not (math.isfinite(g) and g > 0.0):
        raise ValueError("positive finite gravity required")
    sin_hi = 1.0 if c < 0.0 else FULL.up(math.sqrt(max(
        0.0, FULL.up(1.0 - FULL.down(c * c)))))
    tangent = FULL.up(g * sin_hi)
    axial = FULL.up(g * max(0.0, FULL.up(1.0 - c)))
    decay = FULL.up(1.0 - a)
    return {
        "post_prediction_true_gravity_cosine_lower": c,
        "gravity_tangent_residual_upper_mps2": tangent,
        "gravity_axial_residual_upper_mps2": axial,
        "ou_uncancelled_factor_upper": decay,
        "tangent_gravity_decay_remainder_upper_mps2": FULL.up(decay * tangent),
        "axial_gravity_decay_remainder_upper_mps2": FULL.up(decay * axial),
    }


def _split_signed_residual_components(*, row: dict, parent: dict,
                                      alpha: Interval, gravity: float,
                                      gravity_detail: dict):
    rt = Interval.outward_bounds(*map(float, row["first_tangent_residual_magnitude_mps2"]))
    rz = Interval.outward_bounds(*map(float, row["first_axial_residual_mps2"]))
    fyz = parent.get("sample1_force_components_yz_mps2", [])
    if len(fyz) != 2:
        raise RuntimeError("V21 sample-1 force components missing")
    fy, fz = V27.V22._I(fyz[0]), V27.V22._I(fyz[1])
    rho = float(row["sample1_full_residual_norm_upper_mps2"])
    rho_x = min(rho, float(row["sample1_combined_source_x_residual_upper_mps2"]))
    bias = float(row["bias_difference_upper_mps2"])
    geom = float(row["rotation_mismatch_residual_upper_mps2"])
    decay_t = float(gravity_detail["tangent_gravity_decay_remainder_upper_mps2"])
    decay_z = float(gravity_detail["axial_gravity_decay_remainder_upper_mps2"])
    vals = (rho, rho_x, bias, geom, decay_t, decay_z)
    if not all(math.isfinite(x) and x >= 0.0 for x in vals):
        raise RuntimeError("invalid split signed-residual bounds")
    rem_t_hi = FULL.up(decay_t + FULL.up(bias + geom))
    rem_z_hi = FULL.up(decay_z + FULL.up(bias + geom))
    rem_t = Interval.outward_bounds(-rem_t_hi, rem_t_hi)
    rem_z = Interval.outward_bounds(-rem_z_hi, rem_z_hi)
    rx = Interval.outward_bounds(-rho_x, rho_x)
    ry = alpha * rt + fy + rem_t
    rz1 = alpha * rz - (fz - FULL.I(float(gravity))) + rem_z
    return [rx, ry, rz1], {
        "first_tangent_residual_mps2": rt.as_list(),
        "first_axial_residual_mps2": rz.as_list(),
        "sample1_force_y_mps2": fy.as_list(),
        "sample1_force_z_mps2": fz.as_list(),
        "combined_x_residual_abs_upper_mps2": rho_x,
        "bias_difference_remainder_upper_mps2": bias,
        "transport_series_remainder_upper_mps2": geom,
        **gravity_detail,
        "tangent_component_remainder_upper_mps2": rem_t_hi,
        "axial_component_remainder_upper_mps2": rem_z_hi,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6,
          current_component_pieces: int = 4) -> dict:
    path = Path(domain_path).resolve()
    dom = json.loads(path.read_text(encoding="utf-8"))
    V10 = V27.V23.V22.V21B.V21.V12D.V11.V10
    FIRST = V10.FIRST
    first = FIRST.build(path, source_pieces=source_pieces)
    failures = [f"first: {x}" for x in FIRST.validate(first)]
    try:
        src, phase = V10.RG._source_phase_children(source_pieces)[source_cell_index]
        if phase != "due":
            raise RuntimeError("V28 focused witness requires first due source cell")
        F, _Q, _ = V10.FULL._transition_and_Q(src, dom)
        alpha = F[15][15]
        gravity_detail = _gravity_component_decay_bounds(
            cosine_lower=float(first["post_prediction_true_gravity_cosine_lower"]),
            alpha_lower=float(alpha.lo),
            gravity=float(dom["startup"]["gravity_mps2"]))
    except Exception as exc:
        failures.append(f"gravity component bounds: {exc}")
        gravity_detail = {
            "tangent_gravity_decay_remainder_upper_mps2": math.inf,
            "axial_gravity_decay_remainder_upper_mps2": math.inf,
        }

    original = V27._signed_residual_components
    def refined(*, row: dict, parent: dict, alpha: Interval, gravity: float):
        return _split_signed_residual_components(
            row=row, parent=parent, alpha=alpha, gravity=gravity,
            gravity_detail=gravity_detail)
    V27._signed_residual_components = refined
    try:
        parent = V27.build(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces,
            current_component_pieces=current_component_pieces)
    finally:
        V27._signed_residual_components = original

    failures += [f"V27: {x}" for x in V27.validate(parent)]
    if parent.get("P5_SAMPLE1_SIGNED_POST_FIRST_AW_COMPONENTS_V27") != "PASS":
        failures.append("V27 signed-component prerequisite did not pass")

    out = dict(parent)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_SPLIT_GRAVITY_SIGNED_COMPONENTS_V28",
        "V27_signed_post_first_aw_parent_retained": True,
        "first_accel_source_certificate_revalidated": True,
        "post_prediction_gravity_cosine_source_bound_used": True,
        "tangent_and_axial_gravity_decay_split": True,
        "shared_tangent_gravity_decay_for_axial_channel_retired": True,
        "split_gravity_detail": gravity_detail,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_SPLIT_GRAVITY_SIGNED_COMPONENTS_V28": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_V28_SPLIT_GRAVITY_SIGNED_COMPONENTS_OVER_ALL_V23_CURRENT_SUBBOXES"
            if parent.get("first_open_subbox_closed_inside_q8") is True and not failures
            else "DERIVE_DIRECTIONAL_V12D_PSD_S_CORRECTION_REMAINDER_AT_FIRST_OPEN_SUBBOX"),
        "failures": list(dict.fromkeys(failures)),
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_SPLIT_GRAVITY_SIGNED_COMPONENTS_V28":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V27_signed_post_first_aw_parent_retained",
        "first_accel_source_certificate_revalidated",
        "post_prediction_gravity_cosine_source_bound_used",
        "tangent_and_axial_gravity_decay_split",
        "shared_tangent_gravity_decay_for_axial_channel_retired",
        "V23_first_open_subbox_retained",
        "V10_exact_first_update_OU_cancellation_used",
        "signed_tangent_axial_first_residual_cell_retained",
        "V21_signed_one_plus_two_gain_components_used",
        "V12D_correction_perturbation_retained_as_single_ball",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "q8_composed_here", "q8_word_promoted_here", "whole_word_promoted_here",
        "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    gd = d.get("split_gravity_detail", {})
    for k in ("tangent_gravity_decay_remainder_upper_mps2",
              "axial_gravity_decay_remainder_upper_mps2"):
        x = gd.get(k)
        if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or float(x) < 0.0:
            f.append(f"invalid {k}")
    if (isinstance(gd.get("post_prediction_true_gravity_cosine_lower"), (int, float))
            and float(gd["post_prediction_true_gravity_cosine_lower"]) >= 0.0
            and float(gd.get("axial_gravity_decay_remainder_upper_mps2", math.inf))
            > float(gd.get("tangent_gravity_decay_remainder_upper_mps2", -math.inf))):
        f.append("nonnegative-cosine axial gravity remainder exceeded tangent remainder")
    if d.get("P5_SAMPLE1_SPLIT_GRAVITY_SIGNED_COMPONENTS_V28") not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid V28 status")
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
    vf = validate(d); d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_SPLIT_GRAVITY_SIGNED_COMPONENTS_V28"],
        "gravity": d.get("split_gravity_detail"),
        "residual_box": d.get("sample1_signed_residual_box_mps2"),
        "nominal_correction_box": d.get("nominal_signed_correction_box_rad"),
        "source_radial_upper": d.get("source_correlated_radial_upper_rad"),
        "q_current": d.get("V23_first_open_current_q_upper"),
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
