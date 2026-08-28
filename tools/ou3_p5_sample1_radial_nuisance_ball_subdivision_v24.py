#!/usr/bin/env python3
"""V24: preserve sample-1 nuisance and PSD/S correction remainders as radial balls.

V23 keeps the current Cayley subbox tied to its exact nonlinear residual and
signed correction, but its component correction enclosure must still widen a
Euclidean nuisance/remainder ball independently in x, y, and z.  That box is
needed for signed component intersection, yet its Euclidean norm is not the
sharp radial certificate.

For the exact V10 one-plus-two gain,

    ||K_theta|| = max(||[K_xy,K_xz]||, ||[g_y,g_z]||).

The physical rotated-a_w plus accelerometer-bias residual has certified norm
rho_n, and V12D supplies an independent correction perturbation norm eta.  Thus
for the exact-rotation nominal correction d_nom(c),

    ||d|| <= ||d_nom(c)|| + ||K_theta|| rho_n + eta.

V24 retains V23's componentwise box for sign/feasibility and intersects this
independent radial-ball upper with the V13E radial parent before the same V16,
V15, and V18 q8 checks.  No source domain, estimator, six-radian shipping limit,
or theorem-promotion state is changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_current_exact_residual_subdivision_v23 as V23
import ou3_p5_sample1_exact_nonlinear_residual_v22 as V22
import ou3_p5_sample1_signed_cayley_q8_v14 as V14
import ou3_p5_sample1_signed_cayley_q8_v14d as V14D
import ou3_p5_sample1_signed_cayley_q8_v15 as V15
import ou3_p5_sample1_signed_cayley_q8_v16 as V16
import ou3_p5_sample1_signed_cayley_q8_v18 as V18

DEFAULT_DOMAIN = V23.DEFAULT_DOMAIN
SCHEMA = 2400
FULL = V14.FULL
Q_TARGET = V14.Q_TARGET


def _gain_operator_norm(parent: dict) -> float:
    gain = parent["gain_detail"]
    gy, gz = [V22._I(x) for x in gain["perpendicular_gain_components"]]
    kxy, kxz = [V22._I(x) for x in gain["parallel_gain_components"]]
    return max(V22._norm2_upper(gy, gz), V22._norm2_upper(kxy, kxz))


def _nominal_exact_correction(c, parent: dict):
    fyz = parent["sample1_force_components_yz_mps2"]
    force = [FULL.I(0.0), V22._I(fyz[0]), V22._I(fyz[1])]
    yR = V22.exact_rotation_residual(c, force)
    gain = parent["gain_detail"]
    gy, gz = [V22._I(x) for x in gain["perpendicular_gain_components"]]
    kxy, kxz = [V22._I(x) for x in gain["parallel_gain_components"]]
    return [
        kxy * yR[1] + kxz * yR[2],
        gy * yR[0],
        gz * yR[0],
    ], yR


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6,
          current_component_pieces: int = 4) -> dict:
    parent = V23.build(
        Path(domain_path).resolve(), source_pieces=source_pieces,
        source_cell_index=source_cell_index, p_pieces=p_pieces,
        tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
        residual_x_pieces=residual_x_pieces, parallel_pieces=parallel_pieces,
        current_component_pieces=current_component_pieces)
    failures = [f"V23: {x}" for x in V23.validate(parent)]
    if parent.get("P5_SAMPLE1_CURRENT_EXACT_RESIDUAL_SUBDIVISION_V23") != "PASS":
        failures.append("V23 current/exact-residual prerequisite did not pass")

    q_parent = float(parent["sample1_current_cayley_norm_upper"])
    c_parent = [V22._I(x) for x in parent["sample1_current_component_box"]]
    baseline_box = [V22._I(x) for x in parent["baseline_correction_box_rad"]]
    baseline_lo = float(parent["baseline_correction_radial_lower_rad"])
    baseline_hi = float(parent["baseline_correction_radial_upper_rad"])
    nuisance = float(parent["combined_rotated_aw_plus_bias_nuisance_norm_upper_mps2"])
    corr_perturb = float(parent["V12D_correction_perturbation_norm_upper_rad"])
    kn = _gain_operator_norm(parent)
    nuisance_correction_ball = FULL.up(kn * nuisance)
    total_correction_ball = FULL.up(nuisance_correction_ball + corr_perturb)
    subboxes = V23._current_subboxes(c_parent, current_component_pieces)

    q_rejected = component_rejected = radial_rejected = closed = open_count = 0
    radial_strict = 0
    first_open = worst_open = first_closed = None
    min_radial_hi = math.inf; max_radial_hi = 0.0
    min_qbest = math.inf; max_qbest = 0.0

    for idx, raw_c in enumerate(subboxes):
        c = V23._clip_box_to_q_ball(raw_c, q_parent)
        if c is None:
            q_rejected += 1
            continue
        q_sub = min(q_parent, V14.CAYLEY1._norm_upper(c))
        source_box, _yr_parent, _n = V23._source_correction_box(c, parent)
        joint = V22._intersect_boxes(baseline_box, source_box)
        if joint is None:
            component_rejected += 1
            continue

        nominal, yR = _nominal_exact_correction(c, parent)
        nominal_hi = V14.CAYLEY1._norm_upper(nominal)
        nominal_lo = V14.CAYLEY2._norm_lower(nominal)
        source_radial_hi = FULL.up(nominal_hi + total_correction_ball)
        source_radial_lo = max(0.0, FULL.down(nominal_lo - total_correction_ball))
        box_hi = V14.CAYLEY1._norm_upper(joint)
        box_lo = V14.CAYLEY2._norm_lower(joint)
        parent_radial_hi = min(baseline_hi, box_hi)
        radial_hi = min(parent_radial_hi, source_radial_hi)
        radial_lo = max(baseline_lo, box_lo, source_radial_lo)
        if radial_lo > radial_hi:
            radial_rejected += 1
            continue
        radial_strict += int(radial_hi < parent_radial_hi)
        min_radial_hi = min(min_radial_hi, radial_hi)
        max_radial_hi = max(max_radial_hi, radial_hi)

        geo = V15._geodesic_q_and_scalar_lower(q_sub, radial_lo, radial_hi)
        geo_q = math.inf if geo is None else float(geo[0])
        wd, vd, branches, narrowed = V16.axis_cone_normalized_shipping_quaternion(
            joint, radial_lower=radial_lo, radial_upper=radial_hi,
            parent=V14D.radial_sinc_normalized_shipping_quaternion)
        cx_min = V14._minimum_abs(c[0])
        yz2 = max(0.0, FULL.up(q_sub * q_sub) - FULL.down(cx_min * cx_min))
        cyz = min(FULL.up(math.sqrt(yz2)), V18._yz_norm_upper(c[1], c[2]))
        chart = {"cx": c[0], "cy": c[1], "cz": c[2], "cyz_norm_upper": cyz}
        parent_W = FULL.I(2.0) * wd - V14.CAYLEY1._dot(vd, c)
        W, _yb, _yj = V18._support_product_scalar(parent_W, wd, vd, chart)
        product_w, product_q = V14._qplus_from_product_scalar(q_sub, W)
        qbest = min(geo_q, product_q)
        is_closed = ((math.isfinite(geo_q) and geo_q < Q_TARGET)
                     or (math.isfinite(product_q) and product_q < Q_TARGET
                         and product_w > 0.0))
        row = {
            "subbox_index": idx,
            "q_ball_projected_current_component_box": [x.as_list() for x in c],
            "current_q_upper": q_sub,
            "exact_rotation_residual_box_mps2": [x.as_list() for x in yR],
            "nominal_exact_correction_box_rad": [x.as_list() for x in nominal],
            "component_joint_correction_box_rad": [x.as_list() for x in joint],
            "gain_operator_norm_upper": kn,
            "physical_nuisance_correction_ball_upper_rad": nuisance_correction_ball,
            "V12D_correction_perturbation_ball_upper_rad": corr_perturb,
            "total_additive_correction_ball_upper_rad": total_correction_ball,
            "source_radial_upper_rad": source_radial_hi,
            "correction_radial_lower_rad": radial_lo,
            "correction_radial_upper_rad": radial_hi,
            "axis_cone_narrowed": narrowed,
            "quaternion_branches": branches,
            "geodesic_q_upper": geo_q,
            "product_abs_W_lower": product_w,
            "product_q_upper": product_q,
            "best_q_upper": qbest,
            "closed_inside_q8": is_closed,
        }
        if math.isfinite(qbest):
            min_qbest = min(min_qbest, qbest)
            if math.isfinite(max_qbest):
                max_qbest = max(max_qbest, qbest)
        else:
            max_qbest = math.inf
        if is_closed:
            closed += 1
            if first_closed is None: first_closed = row
        else:
            open_count += 1
            if first_open is None: first_open = row
            if worst_open is None or qbest > worst_open["best_q_upper"]:
                worst_open = row

    candidate = len(subboxes)
    compatible = candidate - q_rejected - component_rejected - radial_rejected
    accounted = candidate == q_rejected + component_rejected + radial_rejected + closed + open_count
    focused_closed = bool(not failures and accounted and open_count == 0)
    return {
        **parent,
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_RADIAL_NUISANCE_BALL_SUBDIVISION_V24",
        "V23_current_exact_residual_subdivision_parent_retained": True,
        "structured_gain_operator_norm_used_for_physical_nuisance_ball": True,
        "V12D_correction_perturbation_retained_as_radial_ball": True,
        "componentwise_correction_box_retained_for_sign_and_feasibility": True,
        "independent_radial_ball_intersected_before_q8_test": True,
        "gain_operator_norm_upper": kn,
        "physical_nuisance_correction_ball_upper_rad": nuisance_correction_ball,
        "V12D_correction_perturbation_ball_upper_rad": corr_perturb,
        "total_additive_correction_ball_upper_rad": total_correction_ball,
        "candidate_current_subboxes": candidate,
        "q_ball_rejected_current_subboxes": q_rejected,
        "component_incompatible_current_subboxes": component_rejected,
        "radial_incompatible_current_subboxes": radial_rejected,
        "compatible_current_subboxes": compatible,
        "radial_ball_strictly_refined_current_subboxes": radial_strict,
        "closed_current_subboxes": closed,
        "open_current_subboxes": open_count,
        "all_candidate_current_subboxes_accounted": accounted,
        "minimum_compatible_radial_upper_rad": None if math.isinf(min_radial_hi) else min_radial_hi,
        "maximum_compatible_radial_upper_rad": max_radial_hi,
        "minimum_best_q_upper": None if math.isinf(min_qbest) else min_qbest,
        "maximum_best_q_upper": max_qbest,
        "first_closed_current_subbox": first_closed,
        "first_open_current_subbox": first_open,
        "worst_open_current_subbox": worst_open,
        "focused_first_witness_signed_subcell_closed_by_radial_ball": focused_closed,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_RADIAL_NUISANCE_BALL_SUBDIVISION_V24": "PASS" if not failures else "NOT_ESTABLISHED",
        "next_obligation": (
            "LIFT_V24_RADIAL_NUISANCE_BALL_ROUTE_OVER_FIRST_BASE_ROW"
            if focused_closed else
            "INCREASE_CURRENT_SUBDIVISION_OR_DERIVE_SOURCE_CORRELATED_AW_COMPONENTS_AT_FIRST_OPEN_SUBBOX"),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA: f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_RADIAL_NUISANCE_BALL_SUBDIVISION_V24":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V23_current_exact_residual_subdivision_parent_retained",
        "structured_gain_operator_norm_used_for_physical_nuisance_ball",
        "V12D_correction_perturbation_retained_as_radial_ball",
        "componentwise_correction_box_retained_for_sign_and_feasibility",
        "independent_radial_ball_intersected_before_q8_test",
        "all_candidate_current_subboxes_accounted",
    ):
        if d.get(k) is not True: f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed", "deployed_correction_limit_increased",
              "q8_composed_here", "q8_word_promoted_here", "whole_word_promoted_here",
              "N_H_words_set_here"):
        if d.get(k) is not False: f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET: f.append("q target changed")
    candidate = int(d.get("candidate_current_subboxes", -1))
    qrej = int(d.get("q_ball_rejected_current_subboxes", -1))
    crej = int(d.get("component_incompatible_current_subboxes", -1))
    rrej = int(d.get("radial_incompatible_current_subboxes", -1))
    comp = int(d.get("compatible_current_subboxes", -1))
    closed = int(d.get("closed_current_subboxes", -1)); opened = int(d.get("open_current_subboxes", -1))
    if candidate < 1 or candidate != qrej + crej + rrej + closed + opened:
        f.append("invalid V24 current-subbox accounting")
    if comp != closed + opened: f.append("invalid V24 compatible-subbox accounting")
    if d.get("focused_first_witness_signed_subcell_closed_by_radial_ball") is True and opened != 0:
        f.append("V24 claims focused closure with open subboxes")
    if d.get("P5_SAMPLE1_RADIAL_NUISANCE_BALL_SUBDIVISION_V24") not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid V24 status")
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
    d = build(x.domain, source_pieces=x.source_pieces, source_cell_index=x.source_cell_index,
              p_pieces=x.p_pieces, tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces,
              residual_x_pieces=x.residual_x_pieces, parallel_pieces=x.parallel_pieces,
              current_component_pieces=x.current_component_pieces)
    vf = validate(d); d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_RADIAL_NUISANCE_BALL_SUBDIVISION_V24"],
        "gain_norm": d["gain_operator_norm_upper"],
        "nuisance_correction_ball": d["physical_nuisance_correction_ball_upper_rad"],
        "V12D_correction_ball": d["V12D_correction_perturbation_ball_upper_rad"],
        "radial_refined": d["radial_ball_strictly_refined_current_subboxes"],
        "candidate": d["candidate_current_subboxes"],
        "compatible": d["compatible_current_subboxes"],
        "closed": d["closed_current_subboxes"],
        "open": d["open_current_subboxes"],
        "min_radial_upper": d["minimum_compatible_radial_upper_rad"],
        "max_radial_upper": d["maximum_compatible_radial_upper_rad"],
        "min_best_q": d["minimum_best_q_upper"],
        "max_best_q": d["maximum_best_q_upper"],
        "focused_closed": d["focused_first_witness_signed_subcell_closed_by_radial_ball"],
        "first_open": d["first_open_current_subbox"],
        "worst_open": d["worst_open_current_subbox"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__": raise SystemExit(main())
