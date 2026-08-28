#!/usr/bin/env python3
"""V31: lift the V30 theta-x gain refinement over every V23 current subbox.

V30 closes V23's first open current-Cayley subbox by replacing only the
attitude-x row of V12D's gain-perturbation ball with the rowwise resolvent

    ||Delta K_x|| <= (dC_theta + k_parallel dS) ||S'^{-1}||.

The signed post-first residual enclosure from V27/V28 is source generated and
does not depend on V23's artificial Cartesian current subdivision.  Therefore
it can be intersected with the current-dependent exact V23 residual enclosure
on every subbox.  V31 performs that lift over all ``pieces^3`` current boxes,
retaining

  * V23's q-ball projection and exact nonlinear current-dependent residual;
  * V28's split tangent/axial gravity source enclosure;
  * V29's Euclidean yz and total radial perturbation parents;
  * V30's row-specific theta-x gain perturbation;
  * the V13E signed correction parent, V16 axis cone, V15 geodesic route, and
    V18 signed-product support.

Empty intersections discard only the corresponding artificial current subbox.
No source assumption, estimator setting, six-radian correction limit, q<8
target, or theorem-promotion state is changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_theta_x_gain_perturbation_v30 as V30

DEFAULT_DOMAIN = V30.DEFAULT_DOMAIN
SCHEMA = 3100
V29 = V30.V29
V28 = V29.V28
V27 = V28.V27
V23 = V27.V23
V22 = V23.V22
V14 = V27.V14
V14D = V27.V14D
V15 = V27.V15
V16 = V27.V16
V18 = V27.V18
FULL = V30.FULL
Q_TARGET = V30.Q_TARGET
WITNESS = V30.WITNESS


def _I(x):
    return Interval.outward_bounds(*map(float, x))


def _refined_caps(*, base: dict, vr: dict, row_detail: dict) -> dict:
    kperp = float(base["Ktheta_perpendicular_block_upper"])
    kpar = float(base["Ktheta_parallel_block_upper"])
    drho = float(vr["total_residual_perturbation_upper_mps2"])
    dk = float(vr["sample1_attitude_gain_operator_perturbation_upper"])
    rho = float(base["sample1_full_residual_norm_upper_mps2"])
    parent = V29._directional_perturbation_caps(
        k_perp=kperp, k_parallel=kpar, drho=drho, dk=dk, rho=rho)
    dkx = float(row_detail["theta_x_gain_perturbation_operator_upper"])
    rho_plus = FULL.up(rho + drho)
    gain_x = FULL.up(dkx * rho_plus)
    ex = FULL.up(FULL.up(kpar * drho) + gain_x)
    if ex > FULL.up(float(parent["x_correction_perturbation_abs_upper_rad"])):
        raise RuntimeError("V31 theta-x cap exceeded V29 parent")
    out = dict(parent)
    out.update({
        "theta_x_gain_perturbation_ball_upper_rad": gain_x,
        "V29_full_gain_perturbation_ball_upper_rad": parent["gain_perturbation_ball_upper_rad"],
        "x_correction_perturbation_abs_upper_rad": ex,
    })
    return out


def _source_directional_enclosure(*, path: Path, dom: dict, p22: dict,
                                  base: dict, vr: dict, first: dict,
                                  src: dict) -> dict:
    F, _Q, _ = V23.V22.V21B.V21.V12D.V11.V10.FULL._transition_and_Q(src, dom)
    alpha = F[15][15]
    gravity = float(dom["startup"]["gravity_mps2"])
    gravity_detail = V28._gravity_component_decay_bounds(
        cosine_lower=float(first["post_prediction_true_gravity_cosine_lower"]),
        alpha_lower=float(alpha.lo), gravity=gravity)
    residual, residual_detail = V28._split_signed_residual_components(
        row=base, parent=p22, alpha=alpha, gravity=gravity,
        gravity_detail=gravity_detail)
    nominal = V27._nominal_correction(residual, p22)
    row_detail = V30._theta_x_gain_perturbation_upper(vr, base)
    caps = _refined_caps(base=base, vr=vr, row_detail=row_detail)

    ex = float(caps["x_correction_perturbation_abs_upper_rad"])
    eyz = float(caps["yz_correction_perturbation_norm_upper_rad"])
    eall = float(caps["total_correction_perturbation_norm_upper_rad"])
    source_box = [
        nominal[0] + Interval.outward_bounds(-ex, ex),
        nominal[1] + Interval.outward_bounds(-eyz, eyz),
        nominal[2] + Interval.outward_bounds(-eyz, eyz),
    ]
    nominal_yz = V18._yz_norm_upper(nominal[1], nominal[2])
    yz_source_hi = FULL.up(nominal_yz + eyz)
    nominal_hi = V14.CAYLEY1._norm_upper(nominal)
    nominal_lo = V14.CAYLEY2._norm_lower(nominal)
    source_radial_hi = min(
        FULL.up(nominal_hi + eall),
        FULL.up(float(base["combined_directional_correction_norm_upper_rad"]) + eall),
    )
    source_radial_lo = max(0.0, FULL.down(nominal_lo - eall))
    return {
        "residual": residual,
        "residual_detail": residual_detail,
        "nominal": nominal,
        "row_detail": row_detail,
        "caps": caps,
        "source_box": source_box,
        "yz_source_hi": yz_source_hi,
        "source_radial_hi": source_radial_hi,
        "source_radial_lo": source_radial_lo,
        "gravity_detail": gravity_detail,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6,
          current_component_pieces: int = 4) -> dict:
    path = Path(domain_path).resolve()
    dom = json.loads(path.read_text(encoding="utf-8"))
    V12D = V23.V22.V21B.V21.V12D
    V10 = V12D.V11.V10
    FIRST = V10.FIRST

    p22 = V22.build(
        path, source_pieces=source_pieces, source_cell_index=source_cell_index,
        p_pieces=p_pieces, tangent_pieces=tangent_pieces,
        axial_pieces=axial_pieces, residual_x_pieces=residual_x_pieces,
        parallel_pieces=parallel_pieces)
    v12 = V12D.build(path, source_pieces=source_pieces,
                     source_cell_index=source_cell_index, p_pieces=p_pieces,
                     tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    core = V10.build(path, source_pieces=source_pieces,
                     source_cell_index=source_cell_index, p_pieces=p_pieces,
                     tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    first = FIRST.build(path, source_pieces=source_pieces)

    failures = [f"V22: {x}" for x in V22.validate(p22)]
    failures += [f"V12D: {x}" for x in V12D.validate(v12)]
    failures += [f"V10: {x}" for x in V10.validate(core)]
    failures += [f"first: {x}" for x in FIRST.validate(first)]
    if p22.get("P5_SAMPLE1_EXACT_NONLINEAR_RESIDUAL_WITNESS_V22") != "PASS":
        failures.append("V22 exact nonlinear residual prerequisite did not pass")

    try:
        vr = V30._witness_row(v12)
        base = V30._witness_row(core)
        src, phase = V10.RG._source_phase_children(source_pieces)[source_cell_index]
        if phase != "due":
            raise RuntimeError("V31 focused lift requires first due source cell")
        source = _source_directional_enclosure(
            path=path, dom=dom, p22=p22, base=base, vr=vr,
            first=first, src=src)
        old_eta = float(p22["V12D_correction_perturbation_norm_upper_rad"])
        if float(source["caps"]["total_correction_perturbation_norm_upper_rad"]) > FULL.up(old_eta):
            raise RuntimeError("V31 directional source enclosure exceeded V12D radial parent")
    except Exception as exc:
        failures.append(f"V31 source directional enclosure: {exc}")
        source = None; old_eta = math.inf

    q_parent = float(p22.get("sample1_current_cayley_norm_upper", math.inf))
    c_parent = [_I(x) for x in p22.get("sample1_current_component_box", [])]
    baseline_box = [_I(x) for x in p22.get("baseline_correction_box_rad", [])]
    if len(c_parent) != 3 or len(baseline_box) != 3:
        failures.append("V31 parent current/correction box missing")
        subboxes = []
    else:
        subboxes = V23._current_subboxes(c_parent, current_component_pieces)
    baseline_lo = float(p22.get("baseline_correction_radial_lower_rad", 0.0))
    baseline_hi = float(p22.get("baseline_correction_radial_upper_rad", math.inf))

    q_rejected = current_source_rejected = directional_rejected = 0
    yz_rejected = radial_rejected = closed = open_count = evaluated = 0
    first_open = worst_open = first_closed = None
    min_best = math.inf; max_best = 0.0
    min_radial_hi = math.inf; max_radial_hi = 0.0

    if source is not None:
        nominal = source["nominal"]
        source_box = source["source_box"]
        yz_source_hi = float(source["yz_source_hi"])
        source_radial_hi = float(source["source_radial_hi"])
        source_radial_lo = float(source["source_radial_lo"])

        for idx, raw_c in enumerate(subboxes):
            c = V23._clip_box_to_q_ball(raw_c, q_parent)
            if c is None:
                q_rejected += 1
                continue
            evaluated += 1
            q_sub = min(q_parent, V14.CAYLEY1._norm_upper(c))
            v23_source, yR, nuisance = V23._source_correction_box(c, p22)
            joint0 = V29._intersect_boxes(baseline_box, v23_source)
            if joint0 is None:
                current_source_rejected += 1
                continue
            joint = V29._intersect_boxes(joint0, source_box)
            if joint is None:
                directional_rejected += 1
                continue
            yz = V29._clip_yz_to_radius(joint[1], joint[2], yz_source_hi)
            if yz is None:
                yz_rejected += 1
                continue
            joint[1], joint[2] = yz

            radial_hi = min(baseline_hi, source_radial_hi,
                            V14.CAYLEY1._norm_upper(joint))
            radial_lo = max(baseline_lo, source_radial_lo,
                            V14.CAYLEY2._norm_lower(joint))
            if radial_lo > radial_hi:
                radial_rejected += 1
                continue
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
                "current_source_nuisance_norm_upper_mps2": nuisance,
                "V23_current_source_correction_box_rad": [x.as_list() for x in v23_source],
                "V23_joint_correction_box_rad": [x.as_list() for x in joint0],
                "V31_joint_directional_correction_box_rad": [x.as_list() for x in joint],
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
                min_best = min(min_best, qbest)
                if math.isfinite(max_best):
                    max_best = max(max_best, qbest)
            else:
                max_best = math.inf
            if is_closed:
                closed += 1
                if first_closed is None:
                    first_closed = row
            else:
                open_count += 1
                if first_open is None:
                    first_open = row
                if worst_open is None or qbest > worst_open["best_q_upper"]:
                    worst_open = row

    compatible = closed + open_count
    accounted = (len(subboxes) == q_rejected + current_source_rejected
                 + directional_rejected + yz_rejected + radial_rejected
                 + closed + open_count)
    witness_closed = bool(subboxes and not failures and accounted and open_count == 0)
    status = "PASS" if not failures else "NOT_ESTABLISHED"
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_V30_CURRENT_SUBBOX_LIFT_V31",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V22_exact_current_residual_parent_retained": True,
        "V28_split_gravity_signed_source_enclosure_retained": True,
        "V29_yz_and_radial_perturbation_parents_retained": True,
        "V30_theta_x_row_resolvent_retained": True,
        "V23_current_partition_and_q_ball_projection_retained": True,
        "current_dependent_and_source_directional_correction_enclosures_intersected": True,
        "V16_axis_cone_V15_geodesic_V18_yz_support_retained": True,
        "witness_source_row": list(WITNESS),
        "theta_x_gain_perturbation_detail": None if source is None else source["row_detail"],
        "directional_perturbation_detail": None if source is None else source["caps"],
        "split_gravity_detail": None if source is None else source["gravity_detail"],
        "sample1_signed_residual_box_mps2": None if source is None else [x.as_list() for x in source["residual"]],
        "nominal_signed_correction_box_rad": None if source is None else [x.as_list() for x in source["nominal"]],
        "source_directional_correction_box_rad": None if source is None else [x.as_list() for x in source["source_box"]],
        "source_directional_yz_norm_upper_rad": None if source is None else source["yz_source_hi"],
        "source_directional_radial_lower_rad": None if source is None else source["source_radial_lo"],
        "source_directional_radial_upper_rad": None if source is None else source["source_radial_hi"],
        "previous_isotropic_V12D_correction_perturbation_upper_rad": old_eta,
        "current_component_pieces": int(current_component_pieces),
        "candidate_current_subboxes": len(subboxes),
        "q_ball_rejected_current_subboxes": q_rejected,
        "evaluated_q_compatible_current_subboxes": evaluated,
        "current_source_rejected_current_subboxes": current_source_rejected,
        "directional_source_rejected_current_subboxes": directional_rejected,
        "directional_yz_rejected_current_subboxes": yz_rejected,
        "radial_constraint_rejected_current_subboxes": radial_rejected,
        "compatible_current_subboxes": compatible,
        "closed_current_subboxes": closed,
        "open_current_subboxes": open_count,
        "all_candidate_current_subboxes_accounted": accounted,
        "minimum_compatible_radial_upper_rad": None if math.isinf(min_radial_hi) else min_radial_hi,
        "maximum_compatible_radial_upper_rad": max_radial_hi,
        "minimum_best_q_upper": None if math.isinf(min_best) else min_best,
        "maximum_best_q_upper": max_best,
        "first_closed_current_subbox": first_closed,
        "first_open_current_subbox": first_open,
        "worst_open_current_subbox": worst_open,
        "focused_first_witness_signed_subcell_closed_by_V30_lift": witness_closed,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_V30_CURRENT_SUBBOX_LIFT_V31": status,
        "next_obligation": (
            "LIFT_V31_THETA_X_GAIN_CURRENT_SUBBOX_CLOSURE_INTO_FULL_V18B_Q8_COVER"
            if witness_closed else
            "REFINE_FIRST_REMAINING_V31_CURRENT_SUBBOX_WITH_EXACT_THETA_X_DELTAC_ROW"),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_V30_CURRENT_SUBBOX_LIFT_V31":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V22_exact_current_residual_parent_retained",
        "V28_split_gravity_signed_source_enclosure_retained",
        "V29_yz_and_radial_perturbation_parents_retained",
        "V30_theta_x_row_resolvent_retained",
        "V23_current_partition_and_q_ball_projection_retained",
        "current_dependent_and_source_directional_correction_enclosures_intersected",
        "V16_axis_cone_V15_geodesic_V18_yz_support_retained",
        "all_candidate_current_subboxes_accounted",
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
    n = int(d.get("candidate_current_subboxes", -1))
    pieces = int(d.get("current_component_pieces", 0))
    if pieces < 2 or n != pieces ** 3:
        f.append("current subbox count mismatch")
    if int(d.get("closed_current_subboxes", -1)) < 0 or int(d.get("open_current_subboxes", -1)) < 0:
        f.append("invalid V31 closure counts")
    rd = d.get("theta_x_gain_perturbation_detail") or {}
    dx = float(rd.get("theta_x_gain_perturbation_operator_upper", math.inf))
    dp = float(rd.get("V12D_full_attitude_gain_perturbation_operator_upper", -math.inf))
    if not (math.isfinite(dx) and 0.0 <= dx <= FULL.up(dp)):
        f.append("invalid V30 theta-x gain refinement in lift")
    caps = d.get("directional_perturbation_detail") or {}
    ex = float(caps.get("x_correction_perturbation_abs_upper_rad", math.inf))
    eyz = float(caps.get("yz_correction_perturbation_norm_upper_rad", math.inf))
    eall = float(caps.get("total_correction_perturbation_norm_upper_rad", math.inf))
    old = float(d.get("previous_isotropic_V12D_correction_perturbation_upper_rad", math.inf))
    if not all(math.isfinite(x) and x >= 0.0 for x in (ex, eyz, eall, old)):
        f.append("invalid V31 directional perturbation caps")
    elif eall > FULL.up(old):
        f.append("V31 radial perturbation exceeds V12D parent")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    if d.get("P5_SAMPLE1_V30_CURRENT_SUBBOX_LIFT_V31") not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid V31 status")
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
        "status": d["P5_SAMPLE1_V30_CURRENT_SUBBOX_LIFT_V31"],
        "candidate": d.get("candidate_current_subboxes"),
        "q_rejected": d.get("q_ball_rejected_current_subboxes"),
        "current_source_rejected": d.get("current_source_rejected_current_subboxes"),
        "directional_source_rejected": d.get("directional_source_rejected_current_subboxes"),
        "yz_rejected": d.get("directional_yz_rejected_current_subboxes"),
        "radial_rejected": d.get("radial_constraint_rejected_current_subboxes"),
        "closed": d.get("closed_current_subboxes"),
        "open": d.get("open_current_subboxes"),
        "minimum_best_q": d.get("minimum_best_q_upper"),
        "maximum_best_q": d.get("maximum_best_q_upper"),
        "first_open": d.get("first_open_current_subbox"),
        "worst_open": d.get("worst_open_current_subbox"),
        "witness_closed": d.get("focused_first_witness_signed_subcell_closed_by_V30_lift"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
