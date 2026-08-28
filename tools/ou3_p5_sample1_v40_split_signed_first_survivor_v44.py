#!/usr/bin/env python3
"""V44: test V28/V31 split-signed source geometry on V41's first survivor.

V41 installs V40's exact Joseph-component transport in the complete source-cell-0
sample-1 q<8 cover.  Its first remaining signed-Cayley cell is source row
(p,t,a)=(0,0,23).  V42 final-current bisection and V43's exact nonlinear
current residual do not narrow that cell: both retain the same q=8.3445...
first survivor.

The focused V27/V28/V31 chain contains a different dependency that V43 did not
lift globally.  In the canonical un-Rx source gauge the first useful residual
has signed tangent/axial coordinates and the one-step OU cancellation gives
componentwise sample-1 residual bounds.  V28 further separates tangent and
axial gravity-decay remainders, while V31 keeps distinct x, yz and total
V12D correction perturbation radii.

This producer applies that construction to the actual V41 first-survivor row
under the V40 parent.  It reconstructs the row's 6x6 V14 correction subcells,
checks the V18B/V16/V15 baseline, intersects each correction subcell with the
V28 signed source correction box plus V31 directional perturbation caps, then
reruns the same axis-cone, geodesic and current-yz product checks.

This is deliberately a focused diagnostic before another 461376-cell run.
It changes no estimator setting, theorem domain, source branch, six-radian
shipping correction limit, q<8 target, source-word language, promotion flag,
or N_H.  A successful diagnostic only authorizes a global V45 lift.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_first_psd_exact_joseph_components_v40 as V40
import ou3_p5_sample1_signed_cayley_q8_v14 as V14
import ou3_p5_sample1_signed_cayley_q8_v14d as V14D
import ou3_p5_sample1_signed_cayley_q8_v15 as V15
import ou3_p5_sample1_signed_cayley_q8_v16 as V16
import ou3_p5_sample1_signed_cayley_q8_v18 as V18
import ou3_p5_sample1_signed_cayley_q8_v18b as V18B
import ou3_p5_sample1_effective_input_correction_v21 as V21
import ou3_p5_sample1_split_gravity_signed_components_v28 as V28
import ou3_p5_sample1_theta_x_gain_perturbation_v30 as V30
import ou3_p5_sample1_v30_current_subbox_lift_v31 as V31

DEFAULT_DOMAIN = V40.DEFAULT_DOMAIN
SCHEMA = 4400
Q_TARGET = 8.0
FULL = V14.FULL
V12D = V40.V12D
WITNESS = (0, 0, 23)
V41_REFERENCE_FIRST_Q = 8.344528951460543


def _intersect_boxes(a, b):
    if len(a) != 3 or len(b) != 3:
        raise ValueError("three-component boxes required")
    out = []
    for x, y in zip(a, b):
        lo = max(x.lo, y.lo)
        hi = min(x.hi, y.hi)
        if hi < lo:
            return None
        out.append(Interval(lo, hi))
    return out


def _find_row(rows, witness=WITNESS):
    for row in rows:
        ids = (int(row["p_cell"]), int(row["tangent_residual_cell"]),
               int(row["axial_residual_cell"]))
        if ids == tuple(witness):
            return row
    raise RuntimeError(f"row {tuple(witness)} not found")


def _eval_q(*, q: float, chart: dict, dbox, radial_lo: float,
            radial_hi: float):
    if radial_lo > radial_hi:
        return {
            "closed": True, "incompatible": True,
            "geodesic_q": 0.0, "product_q": 0.0,
            "product_w": math.inf, "best_q": 0.0,
            "axis_narrowed": False,
        }
    geo = V15._geodesic_q_and_scalar_lower(q, radial_lo, radial_hi)
    geo_q = math.inf if geo is None else float(geo[0])
    wd, vd, branches, narrowed = V16.axis_cone_normalized_shipping_quaternion(
        dbox, radial_lower=radial_lo, radial_upper=radial_hi,
        parent=V14D.radial_sinc_normalized_shipping_quaternion)
    cx = chart["cx"]
    cyz = float(chart["cyz_norm_upper"])
    xdot = vd[0] * cx
    vdyz = V18._yz_norm_upper(vd[1], vd[2])
    yzdot = FULL.up(vdyz * cyz)
    parent_W = FULL.I(2.0) * wd - (
        xdot + Interval.outward_bounds(-yzdot, yzdot))
    W, _yb, _yj = V18._support_product_scalar(parent_W, wd, vd, chart)
    product_w, product_q = V14._qplus_from_product_scalar(q, W)
    closed = ((math.isfinite(geo_q) and geo_q < Q_TARGET)
              or (math.isfinite(product_q) and product_q < Q_TARGET
                  and product_w > 0.0))
    return {
        "closed": closed, "incompatible": False,
        "geodesic_q": geo_q, "product_q": float(product_q),
        "product_w": float(product_w),
        "best_q": min(geo_q, float(product_q)),
        "axis_narrowed": bool(narrowed),
        "branches": list(branches),
        "product_W": W.as_list(),
    }


def _row_source_directional(*, path: Path, dom: dict, src: dict,
                            first: dict, base: dict, vr: dict,
                            p, t, Y, r, alpha, qaw) -> dict:
    g = float(dom["startup"]["gravity_mps2"])
    rt = Interval.outward_bounds(
        *map(float, base["first_tangent_residual_magnitude_mps2"]))
    rz = Interval.outward_bounds(*map(float, base["first_axial_residual_mps2"]))
    d0 = Interval.outward_bounds(*map(float, base["first_attitude_correction_rad"]))
    D = FULL.I(g * g) * t + p + r
    a = t * (p + r) / D
    c0 = -(FULL.I(g) * t * p / D)
    b = p * (FULL.I(g * g) * t + r) / D
    bz = p * r / (p + r)
    det_first = t * p * r / D
    fy = -(alpha * (p / D) * rt)
    fz = FULL.I(g) + alpha * (p / (p + r)) * rz
    (gy, gz), (kxy, kxz), gain_detail = V21.V13._signed_gain_components(
        a=a, Y=Y, c0=c0, alpha=alpha, qaw=qaw, b=b, bz=bz,
        det_first=det_first, d=d0, fy=fy, fz=fz, r=r)

    gravity_detail = V28._gravity_component_decay_bounds(
        cosine_lower=float(first["post_prediction_true_gravity_cosine_lower"]),
        alpha_lower=float(alpha.lo), gravity=g)
    residual, residual_detail = V28._split_signed_residual_components(
        row=base,
        parent={"sample1_force_components_yz_mps2": [fy.as_list(), fz.as_list()]},
        alpha=alpha, gravity=g, gravity_detail=gravity_detail)
    rx, ry, rz1 = residual
    nominal = [kxy * ry + kxz * rz1, gy * rx, gz * rx]

    row_detail = V30._theta_x_gain_perturbation_upper(vr, base)
    caps = V31._refined_caps(base=base, vr=vr, row_detail=row_detail)
    ex = float(caps["x_correction_perturbation_abs_upper_rad"])
    eyz = float(caps["yz_correction_perturbation_norm_upper_rad"])
    eall = float(caps["total_correction_perturbation_norm_upper_rad"])
    if not all(math.isfinite(x) and x >= 0.0 for x in (ex, eyz, eall)):
        raise RuntimeError("invalid directional correction caps")
    source_box = [
        nominal[0] + Interval.outward_bounds(-ex, ex),
        nominal[1] + Interval.outward_bounds(-eyz, eyz),
        nominal[2] + Interval.outward_bounds(-eyz, eyz),
    ]
    nominal_yz = V18._yz_norm_upper(nominal[1], nominal[2])
    yz_hi = FULL.up(nominal_yz + eyz)
    radial_hi = min(
        FULL.up(V14.CAYLEY1._norm_upper(nominal) + eall),
        FULL.up(float(base["combined_directional_correction_norm_upper_rad"]) + eall),
    )
    radial_lo = max(0.0, FULL.down(V14.CAYLEY2._norm_lower(nominal) - eall))
    return {
        "force": [FULL.I(0.0), fy, fz],
        "gain_detail": gain_detail,
        "residual": residual,
        "residual_detail": residual_detail,
        "gravity_detail": gravity_detail,
        "nominal": nominal,
        "caps": caps,
        "source_box": source_box,
        "yz_hi": yz_hi,
        "radial_lo": radial_lo,
        "radial_hi": radial_hi,
        "theta_x_gain_detail": row_detail,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    path = Path(domain_path).resolve()
    dom = json.loads(path.read_text(encoding="utf-8"))
    src, phase = V14.RG._source_phase_children(source_pieces)[source_cell_index]
    failures = []
    if phase != "due":
        failures.append("V44 focused witness requires first-due source cell")

    original_psd = V12D._first_psd_perturbation_tangent
    V12D._first_psd_perturbation_tangent = \
        V40._first_psd_perturbation_exact_joseph_components
    try:
        v12 = V12D.build(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    finally:
        V12D._first_psd_perturbation_tangent = original_psd
    failures += [f"V40/V12D: {x}" for x in V12D.validate(v12)]
    if V12D._first_psd_perturbation_tangent is not original_psd:
        failures.append("V40 temporary PSD hook was not restored")

    core = V12D.V11.V10.build(
        path, source_pieces=source_pieces,
        source_cell_index=source_cell_index, p_pieces=p_pieces,
        tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    failures += [f"V10: {x}" for x in V12D.V11.V10.validate(core)]
    first_all = V12D.V11.FIRST.build(path, source_pieces=source_pieces)
    failures += [f"first: {x}" for x in V12D.V11.FIRST.validate(first_all)]
    if len(core.get("rows", [])) != len(v12.get("rows", [])):
        failures.append("V10/V40-V12D row counts differ")

    base = _find_row(core.get("rows", []))
    vr = _find_row(v12.get("rows", []))
    first = dict(first_all)
    fr = first_all["source_cells"][source_cell_index]
    p_all = Interval.outward_bounds(*map(float, fr["P_aw_variance_interval"]))
    pcells = V12D.V11.SUB.parts(p_all.lo, p_all.hi, p_pieces)
    p = pcells[WITNESS[0]]
    h = float(src["dt_s"])
    g = float(dom["startup"]["gravity_mps2"])
    tilt, yaw, eps = V14.RG._attitude_covariance_epsilon(path, h)
    t = Interval.outward_bounds(tilt, FULL.up(tilt + eps))
    Y = Interval.outward_bounds(yaw, FULL.up(yaw + eps))
    vec = V12D.V11.VECTOR.build()
    r = FULL._R_diag(float(
        vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))[0][0]
    F, Q, _ = FULL._transition_and_Q(src, dom)
    alpha = F[15][15]
    qaw = Q[15][15]
    ds = float(v12.get("sample1_S_perturbation_bounds", {}).get(
        "sample1_S_attitude_correction_upper_rad", 0.0))

    original_rotate = V18._rotate_yz_rx_transpose
    V18._rotate_yz_rx_transpose = V18B._rotate_yz_rx_transpose
    try:
        chart = V21._current_component_chart(
            first=first, base=base, vr=vr, dom=dom, src=src,
            sample1_s_angle=ds)
    finally:
        V18._rotate_yz_rx_transpose = original_rotate
    if V18._rotate_yz_rx_transpose is not original_rotate:
        failures.append("V18 signed-angle rotation hook was not restored")

    source = _row_source_directional(
        path=path, dom=dom, src=src, first=first, base=base, vr=vr,
        p=p, t=t, Y=Y, r=r, alpha=alpha, qaw=qaw)

    kperp = float(base["Ktheta_perpendicular_block_upper"])
    kpar = float(base["Ktheta_parallel_block_upper"])
    drho = float(vr["total_residual_perturbation_upper_mps2"])
    dk = float(vr["sample1_attitude_gain_operator_perturbation_upper"])
    rho = float(base["sample1_full_residual_norm_upper_mps2"])
    rho_x = min(rho, float(base["sample1_combined_source_x_residual_upper_mps2"]))
    eta = FULL.up(FULL.up(max(kperp, kpar) * drho)
                  + FULL.up(dk * FULL.up(rho + drho)))
    global_hi = float(vr["V12C_correction_norm_upper_rad"])
    e = Interval.outward_bounds(-eta, eta)

    total = parent_open = refined_closed = incompatible = refined = 0
    first_parent_open = first_parent_refined = worst_parent_open = None
    rx_cells = V12D.V11.SUB.parts(-rho_x, rho_x, residual_x_pieces)
    for rxi, rxc in enumerate(rx_cells):
        rx_min = V14._minimum_abs(rxc)
        rem2 = max(0.0, FULL.up(rho * rho) - FULL.down(rx_min * rx_min))
        ryz_hi = FULL.up(math.sqrt(rem2))
        ux_hi = FULL.up(kpar * ryz_hi)
        u_cells = V12D.V11.SUB.parts(-ux_hi, ux_hi, parallel_pieces)
        for ui, uc in enumerate(u_cells):
            total += 1
            dbox = [uc + e,
                    source["gain_detail"]["perpendicular_gain_components"][0]
                    if False else FULL.I(0.0),
                    FULL.I(0.0)]
            # Reconstruct the exact V14 component parent from the row's signed
            # perpendicular gains; the source directional box below is an
            # independent enclosure of the same correction.
            gain = source["gain_detail"]
            perp = gain["perpendicular_gain_components"]
            gy = Interval.outward_bounds(*map(float, perp[0]))
            gz = Interval.outward_bounds(*map(float, perp[1]))
            dbox = [uc + e, gy * rxc + e, gz * rxc + e]
            box_lo = V14.CAYLEY2._norm_lower(dbox)
            rx_abs = rxc.abs_upper(); u_abs = uc.abs_upper()
            nominal_hi = FULL.up(math.sqrt(FULL.up(
                FULL.up(u_abs * u_abs) + FULL.up((kperp * rx_abs) ** 2))))
            radial_hi = min(global_hi, FULL.up(nominal_hi + eta))
            radial_lo = min(box_lo, radial_hi)
            q = float(chart["q1"])
            parent_eval = _eval_q(
                q=q, chart=chart, dbox=dbox,
                radial_lo=radial_lo, radial_hi=radial_hi)
            if parent_eval["closed"]:
                continue
            parent_open += 1
            prow = {
                "rx_index": rxi, "parallel_index": ui,
                "nominal_rx_subcell_mps2": rxc.as_list(),
                "nominal_parallel_correction_subcell_rad": uc.as_list(),
                "correction_component_box_rad": [x.as_list() for x in dbox],
                "correction_radial_lower_rad": radial_lo,
                "correction_radial_upper_rad": radial_hi,
                "current_q_upper": q,
                **parent_eval,
            }
            if first_parent_open is None:
                first_parent_open = prow
            if (worst_parent_open is None
                    or float(parent_eval["best_q"]) > float(worst_parent_open["best_q"])):
                worst_parent_open = prow

            joint = _intersect_boxes(dbox, source["source_box"])
            if joint is None:
                incompatible += 1
                refined_closed += 1
                rev = {**prow, "source_incompatible": True,
                       "refined_closed": True, "refined_best_q": 0.0}
            else:
                yz = V31.V29._clip_yz_to_radius(
                    joint[1], joint[2], float(source["yz_hi"]))
                if yz is None:
                    incompatible += 1
                    refined_closed += 1
                    rev = {**prow, "source_incompatible": True,
                           "refined_closed": True, "refined_best_q": 0.0}
                else:
                    joint[1], joint[2] = yz
                    rhi = min(radial_hi, float(source["radial_hi"]),
                               V14.CAYLEY1._norm_upper(joint))
                    rlo = max(radial_lo, float(source["radial_lo"]),
                               V14.CAYLEY2._norm_lower(joint))
                    if rlo > rhi:
                        incompatible += 1
                        refined_closed += 1
                        rev = {**prow, "source_incompatible": True,
                               "refined_closed": True, "refined_best_q": 0.0}
                    else:
                        reval = _eval_q(
                            q=q, chart=chart, dbox=joint,
                            radial_lo=rlo, radial_hi=rhi)
                        narrowed = any(
                            joint[i].lo > dbox[i].lo or joint[i].hi < dbox[i].hi
                            for i in range(3)) or rlo > radial_lo or rhi < radial_hi
                        refined += int(narrowed)
                        refined_closed += int(reval["closed"])
                        rev = {
                            **prow,
                            "source_incompatible": False,
                            "source_joint_correction_box_rad": [x.as_list() for x in joint],
                            "source_radial_lower_rad": rlo,
                            "source_radial_upper_rad": rhi,
                            "refined_closed": bool(reval["closed"]),
                            "refined_best_q": float(reval["best_q"]),
                            "refined_product_q": float(reval["product_q"]),
                            "refined_geodesic_q": float(reval["geodesic_q"]),
                            "refined_product_w": float(reval["product_w"]),
                        }
            if first_parent_refined is None:
                first_parent_refined = rev

    ref_match = (first_parent_open is not None
                 and math.isfinite(float(first_parent_open["best_q"]))
                 and abs(float(first_parent_open["best_q"])
                         - V41_REFERENCE_FIRST_Q) <= 2.0e-9)
    if not ref_match:
        failures.append("reconstructed first V41 survivor does not match V41 q reference")

    first_closed = bool(first_parent_refined and first_parent_refined["refined_closed"])
    all_row_closed = parent_open > 0 and refined_closed == parent_open
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_V40_SPLIT_SIGNED_FIRST_SURVIVOR_V44",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V40_exact_Joseph_parent_used": True,
        "V41_first_survivor_row": list(WITNESS),
        "V41_reference_first_q_upper": V41_REFERENCE_FIRST_Q,
        "V41_first_survivor_reconstructed": ref_match,
        "V28_split_tangent_axial_gravity_residual_used": True,
        "V31_directional_x_yz_total_perturbation_caps_used": True,
        "V18B_signed_full_angle_current_chart_used": True,
        "V16_axis_cone_V15_geodesic_V18_yz_support_retained": True,
        "temporary_V40_and_V18B_hooks_restored": (
            V12D._first_psd_perturbation_tangent is original_psd
            and V18._rotate_yz_rx_transpose is original_rotate),
        "evaluated_witness_correction_subcells": total,
        "parent_open_witness_subcells": parent_open,
        "source_directionally_refined_parent_open_subcells": refined,
        "source_incompatible_parent_open_subcells": incompatible,
        "source_closed_parent_open_subcells": refined_closed,
        "first_parent_open_subcell": first_parent_open,
        "first_parent_open_after_V44": first_parent_refined,
        "worst_parent_open_subcell": worst_parent_open,
        "first_V41_survivor_closed_by_V44": first_closed,
        "all_parent_open_subcells_in_witness_row_closed_by_V44": all_row_closed,
        "sample1_current_component_box": [
            chart["cx"].as_list(), chart["cy"].as_list(), chart["cz"].as_list()],
        "sample1_current_q_upper": float(chart["q1"]),
        "V44_signed_residual_box_mps2": [x.as_list() for x in source["residual"]],
        "V44_nominal_signed_correction_box_rad": [x.as_list() for x in source["nominal"]],
        "V44_source_correction_box_rad": [x.as_list() for x in source["source_box"]],
        "V44_source_yz_correction_norm_upper_rad": float(source["yz_hi"]),
        "V44_source_radial_lower_rad": float(source["radial_lo"]),
        "V44_source_radial_upper_rad": float(source["radial_hi"]),
        "V44_directional_perturbation_caps": dict(source["caps"]),
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_V40_SPLIT_SIGNED_FIRST_SURVIVOR_V44": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_V44_SPLIT_SIGNED_DIRECTIONAL_INTERSECTION_OVER_FULL_V41_SOURCE_CELL0_COVER"
            if first_closed and not failures else
            "REFINE_CURRENT_COMPONENT_DEPENDENCY_INSIDE_V44_FIRST_SURVIVOR"
        ),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_V40_SPLIT_SIGNED_FIRST_SURVIVOR_V44":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit", "V40_exact_Joseph_parent_used",
        "V41_first_survivor_reconstructed",
        "V28_split_tangent_axial_gravity_residual_used",
        "V31_directional_x_yz_total_perturbation_caps_used",
        "V18B_signed_full_angle_current_chart_used",
        "V16_axis_cone_V15_geodesic_V18_yz_support_retained",
        "temporary_V40_and_V18B_hooks_restored",
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
    if int(d.get("evaluated_witness_correction_subcells", 0)) != 36:
        f.append("V44 witness did not evaluate the complete 6x6 correction partition")
    if int(d.get("parent_open_witness_subcells", 0)) <= 0:
        f.append("V44 witness has no parent-open subcells")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    if d.get("P5_SAMPLE1_V40_SPLIT_SIGNED_FIRST_SURVIVOR_V44") not in (
            "PASS", "NOT_ESTABLISHED"):
        f.append("invalid V44 status")
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
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(
        x.domain, source_pieces=x.source_pieces,
        source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
        tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces,
        residual_x_pieces=x.residual_x_pieces,
        parallel_pieces=x.parallel_pieces)
    vf = validate(d); d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_V40_SPLIT_SIGNED_FIRST_SURVIVOR_V44"],
        "parent_open": d["parent_open_witness_subcells"],
        "refined": d["source_directionally_refined_parent_open_subcells"],
        "incompatible": d["source_incompatible_parent_open_subcells"],
        "closed": d["source_closed_parent_open_subcells"],
        "first_parent": d["first_parent_open_subcell"],
        "first_refined": d["first_parent_open_after_V44"],
        "first_closed": d["first_V41_survivor_closed_by_V44"],
        "all_row_closed": d["all_parent_open_subcells_in_witness_row_closed_by_V44"],
        "source_residual": d["V44_signed_residual_box_mps2"],
        "source_correction": d["V44_source_correction_box_rad"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
