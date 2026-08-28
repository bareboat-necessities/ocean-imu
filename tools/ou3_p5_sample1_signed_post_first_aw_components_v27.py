#!/usr/bin/env python3
"""V27: signed tangent/axial sample-1 residual at V23's first open q8 witness.

V26 keeps V10's two-block directional norms but still forgets the signed first
residual cell.  In V10's canonical un-Rx gauge the first useful residual is

    r0 = [0, r_t, r_z]

and the ideal first a_w correction is Delta=[0,Delta_t,Delta_z].  Before the
small certified transport/series gauge mismatch, the next residual satisfies

    r1 = alpha (r0-Delta) + (1-alpha) y_R0
         + (b1-alpha b0).

V21 already exposes the same source-row sample-1 force components

    f_y = -alpha Delta_t,
    f_z = g + alpha Delta_z.

Hence the tangent and axial nominal residual components can be evaluated
without introducing an independent physical-a_w ball:

    r1_y = alpha r_t + f_y + remainder,
    r1_z = alpha r_z - (f_z-g) + remainder.

The remainder retains V10's certified gravity-decay, bias-difference and
transport/series mismatch bounds.  The perpendicular component keeps V10's
already certified combined x cap.  V27 propagates these signed component
intervals through V21's signed one-plus-two gain components, adds V12D's PSD/S
correction perturbation as one Euclidean ball, intersects the result with V23's
authoritative first-open correction box, and reruns the unchanged V16/V15/V18
q<8 tests.

This is a focused diagnostic only.  No estimator, source domain, source branch,
6-rad shipping limit, q target, sample-1/word promotion, or N_H_words changes.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_current_exact_residual_subdivision_v23 as V23
import ou3_p5_sample1_exact_nonlinear_residual_v22 as V22
import ou3_p5_sample1_signed_cayley_q8_v14 as V14
import ou3_p5_sample1_signed_cayley_q8_v14d as V14D
import ou3_p5_sample1_signed_cayley_q8_v15 as V15
import ou3_p5_sample1_signed_cayley_q8_v16 as V16
import ou3_p5_sample1_signed_cayley_q8_v18 as V18

DEFAULT_DOMAIN = V23.DEFAULT_DOMAIN
SCHEMA = 2700
FULL = V14.FULL
Q_TARGET = V14.Q_TARGET
WITNESS = (0, 0, 19)


def _v10_witness_row(core: dict) -> dict:
    for row in core.get("rows", []):
        ids = (int(row["p_cell"]), int(row["tangent_residual_cell"]),
               int(row["axial_residual_cell"]))
        if ids == WITNESS:
            return row
    raise RuntimeError("V10 first-q8-witness row not found")


def _signed_residual_components(*, row: dict, parent: dict,
                                alpha: Interval, gravity: float):
    rt = Interval.outward_bounds(*map(float, row["first_tangent_residual_magnitude_mps2"]))
    rz = Interval.outward_bounds(*map(float, row["first_axial_residual_mps2"]))
    fyz = parent.get("sample1_force_components_yz_mps2", [])
    if len(fyz) != 2:
        raise RuntimeError("V21 sample-1 force components missing")
    fy, fz = V22._I(fyz[0]), V22._I(fyz[1])
    rho = float(row["sample1_full_residual_norm_upper_mps2"])
    rho_x = min(rho, float(row["sample1_combined_source_x_residual_upper_mps2"]))
    decay = float(row["ou_decay_times_first_rotational_residual_upper_mps2"])
    bias = float(row["bias_difference_upper_mps2"])
    geom = float(row["rotation_mismatch_residual_upper_mps2"])
    rem_hi = FULL.up(decay + FULL.up(bias + geom))
    rem = Interval.outward_bounds(-rem_hi, rem_hi)
    rx = Interval.outward_bounds(-rho_x, rho_x)
    ry = alpha * rt + fy + rem
    rz1 = alpha * rz - (fz - FULL.I(float(gravity))) + rem
    return [rx, ry, rz1], {
        "first_tangent_residual_mps2": rt.as_list(),
        "first_axial_residual_mps2": rz.as_list(),
        "sample1_force_y_mps2": fy.as_list(),
        "sample1_force_z_mps2": fz.as_list(),
        "combined_x_residual_abs_upper_mps2": rho_x,
        "gravity_decay_remainder_upper_mps2": decay,
        "bias_difference_remainder_upper_mps2": bias,
        "transport_series_remainder_upper_mps2": geom,
        "component_remainder_upper_mps2": rem_hi,
    }


def _nominal_correction(residual, parent: dict):
    gain = parent.get("gain_detail", {})
    perp = gain.get("perpendicular_gain_components", [])
    para = gain.get("parallel_gain_components", [])
    if len(perp) != 2 or len(para) != 2:
        raise RuntimeError("V21 signed gain components missing")
    gy, gz = V22._I(perp[0]), V22._I(perp[1])
    kxy, kxz = V22._I(para[0]), V22._I(para[1])
    rx, ry, rz = residual
    return [kxy * ry + kxz * rz, gy * rx, gz * rx]


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6,
          current_component_pieces: int = 4) -> dict:
    path = Path(domain_path).resolve()
    dom = json.loads(path.read_text(encoding="utf-8"))
    parent = V23.build(
        path, source_pieces=source_pieces, source_cell_index=source_cell_index,
        p_pieces=p_pieces, tangent_pieces=tangent_pieces,
        axial_pieces=axial_pieces, residual_x_pieces=residual_x_pieces,
        parallel_pieces=parallel_pieces,
        current_component_pieces=current_component_pieces)
    failures = [f"V23: {x}" for x in V23.validate(parent)]
    if parent.get("P5_SAMPLE1_CURRENT_EXACT_RESIDUAL_SUBDIVISION_V23") != "PASS":
        failures.append("V23 prerequisite did not pass")
    first = parent.get("first_open_current_subbox")
    if not isinstance(first, dict):
        failures.append("V23 first open current subbox missing")

    V10 = V23.V22.V21B.V21.V12D.V11.V10
    core = V10.build(path, source_pieces=source_pieces,
                     source_cell_index=source_cell_index, p_pieces=p_pieces,
                     tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    failures += [f"V10: {x}" for x in V10.validate(core)]
    try:
        row = _v10_witness_row(core)
        src, phase = V10.RG._source_phase_children(source_pieces)[source_cell_index]
        if phase != "due":
            raise RuntimeError("V27 focused witness requires first due source cell")
        F, _Q, _ = V10.FULL._transition_and_Q(src, dom)
        alpha = F[15][15]
        residual, detail = _signed_residual_components(
            row=row, parent=parent, alpha=alpha,
            gravity=float(dom["startup"]["gravity_mps2"]))
        nominal = _nominal_correction(residual, parent)
        eta = float(parent["V12D_correction_perturbation_norm_upper_rad"])
        if not (math.isfinite(eta) and eta >= 0.0):
            raise RuntimeError("invalid V12D correction perturbation")
        ec = Interval.outward_bounds(-eta, eta)
        source_box = [x + ec for x in nominal]
        radial_from_components = FULL.up(V14.CAYLEY1._norm_upper(nominal) + eta)
        radial_v10 = FULL.up(float(row["combined_directional_correction_norm_upper_rad"]) + eta)
        source_radial_hi = min(radial_from_components, radial_v10)
    except Exception as exc:
        failures.append(f"signed residual: {exc}")
        residual = nominal = source_box = None; detail = {}
        eta = source_radial_hi = math.inf

    closed = False; incompatible = False
    q = geo_q = product_q = product_w = math.inf
    before = joint = None; radial_lo = radial_hi = math.inf
    branches = []; narrowed = False
    if isinstance(first, dict) and source_box is not None:
        q = float(first["current_q_upper"])
        c = [V22._I(x) for x in first["q_ball_projected_current_component_box"]]
        before = [V22._I(x) for x in first["joint_correction_box_rad"]]
        joint = V22._intersect_boxes(before, source_box)
        incompatible = joint is None
        if not incompatible:
            radial_lo = max(float(first["correction_radial_lower_rad"]),
                            V14.CAYLEY2._norm_lower(joint))
            radial_hi = min(float(first["correction_radial_upper_rad"]),
                            source_radial_hi,
                            V14.CAYLEY1._norm_upper(joint))
            if radial_lo > radial_hi:
                incompatible = True
            else:
                geo = V15._geodesic_q_and_scalar_lower(q, radial_lo, radial_hi)
                geo_q = math.inf if geo is None else float(geo[0])
                wd, vd, branches, narrowed = V16.axis_cone_normalized_shipping_quaternion(
                    joint, radial_lower=radial_lo, radial_upper=radial_hi,
                    parent=V14D.radial_sinc_normalized_shipping_quaternion)
                cx_min = V14._minimum_abs(c[0])
                yz2 = max(0.0, FULL.up(q * q) - FULL.down(cx_min * cx_min))
                cyz = min(FULL.up(math.sqrt(yz2)), V18._yz_norm_upper(c[1], c[2]))
                chart = {"cx": c[0], "cy": c[1], "cz": c[2],
                         "cyz_norm_upper": cyz}
                parent_W = FULL.I(2.0) * wd - V14.CAYLEY1._dot(vd, c)
                W, _yb, _yj = V18._support_product_scalar(parent_W, wd, vd, chart)
                product_w, product_q = V14._qplus_from_product_scalar(q, W)
                closed = ((math.isfinite(geo_q) and geo_q < Q_TARGET)
                          or (math.isfinite(product_q) and product_q < Q_TARGET
                              and product_w > 0.0))

    status = "PASS" if not failures else "NOT_ESTABLISHED"
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_SIGNED_POST_FIRST_AW_COMPONENTS_V27",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V23_first_open_subbox_retained": True,
        "V10_exact_first_update_OU_cancellation_used": True,
        "signed_tangent_axial_first_residual_cell_retained": True,
        "V21_signed_one_plus_two_gain_components_used": True,
        "V10_transport_series_and_bias_remainders_retained": True,
        "V12D_correction_perturbation_retained_as_single_ball": True,
        "signed_residual_detail": detail,
        "sample1_signed_residual_box_mps2": None if residual is None else [x.as_list() for x in residual],
        "nominal_signed_correction_box_rad": None if nominal is None else [x.as_list() for x in nominal],
        "V23_first_open_correction_box_rad": None if before is None else [x.as_list() for x in before],
        "joint_signed_correction_box_rad": None if joint is None else [x.as_list() for x in joint],
        "V12D_correction_perturbation_ball_upper_rad": eta,
        "source_correlated_radial_upper_rad": source_radial_hi,
        "source_constraints_incompatible": incompatible,
        "V23_first_open_current_q_upper": q,
        "directional_radial_lower_rad": radial_lo,
        "directional_radial_upper_rad": radial_hi,
        "axis_cone_narrowed": narrowed,
        "quaternion_branches": branches,
        "geodesic_q_upper": geo_q,
        "product_abs_W_lower": product_w,
        "product_q_upper": product_q,
        "first_open_subbox_closed_inside_q8": bool(closed or incompatible),
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_SIGNED_POST_FIRST_AW_COMPONENTS_V27": status,
        "next_obligation": (
            "LIFT_V27_SIGNED_POST_FIRST_AW_COMPONENTS_OVER_ALL_V23_CURRENT_SUBBOXES"
            if (closed or incompatible) and not failures else
            "REFINE_TANGENT_AXIAL_GRAVITY_ROTATION_COMPONENTS_OR_PSD_S_DIRECTION_AT_FIRST_OPEN_SUBBOX"),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA: f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_SIGNED_POST_FIRST_AW_COMPONENTS_V27":
        f.append("qualification mismatch")
    for k in ("source_generated_not_trajectory_fit", "V23_first_open_subbox_retained",
              "V10_exact_first_update_OU_cancellation_used",
              "signed_tangent_axial_first_residual_cell_retained",
              "V21_signed_one_plus_two_gain_components_used",
              "V10_transport_series_and_bias_remainders_retained",
              "V12D_correction_perturbation_retained_as_single_ball"):
        if d.get(k) is not True: f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed", "deployed_correction_limit_increased",
              "q8_composed_here", "q8_word_promoted_here", "whole_word_promoted_here",
              "N_H_words_set_here"):
        if d.get(k) is not False: f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET: f.append("q target changed")
    if d.get("sample1_signed_residual_box_mps2") is None:
        f.append("signed residual box missing")
    if d.get("nominal_signed_correction_box_rad") is None:
        f.append("nominal signed correction box missing")
    if not d.get("source_constraints_incompatible") and d.get("joint_signed_correction_box_rad") is None:
        f.append("compatible V27 witness lost joint correction box")
    if d.get("P5_SAMPLE1_SIGNED_POST_FIRST_AW_COMPONENTS_V27") not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid V27 status")
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
        "status": d["P5_SAMPLE1_SIGNED_POST_FIRST_AW_COMPONENTS_V27"],
        "detail": d.get("signed_residual_detail"),
        "residual_box": d.get("sample1_signed_residual_box_mps2"),
        "nominal_correction_box": d.get("nominal_signed_correction_box_rad"),
        "source_radial_upper": d.get("source_correlated_radial_upper_rad"),
        "q_current": d.get("V23_first_open_current_q_upper"),
        "radial_lower": d.get("directional_radial_lower_rad"),
        "radial_upper": d.get("directional_radial_upper_rad"),
        "geodesic_q": d.get("geodesic_q_upper"),
        "product_W": d.get("product_abs_W_lower"),
        "product_q": d.get("product_q_upper"),
        "incompatible": d.get("source_constraints_incompatible"),
        "closed_q8": d.get("first_open_subbox_closed_inside_q8"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
