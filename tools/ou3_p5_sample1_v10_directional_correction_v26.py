#!/usr/bin/env python3
"""V26: source-correlated V10 directional correction at V23's first open q8 subbox.

V25 shows that the remaining loss is not the scalar attitude-gain norm: even the
certified V10 block norm is too large when the full physical a_w error is kept
as an independent residual ball.  V10 already proves a stronger directional
statement in the same source row.  In its canonical one-plus-two gauge,

    r_x  -> d_yz with gain k_perp,
    r_yz -> d_x  with gain k_parallel,

and its exact first-update/OU cancellation gives a complete residual norm rho
and a much smaller source-correlated perpendicular cap rho_x.  Thus

    |d_x|    <= k_parallel rho,
    ||d_yz|| <= k_perp min(rho_x,rho).

V12D's later PSD/S correction perturbation eta is retained as one Euclidean
correction ball, so its coordinate projections add eta to both directional
caps while the full radial constraint remains d_V10+eta.

V26 applies these three constraints only to V23's authoritative first open
current-Cayley subbox and reruns the unchanged V16 axis-cone, V15 geodesic and
V18 signed-product checks.  This is a focused diagnostic; it does not promote
sample 1, a word, P5, or N_H_words, and changes no estimator/domain/6-rad/q8
contract.
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
SCHEMA = 2600
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


def _directional_caps(row: dict, eta: float) -> dict:
    rho = float(row["sample1_full_residual_norm_upper_mps2"])
    rho_x = min(rho, float(row["sample1_combined_source_x_residual_upper_mps2"]))
    kp = float(row["Ktheta_perpendicular_block_upper"])
    ka = float(row["Ktheta_parallel_block_upper"])
    d10 = float(row["combined_directional_correction_norm_upper_rad"])
    eta = float(eta)
    vals = (rho, rho_x, kp, ka, d10, eta)
    if not all(math.isfinite(x) and x >= 0.0 for x in vals):
        raise RuntimeError("invalid V10/V12D directional bounds")
    return {
        "sample1_full_residual_norm_upper_mps2": rho,
        "sample1_combined_source_x_residual_upper_mps2": rho_x,
        "Ktheta_perpendicular_block_upper": kp,
        "Ktheta_parallel_block_upper": ka,
        "V10_directional_correction_upper_rad": d10,
        "V12D_correction_perturbation_ball_upper_rad": eta,
        "dx_abs_upper_rad": FULL.up(FULL.up(ka * rho) + eta),
        "dyz_norm_upper_rad": FULL.up(FULL.up(kp * rho_x) + eta),
        "radial_upper_rad": FULL.up(d10 + eta),
    }


def _intersect(x: Interval, lo: float, hi: float) -> Interval | None:
    a = max(x.lo, float(lo)); b = min(x.hi, float(hi))
    return None if b < a else Interval(a, b)


def _clip_yz_ball(y: Interval, z: Interval, radius: float):
    r = float(radius)
    if not (math.isfinite(r) and r >= 0.0):
        raise ValueError("finite nonnegative yz radius required")
    yy, zz = y, z
    for _ in range(2):
        ymin = V14._minimum_abs(yy); zmin = V14._minimum_abs(zz)
        if FULL.down(ymin * ymin + zmin * zmin) > FULL.up(r * r):
            return None
        y2 = max(0.0, FULL.up(FULL.up(r * r) - FULL.down(zmin * zmin)))
        yy = _intersect(yy, -FULL.up(math.sqrt(y2)), FULL.up(math.sqrt(y2)))
        if yy is None:
            return None
        ymin = V14._minimum_abs(yy)
        z2 = max(0.0, FULL.up(FULL.up(r * r) - FULL.down(ymin * ymin)))
        zz = _intersect(zz, -FULL.up(math.sqrt(z2)), FULL.up(math.sqrt(z2)))
        if zz is None:
            return None
    return yy, zz


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6,
          current_component_pieces: int = 4) -> dict:
    path = Path(domain_path).resolve()
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
        caps = _directional_caps(
            row, float(parent["V12D_correction_perturbation_norm_upper_rad"]))
    except Exception as exc:
        failures.append(f"directional caps: {exc}")
        caps = {"dx_abs_upper_rad": math.inf,
                "dyz_norm_upper_rad": math.inf,
                "radial_upper_rad": math.inf}

    closed = False; incompatible = False
    q = geo_q = product_q = product_w = math.inf
    before = after = None; radial_lo = radial_hi = math.inf
    branches = []; narrowed = False
    if isinstance(first, dict):
        q = float(first["current_q_upper"])
        c = [V22._I(x) for x in first["q_ball_projected_current_component_box"]]
        before = [V22._I(x) for x in first["joint_correction_box_rad"]]
        x = _intersect(before[0], -float(caps["dx_abs_upper_rad"]),
                       float(caps["dx_abs_upper_rad"]))
        yz = None if x is None else _clip_yz_ball(
            before[1], before[2], float(caps["dyz_norm_upper_rad"]))
        incompatible = x is None or yz is None
        if not incompatible:
            after = [x, yz[0], yz[1]]
            radial_lo = max(float(first["correction_radial_lower_rad"]),
                            V14.CAYLEY2._norm_lower(after))
            radial_hi = min(float(first["correction_radial_upper_rad"]),
                            float(caps["radial_upper_rad"]),
                            V14.CAYLEY1._norm_upper(after))
            if radial_lo > radial_hi:
                incompatible = True
            else:
                geo = V15._geodesic_q_and_scalar_lower(q, radial_lo, radial_hi)
                geo_q = math.inf if geo is None else float(geo[0])
                wd, vd, branches, narrowed = V16.axis_cone_normalized_shipping_quaternion(
                    after, radial_lower=radial_lo, radial_upper=radial_hi,
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
        "qualification": "OU3_P5_SAMPLE1_V10_DIRECTIONAL_FIRST_OPEN_WITNESS_V26",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V23_first_open_subbox_retained": True,
        "V10_combined_perpendicular_residual_identity_revalidated": True,
        "V10_one_plus_two_directional_caps_used": True,
        "V12D_correction_perturbation_retained_as_single_ball": True,
        "V16_axis_cone_and_V18_signed_product_retained": True,
        "V10_directional_source_detail": caps,
        "V23_first_open_current_q_upper": q,
        "V23_first_open_correction_box_rad": None if before is None else [x.as_list() for x in before],
        "directional_correction_box_rad": None if after is None else [x.as_list() for x in after],
        "directional_constraints_incompatible": incompatible,
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
        "P5_SAMPLE1_V10_DIRECTIONAL_FIRST_OPEN_WITNESS_V26": status,
        "next_obligation": (
            "LIFT_V26_DIRECTIONAL_CORRECTION_OVER_ALL_V23_CURRENT_SUBBOXES"
            if (closed or incompatible) and not failures else
            "DERIVE_SIGNED_TANGENT_AXIAL_POST_FIRST_AW_COMPONENTS_AT_FIRST_OPEN_SUBBOX"),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA: f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_V10_DIRECTIONAL_FIRST_OPEN_WITNESS_V26":
        f.append("qualification mismatch")
    for k in ("source_generated_not_trajectory_fit", "V23_first_open_subbox_retained",
              "V10_combined_perpendicular_residual_identity_revalidated",
              "V10_one_plus_two_directional_caps_used",
              "V12D_correction_perturbation_retained_as_single_ball",
              "V16_axis_cone_and_V18_signed_product_retained"):
        if d.get(k) is not True: f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed",
              "deployed_correction_limit_increased", "q8_composed_here",
              "q8_word_promoted_here", "whole_word_promoted_here",
              "N_H_words_set_here"):
        if d.get(k) is not False: f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET: f.append("q target changed")
    detail = d.get("V10_directional_source_detail", {})
    for k in ("dx_abs_upper_rad", "dyz_norm_upper_rad", "radial_upper_rad"):
        x = detail.get(k)
        if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or float(x) < 0.0:
            f.append(f"invalid {k}")
    if not d.get("directional_constraints_incompatible") and d.get("directional_correction_box_rad") is None:
        f.append("compatible V26 witness lost directional correction box")
    if d.get("P5_SAMPLE1_V10_DIRECTIONAL_FIRST_OPEN_WITNESS_V26") not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid V26 status")
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
        "status": d["P5_SAMPLE1_V10_DIRECTIONAL_FIRST_OPEN_WITNESS_V26"],
        "detail": d.get("V10_directional_source_detail"),
        "q_current": d.get("V23_first_open_current_q_upper"),
        "radial_lower": d.get("directional_radial_lower_rad"),
        "radial_upper": d.get("directional_radial_upper_rad"),
        "geodesic_q": d.get("geodesic_q_upper"),
        "product_W": d.get("product_abs_W_lower"),
        "product_q": d.get("product_q_upper"),
        "incompatible": d.get("directional_constraints_incompatible"),
        "closed_q8": d.get("first_open_subbox_closed_inside_q8"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
