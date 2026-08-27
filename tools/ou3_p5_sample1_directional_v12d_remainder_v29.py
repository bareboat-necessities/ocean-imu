#!/usr/bin/env python3
"""V29: preserve V12D correction-perturbation direction at the first q8 witness.

V28 keeps the signed first-residual cell and splits the tangent/axial gravity
remainder, but it still adds V12D's complete PSD/S correction perturbation as
one scalar interval independently to x, y, and z.  V12D already proves the
stronger decomposition

    delta d = K0 delta r + delta K (r + delta r).

For the exact one-plus-two nominal gain, residual perturbations map into
orthogonal correction blocks:

    |delta d_x|       <= k_parallel ||delta r||,
    ||delta d_yz||_2  <= k_perp     ||delta r||.

Only the gain-perturbation term delta K(r+delta r) remains direction-free.
Thus, with drho=||delta r|| and dk=||delta K||,

    eK = dk (rho + drho),
    e_x  = k_parallel drho + eK,
    e_yz = k_perp     drho + eK,
    e_all= max(k_parallel,k_perp) drho + eK.

V29 uses V28's signed nominal correction, intersects x with e_x, retains a
Euclidean yz support e_yz, and retains e_all for the radial certificate.  The
same V23 first-open current box and V16/V15/V18 q<8 checks are then evaluated.
No estimator, source domain, shipping limit, q target, or promotion state is
changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_split_gravity_signed_components_v28 as V28

DEFAULT_DOMAIN = V28.DEFAULT_DOMAIN
SCHEMA = 2900
FULL = V28.FULL
Q_TARGET = V28.Q_TARGET
WITNESS = (0, 0, 19)


def _I(x):
    return Interval.outward_bounds(*map(float, x))


def _witness_row(core: dict) -> dict:
    for row in core.get("rows", []):
        ids = (int(row["p_cell"]), int(row["tangent_residual_cell"]),
               int(row["axial_residual_cell"]))
        if ids == WITNESS:
            return row
    raise RuntimeError("first-q8 V12D/V10 witness row not found")


def _directional_perturbation_caps(*, k_perp: float, k_parallel: float,
                                   drho: float, dk: float, rho: float) -> dict:
    vals = (k_perp, k_parallel, drho, dk, rho)
    if not all(math.isfinite(float(x)) and float(x) >= 0.0 for x in vals):
        raise ValueError("finite nonnegative directional perturbation inputs required")
    gain_ball = FULL.up(dk * FULL.up(rho + drho))
    ex = FULL.up(FULL.up(k_parallel * drho) + gain_ball)
    eyz = FULL.up(FULL.up(k_perp * drho) + gain_ball)
    eall = FULL.up(FULL.up(max(k_perp, k_parallel) * drho) + gain_ball)
    return {
        "gain_perturbation_ball_upper_rad": gain_ball,
        "x_correction_perturbation_abs_upper_rad": ex,
        "yz_correction_perturbation_norm_upper_rad": eyz,
        "total_correction_perturbation_norm_upper_rad": eall,
    }


def _intersect(a: Interval, b: Interval):
    lo = max(a.lo, b.lo); hi = min(a.hi, b.hi)
    return None if hi < lo else Interval(lo, hi)


def _intersect_boxes(a, b):
    out = []
    for x, y in zip(a, b):
        z = _intersect(x, y)
        if z is None:
            return None
        out.append(z)
    return out


def _clip_yz_to_radius(y: Interval, z: Interval, radius: float):
    R = float(radius)
    if not (math.isfinite(R) and R >= 0.0):
        raise ValueError("finite nonnegative yz radius required")
    ymin = V28.V27.V14._minimum_abs(y)
    zmin = V28.V27.V14._minimum_abs(z)
    if FULL.down(ymin * ymin + zmin * zmin) > FULL.up(R * R):
        return None
    ycap = FULL.up(math.sqrt(max(0.0, FULL.up(R * R) - FULL.down(zmin * zmin))))
    zcap = FULL.up(math.sqrt(max(0.0, FULL.up(R * R) - FULL.down(ymin * ymin))))
    yy = _intersect(y, Interval.outward_bounds(-ycap, ycap))
    zz = _intersect(z, Interval.outward_bounds(-zcap, zcap))
    if yy is None or zz is None:
        return None
    return yy, zz


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6,
          current_component_pieces: int = 4) -> dict:
    path = Path(domain_path).resolve()
    parent = V28.build(
        path, source_pieces=source_pieces, source_cell_index=source_cell_index,
        p_pieces=p_pieces, tangent_pieces=tangent_pieces,
        axial_pieces=axial_pieces, residual_x_pieces=residual_x_pieces,
        parallel_pieces=parallel_pieces,
        current_component_pieces=current_component_pieces)
    failures = [f"V28: {x}" for x in V28.validate(parent)]
    if parent.get("P5_SAMPLE1_SPLIT_GRAVITY_SIGNED_COMPONENTS_V28") != "PASS":
        failures.append("V28 prerequisite did not pass")

    V12D = V28.V27.V23.V22.V21B.V21.V12D
    V10 = V12D.V11.V10
    v12 = V12D.build(path, source_pieces=source_pieces,
                     source_cell_index=source_cell_index, p_pieces=p_pieces,
                     tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    core = V10.build(path, source_pieces=source_pieces,
                     source_cell_index=source_cell_index, p_pieces=p_pieces,
                     tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    failures += [f"V12D: {x}" for x in V12D.validate(v12)]
    failures += [f"V10: {x}" for x in V10.validate(core)]

    try:
        vr = _witness_row(v12); base = _witness_row(core)
        kperp = float(base["Ktheta_perpendicular_block_upper"])
        kpar = float(base["Ktheta_parallel_block_upper"])
        drho = float(vr["total_residual_perturbation_upper_mps2"])
        dk = float(vr["sample1_attitude_gain_operator_perturbation_upper"])
        rho = float(base["sample1_full_residual_norm_upper_mps2"])
        caps = _directional_perturbation_caps(
            k_perp=kperp, k_parallel=kpar, drho=drho, dk=dk, rho=rho)
        old_eta = float(parent["V12D_correction_perturbation_ball_upper_rad"])
        if caps["total_correction_perturbation_norm_upper_rad"] > FULL.up(old_eta):
            raise RuntimeError("directional V12D decomposition exceeded parent scalar perturbation")
    except Exception as exc:
        failures.append(f"directional V12D decomposition: {exc}")
        caps = {}; old_eta = math.inf

    V23 = V28.V27.V23
    p23 = V23.build(
        path, source_pieces=source_pieces, source_cell_index=source_cell_index,
        p_pieces=p_pieces, tangent_pieces=tangent_pieces,
        axial_pieces=axial_pieces, residual_x_pieces=residual_x_pieces,
        parallel_pieces=parallel_pieces,
        current_component_pieces=current_component_pieces)
    failures += [f"V23: {x}" for x in V23.validate(p23)]
    first = p23.get("first_open_current_subbox")
    if not isinstance(first, dict):
        failures.append("V23 first open current subbox missing")

    incompatible = False; joint = None
    q = geo_q = product_q = product_w = math.inf
    radial_lo = radial_hi = math.inf; narrowed = False; branches = []
    yz_source_hi = math.inf
    if isinstance(first, dict) and caps and parent.get("nominal_signed_correction_box_rad") is not None:
        nominal = [_I(x) for x in parent["nominal_signed_correction_box_rad"]]
        before = [_I(x) for x in first["joint_correction_box_rad"]]
        ex = float(caps["x_correction_perturbation_abs_upper_rad"])
        eyz = float(caps["yz_correction_perturbation_norm_upper_rad"])
        eall = float(caps["total_correction_perturbation_norm_upper_rad"])
        sx = nominal[0] + Interval.outward_bounds(-ex, ex)
        sy = nominal[1] + Interval.outward_bounds(-eyz, eyz)
        sz = nominal[2] + Interval.outward_bounds(-eyz, eyz)
        joint = _intersect_boxes(before, [sx, sy, sz])
        if joint is None:
            incompatible = True
        else:
            nominal_yz = V28.V27.V18._yz_norm_upper(nominal[1], nominal[2])
            yz_source_hi = FULL.up(nominal_yz + eyz)
            yz = _clip_yz_to_radius(joint[1], joint[2], yz_source_hi)
            if yz is None:
                incompatible = True
            else:
                joint[1], joint[2] = yz
                q = float(first["current_q_upper"])
                c = [_I(x) for x in first["q_ball_projected_current_component_box"]]
                nominal_hi = V28.V27.V14.CAYLEY1._norm_upper(nominal)
                nominal_lo = V28.V27.V14.CAYLEY2._norm_lower(nominal)
                radial_hi = min(
                    float(first["correction_radial_upper_rad"]),
                    V28.V27.V14.CAYLEY1._norm_upper(joint),
                    FULL.up(nominal_hi + eall))
                radial_lo = max(
                    float(first["correction_radial_lower_rad"]),
                    V28.V27.V14.CAYLEY2._norm_lower(joint),
                    max(0.0, FULL.down(nominal_lo - eall)))
                if radial_lo > radial_hi:
                    incompatible = True
                else:
                    geo = V28.V27.V15._geodesic_q_and_scalar_lower(q, radial_lo, radial_hi)
                    geo_q = math.inf if geo is None else float(geo[0])
                    wd, vd, branches, narrowed = V28.V27.V16.axis_cone_normalized_shipping_quaternion(
                        joint, radial_lower=radial_lo, radial_upper=radial_hi,
                        parent=V28.V27.V14D.radial_sinc_normalized_shipping_quaternion)
                    cx_min = V28.V27.V14._minimum_abs(c[0])
                    yz2 = max(0.0, FULL.up(q*q) - FULL.down(cx_min*cx_min))
                    cyz = min(FULL.up(math.sqrt(yz2)),
                               V28.V27.V18._yz_norm_upper(c[1], c[2]))
                    chart = {"cx": c[0], "cy": c[1], "cz": c[2],
                             "cyz_norm_upper": cyz}
                    parent_W = FULL.I(2.0) * wd - V28.V27.V14.CAYLEY1._dot(vd, c)
                    W, _yb, _yj = V28.V27.V18._support_product_scalar(parent_W, wd, vd, chart)
                    product_w, product_q = V28.V27.V14._qplus_from_product_scalar(q, W)

    closed = bool(incompatible or
                  (math.isfinite(geo_q) and geo_q < Q_TARGET) or
                  (math.isfinite(product_q) and product_q < Q_TARGET and product_w > 0.0))
    status = "PASS" if not failures else "NOT_ESTABLISHED"
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_DIRECTIONAL_V12D_REMAINDER_V29",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V28_split_gravity_parent_retained": True,
        "V12D_exact_deltaK_deltaR_decomposition_retained": True,
        "V10_one_plus_two_orthogonal_gain_blocks_used": True,
        "V12D_residual_perturbation_mapped_directionally": True,
        "V12D_gain_perturbation_retained_as_radial_ball": True,
        "V23_first_open_current_subbox_retained": True,
        "directional_perturbation_detail": caps,
        "previous_isotropic_V12D_correction_perturbation_upper_rad": old_eta,
        "directional_yz_source_norm_upper_rad": yz_source_hi,
        "joint_directional_correction_box_rad": None if joint is None else [x.as_list() for x in joint],
        "source_constraints_incompatible": incompatible,
        "current_q_upper": q,
        "directional_radial_lower_rad": radial_lo,
        "directional_radial_upper_rad": radial_hi,
        "axis_cone_narrowed": narrowed,
        "quaternion_branches": branches,
        "geodesic_q_upper": geo_q,
        "product_abs_W_lower": product_w,
        "product_q_upper": product_q,
        "first_open_subbox_closed_inside_q8": closed,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_DIRECTIONAL_V12D_REMAINDER_V29": status,
        "next_obligation": (
            "LIFT_V29_DIRECTIONAL_V12D_REMAINDER_OVER_ALL_V23_CURRENT_SUBBOXES"
            if closed and not failures else
            "REFINE_V12D_GAIN_PERTURBATION_DIRECTION_OR_SOURCE_BIAS_CORRELATION_AT_FIRST_OPEN_SUBBOX"),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA: f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_DIRECTIONAL_V12D_REMAINDER_V29":
        f.append("qualification mismatch")
    for k in ("source_generated_not_trajectory_fit", "V28_split_gravity_parent_retained",
              "V12D_exact_deltaK_deltaR_decomposition_retained",
              "V10_one_plus_two_orthogonal_gain_blocks_used",
              "V12D_residual_perturbation_mapped_directionally",
              "V12D_gain_perturbation_retained_as_radial_ball",
              "V23_first_open_current_subbox_retained"):
        if d.get(k) is not True: f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed", "deployed_correction_limit_increased",
              "q8_composed_here", "q8_word_promoted_here", "whole_word_promoted_here",
              "N_H_words_set_here"):
        if d.get(k) is not False: f.append(f"{k} is not false")
    caps = d.get("directional_perturbation_detail", {})
    old = float(d.get("previous_isotropic_V12D_correction_perturbation_upper_rad", math.inf))
    eall = float(caps.get("total_correction_perturbation_norm_upper_rad", math.inf))
    ex = float(caps.get("x_correction_perturbation_abs_upper_rad", math.inf))
    eyz = float(caps.get("yz_correction_perturbation_norm_upper_rad", math.inf))
    if not all(math.isfinite(x) and x >= 0.0 for x in (old, eall, ex, eyz)):
        f.append("invalid directional perturbation bounds")
    elif eall > FULL.up(old):
        f.append("directional decomposition exceeds V12D parent")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    if not d.get("source_constraints_incompatible") and d.get("joint_directional_correction_box_rad") is None:
        f.append("compatible V29 witness lost joint correction box")
    if d.get("P5_SAMPLE1_DIRECTIONAL_V12D_REMAINDER_V29") not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid V29 status")
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
        "status": d["P5_SAMPLE1_DIRECTIONAL_V12D_REMAINDER_V29"],
        "perturbation": d.get("directional_perturbation_detail"),
        "old_eta": d.get("previous_isotropic_V12D_correction_perturbation_upper_rad"),
        "yz_source": d.get("directional_yz_source_norm_upper_rad"),
        "q_current": d.get("current_q_upper"),
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
