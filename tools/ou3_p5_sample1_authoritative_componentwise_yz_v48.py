#!/usr/bin/env python3
"""V48: componentwise y/z correction perturbation on the V40/V45 parent.

V41's complete source-cell-0 q<8 cover leaves the authoritative first survivor
(p,t,a)=(0,0,23) at q=8.34452895146.  V42 and the later focused current-only
subdivisions show that splitting the final current box again is not the useful
next direction.  V33 also showed that a standalone theta-y/z Delta-C row bound
can exceed the already certified V12D operator parent and therefore must fail
closed.

This stage keeps both lessons.  It executes V45's authoritative V41/V18B chart
capture with V40's exact Joseph-component parent, retains V12D's full Delta-C
parent for theta-y/z, and reuses V34's safe directional first-row Delta-S
refinement for the exact sparse nominal rows

    K_theta,y = [g_y, 0, 0],   K_theta,z = [g_z, 0, 0].

Instead of duplicating one Euclidean yz perturbation radius onto both y and z,
V48 carries separate component bounds

    e_y <= |g_y| drho + ||Delta K_y|| (rho + drho),
    e_z <= |g_z| drho + ||Delta K_z|| (rho + drho),

intersects each with the existing V31 yz-norm parent, and only then reconstructs
an aggregate yz/total radius for radial and axis-cone checks.  The resulting
component box is intersected with V44's source-derived nominal correction and
with V45's authoritative first parent before the unchanged V15/V16/V18 q<8
composition is evaluated.

This is a focused proof refinement, not a theorem promotion.  It changes no
filter setting, source domain, six-radian correction limit, q<8 target, source
language, whole-word criterion, or N_H state.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_directional_innovation_row_lift_v34 as V34
import ou3_p5_sample1_first_psd_exact_joseph_components_v40 as V40
import ou3_p5_sample1_v40_split_signed_first_survivor_v44 as V44
import ou3_p5_sample1_v41_authoritative_split_signed_v45 as V45

DEFAULT_DOMAIN = V45.DEFAULT_DOMAIN
SCHEMA = 4800
Q_TARGET = 8.0
WITNESS = V45.WITNESS
FULL = V44.FULL
V12D = V40.V12D


def _I(pair):
    """Return an outward-rounded interval from a serialized endpoint pair."""
    return Interval.outward_bounds(*map(float, pair))


def _sum_up(*xs: float) -> float:
    """Sum nonnegative scalars with the proof backend's upward rounding."""
    y = 0.0
    for x in xs:
        y = FULL.up(y + float(x))
    return y


def _norm2_up(a: float, b: float) -> float:
    """Return an upward-rounded Euclidean norm of two nonnegative scalars."""
    aa = FULL.up(float(a) * float(a))
    bb = FULL.up(float(b) * float(b))
    return FULL.up(math.sqrt(FULL.up(aa + bb)))


def _norm3_up(a: float, b: float, c: float) -> float:
    """Return an upward-rounded Euclidean norm of three nonnegative scalars."""
    ab = _norm2_up(a, b)
    return _norm2_up(ab, c)


def _box_subset(child, parent) -> bool:
    """Return true when every child interval is contained in its parent."""
    return len(child) == len(parent) and all(
        x.lo >= y.lo and x.hi <= y.hi for x, y in zip(child, parent))


def _componentwise_yz_caps(*, base: dict, vr: dict, ds_detail: dict,
                           parent_caps: dict) -> dict:
    """Build safe separate y/z correction-perturbation radii.

    The nominal theta-y and theta-z gain rows use only residual component x.
    V34 certifies separate gain-row perturbation operators while retaining the
    full V12D Delta-C parent.  Each component therefore gets its own triangle
    bound.  The existing V31 yz and total radii remain intersection parents.
    """
    drho = float(vr["total_residual_perturbation_upper_mps2"])
    rho = float(base["sample1_full_residual_norm_upper_mps2"])
    rho_plus = FULL.up(rho + drho)
    ky = float(ds_detail["nominal_theta_y_gain_row_norm_upper"])
    kz = float(ds_detail["nominal_theta_z_gain_row_norm_upper"])
    dky = float(ds_detail["theta_y_gain_perturbation_intersected_upper"])
    dkz = float(ds_detail["theta_z_gain_perturbation_intersected_upper"])
    vals = (drho, rho, ky, kz, dky, dkz)
    if not all(math.isfinite(x) and x >= 0.0 for x in vals):
        raise ValueError("finite nonnegative componentwise y/z inputs required")

    ry = FULL.up(ky * drho)
    rz = FULL.up(kz * drho)
    gy = FULL.up(dky * rho_plus)
    gz = FULL.up(dkz * rho_plus)
    ey_candidate = _sum_up(ry, gy)
    ez_candidate = _sum_up(rz, gz)

    parent_yz = float(parent_caps["yz_correction_perturbation_norm_upper_rad"])
    parent_total = float(parent_caps["total_correction_perturbation_norm_upper_rad"])
    ex = float(parent_caps["x_correction_perturbation_abs_upper_rad"])
    if not all(math.isfinite(x) and x >= 0.0
               for x in (parent_yz, parent_total, ex)):
        raise ValueError("finite nonnegative V31 parent correction caps required")

    # A vector 2-norm bound is also a valid bound for either component.
    ey = min(parent_yz, ey_candidate)
    ez = min(parent_yz, ez_candidate)
    eyz = min(parent_yz, _norm2_up(ey, ez))
    eall = min(parent_total, _norm3_up(ex, ey, ez))
    return {
        "theta_y_nominal_residual_term_abs_upper_rad": ry,
        "theta_z_nominal_residual_term_abs_upper_rad": rz,
        "theta_y_gain_perturbation_term_abs_upper_rad": gy,
        "theta_z_gain_perturbation_term_abs_upper_rad": gz,
        "theta_y_component_candidate_abs_upper_rad": ey_candidate,
        "theta_z_component_candidate_abs_upper_rad": ez_candidate,
        "theta_y_component_abs_upper_rad": ey,
        "theta_z_component_abs_upper_rad": ez,
        "componentwise_yz_norm_upper_rad": eyz,
        "x_parent_abs_upper_rad": ex,
        "componentwise_total_norm_upper_rad": eall,
        "V31_parent_yz_norm_upper_rad": parent_yz,
        "V31_parent_total_norm_upper_rad": parent_total,
        "theta_y_strictly_below_duplicated_parent_yz": ey < parent_yz,
        "theta_z_strictly_below_duplicated_parent_yz": ez < parent_yz,
    }


def _build_v40_rows(path: Path, *, source_pieces: int, source_cell_index: int,
                    p_pieces: int, tangent_pieces: int,
                    axial_pieces: int) -> tuple[dict, dict, dict, dict, list[str]]:
    """Build V10/V12D rows for WITNESS with the V40 PSD helper installed."""
    failures: list[str] = []
    original = V12D._first_psd_perturbation_tangent
    V12D._first_psd_perturbation_tangent = \
        V40._first_psd_perturbation_exact_joseph_components
    try:
        v12 = V12D.build(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    finally:
        V12D._first_psd_perturbation_tangent = original
    failures += [f"V40/V12D: {x}" for x in V12D.validate(v12)]
    if V12D._first_psd_perturbation_tangent is not original:
        failures.append("V40 PSD helper was not restored")

    V10 = V12D.V11.V10
    core = V10.build(
        path, source_pieces=source_pieces,
        source_cell_index=source_cell_index, p_pieces=p_pieces,
        tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    failures += [f"V10: {x}" for x in V10.validate(core)]
    base = V44._find_row(core.get("rows", []), witness=WITNESS)
    vr = V44._find_row(v12.get("rows", []), witness=WITNESS)
    return core, v12, base, vr, failures


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    """Evaluate the V48 componentwise y/z refinement on the V45 parent."""
    path = Path(domain_path).resolve()
    failures: list[str] = []

    cap = V45._capture_authoritative_chart(
        path, source_pieces=source_pieces,
        source_cell_index=source_cell_index, p_pieces=p_pieces,
        tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
        residual_x_pieces=residual_x_pieces,
        parallel_pieces=parallel_pieces)
    chart_raw = cap["captured"]
    chart = {
        "q1": float(chart_raw["q1"]),
        "cx": _I(chart_raw["cx"]),
        "cy": _I(chart_raw["cy"]),
        "cz": _I(chart_raw["cz"]),
        "cyz_norm_upper": float(chart_raw["cyz_norm_upper"]),
    }
    chart_matches = (
        V45._matches(chart["q1"], V45.V41_Q_CURRENT)
        and V45._matches(chart["cx"].lo, V45.V41_CX[0])
        and V45._matches(chart["cx"].hi, V45.V41_CX[1])
        and V45._matches(chart["cyz_norm_upper"], V45.V41_CYZ))
    if not chart_matches:
        failures.append("authoritative chart does not reproduce V45/V41 witness")
    if cap.get("hooks_restored") is not True:
        failures.append("authoritative chart hooks were not restored")

    # V44's q chart is known not to be authoritative, but its source geometry
    # is independent of those copied q values.  Reuse only that geometry.
    v44 = V44.build(
        path, source_pieces=source_pieces,
        source_cell_index=source_cell_index, p_pieces=p_pieces,
        tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
        residual_x_pieces=residual_x_pieces,
        parallel_pieces=parallel_pieces)
    parent = v44.get("first_parent_open_subcell") or {}
    nominal_raw = v44.get("V44_nominal_signed_correction_box_rad")
    parent_caps = v44.get("V44_directional_perturbation_caps") or {}
    if not parent or not nominal_raw or not parent_caps:
        failures.append("V44 source geometry required by V48 is missing")

    _core, _v12, base, vr, row_failures = _build_v40_rows(
        path, source_pieces=source_pieces,
        source_cell_index=source_cell_index, p_pieces=p_pieces,
        tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    failures += row_failures

    try:
        ds_detail = V34._directional_delta_s_detail(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            base=base, vr=vr)
        caps = _componentwise_yz_caps(
            base=base, vr=vr, ds_detail=ds_detail,
            parent_caps=parent_caps)
    except Exception as exc:
        failures.append(f"V48 componentwise y/z construction: {exc}")
        ds_detail = None
        caps = None

    parent_eval = refined_eval = None
    joint = None
    source_box = None
    source_yz_hi = math.inf
    source_radial_lo = math.inf
    source_radial_hi = math.inf
    subset = False

    if parent and nominal_raw and parent_caps and caps is not None:
        dbox = [_I(x) for x in parent["correction_component_box_rad"]]
        nominal = [_I(x) for x in nominal_raw]
        parent_eval = V44._eval_q(
            q=chart["q1"], chart=chart, dbox=dbox,
            radial_lo=float(parent["correction_radial_lower_rad"]),
            radial_hi=float(parent["correction_radial_upper_rad"]))
        if not V45._matches(float(parent_eval["best_q"]), V45.V41_Q_POST,
                            atol=3.0e-11):
            failures.append("authoritative V48 parent does not reproduce V45 q")

        ex = float(caps["x_parent_abs_upper_rad"])
        ey = float(caps["theta_y_component_abs_upper_rad"])
        ez = float(caps["theta_z_component_abs_upper_rad"])
        eall = float(caps["componentwise_total_norm_upper_rad"])
        source_box = [
            nominal[0] + Interval.outward_bounds(-ex, ex),
            nominal[1] + Interval.outward_bounds(-ey, ey),
            nominal[2] + Interval.outward_bounds(-ez, ez),
        ]
        old_source_box = [_I(x) for x in v44["V44_source_correction_box_rad"]]
        subset = _box_subset(source_box, old_source_box)
        if not subset:
            failures.append("componentwise source box escaped V44 parent source box")

        nominal_yz = V44.V18._yz_norm_upper(nominal[1], nominal[2])
        source_yz_hi = min(
            float(v44["V44_source_yz_correction_norm_upper_rad"]),
            FULL.up(nominal_yz + float(caps["componentwise_yz_norm_upper_rad"])))
        nominal_hi = V44.V14.CAYLEY1._norm_upper(nominal)
        nominal_lo = V44.V14.CAYLEY2._norm_lower(nominal)
        source_radial_hi = min(
            float(v44["V44_source_radial_upper_rad"]),
            FULL.up(nominal_hi + eall))
        source_radial_lo = max(
            float(v44["V44_source_radial_lower_rad"]),
            max(0.0, FULL.down(nominal_lo - eall)))

        joint = V44._intersect_boxes(dbox, source_box)
        incompatible = joint is None
        if joint is not None:
            yz = V44.V31.V29._clip_yz_to_radius(
                joint[1], joint[2], source_yz_hi)
            if yz is None:
                incompatible = True
                joint = None
            else:
                joint[1], joint[2] = yz
                rhi = min(
                    float(parent["correction_radial_upper_rad"]),
                    source_radial_hi,
                    V44.V14.CAYLEY1._norm_upper(joint))
                rlo = max(
                    float(parent["correction_radial_lower_rad"]),
                    source_radial_lo,
                    V44.V14.CAYLEY2._norm_lower(joint))
                if rlo > rhi:
                    incompatible = True
                    joint = None
                else:
                    refined_eval = V44._eval_q(
                        q=chart["q1"], chart=chart, dbox=joint,
                        radial_lo=rlo, radial_hi=rhi)
                    source_radial_lo = rlo
                    source_radial_hi = rhi
        if incompatible:
            refined_eval = {
                "closed": True, "incompatible": True,
                "geodesic_q": 0.0, "product_q": 0.0,
                "product_w": math.inf, "best_q": 0.0,
                "axis_narrowed": False,
            }

    parent_q = math.inf if parent_eval is None else float(parent_eval["best_q"])
    refined_q = math.inf if refined_eval is None else float(refined_eval["best_q"])
    closed = bool(refined_eval and refined_eval.get("closed"))
    strict_component = bool(caps and (
        caps["theta_y_strictly_below_duplicated_parent_yz"]
        or caps["theta_z_strictly_below_duplicated_parent_yz"]))

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_AUTHORITATIVE_COMPONENTWISE_YZ_V48",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V40_exact_Joseph_parent_used": True,
        "V45_authoritative_chart_capture_used": True,
        "V45_authoritative_chart_matches_archived_V41_witness": chart_matches,
        "V44_source_geometry_reused_without_V44_q_values": True,
        "V34_directional_first_row_DeltaS_used": ds_detail is not None,
        "V12D_full_DeltaC_parent_retained_for_theta_yz": bool(
            ds_detail and ds_detail.get(
                "V12D_full_DeltaC_parent_retained_for_theta_yz") is True),
        "componentwise_yz_bounds_used": caps is not None,
        "componentwise_source_box_subset_of_V44_parent": subset,
        "at_least_one_yz_component_strictly_refined": strict_component,
        "V41_first_survivor_row": list(WITNESS),
        "authoritative_current_chart": chart_raw,
        "directional_DeltaS_detail": ds_detail,
        "componentwise_yz_perturbation_detail": caps,
        "authoritative_parent_best_q_upper": parent_q,
        "archived_V41_post_sample1_q_reference": V45.V41_Q_POST,
        "authoritative_componentwise_source_box_rad": (
            None if source_box is None else [x.as_list() for x in source_box]),
        "authoritative_componentwise_joint_box_rad": (
            None if joint is None else [x.as_list() for x in joint]),
        "authoritative_componentwise_yz_norm_upper_rad": source_yz_hi,
        "authoritative_componentwise_radial_lower_rad": source_radial_lo,
        "authoritative_componentwise_radial_upper_rad": source_radial_hi,
        "authoritative_componentwise_best_q_upper": refined_q,
        "authoritative_componentwise_geodesic_q_upper": (
            None if refined_eval is None else float(refined_eval["geodesic_q"])),
        "authoritative_componentwise_product_q_upper": (
            None if refined_eval is None else float(refined_eval["product_q"])),
        "authoritative_componentwise_product_abs_W_lower": (
            None if refined_eval is None else float(refined_eval["product_w"])),
        "first_V41_survivor_closed_by_V48_componentwise_yz": closed,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_established_here": False,
        "P5_SAMPLE1_AUTHORITATIVE_COMPONENTWISE_YZ_V48": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_V48_COMPONENTWISE_YZ_STRUCTURE_OVER_V41_OPEN_CELLS"
            if closed and not failures else
            "REFINE_THETA_YZ_DELTAC_COMPONENT_MATRIX_OR_RESIDUAL_X_STRUCTURE_ON_AUTHORITATIVE_V45_PARENT"
        ),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    """Validate the fail-closed V48 proof-artifact contract."""
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_AUTHORITATIVE_COMPONENTWISE_YZ_V48":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V40_exact_Joseph_parent_used",
        "V45_authoritative_chart_capture_used",
        "V45_authoritative_chart_matches_archived_V41_witness",
        "V44_source_geometry_reused_without_V44_q_values",
        "V34_directional_first_row_DeltaS_used",
        "V12D_full_DeltaC_parent_retained_for_theta_yz",
        "componentwise_yz_bounds_used",
        "componentwise_source_box_subset_of_V44_parent",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed",
        "deployed_correction_limit_increased", "q8_composed_here",
        "q8_word_promoted_here", "whole_word_promoted_here",
        "N_H_words_set_here", "P5_established_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if tuple(d.get("V41_first_survivor_row", ())) != tuple(WITNESS):
        f.append("V41 witness changed")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    if d.get("P5_SAMPLE1_AUTHORITATIVE_COMPONENTWISE_YZ_V48") not in (
            "PASS", "NOT_ESTABLISHED"):
        f.append("invalid V48 status")
    return list(dict.fromkeys(f))


def main() -> int:
    """Run V48 and write its machine-readable proof artifact."""
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
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_AUTHORITATIVE_COMPONENTWISE_YZ_V48"],
        "parent_q": d["authoritative_parent_best_q_upper"],
        "refined_q": d["authoritative_componentwise_best_q_upper"],
        "geodesic_q": d["authoritative_componentwise_geodesic_q_upper"],
        "product_q": d["authoritative_componentwise_product_q_upper"],
        "closed": d["first_V41_survivor_closed_by_V48_componentwise_yz"],
        "strict_component": d["at_least_one_yz_component_strictly_refined"],
        "caps": d["componentwise_yz_perturbation_detail"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
