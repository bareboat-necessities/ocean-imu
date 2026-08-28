#!/usr/bin/env python3
"""V34: directional first-row innovation perturbation for theta-y/z gain rows.

V33 tried to bound the theta-y and theta-z Delta-C rows independently.  That
candidate is not uniformly smaller than V12D's certified full attitude-block
Delta-C parent, so V33 correctly failed closed.  V34 keeps that parent.

The exact nominal one-plus-two attitude gain nevertheless has

    K_theta,y = [g_y, 0, 0],   K_theta,z = [g_z, 0, 0].

Hence the resolvent terms K_theta,y Delta-S and K_theta,z Delta-S depend only on
the first measurement row of Delta-S, not on its full operator norm.

For S(H,P)=H P H^T + R, write H'=H+E and P'=P+D with
||E||<=dH and ||D||<=dP.  For nominal measurement row h_i and row
q_i = h_i P,

  ||e_i^T Delta-S|| <=
      ||h_i|| dP ||H||
    + dH ||P|| ||H||
    + ||q_i|| dH
    + dH dP ||H||
    + ||h_i|| dP dH
    + dH ||P|| dH
    + dH dP dH.

This is the seven-term expansion of
(H+E)(P+D)(H+E)^T - HPH^T, bounded row by row.  V34 intersects that
first-row bound with V12D's existing full ||Delta-S|| parent, keeps V32's exact
theta-x Delta-C refinement, and uses V12D's full Delta-C parent for theta-y/z.
The resulting y/z gain-row bounds are themselves intersected with V12D's full
attitude-gain perturbation parent before the unchanged V31 64-current-subbox
composition.

No estimator setting, source domain, six-radian correction limit, q<8 target,
or theorem-promotion state is changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_mul
import ou3_p5_sample1_exact_theta_x_deltac_lift_v32 as V32
import ou3_p5_sample1_signed_radial_subcells_v13 as V13

DEFAULT_DOMAIN = V32.DEFAULT_DOMAIN
SCHEMA = 3400
V31 = V32.V31
FULL = V32.FULL
Q_TARGET = V32.Q_TARGET
WITNESS = V32.WITNESS


def _sum_up(*xs: float) -> float:
    y = 0.0
    for x in xs:
        y = FULL.up(y + float(x))
    return y


def _row_norm_upper(row) -> float:
    s = 0.0
    for x in row:
        a = x.abs_upper()
        s = FULL.up(s + FULL.up(a * a))
    return FULL.up(math.sqrt(max(0.0, s)))


def _innovation_row_perturbation_upper(*, h_row_norm: float, h_norm: float,
                                       hp_row_norm: float, p_norm: float,
                                       dP: float, dH: float) -> float:
    vals = (h_row_norm, h_norm, hp_row_norm, p_norm, dP, dH)
    if not all(math.isfinite(float(x)) and float(x) >= 0.0 for x in vals):
        raise ValueError("finite nonnegative innovation-row inputs required")
    return _sum_up(
        FULL.up(h_row_norm * dP * h_norm),
        FULL.up(dH * p_norm * h_norm),
        FULL.up(hp_row_norm * dH),
        FULL.up(dH * dP * h_norm),
        FULL.up(h_row_norm * dP * dH),
        FULL.up(dH * p_norm * dH),
        FULL.up(dH * dP * dH),
    )


def _directional_delta_s_detail(path: Path, *, source_pieces: int,
                                source_cell_index: int, p_pieces: int,
                                base: dict, vr: dict) -> dict:
    V12D = V31.V23.V22.V21B.V21.V12D
    V11 = V12D.V11
    V10 = V11.V10
    dom = json.loads(path.read_text(encoding="utf-8"))
    src, phase = V10.RG._source_phase_children(source_pieces)[source_cell_index]
    if phase != "due":
        raise RuntimeError("V34 focused refinement requires first due source cell")

    first = V11.FIRST.build(path, source_pieces=source_pieces)
    fr = first["source_cells"][source_cell_index]
    p_all = Interval.outward_bounds(*map(float, fr["P_aw_variance_interval"]))
    pcells = V11.SUB.parts(p_all.lo, p_all.hi, p_pieces)
    p = pcells[int(base["p_cell"])]

    hstep = float(src["dt_s"])
    g = float(dom["startup"]["gravity_mps2"])
    tilt, yaw, eps = V10.RG._attitude_covariance_epsilon(path, hstep)
    t = Interval.outward_bounds(tilt, FULL.up(tilt + eps))
    Y = Interval.outward_bounds(yaw, FULL.up(yaw + eps))
    vec = V11.VECTOR.build()
    r = FULL._R_diag(
        float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"])
    )[0][0]
    F, Q, _ = FULL._transition_and_Q(src, dom)
    alpha = F[15][15]
    qaw = Q[15][15]

    rt = Interval.outward_bounds(
        *map(float, base["first_tangent_residual_magnitude_mps2"]))
    rz = Interval.outward_bounds(*map(float, base["first_axial_residual_mps2"]))
    d = Interval.outward_bounds(*map(float, base["first_attitude_correction_rad"]))
    D = FULL.I(g * g) * t + p + r
    a = t * (p + r) / D
    c0 = -(FULL.I(g) * t * p / D)
    b = p * (FULL.I(g * g) * t + r) / D
    bz = p * r / (p + r)
    det_first = t * p * r / D
    fy = -(alpha * (p / D) * rt)
    fz = FULL.I(g) + alpha * (p / (p + r)) * rz

    Pn, Hn, _Sn = V11._nominal_sample1_matrices(
        t=t, Y=Y, p=p, r=r, g=g, alpha=alpha, qaw=qaw,
        d=d, fy=fy, fz=fz)
    HP = matrix_mul(Hn, Pn)
    h_row = _row_norm_upper(Hn[0])
    h_norm = V11._op(Hn)
    hp_row = _row_norm_upper(HP[0])
    p_norm = V11._op(Pn)

    dP = float(vr["total_reduced_covariance_perturbation_upper"])
    dH = float(vr["sample1_H_perturbation_upper"])
    parent_dS = float(vr["sample1_innovation_perturbation_upper"])
    parent_dC = float(vr["sample1_attitude_cross_covariance_perturbation_upper"])
    parent_dk = float(vr["sample1_attitude_gain_operator_perturbation_upper"])
    inv = float(vr["actual_innovation_inverse_operator_upper"])

    candidate_dS0 = _innovation_row_perturbation_upper(
        h_row_norm=h_row, h_norm=h_norm, hp_row_norm=hp_row,
        p_norm=p_norm, dP=dP, dH=dH)
    dS0 = min(parent_dS, candidate_dS0)

    (gy, gz), (_kxy, _kxz), _gain = V13._signed_gain_components(
        a=a, Y=Y, c0=c0, alpha=alpha, qaw=qaw, b=b, bz=bz,
        det_first=det_first, d=d, fy=fy, fz=fz, r=r)
    ky = gy.abs_upper()
    kz = gz.abs_upper()

    dKy_candidate = FULL.up(
        FULL.up(parent_dC + FULL.up(ky * dS0)) * inv)
    dKz_candidate = FULL.up(
        FULL.up(parent_dC + FULL.up(kz * dS0)) * inv)
    dKy = min(parent_dk, dKy_candidate)
    dKz = min(parent_dk, dKz_candidate)

    return {
        "first_measurement_H_row_norm_upper": h_row,
        "nominal_H_operator_upper": h_norm,
        "first_measurement_HP_row_norm_upper": hp_row,
        "nominal_P_operator_upper": p_norm,
        "V12D_full_innovation_perturbation_upper": parent_dS,
        "first_measurement_row_DeltaS_candidate_upper": candidate_dS0,
        "first_measurement_row_DeltaS_intersected_upper": dS0,
        "first_measurement_row_DeltaS_strictly_refined": dS0 < parent_dS,
        "V12D_full_DeltaC_operator_upper": parent_dC,
        "V12D_full_attitude_gain_perturbation_operator_upper": parent_dk,
        "actual_innovation_inverse_operator_upper": inv,
        "nominal_theta_y_gain_row_norm_upper": ky,
        "nominal_theta_z_gain_row_norm_upper": kz,
        "theta_y_gain_perturbation_candidate_upper": dKy_candidate,
        "theta_z_gain_perturbation_candidate_upper": dKz_candidate,
        "theta_y_gain_perturbation_intersected_upper": dKy,
        "theta_z_gain_perturbation_intersected_upper": dKz,
        "theta_y_gain_row_exact_scalar_x_channel": True,
        "theta_z_gain_row_exact_scalar_x_channel": True,
        "V12D_full_DeltaC_parent_retained_for_theta_yz": True,
    }


def _directional_caps(*, base: dict, vr: dict, x_detail: dict,
                      ds_detail: dict, parent_fn) -> dict:
    parent = parent_fn(base=base, vr=vr, row_detail=x_detail)
    drho = float(vr["total_residual_perturbation_upper_mps2"])
    rho = float(base["sample1_full_residual_norm_upper_mps2"])
    rho_plus = FULL.up(rho + drho)
    ky = float(ds_detail["nominal_theta_y_gain_row_norm_upper"])
    kz = float(ds_detail["nominal_theta_z_gain_row_norm_upper"])
    dky = float(ds_detail["theta_y_gain_perturbation_intersected_upper"])
    dkz = float(ds_detail["theta_z_gain_perturbation_intersected_upper"])

    ry = FULL.up(ky * drho)
    rz = FULL.up(kz * drho)
    gy = FULL.up(dky * rho_plus)
    gz = FULL.up(dkz * rho_plus)
    row_residual = FULL.up(math.sqrt(FULL.up(
        FULL.up(ry * ry) + FULL.up(rz * rz))))
    row_gain = FULL.up(math.sqrt(FULL.up(
        FULL.up(gy * gy) + FULL.up(gz * gz))))

    kperp = float(base["Ktheta_perpendicular_block_upper"])
    residual_parent = FULL.up(kperp * drho)
    gain_parent = float(parent["gain_perturbation_ball_upper_rad"])
    residual_yz = min(residual_parent, row_residual)
    gain_yz = min(gain_parent, row_gain)
    eyz_candidate = FULL.up(residual_yz + gain_yz)
    eyz = min(float(parent["yz_correction_perturbation_norm_upper_rad"]),
              eyz_candidate)

    ex = float(parent["x_correction_perturbation_abs_upper_rad"])
    eall_candidate = FULL.up(math.sqrt(FULL.up(
        FULL.up(ex * ex) + FULL.up(eyz * eyz))))
    eall = min(float(parent["total_correction_perturbation_norm_upper_rad"]),
               eall_candidate)

    out = dict(parent)
    out.update({
        "theta_y_residual_perturbation_abs_upper_rad": ry,
        "theta_z_residual_perturbation_abs_upper_rad": rz,
        "theta_y_gain_perturbation_ball_upper_rad": gy,
        "theta_z_gain_perturbation_ball_upper_rad": gz,
        "rowwise_yz_residual_perturbation_norm_upper_rad": residual_yz,
        "rowwise_yz_gain_perturbation_norm_upper_rad": gain_yz,
        "yz_correction_perturbation_candidate_norm_upper_rad": eyz_candidate,
        "yz_correction_perturbation_norm_upper_rad": eyz,
        "total_correction_perturbation_candidate_norm_upper_rad": eall_candidate,
        "total_correction_perturbation_norm_upper_rad": eall,
    })
    return out


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6,
          current_component_pieces: int = 4) -> dict:
    path = Path(domain_path).resolve()
    V12D = V31.V23.V22.V21B.V21.V12D
    V10 = V12D.V11.V10

    v12 = V12D.build(
        path, source_pieces=source_pieces,
        source_cell_index=source_cell_index, p_pieces=p_pieces,
        tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    core = V10.build(
        path, source_pieces=source_pieces,
        source_cell_index=source_cell_index, p_pieces=p_pieces,
        tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    failures = [f"V12D: {x}" for x in V12D.validate(v12)]
    failures += [f"V10: {x}" for x in V10.validate(core)]

    try:
        vr = V31.V30._witness_row(v12)
        base = V31.V30._witness_row(core)
        x_detail = V32._exact_theta_x_gain_detail(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            base=base, vr=vr)
        ds_detail = _directional_delta_s_detail(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            base=base, vr=vr)
    except Exception as exc:
        failures.append(f"V34 directional Delta-S construction: {exc}")
        x_detail = ds_detail = None

    original_x = V31.V30._theta_x_gain_perturbation_upper
    original_caps = V31._refined_caps

    def exact_x(vr_row: dict, base_row: dict) -> dict:
        if x_detail is None:
            return original_x(vr_row, base_row)
        return dict(x_detail)

    def refined_caps(*, base: dict, vr: dict, row_detail: dict) -> dict:
        if ds_detail is None:
            return original_caps(base=base, vr=vr, row_detail=row_detail)
        return _directional_caps(
            base=base, vr=vr, x_detail=row_detail,
            ds_detail=ds_detail, parent_fn=original_caps)

    V31.V30._theta_x_gain_perturbation_upper = exact_x
    V31._refined_caps = refined_caps
    try:
        parent = V31.build(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces,
            current_component_pieces=current_component_pieces)
    finally:
        V31._refined_caps = original_caps
        V31.V30._theta_x_gain_perturbation_upper = original_x

    failures += [f"V31: {x}" for x in V31.validate(parent)]
    if parent.get("P5_SAMPLE1_V30_CURRENT_SUBBOX_LIFT_V31") != "PASS":
        failures.append("V31 current-subbox parent did not pass")

    open_count = int(parent.get("open_current_subboxes", -1))
    strict = bool((ds_detail or {}).get(
        "first_measurement_row_DeltaS_strictly_refined", False))
    out = dict(parent)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_DIRECTIONAL_INNOVATION_ROW_LIFT_V34",
        "V31_current_subbox_lift_parent_retained": True,
        "V32_exact_theta_x_construction_retained": True,
        "V33_invalid_theta_yz_DeltaC_route_retired": True,
        "V12D_full_DeltaC_parent_retained_for_theta_yz": True,
        "exact_sparse_theta_yz_gain_rows_used": True,
        "first_measurement_DeltaS_row_bounded_from_exact_nominal_HP": True,
        "first_measurement_DeltaS_row_intersected_with_V12D_parent": True,
        "theta_yz_gain_row_bounds_intersected_with_V12D_parent": True,
        "theta_yz_correction_bounds_intersected_with_V29_parent": True,
        "directional_innovation_detail": ds_detail,
        "theta_x_exact_gain_detail": x_detail,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_DIRECTIONAL_INNOVATION_ROW_LIFT_V34": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_V34_DIRECTIONAL_INNOVATION_REFINEMENT_INTO_FULL_SOURCE_CELL0_Q8_COVER"
            if open_count == 0 and not failures else
            "REFINE_REMAINING_V34_SUBBOX_WITH_BLOCK_RESOLVENT_ACTION_OR_SOURCE_SPLIT"
            if strict and not failures else
            "REFINE_V12D_SAMPLE1_S_PSD_COMPONENT_STRUCTURE_AT_FIRST_OPEN_SUBBOX"
        ),
        "failures": list(dict.fromkeys(failures)),
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_DIRECTIONAL_INNOVATION_ROW_LIFT_V34":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V31_current_subbox_lift_parent_retained",
        "V32_exact_theta_x_construction_retained",
        "V33_invalid_theta_yz_DeltaC_route_retired",
        "V12D_full_DeltaC_parent_retained_for_theta_yz",
        "exact_sparse_theta_yz_gain_rows_used",
        "first_measurement_DeltaS_row_bounded_from_exact_nominal_HP",
        "first_measurement_DeltaS_row_intersected_with_V12D_parent",
        "theta_yz_gain_row_bounds_intersected_with_V12D_parent",
        "theta_yz_correction_bounds_intersected_with_V29_parent",
        "all_candidate_current_subboxes_accounted",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed",
        "deployed_correction_limit_increased", "q8_composed_here",
        "q8_word_promoted_here", "whole_word_promoted_here",
        "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")

    gd = d.get("directional_innovation_detail") or {}
    parent_ds = float(gd.get("V12D_full_innovation_perturbation_upper", math.inf))
    row_ds = float(gd.get("first_measurement_row_DeltaS_intersected_upper", math.inf))
    parent_dk = float(gd.get("V12D_full_attitude_gain_perturbation_operator_upper", math.inf))
    dky = float(gd.get("theta_y_gain_perturbation_intersected_upper", math.inf))
    dkz = float(gd.get("theta_z_gain_perturbation_intersected_upper", math.inf))
    if not (math.isfinite(row_ds) and 0.0 <= row_ds <= FULL.up(parent_ds)):
        f.append("invalid first-measurement Delta-S row refinement")
    if not (math.isfinite(dky) and 0.0 <= dky <= FULL.up(parent_dk)):
        f.append("invalid theta-y gain-row refinement")
    if not (math.isfinite(dkz) and 0.0 <= dkz <= FULL.up(parent_dk)):
        f.append("invalid theta-z gain-row refinement")
    if gd.get("theta_y_gain_row_exact_scalar_x_channel") is not True:
        f.append("theta-y scalar-x gain-row flag missing")
    if gd.get("theta_z_gain_row_exact_scalar_x_channel") is not True:
        f.append("theta-z scalar-x gain-row flag missing")
    if gd.get("V12D_full_DeltaC_parent_retained_for_theta_yz") is not True:
        f.append("theta-y/z full Delta-C parent not retained")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    if d.get("P5_SAMPLE1_DIRECTIONAL_INNOVATION_ROW_LIFT_V34") not in (
            "PASS", "NOT_ESTABLISHED"):
        f.append("invalid V34 status")
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
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_DIRECTIONAL_INNOVATION_ROW_LIFT_V34"],
        "directional_innovation": d.get("directional_innovation_detail"),
        "candidate": d.get("candidate_current_subboxes"),
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
