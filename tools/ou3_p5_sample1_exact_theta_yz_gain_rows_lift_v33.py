#!/usr/bin/env python3
"""V33: exact theta-y/z gain-row perturbation over the V31 current-subbox lift.

The nominal one-plus-two attitude gain has the exact sparse row form

    K_theta = [[0,k_xy,k_xz], [g_y,0,0], [g_z,0,0]].

V32 refines the x row but keeps V29's full attitude-gain perturbation ball for
the y/z correction.  V33 uses the exact nominal y and z covariance rows from
V11's structured sample-1 covariance and applies the same V12D resolvent row by
row:

    dC_i <= dP ||H_theta|| + ||P_theta[i,:]|| dH + dP dH + dP,
    dK_i <= (dC_i + ||K_theta[i,:]|| dS) ||S'^{-1}||.

The two rowwise gain effects are combined by Euclidean support and intersected
with V29's existing full gain-perturbation parent.  The nominal residual
perturbation in the y/z correction likewise uses |g_y| and |g_z| before
intersecting with the V10 perpendicular-block parent.  X uses V32's exact
``[a,0,0]`` covariance-row refinement.  The expensive V31 numerical parent is
built only once after all three row bounds are prepared.

All source geometry, exact current residuals, radial parents, Cayley semantics,
and promotion guards remain unchanged.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_exact_theta_x_deltac_lift_v32 as V32
import ou3_p5_sample1_signed_radial_subcells_v13 as V13

DEFAULT_DOMAIN = V32.DEFAULT_DOMAIN
SCHEMA = 3300
V31 = V32.V31
FULL = V32.FULL
Q_TARGET = V32.Q_TARGET
WITNESS = V32.WITNESS


def _row_norm_upper(row) -> float:
    s = 0.0
    for x in row:
        a = x.abs_upper()
        s = FULL.up(s + FULL.up(a * a))
    return FULL.up(math.sqrt(max(0.0, s)))


def _exact_theta_yz_gain_detail(path: Path, *, source_pieces: int,
                                source_cell_index: int, p_pieces: int,
                                base: dict, vr: dict) -> dict:
    V12D = V31.V23.V22.V21B.V21.V12D
    V11 = V12D.V11
    V10 = V11.V10
    dom = json.loads(path.read_text(encoding="utf-8"))
    src, phase = V10.RG._source_phase_children(source_pieces)[source_cell_index]
    if phase != "due":
        raise RuntimeError("V33 focused refinement requires first due source cell")

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
    r = FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))[0][0]
    F, Q, _ = FULL._transition_and_Q(src, dom)
    alpha = F[15][15]; qaw = Q[15][15]

    rt = Interval.outward_bounds(*map(float, base["first_tangent_residual_magnitude_mps2"]))
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

    Pn, _Hn, _Sn = V11._nominal_sample1_matrices(
        t=t, Y=Y, p=p, r=r, g=g, alpha=alpha, qaw=qaw,
        d=d, fy=fy, fz=fz)
    py = _row_norm_upper(Pn[1][:3])
    pz = _row_norm_upper(Pn[2][:3])
    (gy, gz), (_kxy, _kxz), _gain = V13._signed_gain_components(
        a=a, Y=Y, c0=c0, alpha=alpha, qaw=qaw, b=b, bz=bz,
        det_first=det_first, d=d, fy=fy, fz=fz, r=r)
    ky = gy.abs_upper(); kz = gz.abs_upper()

    dP = float(vr["total_reduced_covariance_perturbation_upper"])
    dH = float(vr["sample1_H_perturbation_upper"])
    htheta = float(vr["sample1_Htheta_operator_upper"])
    dS = float(vr["sample1_innovation_perturbation_upper"])
    inv = float(vr["actual_innovation_inverse_operator_upper"])
    parent_dC = float(vr["sample1_attitude_cross_covariance_perturbation_upper"])
    parent_dk = float(vr["sample1_attitude_gain_operator_perturbation_upper"])

    dCy = V32._theta_x_deltac_upper(
        dP=dP, dH=dH, htheta_norm=htheta, a_row_norm=py)
    dCz = V32._theta_x_deltac_upper(
        dP=dP, dH=dH, htheta_norm=htheta, a_row_norm=pz)
    if dCy > FULL.up(parent_dC) or dCz > FULL.up(parent_dC):
        raise RuntimeError("theta-y/z Delta-C row exceeded V12D parent")
    dKy = FULL.up(FULL.up(dCy + FULL.up(ky * dS)) * inv)
    dKz = FULL.up(FULL.up(dCz + FULL.up(kz * dS)) * inv)
    if dKy > FULL.up(parent_dk) or dKz > FULL.up(parent_dk):
        raise RuntimeError("theta-y/z gain row exceeded V12D parent")

    return {
        "theta_y_nominal_attitude_covariance_row_norm_upper": py,
        "theta_z_nominal_attitude_covariance_row_norm_upper": pz,
        "nominal_theta_y_gain_row_norm_upper": ky,
        "nominal_theta_z_gain_row_norm_upper": kz,
        "theta_y_DeltaC_operator_upper": dCy,
        "theta_z_DeltaC_operator_upper": dCz,
        "V12D_full_DeltaC_operator_upper": parent_dC,
        "theta_y_gain_perturbation_operator_upper": dKy,
        "theta_z_gain_perturbation_operator_upper": dKz,
        "V12D_full_attitude_gain_perturbation_operator_upper": parent_dk,
        "theta_y_gain_row_exact_scalar_x_channel": True,
        "theta_z_gain_row_exact_scalar_x_channel": True,
    }


def _yz_refined_caps(*, base: dict, vr: dict, x_detail: dict,
                     yz_detail: dict, parent_fn) -> dict:
    parent = parent_fn(base=base, vr=vr, row_detail=x_detail)
    drho = float(vr["total_residual_perturbation_upper_mps2"])
    rho = float(base["sample1_full_residual_norm_upper_mps2"])
    rho_plus = FULL.up(rho + drho)
    ky = float(yz_detail["nominal_theta_y_gain_row_norm_upper"])
    kz = float(yz_detail["nominal_theta_z_gain_row_norm_upper"])
    dky = float(yz_detail["theta_y_gain_perturbation_operator_upper"])
    dkz = float(yz_detail["theta_z_gain_perturbation_operator_upper"])

    ry = FULL.up(ky * drho); rz = FULL.up(kz * drho)
    gy = FULL.up(dky * rho_plus); gz = FULL.up(dkz * rho_plus)
    row_residual = FULL.up(math.sqrt(FULL.up(FULL.up(ry * ry) + FULL.up(rz * rz))))
    row_gain = FULL.up(math.sqrt(FULL.up(FULL.up(gy * gy) + FULL.up(gz * gz))))
    parent_gain = float(parent["gain_perturbation_ball_upper_rad"])
    gain_yz = min(parent_gain, row_gain)
    kperp = float(base["Ktheta_perpendicular_block_upper"])
    residual_yz = min(FULL.up(kperp * drho), row_residual)
    eyz = min(float(parent["yz_correction_perturbation_norm_upper_rad"]),
              FULL.up(residual_yz + gain_yz))
    ex = float(parent["x_correction_perturbation_abs_upper_rad"])
    eall = min(float(parent["total_correction_perturbation_norm_upper_rad"]),
               FULL.up(math.sqrt(FULL.up(FULL.up(ex * ex) + FULL.up(eyz * eyz)))))

    out = dict(parent)
    out.update({
        "theta_y_residual_perturbation_abs_upper_rad": ry,
        "theta_z_residual_perturbation_abs_upper_rad": rz,
        "theta_y_gain_perturbation_ball_upper_rad": gy,
        "theta_z_gain_perturbation_ball_upper_rad": gz,
        "rowwise_yz_gain_perturbation_norm_upper_rad": gain_yz,
        "rowwise_yz_residual_perturbation_norm_upper_rad": residual_yz,
        "yz_correction_perturbation_norm_upper_rad": eyz,
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
    v12 = V12D.build(path, source_pieces=source_pieces,
                     source_cell_index=source_cell_index, p_pieces=p_pieces,
                     tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    core = V10.build(path, source_pieces=source_pieces,
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
        old_x = V31.V30._theta_x_gain_perturbation_upper(vr, base)
        if float(x_detail["theta_x_gain_perturbation_operator_upper"]) > FULL.up(
                float(old_x["theta_x_gain_perturbation_operator_upper"])):
            raise RuntimeError("V33 exact theta-x bound exceeded V30 parent")
        yz_detail = _exact_theta_yz_gain_detail(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            base=base, vr=vr)
    except Exception as exc:
        failures.append(f"V33 row refinement: {exc}")
        x_detail = yz_detail = None

    original_x = V31.V30._theta_x_gain_perturbation_upper
    original_caps = V31._refined_caps
    def exact_x(vr_row: dict, base_row: dict) -> dict:
        if x_detail is None:
            return original_x(vr_row, base_row)
        return dict(x_detail)
    def refined_caps(*, base: dict, vr: dict, row_detail: dict) -> dict:
        if yz_detail is None:
            return original_caps(base=base, vr=vr, row_detail=row_detail)
        return _yz_refined_caps(base=base, vr=vr, x_detail=row_detail,
                                yz_detail=yz_detail, parent_fn=original_caps)
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

    out = dict(parent)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_EXACT_THETA_YZ_GAIN_ROWS_LIFT_V33",
        "V31_current_subbox_lift_parent_retained": True,
        "V32_exact_theta_x_construction_retained": True,
        "exact_sparse_theta_gain_row_structure_used": True,
        "theta_yz_nominal_covariance_rows_bounded_separately": True,
        "theta_yz_gain_perturbation_rows_bounded_separately": True,
        "theta_yz_rowwise_bounds_intersect_full_operator_parent": True,
        "theta_x_exact_gain_detail": x_detail,
        "theta_yz_exact_gain_detail": yz_detail,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_EXACT_THETA_YZ_GAIN_ROWS_LIFT_V33": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_V33_ROW_DIRECTIONAL_REFINEMENT_INTO_FULL_SOURCE_CELL0_Q8_COVER"
            if parent.get("focused_first_witness_signed_subcell_closed_by_V30_lift") is True
            else "REFINE_FIRST_REMAINING_V33_SUBBOX_WITH_DIRECTIONAL_INNOVATION_DELTA_S_STRUCTURE"),
        "failures": list(dict.fromkeys(failures)),
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_EXACT_THETA_YZ_GAIN_ROWS_LIFT_V33":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V31_current_subbox_lift_parent_retained",
        "V32_exact_theta_x_construction_retained",
        "exact_sparse_theta_gain_row_structure_used",
        "theta_yz_nominal_covariance_rows_bounded_separately",
        "theta_yz_gain_perturbation_rows_bounded_separately",
        "theta_yz_rowwise_bounds_intersect_full_operator_parent",
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
    xd = d.get("theta_x_exact_gain_detail") or {}
    xdk = float(xd.get("theta_x_gain_perturbation_operator_upper", math.inf))
    xparent = float(xd.get("V12D_full_attitude_gain_perturbation_operator_upper", -math.inf))
    if not (math.isfinite(xdk) and 0.0 <= xdk <= FULL.up(xparent)):
        f.append("invalid exact theta-x gain refinement")
    gd = d.get("theta_yz_exact_gain_detail") or {}
    parent = float(gd.get("V12D_full_attitude_gain_perturbation_operator_upper", -math.inf))
    for key in ("theta_y_gain_perturbation_operator_upper",
                "theta_z_gain_perturbation_operator_upper"):
        x = float(gd.get(key, math.inf))
        if not (math.isfinite(x) and 0.0 <= x <= FULL.up(parent)):
            f.append(f"invalid {key}")
    for key in ("theta_y_gain_row_exact_scalar_x_channel",
                "theta_z_gain_row_exact_scalar_x_channel"):
        if gd.get(key) is not True:
            f.append(f"{key} is not true")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    if d.get("P5_SAMPLE1_EXACT_THETA_YZ_GAIN_ROWS_LIFT_V33") not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid V33 status")
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
        "status": d["P5_SAMPLE1_EXACT_THETA_YZ_GAIN_ROWS_LIFT_V33"],
        "theta_x_detail": d.get("theta_x_exact_gain_detail"),
        "theta_yz_detail": d.get("theta_yz_exact_gain_detail"),
        "directional_caps": d.get("directional_perturbation_detail"),
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
