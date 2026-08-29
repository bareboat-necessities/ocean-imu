#!/usr/bin/env python3
"""V50: correlated theta-z projected Delta-P H_theta refinement.

V49 identifies ``||Delta P|| ||H_theta||`` as the dominant theta-z Delta-C
term.  A naive entrywise replacement is not automatically tighter: treating
DeltaP_zx, DeltaP_zy and DeltaP_zz as independent can destroy the shared
matrix correlation and exceed the scalar operator parent.

V50 therefore carries the structured first-PSD attitude remainder through the
actual matrix algebra before taking a norm.  V36/V40 establish that the
omitted first-attitude remainder is

    O = a E01 + b E02 + c E12,   |a|,|b|,|c| <= eps/2,

with symmetric basis matrices Eij.  The nominal first Joseph attitude action is
B_theta=diag(beta,beta,1).  With the already-certified nominal reset/gauge map
L_theta, each shared basis variable is transported as

    M_ij = L_theta B_theta Eij B_theta L_theta^T.

For the actual sample-1 accelerometer H_theta=-[f]_x, f_x=0, V50 evaluates the
specific theta-z projected combinations

    -f_z M_zy + f_y M_zz,
     f_z M_zx,
    -f_y M_zx

for each basis and sums the three basis-vector norms with outward rounding.
This preserves the correlation within each matrix basis.  V40's gain-error
cross/quadratic transport, reset-direction perturbation and next-process
remainder are still charged in operator norm, as is the accepted sample-1
S=0 covariance branch.  The resulting candidate is intersected with the V49
scalar parent, so the refinement is fail-safe and never weakens the existing
certificate.

This stage only proves or rejects the projected Delta-P refinement.  It does
not alter the estimator, source domain, six-radian correction limit, q<8
target, whole-word language, P5 status or N_H_words.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_mul, matrix_transpose
import ou3_p5_sample1_authoritative_componentwise_yz_v48 as V48

DEFAULT_DOMAIN = V48.DEFAULT_DOMAIN
SCHEMA = 5001
WITNESS = V48.WITNESS
FULL = V48.FULL
V12D = V48.V12D
Q_TARGET = V48.Q_TARGET


def _sum_up(*xs: float) -> float:
    y = 0.0
    for x in xs:
        y = FULL.up(y + float(x))
    return y


def _mul_up(*xs: float) -> float:
    y = 1.0
    for x in xs:
        y = FULL.up(y * float(x))
    return y


def _norm3_up(a: float, b: float, c: float) -> float:
    s = FULL.up(FULL.up(a * a) + FULL.up(b * b))
    s = FULL.up(s + FULL.up(c * c))
    return FULL.up(math.sqrt(max(0.0, s)))


def _basis(i: int, j: int):
    """Return the symmetric 3x3 E_ij basis matrix."""
    if not (0 <= i < j < 3):
        raise ValueError("basis requires 0 <= i < j < 3")
    z = FULL.I(0.0)
    o = FULL.I(1.0)
    M = [[z for _ in range(3)] for _ in range(3)]
    M[i][j] = o
    M[j][i] = o
    return M


def _project_theta_z_row(M, *, fy_abs: float, fz_abs: float) -> dict:
    """Project row z of M through the actual f_x=0 H_theta structure."""
    zx = M[2][0].abs_upper()
    zy = M[2][1].abs_upper()
    zz = M[2][2].abs_upper()
    c0 = _sum_up(_mul_up(fz_abs, zy), _mul_up(fy_abs, zz))
    c1 = _mul_up(fz_abs, zx)
    c2 = _mul_up(fy_abs, zx)
    return {
        "minus_fz_Mzy_plus_fy_Mzz_abs_upper": c0,
        "fz_Mzx_abs_upper": c1,
        "minus_fy_Mzx_abs_upper": c2,
        "vector_norm_upper": _norm3_up(c0, c1, c2),
    }


def _component_detail(path: Path, *, source_pieces: int,
                      source_cell_index: int, p_pieces: int,
                      base: dict, vr: dict) -> dict:
    """Evaluate the correlated V50 projection on the authoritative witness."""
    V11 = V12D.V11
    V10 = V11.V10
    dom = json.loads(path.read_text(encoding="utf-8"))
    src, phase = V10.RG._source_phase_children(source_pieces)[source_cell_index]
    if phase != "due":
        raise RuntimeError("V50 requires the authoritative first-due source cell")

    first = V11.FIRST.build(path, source_pieces=source_pieces)
    fr = first["source_cells"][source_cell_index]
    p_all = Interval.outward_bounds(*map(float, fr["P_aw_variance_interval"]))
    p = V11.SUB.parts(p_all.lo, p_all.hi, p_pieces)[int(base["p_cell"])]

    h = float(src["dt_s"])
    g = float(dom["startup"]["gravity_mps2"])
    tilt, yaw, eps = V10.RG._attitude_covariance_epsilon(path, h)
    t = Interval.outward_bounds(tilt, FULL.up(tilt + eps))
    vec = V11.VECTOR.build()
    r = FULL._R_diag(float(
        vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))[0][0]
    F, _Q, _ = FULL._transition_and_Q(src, dom)
    alpha = F[15][15]

    rt = Interval.outward_bounds(
        *map(float, base["first_tangent_residual_magnitude_mps2"]))
    rz = Interval.outward_bounds(*map(float, base["first_axial_residual_mps2"]))
    d = Interval.outward_bounds(*map(float, base["first_attitude_correction_rad"]))
    D = FULL.I(g * g) * t + p + r
    beta = (p + r) / D
    fy = -(alpha * (p / D) * rt)
    fz = FULL.I(g) + alpha * (p / (p + r)) * rz
    fy_abs = fy.abs_upper()
    fz_abs = fz.abs_upper()

    L, _Rx = V11.V4._Ltheta(d)
    z = FULL.I(0.0)
    B = [[beta, z, z], [z, beta, z], [z, z, FULL.I(1.0)]]
    A = matrix_mul(L, B)

    offdiag = float(vr["first_PSD_offdiagonal_entry_abs_upper"])
    if not (math.isfinite(offdiag) and offdiag >= 0.0):
        raise RuntimeError("invalid V40 offdiagonal component bound")

    basis_rows = []
    nominal_projected = 0.0
    for i, j in ((0, 1), (0, 2), (1, 2)):
        M = matrix_mul(matrix_mul(A, _basis(i, j)), matrix_transpose(A))
        projected = _project_theta_z_row(M, fy_abs=fy_abs, fz_abs=fz_abs)
        contribution = _mul_up(offdiag, float(projected["vector_norm_upper"]))
        nominal_projected = FULL.up(nominal_projected + contribution)
        basis_rows.append({
            "basis": [i, j],
            "unit_basis_projected_theta_z": projected,
            "offdiagonal_amplitude_upper": offdiag,
            "projected_contribution_upper": contribution,
        })

    cross = float(vr["first_PSD_Joseph_cross_transport_upper"])
    quadratic = float(vr["first_PSD_Joseph_quadratic_transport_upper"])
    reset_direction = float(vr["reset_gauge_transform_perturbation_upper"])
    dhi = max(0.0, d.hi)
    Tnom = FULL.up(math.sqrt(FULL.up(1.0 + FULL.up(0.25 * dhi * dhi))))
    unknown_psd = _sum_up(
        _mul_up(Tnom, Tnom, _sum_up(cross, quadratic)),
        reset_direction,
        float(eps),
    )
    s_parent = float(vr["S_reduced_covariance_perturbation_upper"])
    total_parent = float(vr["total_reduced_covariance_perturbation_upper"])
    htheta = float(vr["sample1_Htheta_operator_upper"])
    if not all(math.isfinite(x) and x >= 0.0 for x in
               (unknown_psd, s_parent, total_parent, htheta)):
        raise RuntimeError("invalid V40/V12D scalar remainder parent")

    scalar_projected_parent = _mul_up(total_parent, htheta)
    operator_remainder = _sum_up(unknown_psd, s_parent)
    remainder_projected = _mul_up(operator_remainder, htheta)
    correlated_candidate = _sum_up(nominal_projected, remainder_projected)
    intersected = min(scalar_projected_parent, correlated_candidate)

    return {
        "first_PSD_offdiagonal_entry_abs_upper": offdiag,
        "beta_interval": beta.as_list(),
        "sample1_force_y_abs_upper_mps2": fy_abs,
        "sample1_force_z_abs_upper_mps2": fz_abs,
        "basis_projection_detail": basis_rows,
        "nominal_correlated_PSD_projected_upper": nominal_projected,
        "V40_cross_transport_operator_upper": cross,
        "V40_quadratic_transport_operator_upper": quadratic,
        "V40_reset_direction_operator_upper": reset_direction,
        "next_process_PSD_remainder_upper": float(eps),
        "V40_unknown_PSD_operator_remainder_upper": unknown_psd,
        "sample1_S_operator_parent_upper": s_parent,
        "combined_unknown_and_S_operator_remainder_upper": operator_remainder,
        "sample1_Htheta_operator_upper": htheta,
        "V12D_total_DeltaP_operator_parent_upper": total_parent,
        "scalar_projected_DeltaP_Htheta_parent_upper": scalar_projected_parent,
        "correlated_component_matrix_projected_candidate_upper": correlated_candidate,
        "theta_z_projected_DeltaP_Htheta_intersected_upper": intersected,
        "candidate_over_scalar_parent_ratio": (
            math.inf if scalar_projected_parent == 0.0
            else correlated_candidate / scalar_projected_parent),
        "intersected_over_scalar_parent_ratio": (
            0.0 if scalar_projected_parent == 0.0
            else intersected / scalar_projected_parent),
        "strict_projected_refinement": intersected < scalar_projected_parent,
        "specific_theta_z_accelerometer_combinations_used": True,
        "shared_offdiagonal_basis_correlations_preserved": True,
        "V40_unknown_terms_retained_as_operator_parent": True,
        "sample1_S_branch_retained_as_operator_parent": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24) -> dict:
    """Build the focused authoritative V50 witness artifact."""
    path = Path(domain_path).resolve()
    failures: list[str] = []
    detail = None
    try:
        _core, _v12, base, vr, row_failures = V48._build_v40_rows(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
        failures += row_failures
        ids = (int(base.get("p_cell", -1)),
               int(base.get("tangent_residual_cell", -1)),
               int(base.get("axial_residual_cell", -1)))
        if ids != tuple(WITNESS):
            failures.append("V50 authoritative witness changed")
        detail = _component_detail(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            base=base, vr=vr)
    except Exception as exc:
        failures.append(f"V50 correlated projected Delta-P construction: {exc}")

    strict = bool(detail and detail.get("strict_projected_refinement"))
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_AUTHORITATIVE_CORRELATED_THETA_Z_PROJECTED_DELTAP_V50",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V40_exact_Joseph_parent_used": True,
        "V49_scalar_projected_parent_replaced_only_by_intersection": True,
        "shared_matrix_basis_correlations_preserved": True,
        "V41_first_survivor_row": list(WITNESS),
        "theta_z_componentwise_detail": detail,
        "strict_projected_refinement": strict,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_established_here": False,
        "P5_SAMPLE1_AUTHORITATIVE_THETA_Z_PROJECTED_DELTAP_V50": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "INJECT_V50_CORRELATED_THETA_Z_PROJECTED_DELTAP_INTO_AUTHORITATIVE_V48_Q_COMPOSITION"
            if strict and not failures else
            "REFINE_SAMPLE1_S_AND_V40_UNKNOWN_PSD_COMPONENT_CORRELATIONS_BEFORE_THETA_Z_PROJECTION"
        ),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    """Validate the V50 fail-closed proof contract."""
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != \
            "OU3_P5_SAMPLE1_AUTHORITATIVE_CORRELATED_THETA_Z_PROJECTED_DELTAP_V50":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V40_exact_Joseph_parent_used",
        "V49_scalar_projected_parent_replaced_only_by_intersection",
        "shared_matrix_basis_correlations_preserved",
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
        f.append("authoritative witness changed")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")

    z = d.get("theta_z_componentwise_detail") or {}
    parent = float(z.get("scalar_projected_DeltaP_Htheta_parent_upper", math.inf))
    candidate = float(z.get("correlated_component_matrix_projected_candidate_upper", math.inf))
    intersected = float(z.get("theta_z_projected_DeltaP_Htheta_intersected_upper", math.inf))
    if not all(math.isfinite(x) and x >= 0.0 for x in
               (parent, candidate, intersected)):
        f.append("nonfinite V50 projected bound")
    elif intersected > FULL.up(parent) or intersected > FULL.up(candidate):
        f.append("V50 projected intersection escaped a parent")
    for k in (
        "specific_theta_z_accelerometer_combinations_used",
        "shared_offdiagonal_basis_correlations_preserved",
        "V40_unknown_terms_retained_as_operator_parent",
        "sample1_S_branch_retained_as_operator_parent",
    ):
        if z.get(k) is not True:
            f.append(f"theta-z detail {k} is not true")
    if d.get("P5_SAMPLE1_AUTHORITATIVE_THETA_Z_PROJECTED_DELTAP_V50") not in (
            "PASS", "NOT_ESTABLISHED"):
        f.append("invalid V50 status")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=24)
    ap.add_argument("--tangent-pieces", type=int, default=24)
    ap.add_argument("--axial-pieces", type=int, default=24)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(
        x.domain, source_pieces=x.source_pieces,
        source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
        tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    z = d.get("theta_z_componentwise_detail") or {}
    print(json.dumps({
        "status": d["P5_SAMPLE1_AUTHORITATIVE_THETA_Z_PROJECTED_DELTAP_V50"],
        "scalar_parent": z.get("scalar_projected_DeltaP_Htheta_parent_upper"),
        "candidate": z.get("correlated_component_matrix_projected_candidate_upper"),
        "intersected": z.get("theta_z_projected_DeltaP_Htheta_intersected_upper"),
        "candidate_ratio": z.get("candidate_over_scalar_parent_ratio"),
        "intersected_ratio": z.get("intersected_over_scalar_parent_ratio"),
        "nominal_structured": z.get("nominal_correlated_PSD_projected_upper"),
        "unknown_psd_operator": z.get("V40_unknown_PSD_operator_remainder_upper"),
        "S_operator": z.get("sample1_S_operator_parent_upper"),
        "strict": d.get("strict_projected_refinement"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
