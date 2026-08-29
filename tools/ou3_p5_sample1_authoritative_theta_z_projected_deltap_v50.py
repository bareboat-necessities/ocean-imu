#!/usr/bin/env python3
"""V50: componentwise theta-z projected Delta-P H_theta on V40/V45.

V49 shows that the remaining authoritative V48 survivor is limited by the
theta-z Delta-C row and, inside that row, by the scalar product

    ||Delta P|| ||H_theta||.

That scalar product discards two pieces of structure already established by
V36/V40.  The omitted first-attitude PSD remainder is a symmetric
zero-diagonal matrix

    O = [[0,a,b],[a,0,c],[b,c,0]],   |a|,|b|,|c| <= eps/2,

and the attitude part of the nominal first Joseph transport has the exact
diagonal action

    B_theta = diag(beta,beta,1),   beta=(p+R_a)/(g^2 t+p+R_a).

Hence the attitude block of B_theta O B_theta^T has the entrywise enclosure

    [[0, beta^2 e, beta e],
     [beta^2 e, 0, beta e],
     [beta e, beta e, 0]],        e=eps/2.

V50 transports that 3x3 component matrix through the already-used nominal
first reset/gauge map L_theta.  V40's gain-perturbation cross/quadratic
transport, reset-direction perturbation, and next-process remainder are kept
as scalar operator parents and added to each component.  The sample-1 accepted
S=0 covariance branch is likewise retained through its existing certified
operator parent.  Every resulting entry is intersected with V40/V12D's scalar
operator parent.

For the actual sample-1 accelerometer H_theta=-[f]_x with f_x=0, the theta-z
projected covariance row is then bounded from the specific combinations

    -f_z Delta P_zy + f_y Delta P_zz,
     f_z Delta P_zx,
    -f_y Delta P_zx,

instead of by ||Delta P|| ||H_theta||.  The remaining V49 Delta-C terms and
V34 directional Delta-S term are unchanged.  The refined theta-z gain row is
intersected with both the V34 row parent and V12D full gain parent, then injected
only into V48's componentwise y/z correction construction.  V48's authoritative
V45/V41 current chart and q composition remain unchanged.

This is a focused proof refinement.  It does not change the estimator, source
domain, six-radian correction limit, q<8 target, source language, whole-word
criterion, or N_H_words.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_authoritative_componentwise_yz_v48 as V48

DEFAULT_DOMAIN = V48.DEFAULT_DOMAIN
SCHEMA = 5000
Q_TARGET = V48.Q_TARGET
WITNESS = V48.WITNESS
FULL = V48.FULL
V12D = V48.V12D


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


def _row_norm_upper(row) -> float:
    return _norm3_up(*(x.abs_upper() for x in row))


def _nominal_first_psd_attitude_component_matrix(*, beta_upper: float,
                                                 offdiag_upper: float):
    """Absolute entrywise bound for B_theta O B_theta^T."""
    vals = (beta_upper, offdiag_upper)
    if not all(math.isfinite(float(x)) and float(x) >= 0.0 for x in vals):
        raise ValueError("finite nonnegative nominal PSD component inputs required")
    b = float(beta_upper)
    e = float(offdiag_upper)
    b2e = _mul_up(b, b, e)
    be = _mul_up(b, e)
    return [
        [0.0, b2e, be],
        [b2e, 0.0, be],
        [be, be, 0.0],
    ]


def _abs_congruence_upper(L, M):
    """Entrywise |L M L^T| upper bound for nonnegative component matrix M."""
    if len(L) != 3 or any(len(r) != 3 for r in L):
        raise ValueError("L must be 3x3")
    if len(M) != 3 or any(len(r) != 3 for r in M):
        raise ValueError("M must be 3x3")
    out = [[0.0 for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            s = 0.0
            for a in range(3):
                la = L[i][a].abs_upper()
                for b in range(3):
                    m = float(M[a][b])
                    if m < 0.0 or not math.isfinite(m):
                        raise ValueError("component matrix must be finite and nonnegative")
                    lb = L[j][b].abs_upper()
                    s = FULL.up(s + _mul_up(la, m, lb))
            out[i][j] = s
    return out


def _theta_z_projected_upper(*, dP_zx: float, dP_zy: float, dP_zz: float,
                             fy_abs: float, fz_abs: float) -> dict:
    """Bound the three components of e_z^T Delta-P H_theta^T."""
    vals = (dP_zx, dP_zy, dP_zz, fy_abs, fz_abs)
    if not all(math.isfinite(float(x)) and float(x) >= 0.0 for x in vals):
        raise ValueError("finite nonnegative theta-z projected inputs required")
    c0 = _sum_up(_mul_up(fz_abs, dP_zy), _mul_up(fy_abs, dP_zz))
    c1 = _mul_up(fz_abs, dP_zx)
    c2 = _mul_up(fy_abs, dP_zx)
    n = _norm3_up(c0, c1, c2)
    return {
        "minus_fz_DeltaP_zy_plus_fy_DeltaP_zz_abs_upper": c0,
        "fz_DeltaP_zx_abs_upper": c1,
        "minus_fy_DeltaP_zx_abs_upper": c2,
        "theta_z_projected_DeltaP_Htheta_componentwise_upper": n,
    }


def _authoritative_theta_z_detail(path: Path, *, source_pieces: int,
                                  source_cell_index: int, p_pieces: int,
                                  base: dict, vr: dict,
                                  ds_detail: dict) -> dict:
    """Construct the V50 theta-z Delta-C/gain refinement at WITNESS."""
    V11 = V12D.V11
    V10 = V11.V10
    dom = json.loads(path.read_text(encoding="utf-8"))
    src, phase = V10.RG._source_phase_children(source_pieces)[source_cell_index]
    if phase != "due":
        raise RuntimeError("V50 focused refinement requires first due source cell")

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
    beta = (p + r) / D
    beta_hi = beta.abs_upper()
    offdiag = float(vr["first_PSD_offdiagonal_entry_abs_upper"])
    nominal_first = _nominal_first_psd_attitude_component_matrix(
        beta_upper=beta_hi, offdiag_upper=offdiag)

    # V11/V4's exact nominal shipping-reset plus corrected-body gauge.
    Ltheta, _Rx = V11.V4._Ltheta(d)
    nominal_after_reset = _abs_congruence_upper(Ltheta, nominal_first)

    cross = float(vr["first_PSD_Joseph_cross_transport_upper"])
    quadratic = float(vr["first_PSD_Joseph_quadratic_transport_upper"])
    reset_direction = float(vr["reset_gauge_transform_perturbation_upper"])
    dhi = max(0.0, d.hi)
    Tnom = FULL.up(math.sqrt(FULL.up(1.0 + FULL.up(0.25 * dhi * dhi))))
    unknown_after_reset = _sum_up(
        _mul_up(Tnom, Tnom, _sum_up(cross, quadratic)),
        reset_direction,
        float(eps),
    )

    psd_parent = float(vr["PSD_reduced_covariance_perturbation_upper"])
    s_parent = float(vr["S_reduced_covariance_perturbation_upper"])
    total_parent = float(vr["total_reduced_covariance_perturbation_upper"])
    if min(psd_parent, s_parent, total_parent) < 0.0:
        raise RuntimeError("negative V40/V12D covariance perturbation parent")

    # Only the theta-z row is needed.  Unknown operator remainders and the
    # accepted S=0 branch remain scalar parents; the structured Joseph term is
    # carried entry by entry.
    entry = {}
    for name, j in (("zx", 0), ("zy", 1), ("zz", 2)):
        psd_candidate = _sum_up(nominal_after_reset[2][j], unknown_after_reset)
        psd_component = min(psd_parent, psd_candidate)
        total_candidate = _sum_up(psd_component, s_parent)
        total_component = min(total_parent, total_candidate)
        entry[name] = {
            "nominal_first_Joseph_component_upper": nominal_first[2][j],
            "nominal_after_reset_component_upper": nominal_after_reset[2][j],
            "V40_unknown_PSD_operator_remainder_upper": unknown_after_reset,
            "PSD_component_candidate_upper": psd_candidate,
            "PSD_component_intersected_upper": psd_component,
            "sample1_S_operator_parent_upper": s_parent,
            "total_component_candidate_upper": total_candidate,
            "total_component_intersected_upper": total_component,
        }

    fy = -(alpha * (p / D) * rt)
    fz = FULL.I(g) + alpha * (p / (p + r)) * rz
    fy_abs = fy.abs_upper()
    fz_abs = fz.abs_upper()
    projected = _theta_z_projected_upper(
        dP_zx=float(entry["zx"]["total_component_intersected_upper"]),
        dP_zy=float(entry["zy"]["total_component_intersected_upper"]),
        dP_zz=float(entry["zz"]["total_component_intersected_upper"]),
        fy_abs=fy_abs, fz_abs=fz_abs)

    dP = total_parent
    dH = float(vr["sample1_H_perturbation_upper"])
    htheta = float(vr["sample1_Htheta_operator_upper"])
    pz = _row_norm_upper(
        V11._nominal_sample1_matrices(
            t=t, Y=Y, p=p, r=r, g=g, alpha=alpha, qaw=qaw,
            d=d, fy=fy, fz=fz)[0][2][:3])

    projected_parent = _mul_up(dP, htheta)
    nominal_dH = _mul_up(pz, dH)
    mixed = _mul_up(dP, dH)
    theta_aw = dP
    projected_refined = float(
        projected["theta_z_projected_DeltaP_Htheta_componentwise_upper"])
    candidate_dC = _sum_up(projected_refined, nominal_dH, mixed, theta_aw)

    parent_dC = float(vr["sample1_attitude_cross_covariance_perturbation_upper"])
    dC = min(parent_dC, candidate_dC)

    inv = float(ds_detail["actual_innovation_inverse_operator_upper"])
    kz = float(ds_detail["nominal_theta_z_gain_row_norm_upper"])
    dS0 = float(ds_detail["first_measurement_row_DeltaS_intersected_upper"])
    parent_dK = float(vr["sample1_attitude_gain_operator_perturbation_upper"])
    v34_dK = float(ds_detail["theta_z_gain_perturbation_intersected_upper"])
    dK_candidate = _mul_up(_sum_up(dC, _mul_up(kz, dS0)), inv)
    dK = min(parent_dK, v34_dK, dK_candidate)

    if projected_refined > FULL.up(projected_parent):
        raise RuntimeError("componentwise theta-z projection exceeded scalar parent")
    for name in ("zx", "zy", "zz"):
        if float(entry[name]["total_component_intersected_upper"]) > FULL.up(total_parent):
            raise RuntimeError(f"theta-z {name} component escaped V12D operator parent")
    if dC > FULL.up(parent_dC):
        raise RuntimeError("theta-z Delta-C refinement escaped V12D parent")
    if dK > FULL.up(v34_dK) or dK > FULL.up(parent_dK):
        raise RuntimeError("theta-z gain refinement escaped certified parents")

    out = {
        "beta_first_attitude_transport_upper": beta_hi,
        "first_PSD_offdiagonal_entry_abs_upper": offdiag,
        "nominal_first_PSD_attitude_component_matrix_upper": nominal_first,
        "nominal_reset_Ltheta_interval": [
            [x.as_list() for x in row] for row in Ltheta],
        "nominal_after_reset_PSD_attitude_component_matrix_upper":
            nominal_after_reset,
        "V40_cross_transport_upper": cross,
        "V40_quadratic_transport_upper": quadratic,
        "V40_reset_direction_operator_upper": reset_direction,
        "next_process_PSD_remainder_upper": float(eps),
        "V40_unknown_PSD_operator_remainder_after_reset_upper":
            unknown_after_reset,
        "V40_PSD_operator_parent_upper": psd_parent,
        "sample1_S_operator_parent_upper": s_parent,
        "V12D_total_DeltaP_operator_parent_upper": total_parent,
        "theta_z_DeltaP_component_entries": entry,
        "sample1_force_y_abs_upper_mps2": fy_abs,
        "sample1_force_z_abs_upper_mps2": fz_abs,
        **projected,
        "scalar_projected_DeltaP_Htheta_parent_upper": projected_parent,
        "projected_componentwise_over_scalar_parent_ratio": (
            0.0 if projected_parent == 0.0
            else projected_refined / projected_parent),
        "theta_z_nominal_attitude_covariance_row_norm_upper": pz,
        "nominal_Ptheta_row_DeltaH_upper": nominal_dH,
        "mixed_DeltaP_DeltaH_upper": mixed,
        "theta_aw_cross_block_parent_upper": theta_aw,
        "theta_z_DeltaC_componentwise_candidate_upper": candidate_dC,
        "V12D_full_DeltaC_parent_upper": parent_dC,
        "theta_z_DeltaC_intersected_upper": dC,
        "theta_z_DeltaC_strictly_refined": dC < parent_dC,
        "directional_DeltaS_row_upper": dS0,
        "nominal_theta_z_gain_row_norm_upper": kz,
        "actual_innovation_inverse_operator_upper": inv,
        "theta_z_gain_perturbation_componentwise_candidate_upper": dK_candidate,
        "V34_theta_z_gain_perturbation_parent_upper": v34_dK,
        "V12D_full_attitude_gain_perturbation_parent_upper": parent_dK,
        "theta_z_gain_perturbation_intersected_upper": dK,
        "theta_z_gain_strictly_refined_vs_V34": dK < v34_dK,
        "specific_theta_z_projected_combinations_used": True,
        "V40_unknown_terms_retained_as_operator_parents": True,
        "sample1_S_covariance_branch_retained_as_operator_parent": True,
    }
    return out


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    path = Path(domain_path).resolve()
    failures: list[str] = []

    # V48 already constructs the authoritative V40 witness once.  Refine it
    # exactly where V48 asks for y/z correction caps instead of rebuilding the
    # same 24^3 witness a second time.
    original_caps = V48._componentwise_yz_caps
    injected = {"calls": 0, "last_detail": None}

    def v50_caps(*, base: dict, vr: dict, ds_detail: dict,
                 parent_caps: dict) -> dict:
        local_detail = _authoritative_theta_z_detail(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            base=base, vr=vr, ds_detail=ds_detail)
        injected["calls"] += 1
        injected["last_detail"] = local_detail
        ds2 = dict(ds_detail)
        ds2["theta_z_gain_perturbation_intersected_upper"] = float(
            local_detail["theta_z_gain_perturbation_intersected_upper"])
        out = original_caps(
            base=base, vr=vr, ds_detail=ds2, parent_caps=parent_caps)
        out["V50_theta_z_gain_refinement_injected"] = True
        out["V50_theta_z_gain_parent_upper"] = float(
            ds_detail["theta_z_gain_perturbation_intersected_upper"])
        out["V50_theta_z_gain_intersected_upper"] = float(
            local_detail["theta_z_gain_perturbation_intersected_upper"])
        return out

    V48._componentwise_yz_caps = v50_caps
    try:
        parent = V48.build(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces)
    except Exception as exc:
        failures.append(f"V48 under V50 theta-z injection: {exc}")
        parent = {}
    finally:
        V48._componentwise_yz_caps = original_caps

    failures += [f"V48: {x}" for x in V48.validate(parent)] if parent else []
    if parent and parent.get("P5_SAMPLE1_AUTHORITATIVE_COMPONENTWISE_YZ_V48") != "PASS":
        failures.append("V48 authoritative parent did not pass under V50 injection")
    if injected["calls"] <= 0:
        failures.append("V50 theta-z componentwise cap was not injected into V48")
    if V48._componentwise_yz_caps is not original_caps:
        failures.append("V48 componentwise cap helper was not restored")

    detail = injected["last_detail"]
    caps = parent.get("componentwise_yz_perturbation_detail") or {}
    z_gain_refined = bool(detail and detail.get(
        "theta_z_gain_strictly_refined_vs_V34"))
    z_corr_refined = bool(caps and (
        float(caps.get("theta_z_component_abs_upper_rad", math.inf))
        < float(caps.get("V31_parent_yz_norm_upper_rad", -math.inf))))
    closed = bool(parent.get("first_V41_survivor_closed_by_V48_componentwise_yz"))
    q = float(parent.get("authoritative_componentwise_best_q_upper", math.inf))
    if not math.isfinite(q):
        failures.append("V50 authoritative q is not finite")

    out = dict(parent)
    out.update({
        "schema": SCHEMA,
        "qualification":
            "OU3_P5_SAMPLE1_AUTHORITATIVE_THETA_Z_PROJECTED_DELTAP_V50",
        "V48_authoritative_componentwise_yz_parent_retained": True,
        "V40_exact_Joseph_component_parent_retained": True,
        "V45_V41_authoritative_current_chart_retained": True,
        "theta_z_projected_DeltaP_component_matrix_used": True,
        "theta_z_projected_DeltaP_specific_accelerometer_combinations_used": True,
        "V40_unknown_PSD_terms_retained_as_operator_parents": True,
        "sample1_S_covariance_branch_retained_as_operator_parent": True,
        "theta_z_componentwise_detail": detail,
        "V50_componentwise_cap_injection_calls": int(injected["calls"]),
        "V48_componentwise_cap_helper_restored": (
            V48._componentwise_yz_caps is original_caps),
        "theta_z_gain_strictly_refined_vs_V34": z_gain_refined,
        "theta_z_correction_component_strictly_refined": z_corr_refined,
        "authoritative_V50_best_q_upper": q,
        "first_V41_survivor_closed_by_V50_theta_z_projected_DeltaP": closed,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_AUTHORITATIVE_THETA_Z_PROJECTED_DELTAP_V50": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "PROMOTE_V50_CLOSED_SAMPLE1_CELL_THROUGH_FULL_SOURCE_CELL0_COVER"
            if closed and not failures else
            "REFINE_THETA_Z_REMAINING_V50_DELTAC_OR_RESIDUAL_TERM_ON_AUTHORITATIVE_PARENT"
            if z_gain_refined and not failures else
            "REFINE_V40_UNKNOWN_PSD_OR_SAMPLE1_S_COMPONENT_MATRIX_FOR_THETA_Z"),
        "failures": list(dict.fromkeys(failures)),
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != \
            "OU3_P5_SAMPLE1_AUTHORITATIVE_THETA_Z_PROJECTED_DELTAP_V50":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V48_authoritative_componentwise_yz_parent_retained",
        "V40_exact_Joseph_component_parent_retained",
        "V45_V41_authoritative_current_chart_retained",
        "theta_z_projected_DeltaP_component_matrix_used",
        "theta_z_projected_DeltaP_specific_accelerometer_combinations_used",
        "V40_unknown_PSD_terms_retained_as_operator_parents",
        "sample1_S_covariance_branch_retained_as_operator_parent",
        "V48_componentwise_cap_helper_restored",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed",
        "deployed_correction_limit_increased", "q8_word_promoted_here",
        "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if list(d.get("V41_first_survivor_row", [])) != list(WITNESS):
        f.append("authoritative V41 survivor changed")
    if int(d.get("V50_componentwise_cap_injection_calls", 0)) <= 0:
        f.append("V50 cap injection was not exercised")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")

    detail = d.get("theta_z_componentwise_detail") or {}
    projected = float(detail.get(
        "theta_z_projected_DeltaP_Htheta_componentwise_upper", math.inf))
    projected_parent = float(detail.get(
        "scalar_projected_DeltaP_Htheta_parent_upper", -math.inf))
    if not (math.isfinite(projected) and 0.0 <= projected
            <= FULL.up(projected_parent)):
        f.append("invalid componentwise theta-z projected Delta-P refinement")
    dC = float(detail.get("theta_z_DeltaC_intersected_upper", math.inf))
    dC_parent = float(detail.get("V12D_full_DeltaC_parent_upper", -math.inf))
    if not (math.isfinite(dC) and 0.0 <= dC <= FULL.up(dC_parent)):
        f.append("invalid theta-z Delta-C intersection")
    dK = float(detail.get("theta_z_gain_perturbation_intersected_upper", math.inf))
    dK_parent = float(detail.get(
        "V34_theta_z_gain_perturbation_parent_upper", -math.inf))
    if not (math.isfinite(dK) and 0.0 <= dK <= FULL.up(dK_parent)):
        f.append("invalid theta-z gain intersection")

    q = float(d.get("authoritative_V50_best_q_upper", math.inf))
    if not math.isfinite(q):
        f.append("nonfinite V50 authoritative q")
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
    z = d.get("theta_z_componentwise_detail") or {}
    caps = d.get("componentwise_yz_perturbation_detail") or {}
    print(json.dumps({
        "status": d["P5_SAMPLE1_AUTHORITATIVE_THETA_Z_PROJECTED_DELTAP_V50"],
        "projected_componentwise":
            z.get("theta_z_projected_DeltaP_Htheta_componentwise_upper"),
        "projected_parent":
            z.get("scalar_projected_DeltaP_Htheta_parent_upper"),
        "projected_ratio":
            z.get("projected_componentwise_over_scalar_parent_ratio"),
        "dC": z.get("theta_z_DeltaC_intersected_upper"),
        "dC_parent": z.get("V12D_full_DeltaC_parent_upper"),
        "dKz": z.get("theta_z_gain_perturbation_intersected_upper"),
        "dKz_parent": z.get("V34_theta_z_gain_perturbation_parent_upper"),
        "ez": caps.get("theta_z_component_abs_upper_rad"),
        "q": d.get("authoritative_V50_best_q_upper"),
        "closed":
            d.get("first_V41_survivor_closed_by_V50_theta_z_projected_DeltaP"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
