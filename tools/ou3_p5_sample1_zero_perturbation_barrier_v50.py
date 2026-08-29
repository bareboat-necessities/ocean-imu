#!/usr/bin/env python3
"""V50: attitude-supported Delta-S refinement and zero-perturbation barrier.

V29 through V49 all attack the same object: the perturbation caps that widen the
nominal sample-1 correction box away from V10's exact directional correction.
V48 splits those caps componentwise, V49 decomposes the theta-y/z Delta-C
candidate.  Neither moved the authoritative V41/V45 first survivor, which is
still at q=8.3445...

V50 first supplies the strongest remaining refinement on that route, and then
tests whether the route can close the witness at all.

Part 1 - attitude-supported first-row Delta-S.
----------------------------------------------
In the sample-1 body gauge the shipping accelerometer Jacobian is exactly

    H = [ -[f]_x | I ],

because the source block `J_aw = R_wb()` is orthogonal and the gauge places it on
the identity while `J_att = -skew(f_cog_b)` carries the whole modelled
dependence.  The Jacobian perturbation is therefore supported on the attitude
columns,

    E = Delta H = [ E_theta | 0 ].

That is exactly the fact which lets the certified V12C/V12D parent write

    dC = dP ||H_theta|| + ||P_theta|| dH + dP dH + dP

with no `P_theta,aw Delta H_aw` term, so V50 reproduces that parent bit-for-bit
from those four terms before reusing the same support anywhere else.

V34 bounds the first measurement row of Delta-S by seven terms.  The three that
carry E next to a nominal factor each collapse onto the attitude block:

    e_i^T H P E^T = (h_i P)_theta E_theta^T   -> ||(h_i P)_theta|| dH,
    e_i^T E P H^T = u^T (P H^T)_theta         -> dH ||(P H^T)_theta||,
    e_i^T E P E^T = u^T P_theta,theta E_theta^T -> dH ||P_theta,theta|| dH,

replacing V34's `||h_i P||`, `||P|| ||H||` and `||P||` factors.  The four terms
carrying Delta P are unchanged.  Every refined term is taken as a minimum
against its own V34 parent term and the total is still intersected with V12D's
full `||Delta S||` parent, so the refinement is monotone and fail-closed by
construction.

Part 2 - what the remaining perturbation budget is made of.
----------------------------------------------------------
After the refinement the dominant Delta-S term is `||h_i|| dP ||H||` and the
dominant Delta-C term is `dP ||H_theta||`.  Both are proportional to the single
certified sample-1 reduced covariance perturbation `dP`, so V50 decomposes `dP`
into the four constituents V40/V12C add together, reproducing the certified
parent exactly.

Part 3 - the barrier.
---------------------
V50 then reruns the authoritative V48 composition with every componentwise
correction-perturbation cap forced to zero.  This is not a proof claim about the
filter: it is an audit of the proof route.  If the resulting q is still at or
above the q<8 target, then no refinement of the perturbation bounds that
composition consumes - V34's Delta-S, V48's componentwise split, V49's Delta-C
terms, or Part 1 above - can close the authoritative first survivor.  The
verdict is taken from the actual composed q, covering the geodesic and the
product branch alike, not from the geodesic branch alone.

The composed radial bound is additionally intersected with V44's parent open
subcell, whose own radius already carries the V29/V31 caps.  V50 therefore also
reports the geodesic branch evaluated with the entire certified perturbation
budget subtracted from that parent radius too.  That second number is a route
counterfactual rather than an enclosure and is labelled a diagnostic, but it
answers the same question for the part of the radial bound the zeroed run
cannot reach.

This changes no filter setting, source domain, six-radian correction limit, q<8
target, source language, whole-word criterion, or `N_H` state.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_mul, matrix_transpose
import ou3_p5_sample1_authoritative_componentwise_yz_v48 as V48
import ou3_p5_sample1_directional_innovation_row_lift_v34 as V34
import ou3_p5_sample1_signed_radial_subcells_v13 as V13

DEFAULT_DOMAIN = V48.DEFAULT_DOMAIN
SCHEMA = 5000
Q_TARGET = V48.Q_TARGET
WITNESS = V48.WITNESS
FULL = V48.FULL
V12D = V48.V12D
V11 = V12D.V11
V10 = V11.V10

_ZEROED_CAP_KEYS = (
    "x_parent_abs_upper_rad",
    "theta_y_component_abs_upper_rad",
    "theta_z_component_abs_upper_rad",
    "componentwise_yz_norm_upper_rad",
    "componentwise_total_norm_upper_rad",
)


def _sum_up(*xs: float) -> float:
    """Sum nonnegative scalars with the proof backend's upward rounding."""
    y = 0.0
    for x in xs:
        y = FULL.up(y + float(x))
    return y


def _finite_nonneg(*xs: float) -> bool:
    return all(math.isfinite(float(x)) and float(x) >= 0.0 for x in xs)


def _row_terms(*, h_row_norm: float, h_norm: float, hp_row_norm: float,
               hp_row_theta_norm: float, p_norm: float, p_theta_norm: float,
               c_theta_norm: float, dP: float, dH: float) -> dict:
    """Return V34's seven Delta-S row terms and their attitude-supported form.

    Terms carrying `Delta P` next to a nominal factor are kept exactly as V34
    states them.  The three carrying `Delta H` collapse onto attitude-restricted
    nominal factors.  Each refined term is a minimum against its own parent
    term, so no term can grow.
    """
    vals = (h_row_norm, h_norm, hp_row_norm, hp_row_theta_norm, p_norm,
            p_theta_norm, c_theta_norm, dP, dH)
    if not _finite_nonneg(*vals):
        raise ValueError("finite nonnegative Delta-S row inputs required")

    # Each expression matches V34's own rounding form term by term, so the
    # parent total reproduces V34's row candidate exactly.
    parent = {
        "hi_dP_H": FULL.up(h_row_norm * dP * h_norm),
        "dH_P_H": FULL.up(dH * p_norm * h_norm),
        "hiP_dH": FULL.up(hp_row_norm * dH),
        "dH_dP_H": FULL.up(dH * dP * h_norm),
        "hi_dP_dH": FULL.up(h_row_norm * dP * dH),
        "dH_P_dH": FULL.up(dH * p_norm * dH),
        "dH_dP_dH": FULL.up(dH * dP * dH),
    }
    supported = dict(parent)
    supported["dH_P_H"] = FULL.up(dH * c_theta_norm)
    supported["hiP_dH"] = FULL.up(hp_row_theta_norm * dH)
    supported["dH_P_dH"] = FULL.up(dH * p_theta_norm * dH)

    order = ("hi_dP_H", "dH_P_H", "hiP_dH", "dH_dP_H", "hi_dP_dH",
             "dH_P_dH", "dH_dP_dH")
    refined = {k: min(parent[k], supported[k]) for k in order}
    parent_total = _sum_up(*(parent[k] for k in order))
    refined_total = _sum_up(*(refined[k] for k in order))
    dominant = max(order, key=lambda k: refined[k])
    return {
        "term_order": list(order),
        "V34_parent_terms": parent,
        "attitude_supported_terms": supported,
        "refined_terms": refined,
        "V34_parent_row_candidate_upper": parent_total,
        "attitude_supported_row_candidate_upper": refined_total,
        "row_candidate_strictly_refined": refined_total < parent_total,
        "refined_dominant_term": dominant,
        "refined_dominant_fraction": (
            0.0 if refined_total == 0.0 else refined[dominant] / refined_total),
        "attitude_supported_terms_never_exceed_parent": all(
            refined[k] <= parent[k] for k in order),
    }


def _reduced_covariance_terms(*, vr: dict, dhi: float, eps: float) -> dict:
    """Decompose the certified sample-1 reduced covariance perturbation `dP`.

    V40 forms the PSD contribution as `T_nom^2 dP_plus + dirterm + eps` and
    V12C adds the S-update contribution.  V50 recomputes that sum in V40's own
    expression order and requires it to reproduce the stored parent exactly, so
    the reported fractions decompose the certified quantity instead of
    estimating it.
    """
    dPplus = float(vr["first_posterior_covariance_perturbation_upper"])
    dirterm = float(vr["reset_gauge_transform_perturbation_upper"])
    psd = float(vr["sample1_reduced_covariance_PSD_perturbation_upper"])
    s_part = float(vr["S_reduced_covariance_perturbation_upper"])
    dP = float(vr["total_reduced_covariance_perturbation_upper"])
    if not _finite_nonneg(dPplus, dirterm, psd, s_part, dP, dhi, eps):
        raise ValueError("finite nonnegative reduced-covariance inputs required")

    tnom = FULL.up(math.sqrt(FULL.up(1.0 + FULL.up(0.25 * dhi * dhi))))
    transported = FULL.up(tnom * tnom * dPplus)
    psd_check = FULL.up(FULL.up(transported + dirterm) + eps)
    total_check = FULL.up(psd_check + s_part)
    terms = {
        "transported_first_posterior_upper": transported,
        "reset_gauge_direction_upper": dirterm,
        "sample1_prediction_attitude_epsilon_upper": eps,
        "S_update_reduced_covariance_upper": s_part,
    }
    dominant = max(terms, key=lambda k: terms[k])
    return {
        "first_posterior_covariance_perturbation_upper": dPplus,
        "reset_congruence_T_nom_upper": tnom,
        **terms,
        "reconstructed_PSD_reduced_covariance_upper": psd_check,
        "certified_PSD_reduced_covariance_upper": psd,
        "reconstructed_total_reduced_covariance_upper": total_check,
        "certified_total_reduced_covariance_upper": dP,
        "PSD_decomposition_reproduces_certified_parent": psd_check == psd,
        "total_decomposition_reproduces_certified_parent": total_check == dP,
        "dominant_term": dominant,
        "dominant_fraction": 0.0 if dP == 0.0 else terms[dominant] / dP,
    }


_NEXT_BY_DP_TERM = {
    "transported_first_posterior_upper":
        "REFINE_TRANSPORTED_FIRST_POSTERIOR_PSD_REMAINDER_COMPONENT_MATRIX_ON_AUTHORITATIVE_V40_PARENT",
    "reset_gauge_direction_upper":
        "REFINE_RESET_GAUGE_DIRECTION_COVARIANCE_TERM_ON_AUTHORITATIVE_V40_PARENT",
    "sample1_prediction_attitude_epsilon_upper":
        "REFINE_SAMPLE1_PREDICTION_ATTITUDE_COVARIANCE_EPSILON_ON_AUTHORITATIVE_V40_PARENT",
    "S_update_reduced_covariance_upper":
        "REFINE_SAMPLE1_S_UPDATE_REDUCED_COVARIANCE_PERTURBATION_ON_AUTHORITATIVE_V40_PARENT",
}

_BARRIER_OBLIGATION = (
    "REFINE_NOMINAL_V10_FIRST_ACCELEROMETER_CORRECTION_MAGNITUDE_OR_"
    "SAMPLE0_CURRENT_CHART_ON_AUTHORITATIVE_V45_PARENT")
_REPAIR_OBLIGATION = "REPAIR_V50_ZERO_PERTURBATION_BARRIER"


def _attitude_supported_detail(path: Path, *, source_pieces: int,
                               source_cell_index: int, p_pieces: int,
                               base: dict, vr: dict) -> dict:
    """Build the refined first-row Delta-S and theta-y/z gain-row bounds."""
    dom = json.loads(path.read_text(encoding="utf-8"))
    src, phase = V10.RG._source_phase_children(source_pieces)[source_cell_index]
    if phase != "due":
        raise RuntimeError("V50 focused refinement requires first due source cell")

    first = V11.FIRST.build(path, source_pieces=source_pieces)
    fr = first["source_cells"][source_cell_index]
    p_all = Interval.outward_bounds(*map(float, fr["P_aw_variance_interval"]))
    p = V11.SUB.parts(p_all.lo, p_all.hi, p_pieces)[int(base["p_cell"])]

    hstep = float(src["dt_s"])
    g = float(dom["startup"]["gravity_mps2"])
    tilt, yaw, eps = V10.RG._attitude_covariance_epsilon(path, hstep)
    t = Interval.outward_bounds(tilt, FULL.up(tilt + eps))
    Y = Interval.outward_bounds(yaw, FULL.up(yaw + eps))
    vec = V11.VECTOR.build()
    r = FULL._R_diag(float(
        vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))[0][0]
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
        t=t, Y=Y, p=p, r=r, g=g, alpha=alpha, qaw=qaw, d=d, fy=fy, fz=fz)
    HP = matrix_mul(Hn, Pn)
    C = matrix_mul(Pn, matrix_transpose(Hn))

    h_row = V34._row_norm_upper(Hn[0])
    h_norm = V11._op(Hn)
    hp_row = V34._row_norm_upper(HP[0])
    hp_row_theta = V34._row_norm_upper(HP[0][:3])
    p_norm = V11._op(Pn)
    p_theta_norm = V11._op([row[:3] for row in Pn[:3]])
    c_theta_norm = V11._op([row[:] for row in C[:3]])
    htheta_norm = V11.V5._norm2_upper(fy.abs_upper(), fz.abs_upper())

    dP = float(vr["total_reduced_covariance_perturbation_upper"])
    dH = float(vr["sample1_H_perturbation_upper"])
    parent_dS = float(vr["sample1_innovation_perturbation_upper"])
    parent_dC = float(vr["sample1_attitude_cross_covariance_perturbation_upper"])
    parent_dk = float(vr["sample1_attitude_gain_operator_perturbation_upper"])
    inv = float(vr["actual_innovation_inverse_operator_upper"])

    # The certified V12C Delta-C parent is exactly the four-term attitude
    # supported expansion.  Reproducing it bit-for-bit binds V50's use of
    # Delta H_aw = 0 to the parent that already relies on it.
    deltac_check = _sum_up(
        FULL.up(dP * htheta_norm),
        FULL.up(p_theta_norm * dH),
        FULL.up(dP * dH),
        dP,
    )

    v34_candidate = V34._innovation_row_perturbation_upper(
        h_row_norm=h_row, h_norm=h_norm, hp_row_norm=hp_row,
        p_norm=p_norm, dP=dP, dH=dH)
    terms = _row_terms(
        h_row_norm=h_row, h_norm=h_norm, hp_row_norm=hp_row,
        hp_row_theta_norm=hp_row_theta, p_norm=p_norm,
        p_theta_norm=p_theta_norm, c_theta_norm=c_theta_norm,
        dP=dP, dH=dH)
    v34_dS0 = min(parent_dS, float(terms["V34_parent_row_candidate_upper"]))
    dS0 = min(v34_dS0, float(terms["attitude_supported_row_candidate_upper"]))

    (gy, gz), (_kxy, _kxz), _gain = V13._signed_gain_components(
        a=a, Y=Y, c0=c0, alpha=alpha, qaw=qaw, b=b, bz=bz,
        det_first=det_first, d=d, fy=fy, fz=fz, r=r)

    rows = {}
    for label, k in (("theta_y", gy.abs_upper()), ("theta_z", gz.abs_upper())):
        parent_term = FULL.up(k * v34_dS0)
        refined_term = FULL.up(k * dS0)
        refined_numerator = FULL.up(parent_dC + refined_term)
        rows[label] = {
            "nominal_gain_row_norm_upper": k,
            "V34_DeltaS_term_upper": parent_term,
            "refined_DeltaS_term_upper": refined_term,
            "refined_numerator_upper": refined_numerator,
            "DeltaS_share_of_refined_numerator": (
                0.0 if refined_numerator == 0.0
                else refined_term / refined_numerator),
            "V34_gain_perturbation_candidate_upper": FULL.up(
                FULL.up(parent_dC + parent_term) * inv),
            "refined_gain_perturbation_candidate_upper": FULL.up(
                refined_numerator * inv),
            "V12D_full_gain_perturbation_parent_upper": parent_dk,
        }
        refined_candidate = float(
            rows[label]["refined_gain_perturbation_candidate_upper"])
        rows[label]["gain_perturbation_intersected_upper"] = min(
            parent_dk, refined_candidate)
        rows[label]["candidate_unpins_from_V12D_parent"] = (
            refined_candidate < parent_dk)
        rows[label]["candidate_over_V12D_parent_ratio"] = (
            math.inf if parent_dk <= 0.0 else refined_candidate / parent_dk)

    return {
        "first_measurement_H_row_norm_upper": h_row,
        "nominal_H_operator_upper": h_norm,
        "first_measurement_HP_row_norm_upper": hp_row,
        "first_measurement_HP_row_attitude_norm_upper": hp_row_theta,
        "nominal_P_operator_upper": p_norm,
        "nominal_P_attitude_block_operator_upper": p_theta_norm,
        "nominal_PH_transpose_attitude_rows_operator_upper": c_theta_norm,
        "nominal_Htheta_operator_upper": htheta_norm,
        "reduced_covariance_perturbation_upper": dP,
        "H_perturbation_upper": dH,
        "V12D_full_innovation_perturbation_upper": parent_dS,
        "V12D_full_DeltaC_operator_upper": parent_dC,
        "reconstructed_four_term_DeltaC_upper": deltac_check,
        "DeltaC_parent_is_exact_attitude_supported_expansion": deltac_check == parent_dC,
        "V34_first_measurement_row_DeltaS_candidate_upper": v34_candidate,
        "V34_first_measurement_row_DeltaS_intersected_upper": v34_dS0,
        "attitude_supported_row_DeltaS_intersected_upper": dS0,
        "row_DeltaS_strictly_refined": dS0 < v34_dS0,
        "row_DeltaS_refinement_ratio": (
            math.inf if dS0 == 0.0 else v34_dS0 / dS0),
        "V34_candidate_reproduced": (
            v34_candidate == float(terms["V34_parent_row_candidate_upper"])),
        "attitude_covariance_epsilon": eps,
        "first_attitude_correction_upper_rad": max(0.0, d.hi),
        **terms,
        "gain_rows": rows,
    }


def _zeroed_caps(caps: dict) -> dict:
    """Return V48's cap dictionary with every perturbation radius set to zero."""
    out = dict(caps)
    for key in _ZEROED_CAP_KEYS:
        if key not in out:
            raise KeyError(f"V48 cap dictionary is missing {key}")
        out[key] = 0.0
    out["theta_y_strictly_below_duplicated_parent_yz"] = True
    out["theta_z_strictly_below_duplicated_parent_yz"] = True
    return out


def _angle_diagnostics(*, q_current: float, radial_lower: float,
                       radial_upper: float, cap_total: float) -> dict:
    """Explain the barrier on the geodesic branch in principal-angle terms.

    The composed q comes from the audited run.  These numbers only say how the
    geodesic branch would move if the entire certified correction-perturbation
    budget were removed from the radial bound as well, which is a route
    counterfactual rather than an enclosure, and is reported as a diagnostic.
    """
    if not _finite_nonneg(q_current, radial_lower, radial_upper, cap_total):
        raise ValueError("finite nonnegative angle diagnostics required")
    phi_c = 2.0 * math.atan(0.5 * q_current)
    phi_target = 2.0 * math.atan(0.5 * Q_TARGET)
    needed = (phi_c + radial_upper) - phi_target
    free_upper = max(radial_lower, radial_upper - cap_total)
    geo = V48.V44.V15._geodesic_q_and_scalar_lower(
        q_current, radial_lower, free_upper)
    free_q = math.inf if geo is None else float(geo[0])
    return {
        "current_chart_principal_angle_rad": phi_c,
        "correction_principal_angle_rad": radial_upper,
        "q_target_principal_angle_rad": phi_target,
        "principal_angle_reduction_needed_rad": needed,
        "principal_angle_reduction_available_from_caps_rad": cap_total,
        "cap_share_of_needed_reduction": (
            math.inf if needed <= 0.0 else cap_total / needed),
        "cap_free_correction_principal_angle_rad": free_upper,
        "cap_free_geodesic_q_diagnostic_upper": free_q,
        "cap_free_geodesic_q_still_at_or_above_target": free_q >= Q_TARGET,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    """Refine the first-row Delta-S and audit the perturbation route barrier."""
    path = Path(domain_path).resolve()
    failures: list[str] = []
    detail = dp_terms = angles = None
    barrier = None
    hooks_restored = False
    next_obligation = _REPAIR_OBLIGATION

    context: dict = {"rows": None, "caps": None, "cap_calls": 0}
    root_rows = V48._build_v40_rows
    root_caps = V48._componentwise_yz_caps

    def rows_hook(*args, **kwargs):
        out = root_rows(*args, **kwargs)
        if context["rows"] is None:
            context["rows"] = out
        return out

    def caps_hook(*, base, vr, ds_detail, parent_caps):
        real = root_caps(base=base, vr=vr, ds_detail=ds_detail,
                         parent_caps=parent_caps)
        context["cap_calls"] += 1
        if context["caps"] is None:
            context["caps"] = real
        return _zeroed_caps(real)

    try:
        V48._build_v40_rows = rows_hook
        V48._componentwise_yz_caps = caps_hook
        try:
            zero = V48.build(
                path, source_pieces=source_pieces,
                source_cell_index=source_cell_index, p_pieces=p_pieces,
                tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
                residual_x_pieces=residual_x_pieces,
                parallel_pieces=parallel_pieces)
        finally:
            V48._componentwise_yz_caps = root_caps
            V48._build_v40_rows = root_rows
        hooks_restored = (V48._build_v40_rows is root_rows
                          and V48._componentwise_yz_caps is root_caps)
        if not hooks_restored:
            failures.append("V50 temporary V48 hooks were not restored")
        if context["cap_calls"] != 1:
            failures.append(
                f"V48 cap helper ran {context['cap_calls']} times, expected 1")

        failures += [f"V48 zero-perturbation run: {x}"
                     for x in zero.get("failures", [])]
        if zero.get("V45_authoritative_chart_matches_archived_V41_witness") is not True:
            failures.append("zero-perturbation run lost the authoritative chart")
        if tuple(zero.get("V41_first_survivor_row", ())) != tuple(WITNESS):
            failures.append("zero-perturbation run changed the V41 witness")
        if zero.get("componentwise_source_box_subset_of_V44_parent") is not True:
            failures.append("zero-perturbation source box escaped the V44 parent")

        real_caps = context["caps"] or {}
        joint = zero.get("authoritative_componentwise_joint_box_rad")
        zero_q = float(zero.get("authoritative_componentwise_best_q_upper",
                                math.inf))
        parent_q = float(zero.get("authoritative_parent_best_q_upper", math.inf))
        compatible = joint is not None
        can_close = compatible and math.isfinite(zero_q) and zero_q < Q_TARGET
        barrier = {
            "authoritative_parent_best_q_upper": parent_q,
            "zero_perturbation_best_q_upper": zero_q,
            "zero_perturbation_geodesic_q_upper": zero.get(
                "authoritative_componentwise_geodesic_q_upper"),
            "zero_perturbation_product_q_upper": zero.get(
                "authoritative_componentwise_product_q_upper"),
            "zero_perturbation_radial_lower_rad": zero.get(
                "authoritative_componentwise_radial_lower_rad"),
            "zero_perturbation_radial_upper_rad": zero.get(
                "authoritative_componentwise_radial_upper_rad"),
            "zero_perturbation_joint_box_nonempty": compatible,
            "zero_perturbation_q_equals_parent_q": zero_q == parent_q,
            "q_target": Q_TARGET,
            "real_componentwise_caps": real_caps,
            "perturbation_route_can_close_authoritative_witness": can_close,
            "barrier_established": not can_close,
        }
        if not compatible:
            failures.append(
                "zero-perturbation correction box is incompatible with the V44 "
                "parent subcell; the barrier audit is inconclusive")

        rows_out = context["rows"]
        if rows_out is None:
            failures.append("V50 never observed the V40 witness rows")
        else:
            _core, _v12, base, vr, row_failures = rows_out
            failures += row_failures
            if (int(base.get("p_cell", -1)),
                int(base.get("tangent_residual_cell", -1)),
                int(base.get("axial_residual_cell", -1))) != WITNESS:
                failures.append("V50 did not reconstruct authoritative witness")

            detail = _attitude_supported_detail(
                path, source_pieces=source_pieces,
                source_cell_index=source_cell_index, p_pieces=p_pieces,
                base=base, vr=vr)
            if not detail["DeltaC_parent_is_exact_attitude_supported_expansion"]:
                failures.append(
                    "V12D Delta-C parent is not the four-term attitude-supported "
                    "expansion; V50 may not assume Delta H_aw = 0")
            if not detail["V34_candidate_reproduced"]:
                failures.append("V50 did not reproduce V34's row Delta-S candidate")
            if not detail["attitude_supported_terms_never_exceed_parent"]:
                failures.append("attitude-supported term exceeded its V34 parent")
            if float(detail["attitude_supported_row_DeltaS_intersected_upper"]) > \
                    float(detail["V34_first_measurement_row_DeltaS_intersected_upper"]):
                failures.append("refined row Delta-S exceeded its V34 parent")

            dp_terms = _reduced_covariance_terms(
                vr=vr, dhi=float(detail["first_attitude_correction_upper_rad"]),
                eps=float(detail["attitude_covariance_epsilon"]))
            if not dp_terms["total_decomposition_reproduces_certified_parent"]:
                failures.append(
                    "reduced covariance decomposition did not reproduce the "
                    "certified V12D/V40 parent")

            angles = _angle_diagnostics(
                q_current=float(
                    (zero.get("authoritative_current_chart") or {}).get("q1", 0.0)),
                radial_lower=float(zero.get(
                    "authoritative_componentwise_radial_lower_rad", 0.0)),
                radial_upper=float(zero.get(
                    "authoritative_componentwise_radial_upper_rad", 0.0)),
                cap_total=float(real_caps.get(
                    "componentwise_total_norm_upper_rad", 0.0)))

            next_obligation = (
                _BARRIER_OBLIGATION if barrier["barrier_established"]
                else _NEXT_BY_DP_TERM[dp_terms["dominant_term"]])
    except Exception as exc:
        V48._componentwise_yz_caps = root_caps
        V48._build_v40_rows = root_rows
        failures.append(f"V50 zero-perturbation barrier: {exc}")
        next_obligation = _REPAIR_OBLIGATION

    rows = (detail or {}).get("gain_rows") or {}
    unpinned = sorted(k for k, v in rows.items()
                      if v.get("candidate_unpins_from_V12D_parent"))
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_ZERO_PERTURBATION_BARRIER_V50",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "attitude_supported_Jacobian_perturbation_used": True,
        "V12D_full_DeltaC_parent_retained": True,
        "V12D_full_DeltaS_parent_retained_as_intersection": True,
        "V34_seven_term_row_expansion_retained": True,
        "temporary_V48_hooks_restored": hooks_restored,
        "zero_perturbation_run_is_route_audit_not_filter_claim": True,
        "failed_V33_row_candidate_promoted": False,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "V41_first_survivor_row": list(WITNESS),
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_established_here": False,
        "attitude_supported_row_detail": detail,
        "reduced_covariance_perturbation_decomposition": dp_terms,
        "zero_perturbation_barrier": barrier,
        "geodesic_principal_angle_diagnostics": angles,
        "theta_gain_rows_unpinned_from_V12D_parent": unpinned,
        "P5_SAMPLE1_ZERO_PERTURBATION_BARRIER_V50": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": next_obligation,
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    """Validate the fail-closed V50 proof-artifact contract."""
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_ZERO_PERTURBATION_BARRIER_V50":
        f.append("qualification mismatch")
    for k in ("source_generated_not_trajectory_fit",
              "attitude_supported_Jacobian_perturbation_used",
              "V12D_full_DeltaC_parent_retained",
              "V12D_full_DeltaS_parent_retained_as_intersection",
              "V34_seven_term_row_expansion_retained",
              "temporary_V48_hooks_restored",
              "zero_perturbation_run_is_route_audit_not_filter_claim"):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed",
              "failed_V33_row_candidate_promoted",
              "deployed_correction_limit_increased", "q8_composed_here",
              "q8_word_promoted_here", "whole_word_promoted_here",
              "N_H_words_set_here", "P5_established_here"):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if tuple(d.get("V41_first_survivor_row", ())) != tuple(WITNESS):
        f.append("V41 witness changed")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")

    detail = d.get("attitude_supported_row_detail") or {}
    parent = float(detail.get(
        "V34_first_measurement_row_DeltaS_intersected_upper", -math.inf))
    refined = float(detail.get(
        "attitude_supported_row_DeltaS_intersected_upper", math.inf))
    if not (math.isfinite(parent) and parent >= 0.0
            and math.isfinite(refined) and 0.0 <= refined <= parent):
        f.append("invalid attitude-supported row Delta-S")
    if detail.get("DeltaC_parent_is_exact_attitude_supported_expansion") is not True:
        f.append("Delta-C parent binding is not established")

    dp = d.get("reduced_covariance_perturbation_decomposition") or {}
    if dp.get("total_decomposition_reproduces_certified_parent") is not True:
        f.append("reduced covariance decomposition is not bound to its parent")

    barrier = d.get("zero_perturbation_barrier") or {}
    zero_q = float(barrier.get("zero_perturbation_best_q_upper", math.nan))
    if not math.isfinite(zero_q) or zero_q < 0.0:
        f.append("invalid zero-perturbation q")
    if barrier.get("zero_perturbation_joint_box_nonempty") is not True:
        f.append("zero-perturbation correction box is empty")
    if barrier.get("barrier_established") is (
            barrier.get("perturbation_route_can_close_authoritative_witness")):
        f.append("inconsistent barrier verdict")

    allowed = set(_NEXT_BY_DP_TERM.values()) | {
        _BARRIER_OBLIGATION, _REPAIR_OBLIGATION}
    if d.get("next_obligation") not in allowed:
        f.append("invalid V50 next obligation")
    if d.get("P5_SAMPLE1_ZERO_PERTURBATION_BARRIER_V50") not in (
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
    d = build(x.domain, source_pieces=x.source_pieces,
              source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
              tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces,
              residual_x_pieces=x.residual_x_pieces,
              parallel_pieces=x.parallel_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    detail = d.get("attitude_supported_row_detail") or {}
    print(json.dumps({
        "status": d["P5_SAMPLE1_ZERO_PERTURBATION_BARRIER_V50"],
        "V34_row_DeltaS": detail.get(
            "V34_first_measurement_row_DeltaS_intersected_upper"),
        "refined_row_DeltaS": detail.get(
            "attitude_supported_row_DeltaS_intersected_upper"),
        "refinement_ratio": detail.get("row_DeltaS_refinement_ratio"),
        "refined_dominant_term": detail.get("refined_dominant_term"),
        "refined_dominant_fraction": detail.get("refined_dominant_fraction"),
        "unpinned": d.get("theta_gain_rows_unpinned_from_V12D_parent"),
        "dP_decomposition": d.get("reduced_covariance_perturbation_decomposition"),
        "barrier": d.get("zero_perturbation_barrier"),
        "angles": d.get("geodesic_principal_angle_diagnostics"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
