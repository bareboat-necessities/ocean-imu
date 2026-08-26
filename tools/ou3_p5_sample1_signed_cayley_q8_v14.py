#!/usr/bin/env python3
"""V14 source-correlated sample-1 signed Cayley composition inside q<8.

V12D closes the complete focused sample-1 gain/remainder family and V13E proves
that every resulting correction has a signed component enclosure and a radial
cell accepted by the winding-aware deployed-quaternion primitive through 9 rad.
The remaining obligation is not correction magnitude; it is exclusion of the
*resulting* Cayley antipode and proof that the post-correction Cayley norm is <8.

A generic pre-correction ball ||c||<=q1 throws away the crucial source
correlation.  In the V7/V10 gravity gauge the first (sample-0) nominal attitude
correction is +x.  Before that correction the exact startup source geometry
proves a gravity-tangent Cayley radius c_t, hence |c_x|<=c_t even though the yaw
component may be much larger.  For an x-axis deployed correction the post-
correction Cayley x coordinate depends only on c_x.  V12D's tiny off-axis
sample-0 correction remainder is retained by composing the complete component
box, not discarded.

After the nominal Rx proof-gauge change, c_x is unchanged.  The next 5 ms
attitude-error transport and the possible sample-1 S attitude injection are
then included as arbitrary small rotations.  This gives, for each V10/V12D
source row, a signed sample-1 c_x interval and a full q1 upper.  The remaining
current-error component obeys

    ||c_yz|| <= sqrt(q1^2 - min|c_x|^2).

For each V13E signed correction subcell, the source-normalized shipping
correction quaternion [w_d,v_d] is enclosed directly.  Its signed product term
is bounded by

    v_d^T c = v_dx c_x + v_d,yz^T c_yz,

using the signed x interval and a Euclidean yz Cauchy bound.  Thus

    W = 2 w_d - v_d^T c

is source-correlated enough to test the actual resulting Cayley antipode.
Because the shipping correction quaternion is unit after source normalization,
left quaternion multiplication preserves norm:

    W^2 + ||V||^2 = 4 + ||c||^2.

Therefore, once |W|>=W_min>0,

    ||c_plus|| <= 2 sqrt((4+q1^2)/W_min^2 - 1).

This producer requires that upper to be strictly below 8 for every signed
subcell.  It changes no filter, source domain, correction limit, or P5 word
promotion state.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, hull
import ou3_p5_deployed_quaternion_cayley_cell as CAYLEY1
import ou3_p5_deployed_quaternion_cayley_cell_v2 as CAYLEY2
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_sample1_prefix_v2 as PREFIX2
import ou3_p5_sample1_signed_radial_subcells_v13 as V13
import ou3_p5_sample1_signed_radial_subcells_v13e as V13E
import ou3_p5_sample1_structured_full_gain_v12d as V12D

DEFAULT_DOMAIN = V12D.DEFAULT_DOMAIN
SCHEMA = 14
Q_TARGET = 8.0
FULL = V12D.FULL
SERIES = 1.0e-2


def _minimum_abs(x: Interval) -> float:
    if x.lo <= 0.0 <= x.hi:
        return 0.0
    return min(abs(x.lo), abs(x.hi))


def _norm2_upper(v) -> float:
    return CAYLEY1._norm_upper(v)


def _normalized_series_part(dbox, hi: float):
    raw = CAYLEY1._series_homogeneous(dbox, hi)
    if raw is None:
        return None
    w, v = raw
    # The source immediately normalizes this polynomial quaternion.  A positive
    # norm interval makes component-wise division a rigorous enclosure.
    s_lo = FULL.down(w.lo * w.lo)
    s_hi = FULL.up(w.hi * w.hi)
    for x in v:
        amin = _minimum_abs(x)
        if amin > 0.0:
            s_lo = FULL.down(s_lo + FULL.down(amin * amin))
        amax = x.abs_upper()
        s_hi = FULL.up(s_hi + FULL.up(amax * amax))
    if not s_lo > 0.0:
        raise RuntimeError("series correction normalization lost positive norm")
    n = Interval(FULL.down(math.sqrt(s_lo)), FULL.up(math.sqrt(s_hi)))
    return w / n, [x / n for x in v]


def _axis_part(dbox, lo: float, hi: float):
    lo = max(SERIES, float(lo)); hi = float(hi)
    if hi < lo:
        return None
    half = Interval(FULL.down(0.5 * lo), FULL.up(0.5 * hi))
    sin_half, cos_half = CAYLEY2._trig_interval(half.lo, half.hi)
    if not half.lo > 0.0:
        raise RuntimeError("axis correction half-angle lost positive lower")
    k = FULL.I(0.5) * sin_half / half
    return cos_half, [k * x for x in dbox]


def _normalized_shipping_quaternion(dbox, *, radial_lower: float, radial_upper: float):
    """Source-normalized correction quaternion for one signed/radial d cell."""
    hi = float(radial_upper); lo = max(0.0, float(radial_lower))
    if not (math.isfinite(lo) and math.isfinite(hi) and 0.0 <= lo <= hi <= 9.0):
        raise RuntimeError("invalid source correction radial cell")
    if hi > _norm2_upper(dbox) + 64.0 * math.ulp(max(1.0, _norm2_upper(dbox))):
        raise RuntimeError("radial upper exceeds signed component box")
    parts = []
    if lo < SERIES:
        ser_hi = min(hi, SERIES)
        if ser_hi > 0.0:
            p = _normalized_series_part(dbox, ser_hi)
            if p is not None:
                parts.append((p[0], p[1], "SERIES_NORMALIZED"))
    if hi >= SERIES:
        p = _axis_part(dbox, max(lo, SERIES), hi)
        if p is not None:
            parts.append((p[0], p[1], "AXIS_ANGLE_UNIT"))
    if not parts:
        # Exact zero correction.
        return FULL.I(1.0), [FULL.I(0.0) for _ in range(3)], ["ZERO"]
    w = hull(*(p[0] for p in parts))
    v = [hull(*(p[1][i] for p in parts)) for i in range(3)]
    return w, v, [p[2] for p in parts]


def _widen_cx_by_unknown_rotation(cx: Interval, q_before: float, angle_upper: float) -> Interval:
    """Cayley x enclosure after arbitrary rotation of bounded geodesic angle."""
    q = float(q_before); a = float(angle_upper)
    if a <= 0.0:
        return cx
    if not (math.isfinite(q) and q >= 0.0 and math.isfinite(a) and 0.0 <= a < math.pi):
        raise RuntimeError("invalid bounded attitude transport")
    ca = FULL.up(2.0 * math.tan(FULL.up(0.5 * a)))
    dot = FULL.up(0.25 * FULL.up(ca * q))
    den_lo = FULL.down(1.0 - dot); den_hi = FULL.up(1.0 + dot)
    if not den_lo > 0.0:
        raise RuntimeError("small attitude transport can reach Cayley antipode")
    # |a_x + 0.5(a x c)_x| <= |a| (1 + |c|/2).
    bump = FULL.up(ca * FULL.up(1.0 + FULL.up(0.5 * q)))
    num = Interval(FULL.down(cx.lo - bump), FULL.up(cx.hi + bump))
    return num / Interval(den_lo, den_hi)


def _sample1_current_chart(*, first: dict, base: dict, vr: dict,
                           dom: dict, src: dict, sample1_s_angle: float):
    qpre = float(first["post_prediction_full_cayley_norm_upper"])
    ctan = float(first["post_prediction_cayley_tangent_norm_upper"])
    d0 = Interval.outward_bounds(*map(float, base["first_attitude_correction_rad"]))
    dd = float(vr["first_offaxis_attitude_correction_upper_rad"])
    e = Interval.outward_bounds(-dd, dd)
    dbox0 = [d0 + e, e, e]
    d0_hi = FULL.up(max(0.0, d0.hi) + dd)
    d0_lo = max(0.0, FULL.down(max(0.0, d0.lo) - dd))

    # Startup source: both gravity-tangent components are jointly <=ctan and
    # the full Cayley vector is <=qpre.  The component box is used only in this
    # first exact product; the independent qpre radius remains authoritative.
    cpre = [Interval.outward_bounds(-ctan, ctan),
            Interval.outward_bounds(-ctan, ctan),
            Interval.outward_bounds(-qpre, qpre)]
    w0, v0, branches0 = _normalized_shipping_quaternion(
        dbox0, radial_lower=d0_lo, radial_upper=d0_hi)
    dot0 = CAYLEY1._dot(v0, cpre)
    W0 = FULL.I(2.0) * w0 - dot0
    if W0.lo <= 0.0 <= W0.hi:
        raise RuntimeError("sample-0 source-correlated product scalar crosses zero")
    cross0 = CAYLEY1._cross(v0, cpre)
    Vx0 = w0 * cpre[0] + FULL.I(2.0) * v0[0] + cross0[0]
    cx = FULL.I(2.0) * Vx0 / W0

    q0 = PREFIX2._post_correction_q_upper(qpre, d0_hi)
    if not math.isfinite(q0):
        raise RuntimeError("sample-0 q upper left Cayley chart")
    # V7's simultaneous Rx(d0)^T proof-coordinate change is about x and leaves
    # the x coordinate of a Cayley vector unchanged.
    transport = float(first["first_prediction_transport_angle_upper_rad"])
    cx = _widen_cx_by_unknown_rotation(cx, q0, transport)
    qpred = RG._q_after_first_prediction(q0, dom, float(src["dt_s"]))
    if not math.isfinite(qpred):
        raise RuntimeError("sample-1 prediction q upper left Cayley chart")

    ds = max(0.0, float(sample1_s_angle))
    cx = _widen_cx_by_unknown_rotation(cx, qpred, ds)
    q1 = PREFIX2._post_correction_q_upper(qpred, ds)
    if not math.isfinite(q1):
        raise RuntimeError("sample-1 S update q upper left Cayley chart")

    cx_min = _minimum_abs(cx)
    yz2 = max(0.0, FULL.up(q1 * q1) - FULL.down(cx_min * cx_min))
    cyz = FULL.up(math.sqrt(yz2))
    return {
        "cx": cx,
        "q1": q1,
        "cyz_norm_upper": cyz,
        "sample0_correction_radial_lower": d0_lo,
        "sample0_correction_radial_upper": d0_hi,
        "sample0_quaternion_branches": branches0,
        "sample0_product_scalar": W0,
    }


def _qplus_from_product_scalar(q1: float, W: Interval) -> tuple[float, float]:
    if W.lo <= 0.0 <= W.hi:
        return 0.0, math.inf
    wmin = min(abs(W.lo), abs(W.hi))
    if not wmin > 0.0:
        return 0.0, math.inf
    num = FULL.up(4.0 + FULL.up(q1 * q1))
    den = FULL.down(wmin * wmin)
    if not den > 0.0:
        return wmin, math.inf
    ratio = FULL.up(num / den)
    r = max(0.0, FULL.up(ratio - 1.0))
    qplus = FULL.up(2.0 * math.sqrt(r))
    return wmin, qplus


def _prereq_failure(failures: list[str]) -> dict:
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V13E_signed_radial_prerequisite_passed": False,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "evaluated_signed_cayley_cells": 0,
        "unclosed_q8_cells": 0,
        "first_unclosed_q8_cell": None,
        "max_post_sample1_cayley_norm_upper": 0.0,
        "minimum_abs_product_scalar_lower": None,
        "signed_cayley_q8_composed_here": False,
        "complete_sample1_branch_closed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14": "NOT_ESTABLISHED",
        "next_obligation": "ESTABLISH_V13E_SIGNED_RADIAL_PREREQUISITE",
        "failures": failures,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6) -> dict:
    path = Path(domain_path).resolve()
    dom = json.loads(path.read_text(encoding="utf-8"))
    v13e = V13E.build(
        path, source_pieces=source_pieces, source_cell_index=source_cell_index,
        p_pieces=p_pieces, tangent_pieces=tangent_pieces,
        axial_pieces=axial_pieces, residual_x_pieces=residual_x_pieces,
        parallel_pieces=parallel_pieces)
    failures = [f"V13E: {x}" for x in V13E.validate(v13e)]
    if v13e.get("P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13E") != "PASS":
        failures.append("V13E prerequisite did not pass")
        return _prereq_failure(failures)

    v12d = V12D.build(
        path, source_pieces=source_pieces, source_cell_index=source_cell_index,
        p_pieces=p_pieces, tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    failures += [f"V12D: {x}" for x in V12D.validate(v12d)]
    core = V12D.V11.V10.build(
        path, source_pieces=source_pieces, source_cell_index=source_cell_index,
        p_pieces=p_pieces, tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    failures += [f"V10: {x}" for x in V12D.V11.V10.validate(core)]
    first_all = V12D.V11.FIRST.build(path, source_pieces=source_pieces)
    failures += [f"first: {x}" for x in V12D.V11.FIRST.validate(first_all)]
    src, phase = RG._source_phase_children(source_pieces)[source_cell_index]
    if phase != "due":
        failures.append("V14 focused family requires first due source cell")
    first = dict(first_all)
    # Global source geometry is top-level; row-specific covariance/residual data
    # remains in first_all['source_cells'] and V10/V12D rows.
    sb = v12d.get("sample1_S_perturbation_bounds", {})
    ds = float(sb.get("sample1_S_attitude_correction_upper_rad", 0.0))
    if not (math.isfinite(ds) and ds >= 0.0):
        failures.append("invalid sample-1 S attitude correction upper")

    fr = first_all["source_cells"][source_cell_index]
    p_all = Interval.outward_bounds(*map(float, fr["P_aw_variance_interval"]))
    pcells = V12D.V11.SUB.parts(p_all.lo, p_all.hi, p_pieces)
    hstep = float(src["dt_s"]); g = float(dom["startup"]["gravity_mps2"])
    tilt, yaw, eps = RG._attitude_covariance_epsilon(path, hstep)
    t = Interval.outward_bounds(tilt, FULL.up(tilt + eps))
    Y = Interval.outward_bounds(yaw, FULL.up(yaw + eps))
    vec = V12D.V11.VECTOR.build()
    r = FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))[0][0]
    F, Q, _ = FULL._transition_and_Q(src, dom)
    alpha = F[15][15]; qaw = Q[15][15]

    if len(core.get("rows", [])) != len(v12d.get("rows", [])):
        failures.append("V10/V12D row counts differ")

    total = unclosed = antipode = radial_not_ready = 0
    max_qplus = 0.0; min_w = math.inf
    first_bad = None; worst = None
    sample_bad = []
    chart_cache = {}

    for base, vr in zip(core.get("rows", []), v12d.get("rows", [])):
        ids = ("p_cell", "tangent_residual_cell", "axial_residual_cell")
        if any(int(base[k]) != int(vr[k]) for k in ids):
            failures.append("V10/V12D row ordering mismatch")
            break
        pi = int(base["p_cell"]); p = pcells[pi]
        rt = Interval.outward_bounds(*map(float, base["first_tangent_residual_magnitude_mps2"]))
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
        (gy, gz), (_kxy, _kxz), _detail = V13._signed_gain_components(
            a=a, Y=Y, c0=c0, alpha=alpha, qaw=qaw, b=b, bz=bz,
            det_first=det_first, d=d0, fy=fy, fz=fz, r=r)

        kperp = float(base["Ktheta_perpendicular_block_upper"])
        kpar = float(base["Ktheta_parallel_block_upper"])
        drho = float(vr["total_residual_perturbation_upper_mps2"])
        dk = float(vr["sample1_attitude_gain_operator_perturbation_upper"])
        rho = float(base["sample1_full_residual_norm_upper_mps2"])
        rho_x = min(rho, float(base["sample1_combined_source_x_residual_upper_mps2"]))
        k0 = max(kperp, kpar)
        eta = FULL.up(FULL.up(k0 * drho) + FULL.up(dk * FULL.up(rho + drho)))
        global_hi = float(vr["V12C_correction_norm_upper_rad"])

        key = (int(base["tangent_residual_cell"]), int(base["p_cell"]), int(base["axial_residual_cell"]))
        try:
            chart = _sample1_current_chart(
                first=first, base=base, vr=vr, dom=dom, src=src,
                sample1_s_angle=ds)
        except Exception as exc:
            failures.append(f"current chart {key}: {exc}")
            break
        chart_cache[key] = chart

        rx_cells = V12D.V11.SUB.parts(-rho_x, rho_x, residual_x_pieces)
        for rxc in rx_cells:
            rx_min = _minimum_abs(rxc)
            rem2 = max(0.0, FULL.up(rho * rho) - FULL.down(rx_min * rx_min))
            ryz_hi = FULL.up(math.sqrt(rem2))
            ux_hi = FULL.up(kpar * ryz_hi)
            u_cells = V12D.V11.SUB.parts(-ux_hi, ux_hi, parallel_pieces)
            for uc in u_cells:
                e = Interval.outward_bounds(-eta, eta)
                dbox = [uc + e, gy * rxc + e, gz * rxc + e]
                box_lo = CAYLEY2._norm_lower(dbox)
                rx_abs = rxc.abs_upper(); u_abs = uc.abs_upper()
                nominal_hi = FULL.up(math.sqrt(FULL.up(
                    FULL.up(u_abs * u_abs) + FULL.up((kperp * rx_abs) ** 2))))
                radial_hi = min(global_hi, FULL.up(nominal_hi + eta))
                radial_lo = min(box_lo, radial_hi)
                ready = radial_hi <= 6.0 or (radial_lo > 0.0 and radial_hi <= 9.0)
                total += 1
                if not ready:
                    radial_not_ready += 1; unclosed += 1
                    row = {"reason": "radial cell not ready", "p_cell": pi,
                           "tangent_residual_cell": int(base["tangent_residual_cell"]),
                           "axial_residual_cell": int(base["axial_residual_cell"]),
                           "correction_component_box_rad": [x.as_list() for x in dbox],
                           "radial_lower_rad": radial_lo, "radial_upper_rad": radial_hi}
                    if first_bad is None: first_bad = row
                    continue

                try:
                    wd, vd, branches = _normalized_shipping_quaternion(
                        dbox, radial_lower=radial_lo, radial_upper=radial_hi)
                except Exception as exc:
                    unclosed += 1
                    row = {"reason": f"correction quaternion: {exc}", "p_cell": pi,
                           "tangent_residual_cell": int(base["tangent_residual_cell"]),
                           "axial_residual_cell": int(base["axial_residual_cell"])}
                    if first_bad is None: first_bad = row
                    continue

                cx = chart["cx"]; q1 = float(chart["q1"]); cyz = float(chart["cyz_norm_upper"])
                xdot = vd[0] * cx
                vdyz = _norm2_upper([vd[1], vd[2]])
                yzdot = FULL.up(vdyz * cyz)
                dot = xdot + Interval.outward_bounds(-yzdot, yzdot)
                W = FULL.I(2.0) * wd - dot
                wmin, qplus = _qplus_from_product_scalar(q1, W)
                closed = math.isfinite(qplus) and qplus < Q_TARGET and wmin > 0.0
                if W.lo <= 0.0 <= W.hi:
                    antipode += 1
                if wmin > 0.0:
                    min_w = min(min_w, wmin)
                if math.isfinite(qplus):
                    max_qplus = max(max_qplus, qplus)
                row = {
                    "p_cell": pi,
                    "tangent_residual_cell": int(base["tangent_residual_cell"]),
                    "axial_residual_cell": int(base["axial_residual_cell"]),
                    "nominal_rx_subcell_mps2": rxc.as_list(),
                    "nominal_parallel_correction_subcell_rad": uc.as_list(),
                    "correction_component_box_rad": [x.as_list() for x in dbox],
                    "correction_radial_lower_rad": radial_lo,
                    "correction_radial_upper_rad": radial_hi,
                    "correction_quaternion_branches": branches,
                    "sample1_current_cx": cx.as_list(),
                    "sample1_current_cayley_norm_upper": q1,
                    "sample1_current_cyz_norm_upper": cyz,
                    "product_scalar": W.as_list(),
                    "abs_product_scalar_lower": wmin,
                    "post_sample1_cayley_norm_upper": qplus,
                    "closed_inside_q8": closed,
                }
                if worst is None or (math.isfinite(qplus) and qplus > worst.get("post_sample1_cayley_norm_upper", -1.0)):
                    worst = row
                if not closed:
                    unclosed += 1
                    if len(sample_bad) < 24: sample_bad.append(row)
                    if first_bad is None: first_bad = row

    ok = not failures and total > 0 and unclosed == 0 and first_bad is None
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V13E_signed_radial_prerequisite_passed": True,
        "V12D_tangent_channel_prerequisite_passed": v12d.get("P5_SAMPLE1_FIRST_PSD_TANGENT_REFINEMENT_V12D") == "PASS",
        "pre_sample0_gravity_tangent_cayley_bound_used": True,
        "sample0_offaxis_correction_remainder_included": True,
        "sample0_shipping_quaternion_normalization_included": True,
        "V7_Rx_proof_gauge_leaves_cx_unchanged": True,
        "one_step_attitude_transport_included": True,
        "sample1_S_attitude_injection_included": True,
        "signed_vdx_cx_plus_yz_cauchy_product_used": True,
        "unit_quaternion_product_norm_identity_used": True,
        "post_q_target_strict": Q_TARGET,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "evaluated_signed_cayley_cells": total,
        "radial_not_ready_cells": radial_not_ready,
        "product_scalar_antipode_cells": antipode,
        "unclosed_q8_cells": unclosed,
        "minimum_abs_product_scalar_lower": None if min_w is math.inf else min_w,
        "max_post_sample1_cayley_norm_upper": max_qplus,
        "first_unclosed_q8_cell": first_bad,
        "worst_q8_cell": worst,
        "sample_unclosed_q8_cells": sample_bad,
        "signed_cayley_q8_composed_here": ok,
        "complete_sample1_branch_closed_here": ok,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14": "PASS" if ok else "NOT_ESTABLISHED",
        "next_obligation": (
            "LIFT_CLOSED_SAMPLE1_PREFIX_TO_ALL_SOURCE_PHASE_CELLS_AND_CONTINUE_SAMPLE2_PREFIX"
            if ok else "REFINE_CURRENT_CX_OR_CORRECTION_RADIAL_DIRECTION_AT_FIRST_Q8_WITNESS"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit", "V13E_signed_radial_prerequisite_passed",
        "pre_sample0_gravity_tangent_cayley_bound_used",
        "sample0_offaxis_correction_remainder_included",
        "sample0_shipping_quaternion_normalization_included",
        "V7_Rx_proof_gauge_leaves_cx_unchanged",
        "one_step_attitude_transport_included", "sample1_S_attitude_injection_included",
        "signed_vdx_cx_plus_yz_cauchy_product_used",
        "unit_quaternion_product_norm_identity_used",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "q8_word_promoted_here", "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    st = d.get("P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14")
    if st == "PASS":
        if d.get("signed_cayley_q8_composed_here") is not True or d.get("complete_sample1_branch_closed_here") is not True:
            f.append("V14 PASS did not close sample1 q8 branch")
        if int(d.get("unclosed_q8_cells", -1)) != 0:
            f.append("V14 PASS retains unclosed q8 cells")
        if not float(d.get("max_post_sample1_cayley_norm_upper", math.inf)) < Q_TARGET:
            f.append("V14 PASS does not satisfy strict q<8")
        if not float(d.get("minimum_abs_product_scalar_lower", 0.0)) > 0.0:
            f.append("V14 PASS lacks strict product-scalar separation")
    elif st == "NOT_ESTABLISHED":
        if d.get("complete_sample1_branch_closed_here") is not False:
            f.append("V14 nonclosure claims sample1 closure")
        if d.get("first_unclosed_q8_cell") is None and not f:
            f.append("V14 nonclosure lacks witness")
    else:
        f.append("invalid V14 status")
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
        residual_x_pieces=x.residual_x_pieces, parallel_pieces=x.parallel_pieces)
    vf = validate(d); d["validation_failures"] = vf
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14"],
        "cells": d["evaluated_signed_cayley_cells"],
        "radial_not_ready": d["radial_not_ready_cells"],
        "antipode_cells": d["product_scalar_antipode_cells"],
        "unclosed": d["unclosed_q8_cells"],
        "min_abs_W": d["minimum_abs_product_scalar_lower"],
        "max_q_plus": d["max_post_sample1_cayley_norm_upper"],
        "first_unclosed": d["first_unclosed_q8_cell"],
        "worst": d["worst_q8_cell"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
