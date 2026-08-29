#!/usr/bin/env python3
"""V51: exact monotone enclosure of the first-block source gain rationals.

V50 proved that no refinement of any correction-perturbation bound can close the
authoritative V41/V45 first survivor `(p,t,a)=(0,0,23)`, and named the remaining
obligation: sharpen the *nominal* geometry - V10's exact first-accelerometer
directional correction magnitude, or the sample-0 current chart.

V51 discharges that obligation at the witness.

The first-accelerometer block of the shipping filter is described by seven
rational functions of the same three source intervals: the attitude variance
`t`, the `a_w` variance `p`, and the accelerometer noise variance `r`.  With
`D = g^2 t + p + r`,

    a       = t (p+r) / D,        b   = p (g^2 t + r) / D,
    c0      = -g t p / D,         bz  = p r / (p+r),
    det     = t p r / D,          k_theta = g t / D,
    k_aw,t  = p / D,              k_z     = p / (p+r).

Every one of these is evaluated by the parent backend as a straight interval
expression, which loses the dependency between the numerator and the
denominator.  The loss is not cosmetic.  For the witness cell the parent encloses

    k_z = p/(p+r)  in  [0.5594923342554586, 1.0537323143362434],

even though `p/(p+r) < 1` holds identically for positive `p, r`.

Each of the eight expressions is monotone in each of `t, p, r` separately, so its
exact range over the parameter box is attained at a corner:

    d/dp [p/(p+r)]      = r/(p+r)^2      > 0,   d/dr < 0;
    d/dp [p/D]          = (g^2 t + r)/D^2 > 0,  d/dt < 0,  d/dr < 0;
    d/dt [g t/D]        = g(p+r)/D^2     > 0,   d/dp < 0,  d/dr < 0;
    a = t u/(g^2 t + u) with u = p+r:  d/dt = u^2/D^2 > 0,  d/du = g^2 t^2/D^2 > 0;
    b = p w/(p + w)     with w = g^2 t + r:  d/dp = w^2/(p+w)^2 > 0, d/dw > 0;
    bz = p r/(p+r):     increasing in p and in r;
    |c0| = g t p/D:     d/dt = g p (p+r)/D^2 > 0, d/dp = g t (g^2 t + r)/D^2 > 0,
                        d/dr < 0;
    det = t p r/D:      increasing in t, in p and in r.

`u` and `w` each range over an exact interval and are independent of the
remaining variable, so the corner in the reduced variables is a genuine corner of
the original box.  V51 therefore evaluates each expression at its extremal corner
with the same outward-rounded backend and intersects the result with the parent
enclosure.  The refinement can never widen a bound and never leaves the parent;
if it ever did, the producer fails closed.

Consequences at the authoritative witness, all reproduced from source rather
than copied:

    k_z upper       1.0537323143362434 -> 0.8001468619320714
    sample-1 |f|    21.395742136954993 -> 18.606777069495593
    sample-1 rho    17.922551201967796 -> 15.229738748335985
    k_perp           0.9753682347137846 -> 0.8468904975139163
    k_parallel       0.09899544770387604 -> 0.09772949400653325
    V10 correction   2.0466720610769817 -> 1.7313776836494923 rad

V50 measured that the geodesic branch needed the correction principal angle to
fall by 0.019476337434169544 rad to reach `q<8`.  This refinement removes
0.3152943774274894 rad, sixteen times that.  Composing the refined correction
with V41's archived sample-0 chart on the same SO(3) triangle gives

    q = 4.8010333986449245 < 8,

against the archived parent `q = 8.344528951460543`.  The authoritative first
survivor closes.

Scope.  V51 evaluates the single authoritative witness cell and composes it with
the archived V41 chart, which is an upper bound computed with the unrefined
gains and therefore a conservative partner for the refined correction.  It does
**not** lift the exact monotone enclosure over the complete V41 source-cell-0
cover, and it does not compose q<8, promote sample 1 or P5, or set `N_H_words`.
That lift is the next obligation.  No filter setting, source domain, six-radian
correction limit, `q<8` target, or source language changes here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_structured_full_gain_v8 as V8
import ou3_p5_sample1_structured_full_gain_v9 as V9
import ou3_p5_sample1_structured_full_gain_v10 as V10
import ou3_p5_sample1_signed_cayley_q8_v15 as V15
import ou3_p5_sample1_v41_authoritative_split_signed_v45 as V45

DEFAULT_DOMAIN = V8.DEFAULT_DOMAIN
SCHEMA = 5100
Q_TARGET = 8.0
WITNESS = (0, 0, 23)

V4 = V8.V4
V5 = V8.V5
V6 = V8.V6
FULL = V8.FULL
SUB = V8.SUB
RG = V8.RG
FIRST = V8.FIRST

V41_Q_CURRENT = V45.V41_Q_CURRENT
V41_Q_POST = V45.V41_Q_POST

#: Bound at import time.  A later stage that installs a refined block into V8
#: therefore cannot silently change what V51 compares its refinement against.
_PARENT_FIRST_BLOCK = V8._first_block_quantities

#: Archived parent values for the authoritative witness cell.  V51 recomputes
#: them from source in parent mode and refuses to continue unless every one is
#: reproduced exactly.
PARENT_WITNESS = {
    "first_tangent_residual_magnitude_mps2": [0.0, 0.5919938745482645],
    "first_axial_residual_mps2": [10.219984799344509, 11.149074326557649],
    "first_attitude_correction_rad": [-5e-324, 0.05872522312813619],
    "sample1_force_norm_upper_mps2": 21.395742136954993,
    "sample1_residual_norm_upper_mps2": 17.922551201967796,
    "combined_x_residual_upper_mps2": 1.0514397697958757,
    "Ktheta_perpendicular_block_upper": 0.9753682347137846,
    "Ktheta_parallel_block_upper": 0.09899544770387604,
    "combined_directional_correction_norm_upper_rad": 2.0466720610769817,
}


def _point(x: float) -> Interval:
    return Interval.point(float(x))


def _monotone(low_corner: Interval, high_corner: Interval,
              parent: Interval, *, name: str) -> Interval:
    """Intersect a corner-evaluated exact range with its parent enclosure.

    ``low_corner`` and ``high_corner`` are the expression evaluated with the
    outward-rounded backend at the corners where it attains its minimum and its
    maximum.  Their outer endpoints therefore bound the exact range.  The
    intersection with ``parent`` keeps a valid enclosure, can only narrow it,
    and fails closed if the two disagree.
    """
    lo = max(parent.lo, low_corner.lo)
    hi = min(parent.hi, high_corner.hi)
    if not (math.isfinite(lo) and math.isfinite(hi)) or lo > hi:
        raise RuntimeError(f"monotone refinement of {name} escaped its parent")
    return Interval(lo, hi)


def _first_block(*, t: Interval, p: Interval, r: Interval, g: float,
                 exact: bool) -> dict:
    """Return the eight first-accelerometer block quantities.

    ``exact=False`` reproduces the parent backend expression by expression.
    ``exact=True`` replaces each by its corner-evaluated exact range,
    intersected with that same parent.
    """
    G = FULL.I(g)
    G2 = FULL.I(g * g)
    parent = _PARENT_FIRST_BLOCK(t=t, p=p, r=r, g=g)
    if not exact:
        return parent

    tl, th = _point(t.lo), _point(t.hi)
    pl, ph = _point(p.lo), _point(p.hi)
    rl, rh = _point(r.lo), _point(r.hi)
    ul, uh = pl + rl, ph + rh          # exact range of u = p + r
    wl, wh = G2 * tl + rl, G2 * th + rh  # exact range of w = g^2 t + r

    mag = _monotone(G * tl * pl / (G2 * tl + pl + rh),
                    G * th * ph / (G2 * th + ph + rl),
                    -parent["c0"], name="|c0|")
    return {
        "a": _monotone(tl * ul / (G2 * tl + ul), th * uh / (G2 * th + uh),
                       parent["a"], name="a"),
        "c0": Interval(-mag.hi, -mag.lo),
        "b": _monotone(pl * wl / (pl + wl), ph * wh / (ph + wh),
                       parent["b"], name="b"),
        "bz": _monotone(pl * rl / (pl + rl), ph * rh / (ph + rh),
                        parent["bz"], name="bz"),
        "det_first": _monotone(tl * pl * rl / (G2 * tl + pl + rl),
                               th * ph * rh / (G2 * th + ph + rh),
                               parent["det_first"], name="det"),
        "ktheta": _monotone(G * tl / (G2 * tl + ph + rh),
                            G * th / (G2 * th + pl + rl),
                            parent["ktheta"], name="k_theta"),
        "kaw_t": _monotone(pl / (G2 * th + pl + rh),
                           ph / (G2 * tl + ph + rl),
                           parent["kaw_t"], name="k_aw_t"),
        "kz": _monotone(pl / (pl + rh), ph / (ph + rl),
                        parent["kz"], name="k_z"),
    }


def _source_context(path: Path, *, source_pieces: int,
                    source_cell_index: int) -> dict:
    """Load the source-derived constants shared by both evaluation modes."""
    V8.FULL3._install_backend()
    dom = json.loads(path.read_text(encoding="utf-8"))
    src, phase = RG._source_phase_children(source_pieces)[source_cell_index]
    if phase != "due":
        raise RuntimeError("V51 witness refinement requires the first due source cell")
    first = FIRST.build(path, source_pieces=source_pieces)
    vec = V8.VECTOR.build()
    fr = first["source_cells"][source_cell_index]

    g = float(dom["startup"]["gravity_mps2"])
    ba = float(dom["startup"]["physical_handoff_coordinate_bounds"][
        "accelerometer_bias_error_norm_upper_mps2"])
    hstep = float(src["dt_s"])
    tilt, yaw, eps = RG._attitude_covariance_epsilon(path, hstep)
    F, Q, _ = FULL._transition_and_Q(src, dom)
    cos0 = float(first["post_prediction_true_gravity_cosine_lower"])
    sin_hi = 1.0 if cos0 < 0.0 else FULL.up(
        math.sqrt(max(0.0, FULL.up(1.0 - FULL.down(cos0 * cos0)))))
    rho0 = float(fr["combined_useful_residual_norm_upper_mps2"])
    aw_pred = float(fr["predicted_aw_error_norm_upper_mps2"])
    yRt = FULL.up(g * sin_hi)
    yRz = FULL.up(g * max(0.0, FULL.up(1.0 - cos0)))
    theta_transport = float(first["first_prediction_transport_angle_upper_rad"])
    return {
        "g": g, "ba": ba, "rho0": rho0, "aw_pred": aw_pred,
        "t": Interval.outward_bounds(tilt, FULL.up(tilt + eps)),
        "Y": Interval.outward_bounds(yaw, FULL.up(yaw + eps)),
        "p_all": Interval.outward_bounds(*map(float, fr["P_aw_variance_interval"])),
        "r": FULL._R_diag(float(vec["configured_measurement_bounds"][
            "acc_measurement_std_mps2"]))[0][0],
        "alpha": F[15][15], "qaw": Q[15][15],
        "yRt": yRt, "yRz": yRz,
        "chord0": V4._gravity_chord_from_cos(cos0),
        "pred_chord": V4._correction_chord_upper(theta_transport),
        "transport_chord": V4._correction_chord_upper(theta_transport),
        "rt_cap": min(rho0, FULL.up(aw_pred + FULL.up(yRt + ba))),
        "rz_cap": min(rho0, FULL.up(aw_pred + FULL.up(yRz + ba))),
    }


def _witness_cell(ctx: dict, *, p_pieces: int, tangent_pieces: int,
                  axial_pieces: int) -> tuple[Interval, Interval, Interval]:
    """Return the (p, tangent residual, axial residual) cells of the witness."""
    pi, ti, zi = WITNESS
    p = SUB.parts(ctx["p_all"].lo, ctx["p_all"].hi, p_pieces)[pi]
    rt0 = SUB.parts(0.0, ctx["rt_cap"], tangent_pieces)[ti]
    rz0 = SUB.parts(-ctx["rz_cap"], ctx["rz_cap"], axial_pieces)[zi]
    child = V6._residual_child(rt0, rz0, ctx["rho0"])
    if child is None:
        raise RuntimeError("authoritative witness residual cell was pruned")
    rt, rz = child
    return p, rt, rz


def _evaluate(ctx: dict, p: Interval, rt: Interval, rz: Interval, *,
              exact: bool) -> dict:
    """Run the shipping V8/V10 witness chain in parent or exact-monotone mode."""
    g = ctx["g"]
    t = ctx["t"]
    Y = ctx["Y"]
    r = ctx["r"]
    alpha = ctx["alpha"]
    qaw = ctx["qaw"]
    ba = ctx["ba"]
    blk = _first_block(t=t, p=p, r=r, g=g, exact=exact)

    d = blk["ktheta"] * rt
    awt = blk["kaw_t"] * rt
    az = blk["kz"] * rz
    fy = -(alpha * awt)
    fz = FULL.I(g) + alpha * az
    fn = V5._norm2_upper(fy.abs_upper(), fz.abs_upper())
    kn, detail = V8._block_gain_bounds(
        blk["a"], Y, blk["c0"], alpha, qaw, blk["b"], blk["bz"],
        blk["det_first"], d, fy, fz, r)

    dhi = max(0.0, d.hi)
    chord = min(2.0, FULL.up(
        ctx["chord0"] + FULL.up(V4._correction_chord_upper(dhi) + ctx["pred_chord"])))
    one_t = FULL.I(1.0) - blk["kaw_t"]
    one_z = FULL.I(1.0) - blk["kz"]
    left_t = (one_t * rt).abs_upper()
    left_z = (one_z * rz).abs_upper()
    post_r = FULL.up(V5._norm2_upper(
        FULL.up(left_t + ctx["yRt"]), FULL.up(left_z + ctx["yRz"])) + ba)
    post_tri = FULL.up(ctx["aw_pred"] + V5._norm2_upper(
        awt.abs_upper(), az.abs_upper()))
    post = min(post_r, post_tri)
    eaw1 = FULL.up(alpha.hi * post)
    rho = FULL.up(FULL.up(fn * chord) + FULL.up(eaw1 + ba))

    comb = V10._combined_x_residual_upper(
        alpha_lo=float(alpha.lo), alpha_hi=float(alpha.hi),
        first_rot_x_upper=ctx["yRt"], bias_upper=ba,
        error_transport_rotation_norm_upper=ctx["transport_chord"],
        series_rotation_mismatch_upper=V10._series_vs_axis_rotation_mismatch_upper(
            float(d.lo), float(d.hi)),
        pre_first_aw_error_norm_upper=ctx["aw_pred"], gravity=g)
    x_residual = float(comb["combined_x_residual_upper_mps2"])
    rho_x = min(rho, x_residual)
    kperp = float(detail["scalar_gain_yz_norm_upper"])
    kpar = float(detail["two_by_two_theta_x_gain_norm_upper"])
    corr = V9._two_block_correction_upper(kperp, kpar, rho_x, rho)
    return {
        "first_block": {k: list(v.as_list()) for k, v in blk.items()},
        "first_tangent_residual_magnitude_mps2": list(rt.as_list()),
        "first_axial_residual_mps2": list(rz.as_list()),
        "first_attitude_correction_rad": list(d.as_list()),
        "attitude_rotation_chord_upper": chord,
        "post_first_aw_left_axial_upper_mps2": left_z,
        "post_prediction_aw_error_norm_upper_mps2": eaw1,
        "sample1_force_norm_upper_mps2": fn,
        "sample1_residual_norm_upper_mps2": rho,
        "combined_x_residual_upper_mps2": x_residual,
        "sample1_combined_source_x_residual_upper_mps2": rho_x,
        "Ktheta_operator_norm_upper": kn,
        "Ktheta_perpendicular_block_upper": kperp,
        "Ktheta_parallel_block_upper": kpar,
        "combined_directional_correction_norm_upper_rad": corr,
    }


def _geodesic_q_upper(correction_upper: float) -> float:
    """Compose V41's archived sample-0 chart with a correction angle bound.

    ``_principal_axis_angle_upper`` returns the interval's upper endpoint for
    any correction radius at or below pi, so evaluating it on the degenerate
    interval at ``correction_upper`` yields the same principal angle as any
    correction interval with that upper bound.
    """
    c = float(correction_upper)
    if not (math.isfinite(c) and V15.SERIES <= c <= math.pi):
        return math.inf
    geo = V15._geodesic_q_and_scalar_lower(V41_Q_CURRENT, c, c)
    return math.inf if geo is None else float(geo[0])


def _provenance(parent: dict) -> list[str]:
    """Require the parent-mode reconstruction to reproduce every archived value."""
    out = []
    for key, want in PARENT_WITNESS.items():
        got = parent.get(key)
        if isinstance(want, list):
            if [float(x) for x in (got or [])] != [float(x) for x in want]:
                out.append(f"parent {key} is {got}, archived {want}")
        elif float(got if got is not None else math.nan) != float(want):
            out.append(f"parent {key} is {got}, archived {want}")
    return out


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24) -> dict:
    """Refine the authoritative witness with exact monotone gain enclosures."""
    path = Path(domain_path).resolve()
    failures: list[str] = []
    parent = refined = comparison = None
    parent_q = refined_q = math.inf
    try:
        ctx = _source_context(path, source_pieces=source_pieces,
                              source_cell_index=source_cell_index)
        p, rt, rz = _witness_cell(ctx, p_pieces=p_pieces,
                                  tangent_pieces=tangent_pieces,
                                  axial_pieces=axial_pieces)
        parent = _evaluate(ctx, p, rt, rz, exact=False)
        refined = _evaluate(ctx, p, rt, rz, exact=True)
        failures += _provenance(parent)

        for key, block in refined["first_block"].items():
            pv = parent["first_block"][key]
            if not (pv[0] <= block[0] and block[1] <= pv[1]):
                failures.append(f"refined {key} escaped its parent enclosure")

        narrowing = {}
        for key in ("sample1_force_norm_upper_mps2",
                    "sample1_residual_norm_upper_mps2",
                    "post_prediction_aw_error_norm_upper_mps2",
                    "post_first_aw_left_axial_upper_mps2",
                    "Ktheta_perpendicular_block_upper",
                    "Ktheta_parallel_block_upper",
                    "combined_directional_correction_norm_upper_rad"):
            pv = float(parent[key])
            rv = float(refined[key])
            if rv > pv:
                failures.append(f"refined {key} exceeded its parent")
            narrowing[key] = {"parent": pv, "refined": rv,
                              "ratio": math.inf if rv == 0.0 else pv / rv}

        corr_parent = float(parent["combined_directional_correction_norm_upper_rad"])
        corr_refined = float(refined["combined_directional_correction_norm_upper_rad"])
        parent_q = _geodesic_q_upper(corr_parent)
        refined_q = _geodesic_q_upper(corr_refined)
        comparison = {
            "narrowing": narrowing,
            "correction_angle_reduction_rad": FULL.up(corr_parent - corr_refined),
            "archived_V41_sample0_chart_q_current": V41_Q_CURRENT,
            "archived_V41_post_sample1_q_reference": V41_Q_POST,
            "parent_geodesic_q_upper": parent_q,
            "refined_geodesic_q_upper": refined_q,
            "q_target": Q_TARGET,
            "authoritative_witness_closed_by_refined_correction": (
                math.isfinite(refined_q) and refined_q < Q_TARGET),
        }
    except Exception as exc:
        failures.append(f"V51 exact monotone source gain: {exc}")

    closed = bool(comparison and comparison[
        "authoritative_witness_closed_by_refined_correction"])
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_EXACT_MONOTONE_SOURCE_GAIN_V51",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "exact_monotone_corner_enclosure_used": True,
        "parent_enclosure_retained_as_intersection": True,
        "archived_parent_witness_reproduced": not _provenance(parent or {}),
        "V41_first_survivor_row": list(WITNESS),
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "full_source_cell0_cover_lifted_here": False,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_established_here": False,
        "parent_witness_chain": parent,
        "exact_monotone_witness_chain": refined,
        "witness_comparison": comparison,
        "authoritative_witness_closed": closed,
        "P5_SAMPLE1_EXACT_MONOTONE_SOURCE_GAIN_V51": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_EXACT_MONOTONE_SOURCE_GAIN_OVER_FULL_V41_SOURCE_CELL0_COVER"
            if closed and not failures else
            "REFINE_NOMINAL_SAMPLE1_RESIDUAL_OR_SAMPLE0_CURRENT_CHART_FURTHER"),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    """Validate the fail-closed V51 proof-artifact contract."""
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_EXACT_MONOTONE_SOURCE_GAIN_V51":
        f.append("qualification mismatch")
    for k in ("source_generated_not_trajectory_fit",
              "exact_monotone_corner_enclosure_used",
              "parent_enclosure_retained_as_intersection",
              "archived_parent_witness_reproduced"):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed",
              "deployed_correction_limit_increased",
              "full_source_cell0_cover_lifted_here", "q8_composed_here",
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

    parent = d.get("parent_witness_chain") or {}
    refined = d.get("exact_monotone_witness_chain") or {}
    if not parent or not refined:
        f.append("witness chains missing")
    else:
        for key, block in (refined.get("first_block") or {}).items():
            pv = (parent.get("first_block") or {}).get(key)
            if not pv or not (pv[0] <= block[0] and block[1] <= pv[1]):
                f.append(f"refined {key} is not inside its parent")
        pc = float(parent.get("combined_directional_correction_norm_upper_rad", -1.0))
        rc = float(refined.get("combined_directional_correction_norm_upper_rad", math.inf))
        if not (math.isfinite(pc) and pc >= 0.0 and math.isfinite(rc)
                and 0.0 <= rc <= pc):
            f.append("invalid refined correction magnitude")

    cmp_ = d.get("witness_comparison") or {}
    rq = float(cmp_.get("refined_geodesic_q_upper", math.nan))
    pq = float(cmp_.get("parent_geodesic_q_upper", math.nan))
    if not (math.isfinite(rq) and math.isfinite(pq) and rq <= pq):
        f.append("invalid geodesic q comparison")
    if bool(cmp_.get("authoritative_witness_closed_by_refined_correction")) != (
            math.isfinite(rq) and rq < Q_TARGET):
        f.append("inconsistent witness closure verdict")
    if bool(d.get("authoritative_witness_closed")) != bool(
            cmp_.get("authoritative_witness_closed_by_refined_correction")):
        f.append("inconsistent reported closure")
    if d.get("P5_SAMPLE1_EXACT_MONOTONE_SOURCE_GAIN_V51") not in (
            "PASS", "NOT_ESTABLISHED"):
        f.append("invalid V51 status")
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
    d = build(x.domain, source_pieces=x.source_pieces,
              source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
              tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_EXACT_MONOTONE_SOURCE_GAIN_V51"],
        "closed": d.get("authoritative_witness_closed"),
        "comparison": d.get("witness_comparison"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
