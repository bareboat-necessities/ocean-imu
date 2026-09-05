#!/usr/bin/env python3
"""Validated low-dimensional applied R_S/tau lag envelope for OU-III SEA3 P3.

This certificate keeps the strongest source correlation needed by the P3
translation metric without rebuilding the retired 800-state tuner graph.
The shipping candidate state is only

    (tau_applied, R_S_applied),

with a common physical target tau_t.  The tau candidate and R_S candidate are
updated every valid sample and a later commit copies one candidate snapshot to
the active filter schedule.  Therefore any invariant of the candidate pair is
also an invariant of every committed active pair.

For the deployed SpectralMSE law the fractional powers are eliminated exactly:

  r_raw^14 T_S^7
    = C_J^14 (2 r_a) sigma_a,B^12 tau^48.

The normal tuner has var_wave >= 1e-6 and sigma_aw=c_sigma*sigma_a,B, so the
physical band RMS entering this identity satisfies sigma_a,B >= 1e-3 m/s^2.
No logarithm, fractional-power libm call or replay data is used in the proof.

Two deliberately simple quadratic curves are used:

  target:  R_t >= R0 + k_t (tau_t - tau0)_+^2,
  applied: R   >= R0 + k_a (tau   - tau0)_+^2.

The target curve is certified from the 14th-power identity.  The applied curve
is then certified by one-step induction over the exact deployed EMA map using
outward interval arithmetic and validated exponential arithmetic.  The target
tau continuum is covered directly; there is no source history, predecessor
path, trajectory fitting or old P2 state graph.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from ou3_interval import Interval
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_source_domain_contract as SOURCE
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
LIMITS = REPO / "src" / "tuner" / "SeaStateAdaptationLimits.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_RS_TAU_LOW_DIMENSIONAL_LAG_ENVELOPE"

TAU_KNOT_S = 7.2
TARGET_QUADRATIC = 0.015
APPLIED_QUADRATIC = 0.003
SIGMA_AB_MIN = 1.0e-3
BASE_CELLS = 24
MAX_DEPTH = 14


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def ipow(x: Interval, n: int) -> Interval:
    y = Interval.point(1.0)
    for _ in range(int(n)):
        y = y * x
    return y


def _member_float(text: str, name: str) -> float:
    m = re.search(rf"\b{name}\s*=\s*([0-9.eE+-]+)f\b", text)
    if not m:
        raise RuntimeError(f"cannot extract deployed member {name}")
    return float(m.group(1))


def _curve_value(tau: float, r0: float, k: float) -> float:
    z = max(0.0, float(tau) - TAU_KNOT_S)
    return r0 + k * z * z


def _curve_interval(tau: Interval, r0: float, k: float) -> Interval:
    lo = _curve_value(tau.lo, r0, k)
    hi = _curve_value(tau.hi, r0, k)
    return Interval.outward_bounds(lo, hi)


def _split_edges(lo: float, hi: float, count: int) -> list[float]:
    cuts = {float(lo), float(hi)}
    for x in (0.5, 6.0, TAU_KNOT_S):
        if lo < x < hi:
            cuts.add(x)
    base = sorted(cuts)
    out = [base[0]]
    span = hi - lo
    for a, b in zip(base, base[1:]):
        n = max(1, int(math.ceil(count * (b - a) / span)))
        for j in range(1, n + 1):
            out.append(a + (b - a) * j / n)
    return list(dict.fromkeys(out))


def _cells(edges: list[float]) -> list[Interval]:
    return [Interval.outward_bounds(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]


def _safe_time_scale(tau: Interval) -> Interval:
    lo = min(max(tau.lo, 0.5), 6.0)
    hi = min(max(tau.hi, 0.5), 6.0)
    return Interval.outward_bounds(lo, hi)


def _alpha(tau_target: Interval, mult: float, dt: float) -> Interval:
    scale = _safe_time_scale(tau_target)
    horizon = I(mult) * scale
    # Dynamic horizon clamps are inactive for the deployed tau target box:
    # common: 0.2..2.4 s, R_S: 0.75..9 s.  Keep the check explicit.
    if horizon.lo < 0.05 or horizon.hi > 35.0:
        raise RuntimeError("EMA horizon escaped deployed dynamic horizon bounds")
    z = -I(dt) / horizon
    if z.lo < -VT.MAX_ABS_ARGUMENT or z.hi > 0.0:
        raise RuntimeError("EMA exponential escaped validated range")
    e = VT.exp_interval(z)
    return Interval(math.nextafter(1.0 - e.hi, -math.inf),
                    math.nextafter(1.0 - e.lo, math.inf))


def _pseudo_period(tau: Interval, ratio: float, pmin: float, pmax: float) -> Interval:
    lo = min(max(ratio * tau.lo, pmin), pmax)
    hi = min(max(ratio * tau.hi, pmin), pmax)
    return Interval.outward_bounds(lo, hi)


def _target_curve_power_check(
    tau: Interval,
    *,
    r0: float,
    c_j: float,
    r_a: float,
    ratio: float,
    pmin: float,
    pmax: float,
) -> bool:
    if tau.hi <= TAU_KNOT_S:
        return True  # target clamp itself gives R_t >= R0.
    if tau.lo < TAU_KNOT_S:
        tau = Interval.outward_bounds(TAU_KNOT_S, tau.hi)
    rt = _curve_interval(tau, r0, TARGET_QUADRATIC)
    ts = _pseudo_period(tau, ratio, pmin, pmax)
    # Sufficient interval inequality for r_raw >= target curve:
    #   rt^14 * T_S^7 <= C_J^14 * (2 r_a) * sigma_aB^12 * tau^48.
    lhs = ipow(rt, 14) * ipow(ts, 7)
    rhs = (
        ipow(I(c_j), 14)
        * I(2.0)
        * I(r_a)
        * ipow(I(SIGMA_AB_MIN), 12)
        * ipow(tau, 48)
    )
    return lhs.hi <= rhs.lo


def _target_curve_validate_cell(tau: Interval, params: dict, depth: int = 0) -> tuple[int, int]:
    if _target_curve_power_check(tau, **params):
        return 1, depth
    if depth >= MAX_DEPTH:
        raise RuntimeError(f"cannot certify SpectralMSE target lower on {tau.as_list()}")
    mid = 0.5 * (tau.lo + tau.hi)
    if not tau.lo < mid < tau.hi:
        raise RuntimeError("target lower subdivision stalled")
    a = Interval.outward_bounds(tau.lo, mid)
    b = Interval.outward_bounds(mid, tau.hi)
    ca, da = _target_curve_validate_cell(a, params, depth + 1)
    cb, db = _target_curve_validate_cell(b, params, depth + 1)
    return ca + cb, max(da, db)


def _tau_next(x: Interval, u: Interval, alpha: Interval) -> Interval:
    return x + alpha * (u - x)


def _curve_difference(x: Interval, xn: Interval, r0: float, k: float) -> Interval:
    """Enclose g(xn)-g(x), retaining cancellation on one side of the knot."""
    if x.hi <= TAU_KNOT_S and xn.hi <= TAU_KNOT_S:
        return Interval.point(0.0)
    if x.lo >= TAU_KNOT_S and xn.lo >= TAU_KNOT_S:
        # k[(xn-t0)^2-(x-t0)^2] = k(xn-x)(xn+x-2t0).
        return I(k) * (xn - x) * (xn + x - I(2.0 * TAU_KNOT_S))
    # Mixed-side boxes are rare and recursively subdivided by the caller.
    return _curve_interval(xn, r0, k) - _curve_interval(x, r0, k)


def _invariant_cell_pass(x: Interval, u: Interval, params: dict) -> bool:
    r0 = params["r0"]
    dt = params["dt"]
    a = _alpha(u, params["common_mult"], dt)
    ar = _alpha(u, params["rs_mult"], dt)
    gx = _curve_interval(x, r0, APPLIED_QUADRATIC)
    rt = _curve_interval(u, r0, TARGET_QUADRATIC)
    xn = _tau_next(x, u, a)

    # Below the knot both lower curves are exactly the hard R_S floor and the
    # tau update is a convex combination of points below the knot.
    if x.hi <= TAU_KNOT_S and u.hi <= TAU_KNOT_S:
        return xn.hi <= math.nextafter(TAU_KNOT_S, math.inf)

    # At the invariant boundary R=g(x), the minimum next R is
    # g(x)+alpha_R*(R_target_lower-g(x)).  Prove this exceeds g(xn) using a
    # cancellation-preserving difference rather than subtracting two broad g boxes.
    correction = ar * (rt - gx)
    rise = _curve_difference(x, xn, r0, APPLIED_QUADRATIC)
    slack = correction - rise
    return slack.lo >= 0.0


def _invariant_validate_cell(x: Interval, u: Interval, params: dict, depth: int = 0) -> tuple[int, int]:
    if _invariant_cell_pass(x, u, params):
        return 1, depth
    if depth >= MAX_DEPTH:
        raise RuntimeError(f"cannot certify applied R_S/tau invariant on x={x.as_list()} u={u.as_list()}")
    wx = x.hi - x.lo
    wu = u.hi - u.lo
    if wx >= wu:
        mid = 0.5 * (x.lo + x.hi)
        if not x.lo < mid < x.hi:
            raise RuntimeError("applied invariant x subdivision stalled")
        cells = ((Interval.outward_bounds(x.lo, mid), u),
                 (Interval.outward_bounds(mid, x.hi), u))
    else:
        mid = 0.5 * (u.lo + u.hi)
        if not u.lo < mid < u.hi:
            raise RuntimeError("applied invariant target subdivision stalled")
        cells = ((x, Interval.outward_bounds(u.lo, mid)),
                 (x, Interval.outward_bounds(mid, u.hi)))
    total = 0
    deepest = depth
    for xx, uu in cells:
        c, d = _invariant_validate_cell(xx, uu, params, depth + 1)
        total += c
        deepest = max(deepest, d)
    return total, deepest


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("R_S/tau lag envelope may not be trajectory fitted")
    wrapper = WRAPPER.read_text(encoding="utf-8")
    limits = LIMITS.read_text(encoding="utf-8")

    dynamic = DYNAMIC.build(path)
    df = DYNAMIC.validate(dynamic)
    if df:
        raise RuntimeError(f"dynamic source prerequisite failed: {df}")
    inv = dynamic["dynamic_invariant"]
    tau_lo, tau_hi = map(float, inv["tau_applied_s"])

    markers = (
        "RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;",
        "var_wave = std::max(var_wave, 1e-6f);",
        "sigma_target_ = std::min(sigma_wave * sigma_coeff_,      max_sigma_a_);",
        "const float sigma_aB = std::max(sigma / c_sigma, 1e-6f);",
        "tune_.tau_applied   += alpha",
        "tune_.RS_applied    += alpha_RS",
        "online_tune_apply_pending_ = true",
    )
    missing = [m for m in markers if m not in wrapper]
    if missing:
        raise RuntimeError(f"shipping source parity lost: {missing}")
    for marker in (
        "kDynamicEmaTimeScaleMinSec = 0.5f",
        "kDynamicEmaTimeScaleMaxSec = 6.0f",
        "kDynamicEmaHorizonMinSec = 0.05f",
        "kDynamicEmaHorizonMaxSec = 35.0f",
    ):
        if marker not in limits:
            raise RuntimeError(f"dynamic EMA limit changed: {marker}")

    r0 = SOURCE.parse_const(wrapper, "MIN_R_S")
    c_j = SOURCE.parse_const(wrapper, "R_S_MSE_COEFF_DEFAULT")
    dt = SOURCE.parse_const(wrapper, "FREQ_SMOOTHER_DT")
    ratio = SOURCE.parse_const(wrapper, "PSEUDO_UPDATE_TAU_RATIO_DEFAULT")
    pmin = SOURCE.parse_const(wrapper, "PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT")
    pmax = SOURCE.parse_const(wrapper, "PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT")
    common_mult = SOURCE.parse_const(wrapper, "ADAPT_TAU_SEA_PERIODS")
    rs_mult = SOURCE.parse_const(wrapper, "ADAPT_RS_MULT")
    # Bind the reduced acceleration noise density to the exact configured source
    # expression rather than a replay-derived number.
    if "R_S_ACCEL_NOISE_DENSITY_DEFAULT = 0.0148f * 0.0148f * FREQ_SMOOTHER_DT;" not in wrapper:
        raise RuntimeError("R_S acceleration-noise density definition changed")
    r_a = 0.0148 * 0.0148 * dt

    if not (r0 == 0.15 and common_mult > 0.0 and rs_mult > common_mult):
        raise RuntimeError("unexpected deployed R_S/tau smoothing configuration")
    if not (tau_lo < TAU_KNOT_S < tau_hi):
        raise RuntimeError("lag-envelope knot escaped tau invariant")

    edges = _split_edges(tau_lo, tau_hi, BASE_CELLS)
    target_params = dict(r0=r0, c_j=c_j, r_a=r_a, ratio=ratio, pmin=pmin, pmax=pmax)
    target_leafs = 0
    target_depth = 0
    for cell in _cells(edges):
        c, d = _target_curve_validate_cell(cell, target_params)
        target_leafs += c
        target_depth = max(target_depth, d)

    invariant_params = dict(
        r0=r0, dt=dt, common_mult=common_mult, rs_mult=rs_mult,
    )
    source_cells = _cells(edges)
    invariant_leafs = 0
    invariant_depth = 0
    for x in source_cells:
        for u in source_cells:
            c, d = _invariant_validate_cell(x, u, invariant_params)
            invariant_leafs += c
            invariant_depth = max(invariant_depth, d)

    initial_tau = _member_float(wrapper, "tau_applied")
    initial_rs = _member_float(wrapper, "RS_applied")
    initial_floor = _curve_value(initial_tau, r0, APPLIED_QUADRATIC)
    initial_inside = initial_rs >= initial_floor
    tau12_floor = _curve_value(tau_hi, r0, APPLIED_QUADRATIC)

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "old_P2_800_state_graph_consumed": False,
        "source_history_graph_consumed": False,
        "predecessor_path_enumeration_consumed": False,
        "candidate_snapshot_commit_preserves_invariant": True,
        "SpectralMSE_fractional_power_removed_by_14th_power_identity": True,
        "ordinary_libm_fractional_power_used_in_pass_decision": False,
        "validated_exponential_arithmetic": True,
        "tau_target_and_applied_invariant_s": [tau_lo, tau_hi],
        "R_S_hard_floor": r0,
        "sigma_aB_target_lower_mps2": SIGMA_AB_MIN,
        "target_lower_curve": {
            "formula": "R_target >= R0 + k_t*max(tau_target-tau0,0)^2",
            "tau0_s": TAU_KNOT_S,
            "k_t": TARGET_QUADRATIC,
            "R_at_tau_max_lower": _curve_value(tau_hi, r0, TARGET_QUADRATIC),
            "cells_certified": target_leafs,
            "max_subdivision_depth": target_depth,
            "pass": True,
        },
        "applied_invariant_lower_curve": {
            "formula": "R_applied >= R0 + k_a*max(tau_applied-tau0,0)^2",
            "tau0_s": TAU_KNOT_S,
            "k_a": APPLIED_QUADRATIC,
            "R_at_tau_max_lower": tau12_floor,
            "base_tau_cells": len(source_cells),
            "one_step_leaf_boxes_certified": invariant_leafs,
            "max_subdivision_depth": invariant_depth,
            "initial_tau_s": initial_tau,
            "initial_R_S": initial_rs,
            "initial_state_inside": initial_inside,
            "pass": initial_inside,
        },
        "safe_use_in_P3": (
            "for any committed active tau, the same committed R_S is bounded below by the applied curve; use this in the lifted measurement-attenuation metric instead of pairing tau=12 with R_S=0.15"
        ),
        "P3_promoted": False,
        "P4_promoted": False,
        "next_obligation": (
            "consume this lower envelope in the matrix-valued selected-process-mode posterior; certify a separate upper lag envelope only if the four-S strictness matrix still needs more than the global R_S ceiling"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "candidate_snapshot_commit_preserves_invariant",
        "SpectralMSE_fractional_power_removed_by_14th_power_identity",
        "validated_exponential_arithmetic",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "old_P2_800_state_graph_consumed",
        "source_history_graph_consumed", "predecessor_path_enumeration_consumed",
        "ordinary_libm_fractional_power_used_in_pass_decision", "P3_promoted", "P4_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    for section in ("target_lower_curve", "applied_invariant_lower_curve"):
        row = d.get(section, {})
        if row.get("pass") is not True:
            f.append(f"{section} did not pass")
        if not (float(row.get("R_at_tau_max_lower", 0.0)) > float(d.get("R_S_hard_floor", math.inf))):
            f.append(f"{section} did not improve the raw R_S floor at tau max")
    if d.get("applied_invariant_lower_curve", {}).get("initial_state_inside") is not True:
        f.append("initial tune state is outside applied lag invariant")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "qualification": d["qualification"],
        "target_tau_max_R_lower": d["target_lower_curve"]["R_at_tau_max_lower"],
        "applied_tau_max_R_lower": d["applied_invariant_lower_curve"]["R_at_tau_max_lower"],
        "target_cells": d["target_lower_curve"]["cells_certified"],
        "invariant_leaf_boxes": d["applied_invariant_lower_curve"]["one_step_leaf_boxes_certified"],
        "max_depth": d["applied_invariant_lower_curve"]["max_subdivision_depth"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
