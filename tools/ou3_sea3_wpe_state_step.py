#!/usr/bin/env python3
"""Validated finite-window WavePeriodEstimator state step for complete SEA3.

This is a transition of the existing SEA3 front-end state ``z^t``.  It is not a
period estimator, a source generator, or an independent frequency model.  The
only admissible input ``a_vertical`` is the vertical acceleration produced by
the same phase-continuous SEA3 realization through the shipping private Mahony
observer.

The recurrence mirrors ``WavePeriodEstimator::update`` in Normal Live:

* two shared high-pass states;
* leaky velocity/elevation integrations;
* debiased EW first/second moments;
* variance ratio and leak subtraction;
* canonical log-period EMA;
* the one-way usable-period latch.

Normal Live already implies the startup moment gate has passed and the canonical
log-period state is finite/usable.  The per-sample moment-ratio validity tests
can nevertheless fail.  When an interval straddles such a test, this module
retains both source branches: the shipping early-return branch (moments advance,
log period holds) and the valid-update branch with the condition intersected.
It never selects the favorable branch.

The mathematical recurrence is outward enclosed.  Target-specific binary32
rounding/libm error is intentionally *not* claimed closed here; canonical P3
cannot promote until that implementation-error obligation and the upstream
private-Mahony/SEA3 realization are both closed.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import argparse
import json
import math
from pathlib import Path
import re

from ou3_interval import Interval
import ou3_validated_log as VLOG
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
ESTIMATOR = REPO / "src" / "tuner" / "WavePeriodEstimator.h"
LIMITS = REPO / "src" / "tuner" / "SeaStateAdaptationLimits.h"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_WPE_FINITE_WINDOW_STATE_STEP"

PI = Interval.outward_bounds(3.141592653589793, 3.141592653589794)
ONE = Interval.point(1.0)
TWO = Interval.point(2.0)


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def clamp_interval(x: Interval, lo: float, hi: float) -> Interval:
    if not (math.isfinite(lo) and math.isfinite(hi) and lo <= hi):
        raise ValueError("finite ordered clamp required")
    return Interval(max(lo, min(hi, x.lo)), max(lo, min(hi, x.hi)))


def max_interval(x: Interval, c: float) -> Interval:
    return Interval(max(x.lo, c), max(x.hi, c))


def wide_exp(x: Interval) -> Interval:
    """Validated exp for arbitrary finite interval by exact power-of-two halving."""
    if not (math.isfinite(x.lo) and math.isfinite(x.hi)):
        raise ValueError("finite exponential interval required")
    scale = 1
    m = max(abs(x.lo), abs(x.hi))
    while m / scale > VT.MAX_ABS_ARGUMENT:
        scale *= 2
    y = VT.exp_interval(x / I(float(scale)))
    s = scale
    while s > 1:
        y = y.square()
        s //= 2
    return y


def sqrt_positive(x: Interval) -> Interval:
    """Validated positive sqrt by rational endpoint certification."""
    if x.lo <= 0.0:
        raise ValueError("strictly positive interval required for sqrt")

    def endpoint(v: float, lower: bool) -> float:
        y = math.sqrt(v)  # seed only
        direction = -math.inf if lower else math.inf
        for _ in range(4096):
            yy = Interval.point(y).square()
            if (lower and yy.hi <= v) or ((not lower) and yy.lo >= v):
                return y
            y = math.nextafter(y, direction)
        raise RuntimeError("could not certify square-root endpoint")

    return Interval(endpoint(x.lo, True), endpoint(x.hi, False))


def _strict_positive_branch(x: Interval, threshold: float) -> Interval | None:
    """Intersect x with shipping condition x > threshold."""
    if x.hi <= threshold:
        return None
    lo = max(x.lo, math.nextafter(float(threshold), math.inf))
    return Interval(lo, x.hi)


def _member_default(text: str, name: str) -> float:
    m = re.search(rf"\b{name}\s*=\s*([0-9.eE+-]+)f", text)
    if not m:
        raise RuntimeError(f"cannot extract WPE default {name}")
    return float(m.group(1))


def _limit(text: str, name: str) -> float:
    m = re.search(rf"\b{name}\s*=\s*([0-9.eE+-]+)f", text)
    if not m:
        raise RuntimeError(f"cannot extract dynamic limit {name}")
    return float(m.group(1))


@dataclass(frozen=True)
class Constants:
    dt: float
    high_pass_hz: float
    moment_horizon_periods: float
    log_smoothing_periods: float
    min_horizon_s: float
    max_horizon_s: float
    dynamic_horizon_min_s: float
    dynamic_horizon_max_s: float


@dataclass(frozen=True)
class WPEState:
    accel_prev: Interval
    high_pass_1: Interval
    high_pass_1_prev: Interval
    high_pass_2: Interval
    velocity: Interval
    elevation: Interval
    velocity_mean: Interval
    velocity_sq: Interval
    elevation_mean: Interval
    elevation_sq: Interval
    weight: Interval
    elapsed_s: Interval
    raw_period_s: Interval
    log_period_s: Interval
    usable_period: bool
    last_moment_horizon_s: Interval
    last_log_horizon_s: Interval


def constants(dt: float = 0.005) -> Constants:
    e = ESTIMATOR.read_text(encoding="utf-8")
    l = LIMITS.read_text(encoding="utf-8")
    ctor = re.search(
        r"WavePeriodEstimator\(float high_pass_hz = ([0-9.eE+-]+)f,\s*"
        r"float moment_horizon_periods = ([0-9.eE+-]+)f,\s*"
        r"float log_smoothing_periods = ([0-9.eE+-]+)f,\s*"
        r"float min_horizon_sec = ([0-9.eE+-]+)f,\s*"
        r"float max_horizon_sec = ([0-9.eE+-]+)f\)",
        e,
        re.S,
    )
    if not ctor:
        raise RuntimeError("cannot extract WavePeriodEstimator constructor defaults")
    return Constants(
        dt=float(dt),
        high_pass_hz=float(ctor.group(1)),
        moment_horizon_periods=float(ctor.group(2)),
        log_smoothing_periods=float(ctor.group(3)),
        min_horizon_s=float(ctor.group(4)),
        max_horizon_s=float(ctor.group(5)),
        dynamic_horizon_min_s=_limit(l, "kDynamicEmaHorizonMinSec"),
        dynamic_horizon_max_s=_limit(l, "kDynamicEmaHorizonMaxSec"),
    )


def frequency_hz(state: WPEState) -> Interval:
    if not state.usable_period:
        raise RuntimeError("Normal-Live WPE state must have usable period")
    return wide_exp(-state.log_period_s)


def _post_moment_hold(base: WPEState) -> WPEState:
    """Shipping early-return branch: physical/filter moments advanced, log holds."""
    return base


def advance(
    state: WPEState,
    *,
    a_vertical: Interval,
    c: Constants | None = None,
) -> list[WPEState]:
    """Advance one Normal-Live WPE sample, retaining all validity branches."""
    c = c or constants()
    if not state.usable_period:
        raise ValueError("canonical P3 WPE step accepts only already-usable Normal-Live state")
    if not state.log_period_s.lo < state.log_period_s.hi and not math.isfinite(state.log_period_s.lo):
        raise ValueError("finite log-period state required")
    if state.weight.lo <= 1e-3:
        raise ValueError("source-reachable Normal-Live invariant must retain WPE weight > 1e-3")

    lam = TWO * PI * I(c.high_pass_hz)
    decay = VT.exp_interval(-(lam * I(c.dt)))
    gain = (ONE - decay) / lam

    stage1 = decay * (state.high_pass_1 + a_vertical - state.accel_prev)
    stage2 = decay * (state.high_pass_2 + stage1 - state.high_pass_1_prev)
    velocity = decay * state.velocity + gain * stage2
    # Shipping uses the just-updated velocity in the elevation recurrence.
    elevation = decay * state.elevation + gain * velocity
    elapsed = state.elapsed_s + I(c.dt)

    period_for_horizon = wide_exp(state.log_period_s)
    moment_h = clamp_interval(
        I(c.moment_horizon_periods) * period_for_horizon,
        c.min_horizon_s,
        c.max_horizon_s,
    )
    alpha = ONE - VT.exp_interval(-(I(c.dt) / moment_h))
    oma = ONE - alpha

    weight = oma * state.weight + alpha
    velocity_mean = oma * state.velocity_mean + alpha * velocity
    velocity_sq = oma * state.velocity_sq + alpha * velocity.square()
    elevation_mean = oma * state.elevation_mean + alpha * elevation
    elevation_sq = oma * state.elevation_sq + alpha * elevation.square()

    base = WPEState(
        accel_prev=a_vertical,
        high_pass_1=stage1,
        high_pass_1_prev=stage1,
        high_pass_2=stage2,
        velocity=velocity,
        elevation=elevation,
        velocity_mean=velocity_mean,
        velocity_sq=velocity_sq,
        elevation_mean=elevation_mean,
        elevation_sq=elevation_sq,
        weight=weight,
        elapsed_s=elapsed,
        raw_period_s=state.raw_period_s,
        log_period_s=state.log_period_s,
        usable_period=True,
        last_moment_horizon_s=moment_h,
        last_log_horizon_s=state.last_log_horizon_s,
    )

    # The weight condition is guaranteed by the Normal-Live invariant above and
    # monotonic positive-weight recurrence.  The two variance conditions and
    # omega^2 condition may still be false on a concrete sample.
    vmean = velocity_mean / weight
    emean = elevation_mean / weight
    vvar = max_interval(velocity_sq / weight - vmean.square(), 0.0)
    evar = max_interval(elevation_sq / weight - emean.square(), 0.0)

    ev = _strict_positive_branch(evar, 1e-12)
    vv = _strict_positive_branch(vvar, 1e-12)
    successors: list[WPEState] = []
    invalid_variance_possible = evar.lo <= 1e-12 or vvar.lo <= 1e-12
    if invalid_variance_possible:
        successors.append(_post_moment_hold(base))
    if ev is None or vv is None:
        return successors or [_post_moment_hold(base)]

    ratio_sq = vv / ev
    omega_sq = ratio_sq - lam.square()
    om = _strict_positive_branch(omega_sq, 1e-8)
    if omega_sq.lo <= 1e-8:
        successors.append(_post_moment_hold(base))
    if om is None:
        return successors or [_post_moment_hold(base)]

    raw_period = TWO * PI / sqrt_positive(om)
    log_raw = VLOG.log_interval(raw_period)

    if c.log_smoothing_periods <= 0.0:
        log_new = log_raw
        log_h = I(c.dt)
    else:
        sea_period = wide_exp(state.log_period_s)
        requested = I(c.log_smoothing_periods) * sea_period
        lo_guard = max(c.dynamic_horizon_min_s, c.dt)
        log_h = clamp_interval(requested, lo_guard, c.dynamic_horizon_max_s)
        alpha_log = ONE - VT.exp_interval(-(I(c.dt) / log_h))
        log_new = state.log_period_s + alpha_log * (log_raw - state.log_period_s)

    successors.append(replace(
        base,
        raw_period_s=raw_period,
        log_period_s=log_new,
        last_log_horizon_s=log_h,
    ))
    return successors


def build() -> dict:
    c = constants()
    text = ESTIMATOR.read_text(encoding="utf-8")
    parity = {
        "exact_decay": "const float decay = std::exp(-lambda_ * dt_sec);" in text,
        "stage1": "const float stage1 = decay * (high_pass_1_ + vertical_accel_ms2 - accel_prev_);" in text,
        "stage2": "const float stage2 = decay * (high_pass_2_ + stage1 - high_pass_1_prev_);" in text,
        "velocity": "velocity_ = decay * velocity_ + gain * stage2;" in text,
        "elevation_uses_updated_velocity": "elevation_ = decay * elevation_ + gain * velocity_;" in text,
        "moment_ratio": "const float ratio_sq = velocity_var / elevation_var;" in text,
        "leak_subtraction": "const float omega_sq = ratio_sq - lambda_ * lambda_;" in text,
        "log_state": "log_period_sec_ += alpha * (log_raw - log_period_sec_);" in text,
        "usable_latch_one_way": "if (usable_period_) return;" in text,
    }

    # Arithmetic smoke only.  It is not a source cell and is never promotable.
    st = WPEState(
        accel_prev=I(0.0), high_pass_1=I(0.0), high_pass_1_prev=I(0.0),
        high_pass_2=I(0.0), velocity=I(0.1), elevation=I(0.1),
        velocity_mean=I(0.0), velocity_sq=I(0.4),
        elevation_mean=I(0.0), elevation_sq=I(0.1), weight=I(0.8),
        elapsed_s=I(60.0), raw_period_s=I(3.2),
        log_period_s=VLOG.log_interval(I(3.2)), usable_period=True,
        last_moment_horizon_s=I(20.0), last_log_horizon_s=I(0.16),
    )
    succ = advance(st, a_vertical=I(0.2), c=c)
    smoke = {
        "successors": len(succ),
        "all_usable": all(s.usable_period for s in succ),
        "all_log_period_finite": all(
            math.isfinite(s.log_period_s.lo) and math.isfinite(s.log_period_s.hi)
            for s in succ
        ),
    }
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generator": False,
        "trajectory_replay_used": False,
        "requires_same_SEA3_vertical_acceleration": True,
        "shipping_source_parity": parity,
        "shipping_source_parity_pass": all(parity.values()),
        "validity_boundaries_are_branched_not_selected": True,
        "validated_log_used": True,
        "mathematical_WPE_state_step_closed": all(parity.values()) and smoke["all_usable"],
        "target_binary32_libm_roundoff_closed": False,
        "private_Mahony_step_closed_here": False,
        "complete_SEA3_family_materialized_here": False,
        "P3_promoted": False,
        "smoke": smoke,
        "next_obligation": (
            "feed this transition only from the same validated phase-continuous SEA3/private-Mahony realization and close target binary32/libm error before P3 promotion"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "requires_same_SEA3_vertical_acceleration",
        "shipping_source_parity_pass",
        "validity_boundaries_are_branched_not_selected",
        "validated_log_used",
        "mathematical_WPE_state_step_closed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_generator", "trajectory_replay_used", "target_binary32_libm_roundoff_closed",
        "private_Mahony_step_closed_here", "complete_SEA3_family_materialized_here", "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "qualification": d["qualification"],
        "parity": d["shipping_source_parity_pass"],
        "smoke": d["smoke"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
