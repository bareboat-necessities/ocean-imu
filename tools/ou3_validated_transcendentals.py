#!/usr/bin/env python3
"""Validated scalar transcendental bounds used by the OU-III proof backend.

The routines here avoid trusting platform ``libm`` for theorem promotion. Input
binary64 endpoints are converted to exact rationals, Taylor series are evaluated
with exact ``Fraction`` arithmetic, and the final rational bounds are converted
back to binary64 with directed ``nextafter`` correction.

The first supported domain is the deployed OU one-step range ``0 <= h/tau <= 1``.
That is enough for the 200 Hz OU transition once the source-domain timing guard
is bound explicitly. Larger arguments are rejected rather than silently losing
validation.
"""
from __future__ import annotations

from fractions import Fraction
import math

from ou3_interval import Interval

MAX_EXP_NEG_X = Fraction(1, 1)
SERIES_TERMS = 48


def _q(x: float) -> Fraction:
    x = float(x)
    if not math.isfinite(x):
        raise ValueError("validated transcendental input must be finite")
    return Fraction.from_float(x)


def _down_q(x: Fraction) -> float:
    y = float(x)
    if Fraction.from_float(y) > x:
        y = math.nextafter(y, -math.inf)
    return math.nextafter(y, -math.inf)


def _up_q(x: Fraction) -> float:
    y = float(x)
    if Fraction.from_float(y) < x:
        y = math.nextafter(y, math.inf)
    return math.nextafter(y, math.inf)


def _exp_neg_q_bounds(x: Fraction) -> tuple[Fraction, Fraction]:
    if x < 0 or x > MAX_EXP_NEG_X:
        raise ValueError(f"validated exp(-x) requires 0 <= x <= 1, got {float(x)}")
    term = Fraction(1, 1)
    total = term
    lower = None
    upper = total
    for k in range(1, SERIES_TERMS + 1):
        term = term * x / k
        total = total - term if k & 1 else total + term
        if k & 1:
            lower = total
        else:
            upper = total
    assert lower is not None and upper is not None
    # Alternating Taylor series with decreasing terms on x<=1: odd partial
    # sums are lower bounds and even partial sums are upper bounds.
    return lower, upper


def exp_neg_scalar_bounds(x: float) -> tuple[float, float]:
    lo, hi = _exp_neg_q_bounds(_q(x))
    return _down_q(lo), _up_q(hi)


def exp_neg(x: Interval) -> Interval:
    if x.lo < 0.0 or x.hi > 1.0:
        raise ValueError("validated exp(-x) interval requires [0,1]")
    lo_q, _ = _exp_neg_q_bounds(_q(x.hi))
    _, hi_q = _exp_neg_q_bounds(_q(x.lo))
    return Interval(_down_q(lo_q), _up_q(hi_q))


def expm1_neg(x: Interval) -> Interval:
    e = exp_neg(x)
    return e - Interval.point(1.0)


def _phi_pa_factor_q_bounds(x: Fraction) -> tuple[Fraction, Fraction]:
    elo, ehi = _exp_neg_q_bounds(x)
    return x - 1 + elo, x - 1 + ehi


def _phi_sa_factor_q_bounds(x: Fraction) -> tuple[Fraction, Fraction]:
    elo, ehi = _exp_neg_q_bounds(x)
    base = x * x / 2 - x + 1
    return base - ehi, base - elo


def ou_discrete_coefficients(h: Interval, tau: Interval) -> dict[str, Interval]:
    """Enclose alpha, expm1(-h/tau), phi_pa and phi_Sa.

    For h>=0 and tau>0, the two cancellation-prone coefficient factors are
    monotone increasing in x=h/tau. We therefore evaluate exact-rational
    endpoint Taylor bounds instead of interval-evaluating ``x + expm1(-x)``.
    """
    if h.lo < 0.0 or tau.lo <= 0.0:
        raise ValueError("OU transition requires h>=0 and tau>0")
    x = h / tau
    if x.lo < 0.0 or x.hi > 1.0:
        raise ValueError(f"validated OU transition requires h/tau <= 1, got {x}")

    alpha = exp_neg(x)
    em1 = alpha - Interval.point(1.0)

    pa_lo_q, _ = _phi_pa_factor_q_bounds(_q(x.lo))
    _, pa_hi_q = _phi_pa_factor_q_bounds(_q(x.hi))
    sa_lo_q, _ = _phi_sa_factor_q_bounds(_q(x.lo))
    _, sa_hi_q = _phi_sa_factor_q_bounds(_q(x.hi))
    pa_factor = Interval(_down_q(pa_lo_q), _up_q(pa_hi_q))
    sa_factor = Interval(_down_q(sa_lo_q), _up_q(sa_hi_q))

    tau2 = tau.square()
    tau3 = tau2 * tau
    return {
        "x": x,
        "alpha": alpha,
        "expm1_neg_x": em1,
        "phi_pa": tau2 * pa_factor,
        "phi_Sa": tau3 * sa_factor,
    }
