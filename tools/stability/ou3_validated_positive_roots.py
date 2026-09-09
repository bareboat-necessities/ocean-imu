#!/usr/bin/env python3
"""Validated positive binary64 root/rational-power endpoint enclosures.

The P4 adaptive source cover needs the deployed SpectralMSE map, including
``sqrt(x)`` and ``x^(6/7)``.  Ordinary libm results are not theorem bounds.
This module uses binary64 only as a search grid: every candidate is converted
exactly to :class:`fractions.Fraction` and accepted only after an exact integer
power comparison proves which side of the real root it lies on.

``math.sqrt``/``math.pow`` may provide an initial search seed, but correctness
does not depend on their rounding.  The loop walks by ``nextafter`` until the
exact inequalities certify adjacent lower/upper binary64 endpoints.
"""
from __future__ import annotations

from fractions import Fraction
import math

from ou3_interval import Interval


def _q(x: float) -> Fraction:
    x=float(x)
    if not (math.isfinite(x) and x >= 0.0):
        raise ValueError("positive-root endpoint must be finite nonnegative")
    return Fraction.from_float(x)


def _bracket_monotone_root(x: float, *, numerator_power: int,
                           root_power: int, seed: float) -> Interval:
    """Enclose y=x^(numerator_power/root_power) by exact comparisons."""
    if not (isinstance(numerator_power,int) and isinstance(root_power,int)
            and numerator_power > 0 and root_power > 0):
        raise ValueError("positive integer powers required")
    x=float(x)
    if x == 0.0:
        return Interval.point(0.0)
    target=_q(x) ** numerator_power
    y=float(seed)
    if not (math.isfinite(y) and y > 0.0):
        y=1.0

    def below_or_equal(v: float) -> bool:
        return Fraction.from_float(v) ** root_power <= target

    # Find one certified lower grid point.
    if below_or_equal(y):
        lo=y
        while True:
            nxt=math.nextafter(lo, math.inf)
            if not math.isfinite(nxt) or not below_or_equal(nxt):
                break
            lo=nxt
    else:
        while not below_or_equal(y):
            nxt=math.nextafter(y, -math.inf)
            if nxt <= 0.0:
                y=0.0
                break
            y=nxt
        lo=y

    if Fraction.from_float(lo) ** root_power == target:
        return Interval.point(lo)

    hi=math.nextafter(lo, math.inf)
    if not math.isfinite(hi):
        raise OverflowError("root enclosure overflow")
    if Fraction.from_float(hi) ** root_power <= target:
        # This should be unreachable after the upward walk, but retain an exact
        # correction loop so correctness is independent of the seed logic.
        while Fraction.from_float(hi) ** root_power <= target:
            lo=hi
            hi=math.nextafter(hi, math.inf)
            if not math.isfinite(hi):
                raise OverflowError("root enclosure overflow")
    return Interval(lo,hi)


def sqrt_point(x: float) -> Interval:
    x=float(x)
    seed=math.sqrt(x) if x >= 0.0 else math.nan
    return _bracket_monotone_root(x,numerator_power=1,root_power=2,seed=seed)


def pow_6_7_point(x: float) -> Interval:
    x=float(x)
    seed=math.pow(x,6.0/7.0) if x > 0.0 else 0.0
    return _bracket_monotone_root(x,numerator_power=6,root_power=7,seed=seed)


def sqrt_interval(x: Interval) -> Interval:
    if not isinstance(x,Interval) or x.lo < 0.0:
        raise ValueError("sqrt interval must be nonnegative")
    return Interval(sqrt_point(x.lo).lo, sqrt_point(x.hi).hi)


def pow_6_7_interval(x: Interval) -> Interval:
    if not isinstance(x,Interval) or x.lo < 0.0:
        raise ValueError("6/7 power interval must be nonnegative")
    return Interval(pow_6_7_point(x.lo).lo, pow_6_7_point(x.hi).hi)


def validate_self_test() -> list[str]:
    failures=[]
    for x in (0.0,1e-12,0.125,1.0,2.0,10.0,1e6):
        s=sqrt_point(x); p=pow_6_7_point(x)
        qx=_q(x)
        if not (Fraction.from_float(s.lo)**2 <= qx <= Fraction.from_float(s.hi)**2):
            failures.append(f"sqrt bracket failed at {x}")
        tgt=qx**6
        if not (Fraction.from_float(p.lo)**7 <= tgt <= Fraction.from_float(p.hi)**7):
            failures.append(f"6/7 bracket failed at {x}")
    return failures
