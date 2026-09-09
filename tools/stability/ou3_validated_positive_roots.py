#!/usr/bin/env python3
"""Validated positive binary64 root/rational-power endpoint enclosures.

The P4 adaptive source cover needs the deployed SpectralMSE map, including
``sqrt(x)``, ``x^(6/7)`` and the cached ``q_eff^(1/14)``.  Ordinary libm results
are not theorem bounds.  This module uses binary64 only as a search grid: every
candidate is converted exactly to :class:`fractions.Fraction` and accepted only
after an exact integer-power comparison proves which side of the real root it
lies on.

``math.sqrt``/``math.pow`` may provide an initial search seed, but correctness
does not depend on their rounding.  The loop walks by ``nextafter`` until exact
inequalities certify adjacent lower/upper binary64 endpoints.
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


def rational_power_point(x: float, numerator_power: int, root_power: int) -> Interval:
    """Enclose ``x**(numerator_power/root_power)`` by exact rational tests."""
    if not (isinstance(numerator_power,int) and isinstance(root_power,int)
            and numerator_power > 0 and root_power > 0):
        raise ValueError("positive integer powers required")
    x=float(x)
    if x == 0.0:
        return Interval.point(0.0)
    target=_q(x) ** numerator_power
    seed=math.pow(x,float(numerator_power)/float(root_power))
    y=float(seed)
    if not (math.isfinite(y) and y > 0.0):
        y=1.0

    def below_or_equal(v: float) -> bool:
        return Fraction.from_float(v) ** root_power <= target

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
    while Fraction.from_float(hi) ** root_power <= target:
        lo=hi
        hi=math.nextafter(hi, math.inf)
        if not math.isfinite(hi):
            raise OverflowError("root enclosure overflow")
    return Interval(lo,hi)


def rational_power_interval(x: Interval, numerator_power: int, root_power: int) -> Interval:
    if not isinstance(x,Interval) or x.lo < 0.0:
        raise ValueError("rational-power interval must be nonnegative")
    return Interval(
        rational_power_point(x.lo,numerator_power,root_power).lo,
        rational_power_point(x.hi,numerator_power,root_power).hi,
    )


def sqrt_point(x: float) -> Interval:
    return rational_power_point(x,1,2)


def pow_6_7_point(x: float) -> Interval:
    return rational_power_point(x,6,7)


def root_14_point(x: float) -> Interval:
    return rational_power_point(x,1,14)


def sqrt_interval(x: Interval) -> Interval:
    return rational_power_interval(x,1,2)


def pow_6_7_interval(x: Interval) -> Interval:
    return rational_power_interval(x,6,7)


def root_14_interval(x: Interval) -> Interval:
    return rational_power_interval(x,1,14)


def validate_self_test() -> list[str]:
    failures=[]
    for x in (0.0,1e-12,0.125,1.0,2.0,10.0,1e6):
        qx=_q(x)
        for p,r,label in ((1,2,"sqrt"),(6,7,"6/7"),(1,14,"1/14")):
            y=rational_power_point(x,p,r)
            target=qx**p
            if not (Fraction.from_float(y.lo)**r <= target <= Fraction.from_float(y.hi)**r):
                failures.append(f"{label} bracket failed at {x}")
    return failures
