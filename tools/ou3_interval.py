#!/usr/bin/env python3
"""Small auditable IEEE-754 outward-rounded interval primitives for OU-III proof work.

The theorem backend must never depend on ordinary round-to-nearest arithmetic for
an enclosure claim.  This module keeps the trusted arithmetic layer deliberately
small: binary64 add/subtract/multiply/divide are evaluated once and widened by
one representable number with ``math.nextafter``.  For the IEEE-754 basic
operations used here, that encloses the exact real result of binary64 inputs.

Transcendental functions are intentionally absent.  They require a separately
validated implementation rather than an unqualified libm call.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

NEG_INF = -math.inf
POS_INF = math.inf


def down(x: float) -> float:
    """Move a finite binary64 result one representable value toward -infinity."""
    x = float(x)
    if math.isnan(x):
        raise ValueError("cannot outward-round NaN")
    return x if x == NEG_INF else math.nextafter(x, NEG_INF)


def up(x: float) -> float:
    """Move a finite binary64 result one representable value toward +infinity."""
    x = float(x)
    if math.isnan(x):
        raise ValueError("cannot outward-round NaN")
    return x if x == POS_INF else math.nextafter(x, POS_INF)


@dataclass(frozen=True)
class Interval:
    lo: float
    hi: float

    def __post_init__(self) -> None:
        lo = float(self.lo)
        hi = float(self.hi)
        if math.isnan(lo) or math.isnan(hi) or lo > hi:
            raise ValueError(f"invalid interval [{lo!r}, {hi!r}]")
        object.__setattr__(self, "lo", lo)
        object.__setattr__(self, "hi", hi)

    @staticmethod
    def point(x: float) -> "Interval":
        x = float(x)
        if not math.isfinite(x):
            raise ValueError("point interval requires a finite value")
        return Interval(x, x)

    @staticmethod
    def outward_bounds(lo: float, hi: float) -> "Interval":
        """Enclose two supplied binary64 endpoint values with explicit slack."""
        lo = float(lo)
        hi = float(hi)
        if not (math.isfinite(lo) and math.isfinite(hi)) or lo > hi:
            raise ValueError(f"invalid finite bounds [{lo!r}, {hi!r}]")
        return Interval(down(lo), up(hi))

    def contains(self, x: float) -> bool:
        x = float(x)
        return self.lo <= x <= self.hi

    def contains_interval(self, other: "Interval") -> bool:
        return self.lo <= other.lo and other.hi <= self.hi

    def width(self) -> float:
        return up(self.hi - self.lo)

    def __add__(self, other: "Interval") -> "Interval":
        return Interval(down(self.lo + other.lo), up(self.hi + other.hi))

    def __sub__(self, other: "Interval") -> "Interval":
        return Interval(down(self.lo - other.hi), up(self.hi - other.lo))

    def __neg__(self) -> "Interval":
        return Interval(-self.hi, -self.lo)

    def __mul__(self, other: "Interval") -> "Interval":
        products = (
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        )
        if any(math.isnan(x) for x in products):
            raise ValueError("indeterminate interval multiplication")
        return Interval(down(min(products)), up(max(products)))

    def reciprocal(self) -> "Interval":
        if self.lo <= 0.0 <= self.hi:
            raise ZeroDivisionError("interval contains zero")
        a = 1.0 / self.lo
        b = 1.0 / self.hi
        return Interval(down(min(a, b)), up(max(a, b)))

    def __truediv__(self, other: "Interval") -> "Interval":
        return self * other.reciprocal()

    def square(self) -> "Interval":
        if self.lo <= 0.0 <= self.hi:
            hi = max(self.lo * self.lo, self.hi * self.hi)
            return Interval(0.0, up(hi))
        a = self.lo * self.lo
        b = self.hi * self.hi
        return Interval(down(min(a, b)), up(max(a, b)))

    def as_list(self) -> list[float]:
        return [self.lo, self.hi]


def hull(*intervals: Interval) -> Interval:
    if not intervals:
        raise ValueError("hull requires at least one interval")
    return Interval(min(x.lo for x in intervals), max(x.hi for x in intervals))
