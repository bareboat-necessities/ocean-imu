#!/usr/bin/env python3
"""Small auditable IEEE-754 outward-rounded interval primitives for OU-III proof work.

The theorem backend must never depend on ordinary round-to-nearest arithmetic for
an enclosure claim.  This module keeps the trusted arithmetic layer deliberately
small: binary64 add/subtract/multiply/divide are evaluated once and widened by
one representable number with ``math.nextafter``.  For the IEEE-754 basic
operations used here, that encloses the exact real result of binary64 inputs.

The matrix layer below is intentionally elementary.  It is built entirely from
the scalar interval operations and uses Gershgorin/absolute-row-sum bounds, not
ordinary floating-point eigenvalue or singular-value routines.  Those bounds
can be conservative, but every claimed matrix inequality remains independently
auditable and outward rounded.

Transcendental functions are intentionally absent.  They require a separately
validated implementation rather than an unqualified libm call.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

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

    def abs_upper(self) -> float:
        """Outward upper bound on ``abs(x)`` for every x in this interval."""
        return up(max(abs(self.lo), abs(self.hi)))

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


IntervalMatrix = list[list[Interval]]


def _shape(A: Sequence[Sequence[Interval]]) -> tuple[int, int]:
    rows = len(A)
    cols = len(A[0]) if rows else 0
    if any(len(row) != cols for row in A):
        raise ValueError("ragged interval matrix")
    return rows, cols


def matrix_point(values: Sequence[Sequence[float]]) -> IntervalMatrix:
    """Convert a finite rectangular numeric matrix to point intervals."""
    rows = [list(row) for row in values]
    if rows:
        cols = len(rows[0])
        if any(len(row) != cols for row in rows):
            raise ValueError("ragged point matrix")
    return [[Interval.point(float(x)) for x in row] for row in rows]


def matrix_identity(n: int) -> IntervalMatrix:
    if n < 0:
        raise ValueError("matrix dimension must be nonnegative")
    return [
        [Interval.point(1.0 if i == j else 0.0) for j in range(n)]
        for i in range(n)
    ]


def matrix_transpose(A: Sequence[Sequence[Interval]]) -> IntervalMatrix:
    rows, cols = _shape(A)
    return [[A[i][j] for i in range(rows)] for j in range(cols)]


def matrix_add(A: Sequence[Sequence[Interval]],
               B: Sequence[Sequence[Interval]]) -> IntervalMatrix:
    sa, sb = _shape(A), _shape(B)
    if sa != sb:
        raise ValueError(f"matrix-add shape mismatch {sa} != {sb}")
    return [[A[i][j] + B[i][j] for j in range(sa[1])] for i in range(sa[0])]


def matrix_sub(A: Sequence[Sequence[Interval]],
               B: Sequence[Sequence[Interval]]) -> IntervalMatrix:
    sa, sb = _shape(A), _shape(B)
    if sa != sb:
        raise ValueError(f"matrix-sub shape mismatch {sa} != {sb}")
    return [[A[i][j] - B[i][j] for j in range(sa[1])] for i in range(sa[0])]


def matrix_mul(A: Sequence[Sequence[Interval]],
               B: Sequence[Sequence[Interval]]) -> IntervalMatrix:
    ra, ca = _shape(A)
    rb, cb = _shape(B)
    if ca != rb:
        raise ValueError(f"matrix-mul shape mismatch {(ra, ca)} x {(rb, cb)}")
    out: IntervalMatrix = []
    for i in range(ra):
        row: list[Interval] = []
        for j in range(cb):
            total = Interval.point(0.0)
            for k in range(ca):
                total = total + A[i][k] * B[k][j]
            row.append(total)
        out.append(row)
    return out


def matrix_abs_row_sum_upper(A: Sequence[Sequence[Interval]]) -> float:
    """Validated upper bound on max absolute row sum, hence on ||A||_inf."""
    rows, cols = _shape(A)
    best = 0.0
    for i in range(rows):
        total = 0.0
        for j in range(cols):
            total = up(total + A[i][j].abs_upper())
        best = max(best, total)
    return up(best)


def matrix_abs_col_sum_upper(A: Sequence[Sequence[Interval]]) -> float:
    """Validated upper bound on max absolute column sum, hence on ||A||_1."""
    rows, cols = _shape(A)
    best = 0.0
    for j in range(cols):
        total = 0.0
        for i in range(rows):
            total = up(total + A[i][j].abs_upper())
        best = max(best, total)
    return up(best)


def matrix_spectral_norm_upper(A: Sequence[Sequence[Interval]]) -> float:
    """Validated, sqrt-free upper bound on spectral norm.

    ``||A||_2 <= sqrt(||A||_1 ||A||_inf) <= max(||A||_1,||A||_inf)``.
    The last expression avoids adding a transcendental square-root operation to
    the trusted core.
    """
    return up(max(matrix_abs_row_sum_upper(A), matrix_abs_col_sum_upper(A)))


def symmetric_gershgorin_lower(A: Sequence[Sequence[Interval]]) -> float:
    """Validated lower bound on every eigenvalue of every symmetric A in box.

    The caller is responsible for the semantic fact that concrete matrices are
    symmetric.  Interval off-diagonal entries may be wider/asymmetric as boxes;
    the absolute radius safely covers either orientation.
    """
    n, m = _shape(A)
    if n != m:
        raise ValueError("Gershgorin eigenvalue bound requires square matrix")
    lower = POS_INF
    for i in range(n):
        radius = 0.0
        for j in range(n):
            if i != j:
                radius = up(radius + max(A[i][j].abs_upper(), A[j][i].abs_upper()))
        disc = down(A[i][i].lo - radius)
        lower = min(lower, disc)
    return lower


def symmetric_gershgorin_upper(A: Sequence[Sequence[Interval]]) -> float:
    """Validated upper bound on every eigenvalue of every symmetric A in box."""
    n, m = _shape(A)
    if n != m:
        raise ValueError("Gershgorin eigenvalue bound requires square matrix")
    upper_bound = NEG_INF
    for i in range(n):
        radius = 0.0
        for j in range(n):
            if i != j:
                radius = up(radius + max(A[i][j].abs_upper(), A[j][i].abs_upper()))
        disc = up(A[i][i].hi + radius)
        upper_bound = max(upper_bound, disc)
    return upper_bound


def symmetric_positive_definite_gershgorin(
    A: Sequence[Sequence[Interval]],
) -> tuple[bool, float]:
    """Certify SPD when the outward Gershgorin lower bound is strictly positive."""
    lower = symmetric_gershgorin_lower(A)
    return lower > 0.0, lower
