#!/usr/bin/env python3
"""High-precision outward-rounded algebraic interval primitives for OU-III P3.

The binary64 interval layer is intentionally tiny and auditable, but a 600-sample
full Joseph word can accumulate enough entrywise wrapping that the *box* around a
positive covariance contains non-covariance matrices.  This module provides a
second arithmetic precision lane for the same algebraic proof object.  It does
not change the source/event family or the theorem inequality.

Only basic algebra is implemented here.  Inputs are converted *exactly* from the
binary64 endpoints already certified by the source/model producers, then every
+,-,*,/ operation is rounded outwards with Python ``decimal`` at a fixed high
precision.  No transcendental is evaluated here.

This module is suitable for the complete H18/A21 P/Psi/Omega word.  It is not a
scalar surrogate, a midpoint calculation, or a post-hoc floating eigensolve.
"""
from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from typing import Sequence

from ou3_interval import Interval, IntervalMatrix

PRECISION = 90


def _exact_float(x: float) -> Decimal:
    return Decimal.from_float(float(x))


def _op1(x: Decimal, *, rounding: str, op) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = PRECISION
        ctx.rounding = rounding
        return +op(x)


def _op2(a: Decimal, b: Decimal, *, rounding: str, op) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = PRECISION
        ctx.rounding = rounding
        return +op(a, b)


def _add_lo(a: Decimal, b: Decimal) -> Decimal:
    return _op2(a, b, rounding=ROUND_FLOOR, op=lambda x, y: x + y)


def _add_hi(a: Decimal, b: Decimal) -> Decimal:
    return _op2(a, b, rounding=ROUND_CEILING, op=lambda x, y: x + y)


def _sub_lo(a: Decimal, b: Decimal) -> Decimal:
    return _op2(a, b, rounding=ROUND_FLOOR, op=lambda x, y: x - y)


def _sub_hi(a: Decimal, b: Decimal) -> Decimal:
    return _op2(a, b, rounding=ROUND_CEILING, op=lambda x, y: x - y)


def _mul_lo(a: Decimal, b: Decimal) -> Decimal:
    return _op2(a, b, rounding=ROUND_FLOOR, op=lambda x, y: x * y)


def _mul_hi(a: Decimal, b: Decimal) -> Decimal:
    return _op2(a, b, rounding=ROUND_CEILING, op=lambda x, y: x * y)


def _div_lo(a: Decimal, b: Decimal) -> Decimal:
    return _op2(a, b, rounding=ROUND_FLOOR, op=lambda x, y: x / y)


def _div_hi(a: Decimal, b: Decimal) -> Decimal:
    return _op2(a, b, rounding=ROUND_CEILING, op=lambda x, y: x / y)


@dataclass(frozen=True)
class DInterval:
    lo: Decimal
    hi: Decimal

    def __post_init__(self) -> None:
        if self.lo.is_nan() or self.hi.is_nan() or self.lo > self.hi:
            raise ValueError(f"invalid Decimal interval [{self.lo}, {self.hi}]")

    @staticmethod
    def zero() -> "DInterval":
        z = Decimal(0)
        return DInterval(z, z)

    @staticmethod
    def one() -> "DInterval":
        o = Decimal(1)
        return DInterval(o, o)

    @staticmethod
    def from_float_point(x: float) -> "DInterval":
        d = _exact_float(x)
        return DInterval(d, d)

    @staticmethod
    def from_binary_interval(x: Interval) -> "DInterval":
        return DInterval(_exact_float(x.lo), _exact_float(x.hi))

    def is_exact_zero(self) -> bool:
        return self.lo == 0 and self.hi == 0

    def contains_zero(self) -> bool:
        return self.lo <= 0 <= self.hi

    def width(self) -> Decimal:
        return _sub_hi(self.hi, self.lo)

    def midpoint_float(self) -> float:
        with localcontext() as ctx:
            ctx.prec = PRECISION
            m = (self.lo + self.hi) / Decimal(2)
        return float(m)

    def __neg__(self) -> "DInterval":
        return DInterval(-self.hi, -self.lo)

    def __add__(self, other: "DInterval") -> "DInterval":
        return DInterval(_add_lo(self.lo, other.lo), _add_hi(self.hi, other.hi))

    def __sub__(self, other: "DInterval") -> "DInterval":
        return DInterval(_sub_lo(self.lo, other.hi), _sub_hi(self.hi, other.lo))

    def __mul__(self, other: "DInterval") -> "DInterval":
        lo_candidates = (
            _mul_lo(self.lo, other.lo),
            _mul_lo(self.lo, other.hi),
            _mul_lo(self.hi, other.lo),
            _mul_lo(self.hi, other.hi),
        )
        hi_candidates = (
            _mul_hi(self.lo, other.lo),
            _mul_hi(self.lo, other.hi),
            _mul_hi(self.hi, other.lo),
            _mul_hi(self.hi, other.hi),
        )
        return DInterval(min(lo_candidates), max(hi_candidates))

    def reciprocal(self) -> "DInterval":
        if self.contains_zero():
            raise ZeroDivisionError("Decimal interval contains zero")
        lo_candidates = (_div_lo(Decimal(1), self.lo), _div_lo(Decimal(1), self.hi))
        hi_candidates = (_div_hi(Decimal(1), self.lo), _div_hi(Decimal(1), self.hi))
        return DInterval(min(lo_candidates), max(hi_candidates))

    def __truediv__(self, other: "DInterval") -> "DInterval":
        return self * other.reciprocal()

    def square(self) -> "DInterval":
        if self.contains_zero():
            a = _mul_hi(self.lo, self.lo)
            b = _mul_hi(self.hi, self.hi)
            return DInterval(Decimal(0), max(a, b))
        lo_candidates = (_mul_lo(self.lo, self.lo), _mul_lo(self.hi, self.hi))
        hi_candidates = (_mul_hi(self.lo, self.lo), _mul_hi(self.hi, self.hi))
        return DInterval(min(lo_candidates), max(hi_candidates))

    def as_strings(self) -> list[str]:
        return [str(self.lo), str(self.hi)]


DMatrix = list[list[DInterval]]


def shape(A: Sequence[Sequence[DInterval]]) -> tuple[int, int]:
    r = len(A)
    c = len(A[0]) if r else 0
    if any(len(row) != c for row in A):
        raise ValueError("ragged Decimal interval matrix")
    return r, c


def from_binary_matrix(A: Sequence[Sequence[Interval]]) -> DMatrix:
    shape_binary = (len(A), len(A[0]) if A else 0)
    if any(len(row) != shape_binary[1] for row in A):
        raise ValueError("ragged binary interval matrix")
    return [[DInterval.from_binary_interval(x) for x in row] for row in A]


def zero_matrix(r: int, c: int) -> DMatrix:
    return [[DInterval.zero() for _ in range(c)] for _ in range(r)]


def identity(n: int) -> DMatrix:
    return [[DInterval.one() if i == j else DInterval.zero() for j in range(n)] for i in range(n)]


def copy_matrix(A: Sequence[Sequence[DInterval]]) -> DMatrix:
    shape(A)
    return [[x for x in row] for row in A]


def transpose(A: Sequence[Sequence[DInterval]]) -> DMatrix:
    r, c = shape(A)
    return [[A[i][j] for i in range(r)] for j in range(c)]


def add(A: Sequence[Sequence[DInterval]], B: Sequence[Sequence[DInterval]]) -> DMatrix:
    sa, sb = shape(A), shape(B)
    if sa != sb:
        raise ValueError(f"matrix add mismatch {sa} != {sb}")
    return [[A[i][j] + B[i][j] for j in range(sa[1])] for i in range(sa[0])]


def sub(A: Sequence[Sequence[DInterval]], B: Sequence[Sequence[DInterval]]) -> DMatrix:
    sa, sb = shape(A), shape(B)
    if sa != sb:
        raise ValueError(f"matrix sub mismatch {sa} != {sb}")
    return [[A[i][j] - B[i][j] for j in range(sa[1])] for i in range(sa[0])]


def mul(A: Sequence[Sequence[DInterval]], B: Sequence[Sequence[DInterval]]) -> DMatrix:
    ra, ca = shape(A)
    rb, cb = shape(B)
    if ca != rb:
        raise ValueError(f"matrix mul mismatch {(ra, ca)} x {(rb, cb)}")
    out = zero_matrix(ra, cb)
    # Sparse-aware ijk order.  F, H and I-KH remain structurally sparse in the
    # complete word; skipping exact zeros is a large speedup without changing
    # any interval operation.
    for i in range(ra):
        for k in range(ca):
            aik = A[i][k]
            if aik.is_exact_zero():
                continue
            for j in range(cb):
                bkj = B[k][j]
                if bkj.is_exact_zero():
                    continue
                out[i][j] = out[i][j] + aik * bkj
    return out


def symmetric_hull(A: Sequence[Sequence[DInterval]]) -> DMatrix:
    n, m = shape(A)
    if n != m:
        raise ValueError("symmetric hull requires square matrix")
    out = copy_matrix(A)
    for i in range(n):
        for j in range(i + 1, n):
            lo = min(A[i][j].lo, A[j][i].lo)
            hi = max(A[i][j].hi, A[j][i].hi)
            x = DInterval(lo, hi)
            out[i][j] = x
            out[j][i] = x
    return out


def scale(A: Sequence[Sequence[DInterval]], c: float) -> DMatrix:
    ci = DInterval.from_float_point(c)
    return [[ci * x for x in row] for row in A]


class DecimalIntervalPivotError(RuntimeError):
    pass


def solve_gauss_jordan(A: Sequence[Sequence[DInterval]], B: Sequence[Sequence[DInterval]]) -> DMatrix:
    n, m = shape(A)
    rb, cb = shape(B)
    if n != m or rb != n:
        raise ValueError("Decimal interval solve shape mismatch")
    aug = [[x for x in A[i]] + [x for x in B[i]] for i in range(n)]
    width = n + cb
    z, o = DInterval.zero(), DInterval.one()
    for k in range(n):
        pivot = aug[k][k]
        if pivot.contains_zero():
            raise DecimalIntervalPivotError(
                f"pivot {k} contains zero: [{pivot.lo}, {pivot.hi}]"
            )
        for j in range(width):
            if j != k:
                aug[k][j] = aug[k][j] / pivot
        aug[k][k] = o
        for i in range(n):
            if i == k:
                continue
            factor = aug[i][k]
            if factor.is_exact_zero():
                aug[i][k] = z
                continue
            for j in range(width):
                if j != k:
                    aug[i][j] = aug[i][j] - factor * aug[k][j]
            aug[i][k] = z
    return [[aug[i][n + j] for j in range(cb)] for i in range(n)]


def inverse(A: Sequence[Sequence[DInterval]]) -> DMatrix:
    n, m = shape(A)
    if n != m:
        raise ValueError("Decimal interval inverse requires square matrix")
    return solve_gauss_jordan(A, identity(n))


def spd_ldlt(A: Sequence[Sequence[DInterval]]) -> tuple[bool, list[DInterval]]:
    n, m = shape(A)
    if n != m:
        raise ValueError("Decimal interval LDLT requires square matrix")
    L = zero_matrix(n, n)
    d = [DInterval.zero() for _ in range(n)]
    for i in range(n):
        L[i][i] = DInterval.one()
    for j in range(n):
        pivot = A[j][j]
        for k in range(j):
            pivot = pivot - L[j][k].square() * d[k]
        d[j] = pivot
        if not pivot.lo > 0:
            return False, d[: j + 1]
        for i in range(j + 1, n):
            num = A[i][j]
            for k in range(j):
                num = num - L[i][k] * L[j][k] * d[k]
            L[i][j] = num / d[j]
    return True, d


def contains_zero_matrix(A: Sequence[Sequence[DInterval]]) -> bool:
    return all(x.contains_zero() for row in A for x in row)


def midpoint_float_matrix(A: Sequence[Sequence[DInterval]]) -> list[list[float]]:
    shape(A)
    return [[x.midpoint_float() for x in row] for row in A]


def max_width(A: Sequence[Sequence[DInterval]]) -> Decimal:
    shape(A)
    best = Decimal(0)
    for row in A:
        for x in row:
            best = max(best, x.width())
    return best
