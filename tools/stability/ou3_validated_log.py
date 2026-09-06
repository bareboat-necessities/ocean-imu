#!/usr/bin/env python3
"""Validated natural-log enclosures for the complete SEA3 source state.

WavePeriodEstimator stores its canonical period in log space.  Canonical P3 may
therefore not call ordinary ``math.log`` and then pretend the rounded result is
an enclosure.  This small trusted layer evaluates log from exact rational input
using

    log(x) = k log(2) + 2 * sum_{j>=0} z^(2j+1)/(2j+1),
    z = (m-1)/(m+1),  x = 2^k m,  1 <= m < 2.

For 1 <= m <= 2, 0 <= z <= 1/3 and every term is nonnegative.  After N terms,
the omitted tail is bounded by

    2 z^(2N+1) / ((2N+1) (1-z^2)).

Binary64 point inputs are converted exactly to ``Fraction`` before range
reduction.  The final rational endpoints are converted back to binary64 with an
explicit directed correction.  No libm logarithm participates in a proof
claim.
"""
from __future__ import annotations

from fractions import Fraction
import math

from ou3_interval import Interval

DEFAULT_TERMS = 40


def _down_fraction(q: Fraction) -> float:
    f = float(q)
    if not math.isfinite(f):
        raise OverflowError(f"rational does not fit finite binary64: {q!r}")
    if Fraction.from_float(f) > q:
        f = math.nextafter(f, -math.inf)
    return f


def _up_fraction(q: Fraction) -> float:
    f = float(q)
    if not math.isfinite(f):
        raise OverflowError(f"rational does not fit finite binary64: {q!r}")
    if Fraction.from_float(f) < q:
        f = math.nextafter(f, math.inf)
    return f


def _unit_log_bounds(q: Fraction, terms: int) -> tuple[Fraction, Fraction]:
    """Bounds log(q) for exact rational 1 <= q <= 2."""
    if not (Fraction(1, 1) <= q <= Fraction(2, 1)):
        raise ValueError("unit logarithm argument must lie in [1,2]")
    if terms < 1:
        raise ValueError("terms must be positive")
    z = (q - 1) / (q + 1)
    z2 = z * z
    power = z
    total = Fraction(0, 1)
    for j in range(terms):
        total += power / Fraction(2 * j + 1, 1)
        power *= z2
    partial = Fraction(2, 1) * total
    first_omitted_degree = 2 * terms + 1
    if z == 0:
        return partial, partial
    tail = (
        Fraction(2, 1)
        * power
        / Fraction(first_omitted_degree, 1)
        / (Fraction(1, 1) - z2)
    )
    return partial, partial + tail


def _log_fraction_bounds(q: Fraction, terms: int) -> tuple[Fraction, Fraction]:
    if q <= 0:
        raise ValueError("logarithm argument must be positive")

    m = q
    k = 0
    while m < 1:
        m *= 2
        k -= 1
    while m >= 2:
        m /= 2
        k += 1

    m_lo, m_hi = _unit_log_bounds(m, terms)
    ln2_lo, ln2_hi = _unit_log_bounds(Fraction(2, 1), terms)
    if k >= 0:
        return m_lo + k * ln2_lo, m_hi + k * ln2_hi
    # Multiplication by a negative integer reverses the log(2) endpoint used
    # for lower/upper enclosure.
    return m_lo + k * ln2_hi, m_hi + k * ln2_lo


def log_point(x: float, terms: int = DEFAULT_TERMS) -> Interval:
    x = float(x)
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError("validated logarithm requires a finite positive point")
    lo, hi = _log_fraction_bounds(Fraction.from_float(x), terms)
    return Interval(_down_fraction(lo), _up_fraction(hi))


def log_interval(x: Interval, terms: int = DEFAULT_TERMS) -> Interval:
    if not (math.isfinite(x.lo) and math.isfinite(x.hi) and x.lo > 0.0):
        raise ValueError("validated logarithm interval must be finite and positive")
    # log is monotone.
    return Interval(log_point(x.lo, terms).lo, log_point(x.hi, terms).hi)
