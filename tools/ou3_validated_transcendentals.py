#!/usr/bin/env python3
"""Validated transcendental enclosures for OU-III proof work.

The deployment proof must not promote bounds obtained from ordinary ``libm``
round-to-nearest calls.  This module supplies the trusted transcendental layer
on top of :mod:`ou3_interval` without using ``math.exp``, ``math.sin`` or
``math.cos`` in a proof calculation.

For |x| <= 1/2, exp(x) is enclosed by an exact-rational Taylor polynomial plus
its Lagrange remainder.  The elementary bound exp(|x|) < 2 on this interval is
sufficient.

For the P4 SO(3) certificate, sin/cos and the Rodrigues kernels

    sinc(x) = sin(x)/x,
    cosc(x) = (1-cos(x))/x^2

are enclosed by exact-rational alternating Taylor polynomials with a rigorous
next-term/Lagrange bound.  Their audited P4 range is |x| <= 4, which covers the
entire geodesic chart [0, pi) and finite corrections used by the certificate.
The removable values sinc(0)=1 and cosc(0)=1/2 are evaluated from their direct
series, so no zero-crossing division is present in the trusted path.

Binary64 inputs are converted exactly to ``Fraction``.  Only after exact
rational lower/upper endpoints have been obtained are they converted back to
binary64, with an explicit nextafter correction if round-to-nearest landed on
the wrong side.
"""
from __future__ import annotations

from fractions import Fraction
import math

from ou3_interval import Interval

MAX_ABS_ARGUMENT = 0.5
MAX_TRIG_ARGUMENT = 4.0
DEFAULT_ORDER = 20
DEFAULT_TRIG_ORDER = 36


def _down_fraction(q: Fraction) -> float:
    """Greatest convenient binary64 lower enclosure of exact rational ``q``."""
    f = float(q)
    if not math.isfinite(f):
        raise OverflowError(f"rational does not fit finite binary64: {q!r}")
    if Fraction.from_float(f) > q:
        f = math.nextafter(f, -math.inf)
    return f


def _up_fraction(q: Fraction) -> float:
    """Smallest convenient binary64 upper enclosure of exact rational ``q``."""
    f = float(q)
    if not math.isfinite(f):
        raise OverflowError(f"rational does not fit finite binary64: {q!r}")
    if Fraction.from_float(f) < q:
        f = math.nextafter(f, math.inf)
    return f


def _check_point(x: float) -> Fraction:
    x = float(x)
    if not math.isfinite(x) or abs(x) > MAX_ABS_ARGUMENT:
        raise ValueError(
            f"validated exponential argument must be finite and |x| <= "
            f"{MAX_ABS_ARGUMENT}, got {x!r}"
        )
    return Fraction.from_float(x)


def _check_trig_point(x: float) -> Fraction:
    x = float(x)
    if not math.isfinite(x) or abs(x) > MAX_TRIG_ARGUMENT:
        raise ValueError(
            f"validated trigonometric argument must be finite and |x| <= "
            f"{MAX_TRIG_ARGUMENT}, got {x!r}"
        )
    return Fraction.from_float(x)


def _remainder_bound(x: Fraction, order: int) -> Fraction:
    if order < 1:
        raise ValueError("Taylor order must be positive")
    return Fraction(2, 1) * abs(x) ** (order + 1) / Fraction(math.factorial(order + 1), 1)


def exp_point(x: float, order: int = DEFAULT_ORDER) -> Interval:
    q = _check_point(x)
    term = Fraction(1, 1)
    total = term
    for n in range(1, order + 1):
        term = term * q / Fraction(n, 1)
        total += term
    rem = _remainder_bound(q, order)
    return Interval(_down_fraction(total - rem), _up_fraction(total + rem))


def expm1_point(x: float, order: int = DEFAULT_ORDER) -> Interval:
    q = _check_point(x)
    term = Fraction(1, 1)
    total = Fraction(0, 1)
    for n in range(1, order + 1):
        term = term * q / Fraction(n, 1)
        total += term
    rem = _remainder_bound(q, order)
    return Interval(_down_fraction(total - rem), _up_fraction(total + rem))


def exp_interval(x: Interval, order: int = DEFAULT_ORDER) -> Interval:
    _check_point(x.lo)
    _check_point(x.hi)
    return Interval(exp_point(x.lo, order).lo, exp_point(x.hi, order).hi)


def expm1_interval(x: Interval, order: int = DEFAULT_ORDER) -> Interval:
    _check_point(x.lo)
    _check_point(x.hi)
    return Interval(expm1_point(x.lo, order).lo, expm1_point(x.hi, order).hi)


def _positive_ou_kernel_point(
    x: float, *, start_order: int, first_sign: int, order: int = DEFAULT_ORDER
) -> Interval:
    q = _check_point(x)
    if q < 0:
        raise ValueError("OU kernel argument must be nonnegative")
    total = Fraction(0, 1)
    for n in range(start_order, order + 1):
        sign = first_sign if (n - start_order) % 2 == 0 else -first_sign
        total += Fraction(sign, 1) * q**n / Fraction(math.factorial(n), 1)
    rem = _remainder_bound(q, order)
    return Interval(_down_fraction(total - rem), _up_fraction(total + rem))


def ou_phi_pa_kernel_point(x: float, order: int = DEFAULT_ORDER) -> Interval:
    return _positive_ou_kernel_point(x, start_order=2, first_sign=1, order=order)


def ou_phi_Sa_kernel_point(x: float, order: int = DEFAULT_ORDER) -> Interval:
    return _positive_ou_kernel_point(x, start_order=3, first_sign=1, order=order)


def ou_phi_pa_kernel_interval(x: Interval, order: int = DEFAULT_ORDER) -> Interval:
    if x.lo < 0.0:
        raise ValueError("OU kernel interval must be nonnegative")
    return Interval(
        ou_phi_pa_kernel_point(x.lo, order).lo,
        ou_phi_pa_kernel_point(x.hi, order).hi,
    )


def ou_phi_Sa_kernel_interval(x: Interval, order: int = DEFAULT_ORDER) -> Interval:
    if x.lo < 0.0:
        raise ValueError("OU kernel interval must be nonnegative")
    return Interval(
        ou_phi_Sa_kernel_point(x.lo, order).lo,
        ou_phi_Sa_kernel_point(x.hi, order).hi,
    )


def _sin_fraction(q: Fraction, order: int) -> tuple[Fraction, Fraction]:
    if order < 3:
        raise ValueError("trigonometric Taylor order must be >= 3")
    total = Fraction(0, 1)
    n = 0
    while 2 * n + 1 <= order:
        k = 2 * n + 1
        total += (-1 if n & 1 else 1) * q**k / Fraction(math.factorial(k), 1)
        n += 1
    next_k = 2 * n + 1
    rem = abs(q) ** next_k / Fraction(math.factorial(next_k), 1)
    return total - rem, total + rem


def _cos_fraction(q: Fraction, order: int) -> tuple[Fraction, Fraction]:
    if order < 2:
        raise ValueError("trigonometric Taylor order must be >= 2")
    total = Fraction(0, 1)
    n = 0
    while 2 * n <= order:
        k = 2 * n
        total += (-1 if n & 1 else 1) * q**k / Fraction(math.factorial(k), 1)
        n += 1
    next_k = 2 * n
    rem = abs(q) ** next_k / Fraction(math.factorial(next_k), 1)
    return total - rem, total + rem


def sin_point(x: float, order: int = DEFAULT_TRIG_ORDER) -> Interval:
    q = _check_trig_point(x)
    lo, hi = _sin_fraction(q, order)
    return Interval(_down_fraction(lo), _up_fraction(hi))


def cos_point(x: float, order: int = DEFAULT_TRIG_ORDER) -> Interval:
    q = _check_trig_point(x)
    lo, hi = _cos_fraction(q, order)
    return Interval(_down_fraction(lo), _up_fraction(hi))


def sinc_point(x: float, order: int = DEFAULT_TRIG_ORDER) -> Interval:
    q = _check_trig_point(x)
    x2 = q * q
    total = Fraction(0, 1)
    n = 0
    while 2 * n <= order:
        total += (-1 if n & 1 else 1) * x2**n / Fraction(math.factorial(2 * n + 1), 1)
        n += 1
    next_term = x2**n / Fraction(math.factorial(2 * n + 1), 1)
    return Interval(_down_fraction(total - abs(next_term)), _up_fraction(total + abs(next_term)))


def cosc_point(x: float, order: int = DEFAULT_TRIG_ORDER) -> Interval:
    q = _check_trig_point(x)
    x2 = q * q
    total = Fraction(0, 1)
    n = 0
    while 2 * n <= order:
        total += (-1 if n & 1 else 1) * x2**n / Fraction(math.factorial(2 * n + 2), 1)
        n += 1
    next_term = x2**n / Fraction(math.factorial(2 * n + 2), 1)
    return Interval(_down_fraction(total - abs(next_term)), _up_fraction(total + abs(next_term)))


def sinc_interval(x: Interval, order: int = DEFAULT_TRIG_ORDER) -> Interval:
    """Enclose sinc(X) for a nonnegative X subset [0,4].

    On [0,pi] sinc is decreasing.  P4 only calls this routine for geodesic or
    correction-norm intervals that are separately certified below pi.  For the
    wider audited numerical range we conservatively include zero and both
    endpoints rather than assert monotonicity beyond pi.
    """
    if x.lo < 0.0 or x.hi > MAX_TRIG_ARGUMENT:
        raise ValueError("sinc interval must lie in [0,4]")
    a, b = sinc_point(x.lo, order), sinc_point(x.hi, order)
    lo = min(a.lo, b.lo, 0.0 if x.hi >= 3.0 else min(a.lo, b.lo))
    hi = max(a.hi, b.hi, 1.0 if x.lo == 0.0 else max(a.hi, b.hi))
    return Interval(lo, hi)


def cosc_interval(x: Interval, order: int = DEFAULT_TRIG_ORDER) -> Interval:
    """Enclose cosc(X) for a nonnegative X subset [0,4].

    cosc is positive and decreasing on [0,pi].  Beyond that P4 does not promote
    a chart, so the conservative hull additionally includes zero and 1/2.
    """
    if x.lo < 0.0 or x.hi > MAX_TRIG_ARGUMENT:
        raise ValueError("cosc interval must lie in [0,4]")
    a, b = cosc_point(x.lo, order), cosc_point(x.hi, order)
    lo = min(a.lo, b.lo, 0.0 if x.hi >= 3.0 else min(a.lo, b.lo))
    hi = max(a.hi, b.hi, 0.5)
    return Interval(lo, hi)
