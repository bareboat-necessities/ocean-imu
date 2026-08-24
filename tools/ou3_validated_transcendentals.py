#!/usr/bin/env python3
"""Validated small-argument transcendental enclosures for OU-III proof work.

The deployment proof must not promote bounds obtained from ordinary ``libm``
round-to-nearest calls.  This module supplies the first transcendental layer on
top of :mod:`ou3_interval` without using ``math.exp`` or ``math.expm1`` in any
proof calculation.

For |x| <= 1/2, exp(x) is enclosed by an exact-rational Taylor polynomial plus
its Lagrange remainder.  The elementary bound exp(|x|) < 2 on this interval is
sufficient: for 0 <= u <= 1/2,

    exp(u) = 1 + u + sum_{n>=2} u^n/n!
           <= 1 + 1/2 + sum_{n>=2} (1/2)^n < 2.

Binary64 inputs are converted exactly to ``Fraction``.  Only after the exact
rational lower/upper endpoints have been obtained are they converted back to
binary64, with an explicit nextafter correction if round-to-nearest landed on
the wrong side.  Thus the enclosure does not rely on platform transcendental
accuracy.

The range is intentionally narrow.  The current nominal 200 Hz OU-III one-step
source box has h/tau <= 0.25.  Larger ranges require a separately audited range
reduction rather than silently extending the trusted core.
"""
from __future__ import annotations

from fractions import Fraction
import math

from ou3_interval import Interval

MAX_ABS_ARGUMENT = 0.5
DEFAULT_ORDER = 20


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
            f"validated transcendental argument must be finite and |x| <= "
            f"{MAX_ABS_ARGUMENT}, got {x!r}"
        )
    return Fraction.from_float(x)


def _remainder_bound(x: Fraction, order: int) -> Fraction:
    if order < 1:
        raise ValueError("Taylor order must be positive")
    # Lagrange remainder <= exp(|x|)|x|^(N+1)/(N+1)! and exp(|x|)<2.
    return (
        Fraction(2, 1)
        * abs(x) ** (order + 1)
        / Fraction(math.factorial(order + 1), 1)
    )


def exp_point(x: float, order: int = DEFAULT_ORDER) -> Interval:
    """Outward enclosure of exp(x) for one binary64 point, without libm exp."""
    q = _check_point(x)
    term = Fraction(1, 1)
    total = term
    for n in range(1, order + 1):
        term = term * q / Fraction(n, 1)
        total += term
    rem = _remainder_bound(q, order)
    return Interval(_down_fraction(total - rem), _up_fraction(total + rem))


def expm1_point(x: float, order: int = DEFAULT_ORDER) -> Interval:
    """Outward enclosure of exp(x)-1 without cancellation or libm expm1."""
    q = _check_point(x)
    term = Fraction(1, 1)
    total = Fraction(0, 1)
    for n in range(1, order + 1):
        term = term * q / Fraction(n, 1)
        total += term
    rem = _remainder_bound(q, order)
    return Interval(_down_fraction(total - rem), _up_fraction(total + rem))


def exp_interval(x: Interval, order: int = DEFAULT_ORDER) -> Interval:
    """Outward enclosure of exp(X), using monotonicity at interval endpoints."""
    _check_point(x.lo)
    _check_point(x.hi)
    lo = exp_point(x.lo, order).lo
    hi = exp_point(x.hi, order).hi
    return Interval(lo, hi)


def expm1_interval(x: Interval, order: int = DEFAULT_ORDER) -> Interval:
    """Outward enclosure of exp(X)-1, using monotonicity at endpoints."""
    _check_point(x.lo)
    _check_point(x.hi)
    lo = expm1_point(x.lo, order).lo
    hi = expm1_point(x.hi, order).hi
    return Interval(lo, hi)


def _positive_ou_kernel_point(
    x: float, *, start_order: int, first_sign: int, order: int = DEFAULT_ORDER
) -> Interval:
    """Enclose a positive OU Taylor tail at one x in [0, 1/2]."""
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
    """Enclose x + exp(-x) - 1 for x in [0, 1/2]."""
    # x + exp(-x) - 1 = x^2/2! - x^3/3! + ...
    return _positive_ou_kernel_point(
        x, start_order=2, first_sign=1, order=order
    )


def ou_phi_Sa_kernel_point(x: float, order: int = DEFAULT_ORDER) -> Interval:
    """Enclose x^2/2 - x + 1 - exp(-x) for x in [0, 1/2]."""
    # = x^3/3! - x^4/4! + ...
    return _positive_ou_kernel_point(
        x, start_order=3, first_sign=1, order=order
    )


def ou_phi_pa_kernel_interval(x: Interval, order: int = DEFAULT_ORDER) -> Interval:
    """Enclose the monotone positive phi_pa kernel over X subset [0,1/2]."""
    if x.lo < 0.0:
        raise ValueError("OU kernel interval must be nonnegative")
    # derivative 1-exp(-x) >= 0.
    return Interval(
        ou_phi_pa_kernel_point(x.lo, order).lo,
        ou_phi_pa_kernel_point(x.hi, order).hi,
    )


def ou_phi_Sa_kernel_interval(x: Interval, order: int = DEFAULT_ORDER) -> Interval:
    """Enclose the monotone positive phi_Sa kernel over X subset [0,1/2]."""
    if x.lo < 0.0:
        raise ValueError("OU kernel interval must be nonnegative")
    # derivative x - 1 + exp(-x) = phi_pa kernel >= 0.
    return Interval(
        ou_phi_Sa_kernel_point(x.lo, order).lo,
        ou_phi_Sa_kernel_point(x.hi, order).hi,
    )
