#!/usr/bin/env python3
"""Small outward binary32 arithmetic layer for the SEA3 shipping front end.

The C++ measurement-only front end is implemented in ``float``.  A binary64
interval around the corresponding real formula is therefore not, by itself, a
certificate for the shipping computation.  These helpers quantize every basic
operation to binary32 and expand one adjacent float at each endpoint.  The
extra ulp makes the enclosure independent of host tie direction while retaining
the actual single-precision scale.

No claim about compiler reassociation/FMA contraction is made here; callers
must either prove the source evaluation order or retain that as an explicit
promotion blocker.
"""
from __future__ import annotations

import math
import struct

from ou3_interval import Interval

F32_MAX_BITS = 0x7F7FFFFF
F32_POS_INF_BITS = 0x7F800000
F32_SIGN = 0x80000000


def f32(x: float) -> float:
    return struct.unpack("<f", struct.pack("<f", float(x)))[0]


def f32_bits(x: float) -> int:
    return struct.unpack("<I", struct.pack("<f", f32(x)))[0]


def f32_from_bits(bits: int) -> float:
    return struct.unpack("<f", struct.pack("<I", int(bits) & 0xFFFFFFFF))[0]


def next_f32_up(x: float) -> float:
    x = f32(x)
    if math.isnan(x) or x == math.inf:
        return x
    if x == 0.0:
        return f32_from_bits(1)
    bits = f32_bits(x)
    bits = bits - 1 if (bits & F32_SIGN) else bits + 1
    return f32_from_bits(bits)


def next_f32_down(x: float) -> float:
    x = f32(x)
    if math.isnan(x) or x == -math.inf:
        return x
    if x == 0.0:
        return f32_from_bits(F32_SIGN | 1)
    bits = f32_bits(x)
    bits = bits + 1 if (bits & F32_SIGN) else bits - 1
    return f32_from_bits(bits)


def round_interval(x: Interval) -> Interval:
    if not (math.isfinite(x.lo) and math.isfinite(x.hi)):
        raise ValueError("binary32 enclosure requires finite real interval")
    lo = f32(x.lo)
    hi = f32(x.hi)
    if not (math.isfinite(lo) and math.isfinite(hi)):
        raise OverflowError("binary32 operation overflowed")
    # One adjacent float is deliberate: pack/unpack gives the nearest host
    # binary32; expanding once encloses either neighboring answer at a tie and
    # all exact reals between the supplied endpoints by monotonic rounding.
    return Interval(float(next_f32_down(lo)), float(next_f32_up(hi)))


def point(x: float) -> Interval:
    y = f32(x)
    if not math.isfinite(y):
        raise ValueError("finite binary32 point required")
    return Interval(float(y), float(y))


def add(a: Interval, b: Interval) -> Interval:
    return round_interval(a + b)


def sub(a: Interval, b: Interval) -> Interval:
    return round_interval(a - b)


def mul(a: Interval, b: Interval) -> Interval:
    return round_interval(a * b)


def neg(a: Interval) -> Interval:
    return Interval(-a.hi, -a.lo)


def square(a: Interval) -> Interval:
    return mul(a, a)


def sum3(a: Interval, b: Interval, c: Interval) -> Interval:
    return add(add(a, b), c)


def dot3(a: tuple[Interval, Interval, Interval],
         b: tuple[Interval, Interval, Interval]) -> Interval:
    return add(add(mul(a[0], b[0]), mul(a[1], b[1])), mul(a[2], b[2]))


def positive_float_bit_bounds(x: Interval) -> tuple[int, int]:
    """Bit-index bounds for all positive finite binary32 values in ``x``."""
    if not (math.isfinite(x.lo) and math.isfinite(x.hi) and x.hi > 0.0):
        raise ValueError("positive finite interval required")
    lo_real = max(x.lo, f32_from_bits(1))
    lo = f32(lo_real)
    if lo < lo_real:
        lo = next_f32_up(lo)
    hi = f32(x.hi)
    if hi > x.hi:
        hi = next_f32_down(hi)
    if not (lo > 0.0 and hi > 0.0 and math.isfinite(lo) and math.isfinite(hi)):
        raise ValueError("interval contains no positive finite binary32 value")
    blo, bhi = f32_bits(lo), f32_bits(hi)
    if blo > bhi:
        raise ValueError("interval contains no positive binary32 lattice point")
    return blo, bhi
