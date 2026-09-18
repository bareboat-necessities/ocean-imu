#!/usr/bin/env python3
"""Validated enclosure of Mahony_AHRS<float>::invSqrt.

For positive binary32 input ``number`` the shipping code performs

    i = 0x5f375a86 - (bits(number) >> 1)
    y = reinterpret_float(i)
    y = y * (1.5f - number * 0.5f * y * y)

The positive binary32 bit encoding is numerically ordered. On any contiguous
input-bit cell, the magic-integer map is monotone in the opposite direction, so
its exact initial range is obtained from the two bit endpoints. These helpers
partition an input interval, propagate each cell through outward binary32
operations, and hull the results. They are independent arithmetic utilities:
physical source ranges and theorem-domain membership are supplied by callers.
"""
from __future__ import annotations

import math

from ou3_interval import Interval, hull
import ou3_binary32_interval as F32

MAGIC = 0x5F375A86


def shipping_point(number: float) -> float:
    x = F32.f32(number)
    bits = F32.f32_bits(x)
    y = F32.f32_from_bits(MAGIC - (bits >> 1))
    t = F32.f32(F32.f32(F32.f32(x * F32.f32(0.5)) * y) * y)
    return F32.f32(y * F32.f32(F32.f32(1.5) - t))


def _cell_inverse_sqrt(start: int, end: int) -> tuple[Interval, Interval]:
    x = Interval(float(F32.f32_from_bits(start)), float(F32.f32_from_bits(end)))
    y_lo_bits = MAGIC - (end >> 1)
    y_hi_bits = MAGIC - (start >> 1)
    y0_lo = F32.f32_from_bits(y_lo_bits)
    y0_hi = F32.f32_from_bits(y_hi_bits)
    if not (0.0 < y0_lo <= y0_hi and math.isfinite(y0_hi)):
        raise RuntimeError("fast-invsqrt initial bit image left positive finite range")
    y0 = Interval(float(y0_lo), float(y0_hi))
    half = F32.point(0.5)
    one_half = F32.point(1.5)
    product = F32.mul(F32.mul(F32.mul(x, half), y0), y0)
    correction = F32.sub(one_half, product)
    return x, F32.mul(y0, correction)


def _bit_cells(number: Interval, max_cells: int):
    if max_cells < 1:
        raise ValueError("max_cells must be positive")
    blo, bhi = F32.positive_float_bit_bounds(number)
    count = bhi - blo + 1
    cells = min(max_cells, count)
    stride = max(1, (count + cells - 1) // cells)
    start = blo
    while start <= bhi:
        end = min(bhi, start + stride - 1)
        yield start, end
        start = end + 1


def enclosure(number: Interval, *, max_cells: int = 512) -> Interval:
    pieces = [_cell_inverse_sqrt(start, end)[1]
              for start, end in _bit_cells(number, max_cells)]
    return hull(*pieces)


def normalized_norm2_enclosure(number: Interval, *, max_cells: int = 2048) -> Interval:
    pieces = []
    for start, end in _bit_cells(number, max_cells):
        x, y = _cell_inverse_sqrt(start, end)
        pieces.append(F32.mul(F32.mul(x, y), y))
    return hull(*pieces)


def all_positive_normal_normalized_norm2_enclosure(
    *, mantissa_cells_per_exponent: int = 16
) -> Interval:
    if mantissa_cells_per_exponent < 1:
        raise ValueError("at least one mantissa cell per exponent is required")
    pieces = []
    mant_count = 1 << 23
    stride = (mant_count + mantissa_cells_per_exponent - 1) // mantissa_cells_per_exponent
    for exponent in range(1, 255):
        exp_base = exponent << 23
        mant_start = 0
        while mant_start < mant_count:
            mant_end = min(mant_count - 1, mant_start + stride - 1)
            start = exp_base | mant_start
            end = exp_base | mant_end
            x, y = _cell_inverse_sqrt(start, end)
            pieces.append(F32.mul(F32.mul(x, y), y))
            mant_start = mant_end + 1
    return hull(*pieces)
