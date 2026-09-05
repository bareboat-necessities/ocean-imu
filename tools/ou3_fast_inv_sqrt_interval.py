#!/usr/bin/env python3
"""Validated enclosure of Mahony_AHRS<float>::invSqrt.

For positive binary32 input ``number`` the shipping code performs

    i = 0x5f375a86 - (bits(number) >> 1)
    y = reinterpret_float(i)
    y = y * (1.5f - number * 0.5f * y * y)

The positive binary32 bit encoding is numerically ordered. On any contiguous
input-bit cell, the magic-integer map is monotone in the opposite direction, so
its exact y0 range is obtained from the two bit endpoints. We partition the
input lattice, propagate each cell through outward binary32 operations, and
hull the results. This encloses the actual source algorithm; no ideal reciprocal
square root or empirical relative-error constant is substituted.

``normalized_norm2_enclosure`` keeps the input and reciprocal-square-root image
in the same bit cell while enclosing ``x * invSqrt(x)^2``. The all-normal
variant partitions every positive normal binary32 exponent and mantissa range,
so the Mahony norm shell is not tied to a fitted physical-amplitude interval.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, hull
import ou3_binary32_interval as F32

REPO = Path(__file__).resolve().parents[1]
MAHONY = REPO / "src" / "ahrs" / "Mahony_AHRS.h"
MAGIC = 0x5F375A86
SCHEMA = 3
QUALIFICATION = "OU3_MAHONY_BINARY32_FAST_INVSQRT_INTERVAL_V3"


def _shipping_point(number: float) -> float:
    x = F32.f32(number)
    bits = F32.f32_bits(x)
    y = F32.f32_from_bits(MAGIC - (bits >> 1))
    t = F32.f32(F32.f32(F32.f32(x * F32.f32(0.5)) * y) * y)
    y = F32.f32(y * F32.f32(F32.f32(1.5) - t))
    return y


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
    pieces: list[Interval] = []
    for start, end in _bit_cells(number, max_cells):
        _, y = _cell_inverse_sqrt(start, end)
        pieces.append(y)
    return hull(*pieces)


def normalized_norm2_enclosure(number: Interval, *, max_cells: int = 2048) -> Interval:
    """Enclose squared norm after multiplying by shipping ``invSqrt``."""
    pieces: list[Interval] = []
    for start, end in _bit_cells(number, max_cells):
        x, y = _cell_inverse_sqrt(start, end)
        pieces.append(F32.mul(F32.mul(x, y), y))
    return hull(*pieces)


def all_positive_normal_normalized_norm2_enclosure(
    *, mantissa_cells_per_exponent: int = 16,
) -> Interval:
    """Enclose normalization over every positive *normal* binary32 input.

    Mahony's physical accelerometer norm and its quaternion norm are normal
    floats on the declared SEA3 path. We enumerate all 254 finite normal
    exponents and partition each 23-bit mantissa range, retaining local x/y
    dependence in every cell.
    """
    if mantissa_cells_per_exponent < 1:
        raise ValueError("at least one mantissa cell per exponent is required")
    pieces: list[Interval] = []
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


def build() -> dict:
    text = MAHONY.read_text(encoding="utf-8")
    parity = {
        "lomont_magic": "0x5f375a86 - (i >> 1)" in text,
        "single_newton_step": "y = y * (1.5f - (number * 0.5f * y * y));" in text,
        "float_memcpy_bit_cast": "std::memcpy(&i, &y, sizeof(i));" in text,
    }
    samples = (1.0e-6, 0.01, 1.0, 9.80665 * 9.80665, 400.0)
    checks = []
    for requested in samples:
        x = F32.f32(requested)
        out = enclosure(F32.point(x), max_cells=32)
        actual = _shipping_point(x)
        checks.append({
            "requested_input": requested,
            "shipping_binary32_input": x,
            "actual": actual,
            "enclosure": out.as_list(),
            "contains": out.contains(actual),
        })

    force_n2 = Interval(5.80665 * 5.80665, 13.80665 * 13.80665)
    q_guard_n2 = Interval(0.01, 16.0)
    force_shell = normalized_norm2_enclosure(force_n2)
    q_shell = normalized_norm2_enclosure(q_guard_n2)
    all_normal_shell = all_positive_normal_normalized_norm2_enclosure()
    shells = {
        "physical_specific_force_norm2_input": force_n2.as_list(),
        "physical_specific_force_post_normalization_norm2": force_shell.as_list(),
        "quaternion_raw_norm2_guard_input": q_guard_n2.as_list(),
        "quaternion_post_normalization_norm2": q_shell.as_list(),
        "all_positive_normal_binary32_post_normalization_norm2": all_normal_shell.as_list(),
        "all_positive_normal_shell_positive_finite": (
            all_normal_shell.lo > 0.0 and math.isfinite(all_normal_shell.hi)
        ),
        "all_positive_normal_norm_below_1p1": all_normal_shell.hi < 1.21,
        "guard_shells_positive_finite": (
            force_shell.lo > 0.0 and q_shell.lo > 0.0
            and math.isfinite(force_shell.hi) and math.isfinite(q_shell.hi)
        ),
        "guard_shells_norm_below_1p1": force_shell.hi < 1.21 and q_shell.hi < 1.21,
    }
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "shipping_source_parity": parity,
        "shipping_source_parity_pass": all(parity.values()),
        "positive_binary32_bit_order_used": True,
        "magic_map_endpoint_monotonicity_used": True,
        "newton_step_binary32_outward": True,
        "cell_correlated_normalized_norm2_enclosure": True,
        "all_positive_normal_binary32_normalization_shell_enclosed": True,
        "self_test_inputs_quantized_at_shipping_float_boundary": True,
        "ideal_inverse_sqrt_substituted": False,
        "sample_checks": checks,
        "sample_checks_pass": all(x["contains"] for x in checks),
        "normalization_shell_checks": shells,
        "normalization_shell_checks_pass": (
            shells["all_positive_normal_shell_positive_finite"]
            and shells["all_positive_normal_norm_below_1p1"]
            and shells["guard_shells_positive_finite"]
            and shells["guard_shells_norm_below_1p1"]
        ),
        "compiler_reassociation_or_FMA_closed": False,
        "P3_promoted": False,
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        failures.append("schema/qualification mismatch")
    for key in (
        "shipping_source_parity_pass", "positive_binary32_bit_order_used",
        "magic_map_endpoint_monotonicity_used", "newton_step_binary32_outward",
        "cell_correlated_normalized_norm2_enclosure",
        "all_positive_normal_binary32_normalization_shell_enclosed",
        "self_test_inputs_quantized_at_shipping_float_boundary", "sample_checks_pass",
        "normalization_shell_checks_pass",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in ("ideal_inverse_sqrt_substituted", "compiler_reassociation_or_FMA_closed", "P3_promoted"):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "parity": d["shipping_source_parity"],
        "checks": d["sample_checks"],
        "normalization_shells": d["normalization_shell_checks"],
        "failures": failures,
    }, indent=2))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
