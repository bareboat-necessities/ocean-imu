#!/usr/bin/env python3
"""Validated enclosure of Mahony_AHRS<float>::invSqrt.

For positive binary32 input ``number`` the shipping code performs

    i = 0x5f375a86 - (bits(number) >> 1)
    y = reinterpret_float(i)
    y = y * (1.5f - number * 0.5f * y * y)

The positive binary32 bit encoding is numerically ordered.  On any contiguous
input-bit cell, the magic-integer map is monotone in the opposite direction, so
its exact y0 range is obtained from the two bit endpoints.  We partition the
input lattice, propagate each cell through outward binary32 operations, and
hull the results.  This encloses the actual source algorithm; no ideal reciprocal
square root or empirical relative-error constant is substituted.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import struct

from ou3_interval import Interval, hull
import ou3_binary32_interval as F32

REPO = Path(__file__).resolve().parents[1]
MAHONY = REPO / "src" / "ahrs" / "Mahony_AHRS.h"
MAGIC = 0x5F375A86
SCHEMA = 1
QUALIFICATION = "OU3_MAHONY_BINARY32_FAST_INVSQRT_INTERVAL"


def _shipping_point(number: float) -> float:
    x = F32.f32(number)
    bits = F32.f32_bits(x)
    y = F32.f32_from_bits(MAGIC - (bits >> 1))
    # Force the source's left-associated binary32 arithmetic.
    t = F32.f32(F32.f32(F32.f32(x * F32.f32(0.5)) * y) * y)
    y = F32.f32(y * F32.f32(F32.f32(1.5) - t))
    return y


def enclosure(number: Interval, *, max_cells: int = 512) -> Interval:
    blo, bhi = F32.positive_float_bit_bounds(number)
    count = bhi - blo + 1
    cells = min(max_cells, count)
    stride = max(1, (count + cells - 1) // cells)
    pieces: list[Interval] = []

    start = blo
    half = F32.point(0.5)
    one_half = F32.point(1.5)
    while start <= bhi:
        end = min(bhi, start + stride - 1)
        x = Interval(float(F32.f32_from_bits(start)), float(F32.f32_from_bits(end)))

        # MAGIC-(bits>>1) decreases as positive input bits increase.
        y_lo_bits = MAGIC - (end >> 1)
        y_hi_bits = MAGIC - (start >> 1)
        y0_lo = F32.f32_from_bits(y_lo_bits)
        y0_hi = F32.f32_from_bits(y_hi_bits)
        if not (0.0 < y0_lo <= y0_hi and math.isfinite(y0_hi)):
            raise RuntimeError("fast-invsqrt initial bit image left positive finite range")
        y0 = Interval(float(y0_lo), float(y0_hi))

        product = F32.mul(F32.mul(F32.mul(x, half), y0), y0)
        correction = F32.sub(one_half, product)
        pieces.append(F32.mul(y0, correction))
        start = end + 1

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
    for x in samples:
        out = enclosure(Interval.outward_bounds(x, x), max_cells=32)
        actual = _shipping_point(x)
        checks.append({"input": x, "actual": actual, "enclosure": out.as_list(), "contains": out.contains(actual)})
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "shipping_source_parity": parity,
        "shipping_source_parity_pass": all(parity.values()),
        "positive_binary32_bit_order_used": True,
        "magic_map_endpoint_monotonicity_used": True,
        "newton_step_binary32_outward": True,
        "ideal_inverse_sqrt_substituted": False,
        "sample_checks": checks,
        "sample_checks_pass": all(x["contains"] for x in checks),
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
        "sample_checks_pass",
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
    print(json.dumps({"parity": d["shipping_source_parity"], "checks": d["sample_checks"], "failures": failures}, indent=2))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
