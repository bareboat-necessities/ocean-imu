#!/usr/bin/env python3
"""Exact implementation-input trace induced by the validated continuum SEA3 member.

The physical source remains the continuum member in
``ou3_sea3_validated_continuum_member``.  This file does not discretize that
source into harmonics or define a replacement source family.

For every 200 Hz source time used by the 60 s startup prehistory and the 601
sample Live word, the continuum member supplies an outward Decimal interval for
the exact vertical acceleration.  We require the entire interval to round to
one and the same IEEE-754 binary32 value.  We do the same after subtracting
standard gravity for the actual accelerometer z measurement at identity
attitude.

Consequently the emitted binary32 values are not approximate source samples:
they are the unique inputs that the unchanged float C++ implementation receives
for this exact legal complete-SEA3 continuum member.  If an enclosure ever
straddles a float rounding boundary, generation fails closed instead of choosing
a representative point.
"""
from __future__ import annotations

import argparse
from decimal import Decimal
import json
from pathlib import Path
import struct

import ou3_sea3_validated_continuum_member as MEMBER

SCHEMA = 1
QUALIFICATION = "OU3_SEA3_VALIDATED_CONTINUUM_FLOAT_TRACE_V1"
G_STD = Decimal("9.80665")


def _f32(x: Decimal) -> float:
    """Round a Decimal to IEEE binary32 using host binary64 only as a staging type."""
    return struct.unpack("!f", struct.pack("!f", float(x)))[0]


def _unique_f32(lo: Decimal, hi: Decimal, label: str) -> float:
    flo = _f32(lo)
    fhi = _f32(hi)
    if struct.pack("!f", flo) != struct.pack("!f", fhi):
        raise RuntimeError(
            f"validated interval crosses a binary32 rounding boundary for {label}: "
            f"[{lo},{hi}] -> [{flo},{fhi}]"
        )
    return flo


def build() -> dict:
    membership = MEMBER.self_check()
    if not membership["driver_norm_certificate"]["driver_norm_strictly_below_one"]:
        raise RuntimeError("validated continuum member is not in the complete-SEA3 hard-driver ball")

    pre = int(MEMBER.PREHISTORY_S / MEMBER.DT)
    total = pre + MEMBER.WINDOW_SAMPLES
    rows: list[dict] = []
    max_acc_width = Decimal(0)
    max_sf_width = Decimal(0)
    max_abs_acc = Decimal(0)

    for k in range(total):
        t = Decimal(k) * MEMBER.DT
        a = MEMBER.acceleration_interval(t)
        sf_lo = a.lo - G_STD
        sf_hi = a.hi - G_STD
        a_f32 = _unique_f32(a.lo, a.hi, f"a_z[{k}]")
        sf_f32 = _unique_f32(sf_lo, sf_hi, f"specific_force_z[{k}]")
        max_acc_width = max(max_acc_width, a.hi - a.lo)
        max_sf_width = max(max_sf_width, sf_hi - sf_lo)
        max_abs_acc = max(max_abs_acc, abs(a.lo), abs(a.hi))
        rows.append({
            "k": k,
            "source_time_s": float(t),
            "phase": "prehistory" if k < pre else "live_word",
            "acceleration_interval_mps2": [str(a.lo), str(a.hi)],
            "acceleration_binary32": a_f32,
            "specific_force_z_interval_mps2": [str(sf_lo), str(sf_hi)],
            "specific_force_z_binary32": sf_f32,
            "binary32_unique": True,
        })

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "continuum_member": "OU3_SEA3_VALIDATED_CONTINUUM_MEMBER_V1",
        "prehistory_samples": pre,
        "word_samples": MEMBER.WINDOW_SAMPLES,
        "total_samples": total,
        "dt_s": float(MEMBER.DT),
        "live_entry_source_time_s": float(MEMBER.PREHISTORY_S),
        "all_intervals_round_to_unique_binary32": True,
        "same_continuum_member_before_and_after_live": True,
        "phase_reset_at_live": False,
        "quadrature_used": False,
        "finite_harmonic_source_used": False,
        "trajectory_replay_used": False,
        "independent_sample_choice_used": False,
        "normal_live_acceleration_cap_mps2": 4.0,
        "max_abs_exact_acceleration_mps2": str(max_abs_acc),
        "inside_normal_live_acceleration_cap": max_abs_acc <= Decimal(4),
        "max_acceleration_interval_width": str(max_acc_width),
        "max_specific_force_interval_width": str(max_sf_width),
        "rows": rows,
        "P3_changed": False,
        "P4_promoted": False,
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("trace detached from complete SEA3")
    for key in (
        "all_intervals_round_to_unique_binary32",
        "same_continuum_member_before_and_after_live",
        "inside_normal_live_acceleration_cap",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "phase_reset_at_live",
        "quadrature_used",
        "finite_harmonic_source_used",
        "trajectory_replay_used",
        "independent_sample_choice_used",
        "P3_changed",
        "P4_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} must remain false")
    if int(d.get("prehistory_samples", -1)) != 12000:
        f.append("prehistory is not exactly 60 s at 200 Hz")
    if int(d.get("word_samples", -1)) != 601:
        f.append("word is not 601 samples")
    if int(d.get("total_samples", -1)) != 12601:
        f.append("total trace length mismatch")
    rows = d.get("rows", [])
    if len(rows) != 12601 or not all(row.get("binary32_unique") is True for row in rows):
        f.append("binary32 source trace incomplete")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--shipping-input", type=Path)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.shipping_input is not None:
        args.shipping_input.parent.mkdir(parents=True, exist_ok=True)
        # Python's repr of a binary32 value is a shortest decimal which parses
        # back to the same binary value in C++ binary32 conversion.
        args.shipping_input.write_text(
            "".join(f"{row['specific_force_z_binary32']!r}\n" for row in d["rows"]),
            encoding="utf-8",
        )
    print(json.dumps({
        "validation_pass": not failures,
        "prehistory_samples": d["prehistory_samples"],
        "word_samples": d["word_samples"],
        "all_unique_binary32": d["all_intervals_round_to_unique_binary32"],
        "max_abs_acceleration_mps2": d["max_abs_exact_acceleration_mps2"],
        "max_interval_width": d["max_acceleration_interval_width"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
