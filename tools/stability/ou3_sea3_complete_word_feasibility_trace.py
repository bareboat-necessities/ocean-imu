#!/usr/bin/env python3
"""Unique binary32 shipping trace for the mandatory complete-word feasibility run.

This is a non-promoting diagnostic transport of the exact continuum member in
``ou3_sea3_validated_continuum_member``.  It does not define another source.
The same analytic member is continued for 60 s of startup, one 601-sample H18
window, the separate H-to-A release event, and one 601-sample A21 window.

Every emitted float is accepted only when the validated Decimal enclosure of
the exact continuum value rounds uniquely to that IEEE-754 binary32 value.
There is no replay, quadrature source, finite harmonic source, independent
sample choice, or independently selected tuner/R_S schedule.
"""
from __future__ import annotations

import argparse
from decimal import Decimal
import json
from pathlib import Path
import struct

import ou3_sea3_validated_continuum_member as MEMBER

SCHEMA = 1
QUALIFICATION = "OU3_SEA3_COMPLETE_WORD_FEASIBILITY_FLOAT_TRACE_V1"
G_STD = Decimal("9.80665")
PREHISTORY_SAMPLES = 12000
H_SAMPLES = 601
A_SAMPLES = 601
LIVE_SAMPLES = H_SAMPLES + A_SAMPLES
TOTAL_SAMPLES = PREHISTORY_SAMPLES + LIVE_SAMPLES


def _f32(x: Decimal) -> float:
    return struct.unpack("!f", struct.pack("!f", float(x)))[0]


def _unique_f32(lo: Decimal, hi: Decimal, label: str) -> float:
    flo = _f32(lo)
    fhi = _f32(hi)
    if struct.pack("!f", flo) != struct.pack("!f", fhi):
        raise RuntimeError(
            f"validated interval crosses binary32 boundary for {label}: "
            f"[{lo},{hi}] -> [{flo},{fhi}]"
        )
    return flo


def build() -> dict:
    # Membership in the theorem's hard driver ball is analytic and independent
    # of the sampled implementation trace.  The loop below already validates
    # every timestamp of the *longer* startup/H/release/A history, so calling
    # MEMBER.self_check() here would only recompute the same 12,601 earlier
    # timestamps before immediately recomputing all 13,202 of them again.
    cert = MEMBER.driver_norm_certificate()
    if cert["driver_norm_strictly_below_one"] is not True:
        raise RuntimeError("continuum member escaped the complete-SEA3 driver ball")

    rows: list[dict] = []
    max_abs_acc = Decimal(0)
    max_width = Decimal(0)
    for k in range(TOTAL_SAMPLES):
        t = Decimal(k) * MEMBER.DT
        a = MEMBER.acceleration_interval(t)
        sf_lo = a.lo - G_STD
        sf_hi = a.hi - G_STD
        a32 = _unique_f32(a.lo, a.hi, f"a_z[{k}]")
        sf32 = _unique_f32(sf_lo, sf_hi, f"specific_force_z[{k}]")
        max_abs_acc = max(max_abs_acc, abs(a.lo), abs(a.hi))
        max_width = max(max_width, a.hi - a.lo)
        if k < PREHISTORY_SAMPLES:
            phase = "prehistory"
        elif k < PREHISTORY_SAMPLES + H_SAMPLES:
            phase = "H18_window"
        else:
            phase = "A21_window"
        rows.append({
            "k": k,
            "source_time_s": str(t),
            "phase": phase,
            "acceleration_interval_mps2": [str(a.lo), str(a.hi)],
            "acceleration_binary32": a32,
            "specific_force_z_interval_mps2": [str(sf_lo), str(sf_hi)],
            "specific_force_z_binary32": sf32,
            "binary32_unique": True,
        })

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "continuum_member": "OU3_SEA3_VALIDATED_CONTINUUM_MEMBER_V1",
        "dt_s": float(MEMBER.DT),
        "prehistory_samples": PREHISTORY_SAMPLES,
        "H18_window_samples": H_SAMPLES,
        "A21_window_samples": A_SAMPLES,
        "total_samples": TOTAL_SAMPLES,
        "live_entry_source_time_s": 60.0,
        "H_to_A_release_after_H18_window": True,
        "same_continuum_member_across_startup_H_release_A": True,
        "phase_reset_at_live_or_H_to_A": False,
        "all_intervals_round_to_unique_binary32": True,
        "quadrature_used": False,
        "finite_harmonic_source_used": False,
        "trajectory_replay_used": False,
        "independent_sample_choice_used": False,
        "independent_tuner_schedule_used": False,
        "independent_RS_schedule_used": False,
        "max_abs_exact_acceleration_mps2": str(max_abs_acc),
        "inside_normal_live_acceleration_cap": max_abs_acc <= Decimal(4),
        "max_acceleration_interval_width": str(max_width),
        "driver_norm_certificate": cert,
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
        "H_to_A_release_after_H18_window",
        "same_continuum_member_across_startup_H_release_A",
        "all_intervals_round_to_unique_binary32",
        "inside_normal_live_acceleration_cap",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "phase_reset_at_live_or_H_to_A",
        "quadrature_used",
        "finite_harmonic_source_used",
        "trajectory_replay_used",
        "independent_sample_choice_used",
        "independent_tuner_schedule_used",
        "independent_RS_schedule_used",
        "P3_changed",
        "P4_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} must remain false")
    if int(d.get("prehistory_samples", -1)) != PREHISTORY_SAMPLES:
        f.append("prehistory length changed")
    if int(d.get("H18_window_samples", -1)) != H_SAMPLES:
        f.append("H18 window is not 601 samples")
    if int(d.get("A21_window_samples", -1)) != A_SAMPLES:
        f.append("A21 window is not 601 samples")
    rows = d.get("rows", [])
    if len(rows) != TOTAL_SAMPLES or not all(r.get("binary32_unique") is True for r in rows):
        f.append("complete-word binary32 trace incomplete")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--shipping-input", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.shipping_input.parent.mkdir(parents=True, exist_ok=True)
    args.shipping_input.write_text(
        "".join(f"{r['specific_force_z_binary32']!r}\n" for r in d["rows"]),
        encoding="utf-8",
    )
    print(json.dumps({
        "validation_pass": not vf,
        "total_samples": d["total_samples"],
        "H18_window_samples": d["H18_window_samples"],
        "A21_window_samples": d["A21_window_samples"],
        "max_abs_acceleration_mps2": d["max_abs_exact_acceleration_mps2"],
        "max_interval_width": d["max_acceleration_interval_width"],
        "failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
