#!/usr/bin/env python3
"""Correct the sample-1 prefix SO(3)->Cayley composition helper wiring.

The V1 sample-1 diagnostic reached the first sample-1 S update and established
an ~8e-12 rad correction, but then called angle-conversion helpers that are
implemented locally in the preceding post-reset producer rather than exported
as functions.  This wrapper installs the same outward-rounded formulas used by
that producer before evaluating the complete V1 source family.  No numerical
bound, source domain, or deployment gate changes.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_prefix as BASE

DEFAULT_DOMAIN = BASE.DEFAULT_DOMAIN
SCHEMA = 2


def _post_correction_q_upper(q: float, d: float) -> float:
    if not (math.isfinite(q) and q >= 0.0 and math.isfinite(d) and d >= 0.0):
        return math.inf
    theta = BASE.FULL.up(2.0 * math.atan(0.5 * q))
    total = BASE.FULL.up(theta + d)
    if total >= math.pi:
        return math.inf
    return BASE.FULL.up(2.0 * math.tan(0.5 * total))


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2) -> dict:
    old = BASE._post_correction_q_upper
    BASE._post_correction_q_upper = _post_correction_q_upper
    try:
        out = dict(BASE.build(Path(domain_path).resolve(), source_pieces=source_pieces))
    finally:
        BASE._post_correction_q_upper = old
    out["schema"] = SCHEMA
    out["qualification"] = "OU3_P5_SAMPLE1_SOURCE_CORRELATED_S_ACCEL_PREFIX_DIAGNOSTIC_V2"
    out["angle_conversion_wiring_fixed"] = True
    out["post_correction_composition_formula"] = "theta=2*atan(q/2); theta_plus=theta+|d|; q_plus=2*tan(theta_plus/2)"
    out["same_SO3_triangle_cayley_formula_as_first_post_reset"] = True
    return out


def validate(d: dict) -> list[str]:
    base = dict(d)
    base["schema"] = BASE.SCHEMA
    failures = BASE.validate(base)
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("angle_conversion_wiring_fixed") is not True:
        failures.append("sample-1 angle conversion wiring is not fixed")
    if d.get("same_SO3_triangle_cayley_formula_as_first_post_reset") is not True:
        failures.append("sample-1 post-correction composition does not match first post-reset formula")
    if d.get("post_correction_composition_formula") != "theta=2*atan(q/2); theta_plus=theta+|d|; q_plus=2*tan(theta_plus/2)":
        failures.append("sample-1 post-correction formula mismatch")
    witness = d.get("first_failure")
    if witness is not None and "AttributeError" in str(witness.get("reason", "")):
        failures.append("sample-1 diagnostic still failed on helper wiring")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve(), source_pieces=args.source_pieces)
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_SAMPLE1_S_ACCEL_PREFIX_CERTIFICATE"],
        "paths": out["evaluated_sample1_paths"],
        "phase_counts": out["sample1_phase_counts"],
        "inverse_backends": out["inverse_backend_counts"],
        "max_S_d": out["max_sample1_S_correction_norm_upper_rad"],
        "max_q_after_S": out["max_sample1_q_after_S_upper"],
        "max_acc_residual": out["max_sample1_acc_residual_norm_upper_mps2"],
        "max_acc_d": out["max_sample1_acc_correction_norm_upper_rad"],
        "first_failure": out["first_failure"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
