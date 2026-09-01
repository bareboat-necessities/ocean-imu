#!/usr/bin/env python3
"""V3 signed 30 deg first-accelerometer sector with shared-force-magnitude gain.

V2 preserved the shared ``X/(X+lambda)`` dependency inside the ``KH`` ratios but
still obtained the accelerometer gain rows ``K`` from an ordinary interval
quotient in which the specific-force magnitude ``m`` appears in the numerator
and in the denominator.  On the audited live force cells that costs up to a
factor of ``1.78``, and it is what drove the first-accelerometer correction
family past the ``3`` rad monotone Cayley chart: V2 aborted after 13 of 40960
children with ``max_d = 2.6482`` rad and a composition denominator that could
reach zero.

This stage swaps in ``ou3_p4_shared_force_gain``.  Nothing else changes -- same
filter, same declared 30 deg candidate, same 0.3 g startup ``a_w`` envelope,
same source/alignment/force/tangent children, same signed composition -- and
the resulting bounds are pointwise no larger than V2's.  With the slack removed
the complete 30 deg family becomes evaluable for the first time: every child
composes, the correction family stays inside the monotone Cayley chart, and the
certificate reports a finite worst post-update norm instead of an abort.

That worst norm is still far outside the operation-matched outer sector, so the
certificate remains ``NOT_ESTABLISHED``.  The residual gap is reported here as
a measured quantity rather than as a broken enclosure, and
``ou3_p4_first_accel_nuisance_floor.py`` accounts for where it comes from.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p4_30deg_signed_first_accel_sector_v2 as V2
import ou3_p4_candidate_first_accel_range_v3 as RANGE
import ou3_p4_shared_force_gain as SHARED

DEFAULT_DOMAIN = V2.DEFAULT_DOMAIN
SCHEMA = 3


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2,
          alignment_pieces: int = 16, force_magnitude_pieces: int = 4,
          tangent_pieces: int = 32) -> dict:
    old = RANGE._tangent_structured_gain_bounds
    try:
        RANGE._tangent_structured_gain_bounds = SHARED.shared_force_structured_gain_bounds
        out = dict(V2.build(
            Path(domain_path).resolve(),
            source_pieces=source_pieces,
            alignment_pieces=alignment_pieces,
            force_magnitude_pieces=force_magnitude_pieces,
            tangent_pieces=tangent_pieces,
        ))
    finally:
        RANGE._tangent_structured_gain_bounds = old

    out["schema"] = SCHEMA
    out["qualification"] = "OU3_P4_30DEG_SIGNED_FIRST_ACCEL_SHARED_FORCE_GAIN"
    out["shared_force_magnitude_dependency_preserved"] = True
    out["naive_interval_force_magnitude_gain_used"] = False
    out["first_accelerometer_family_completely_evaluated"] = out["first_failure"] is None
    return out


def validate(d: dict) -> list[str]:
    base = dict(d)
    base["schema"] = V2.SCHEMA
    failures = V2.validate(base)
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("shared_force_magnitude_dependency_preserved") is not True:
        failures.append("shared force-magnitude dependency is not preserved")
    if d.get("naive_interval_force_magnitude_gain_used") is not False:
        failures.append("naive interval force-magnitude gain is still active")
    if d.get("first_accelerometer_family_completely_evaluated") is not True:
        failures.append("shared-force gain did not make the 30deg family evaluable")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--alignment-pieces", type=int, default=16)
    ap.add_argument("--force-magnitude-pieces", type=int, default=4)
    ap.add_argument("--tangent-pieces", type=int, default=32)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain.resolve(), source_pieces=a.source_pieces,
              alignment_pieces=a.alignment_pieces,
              force_magnitude_pieces=a.force_magnitude_pieces,
              tangent_pieces=a.tangent_pieces)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P4_30DEG_SIGNED_FIRST_ACCEL_SECTOR_CERTIFICATE"],
        "evaluable": d["first_accelerometer_family_completely_evaluated"],
        "evaluated": d["evaluated_children"],
        "q_pre": d["post_prediction_q_upper"],
        "q_outer": d["operation_matched_outer_q_upper"],
        "max_d": d["max_correction_norm_upper_rad"],
        "max_qplus": d["max_accepted_or_rejected_post_update_q_upper"],
        "min_den": d["minimum_signed_composition_denominator_lower"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
