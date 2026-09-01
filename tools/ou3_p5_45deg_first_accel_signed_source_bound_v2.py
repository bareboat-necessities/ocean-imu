#!/usr/bin/env python3
"""Schema-binding wrapper for the signed 45 deg first-accelerometer bound.

The V1 numerical producer was written against an early draft name for the P4
candidate table.  The retained entrance producer publishes that table under
P4_complete_word_search.candidate_rows.  This wrapper supplies that exact alias
while V1 executes, then restores the imported producer.  No numerical bound or
proof assumption is changed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p5_45deg_first_accel_signed_source_bound as V1

DEFAULT_DOMAIN = V1.DEFAULT_DOMAIN
SCHEMA = 2


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2,
          tangent_cells: int = V1.DEFAULT_TANGENT_CELLS) -> dict:
    original = V1.ENTRANCE.build

    def bound_entrance(path):
        d = original(path)
        out = dict(d)
        out["P4_complete_word_candidate_sectors"] = list(
            d["P4_complete_word_search"]["candidate_rows"]
        )
        return out

    V1.ENTRANCE.build = bound_entrance
    try:
        out = dict(V1.build(Path(domain_path).resolve(), source_pieces=source_pieces,
                            tangent_cells=tangent_cells))
    finally:
        V1.ENTRANCE.build = original
    out["schema"] = SCHEMA
    out["qualification"] = "OU3_P5_45DEG_FIRST_ACCEL_SIGNED_SOURCE_CORRELATED_BOUND_SCHEMA_BOUND"
    out["entrance_candidate_table_path"] = "P4_complete_word_search.candidate_rows"
    out["entrance_schema_alias_changes_numerics"] = False
    return out


def validate(d: dict) -> list[str]:
    base = dict(d)
    base["schema"] = V1.SCHEMA
    failures = V1.validate(base)
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("entrance_candidate_table_path") != "P4_complete_word_search.candidate_rows":
        failures.append("entrance candidate table path mismatch")
    if d.get("entrance_schema_alias_changes_numerics") is not False:
        failures.append("schema alias changes numerics")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--tangent-cells", type=int, default=V1.DEFAULT_TANGENT_CELLS)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve(), source_pieces=args.source_pieces,
                tangent_cells=args.tangent_cells)
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_45DEG_FIRST_ACCEL_SIGNED_SOURCE_BOUND"],
        "q_pre": out["pre_update_q_upper"],
        "q_scalar": out["sign_agnostic_scalar_post_update_q_upper"],
        "q_signed": out["signed_source_correlated_post_update_q_upper"],
        "improvement_factor": out["q_upper_improvement_factor"],
        "max_d": out["max_signed_decomposition_correction_norm_upper_rad"],
        "min_den": out["minimum_signed_composition_denominator_lower"],
        "returned_to_30deg": out["returned_to_30deg_P4_sector_here"],
        "validation_failures": vf,
        "next": out["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
