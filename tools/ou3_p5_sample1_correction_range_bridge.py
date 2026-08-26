#!/usr/bin/env python3
"""Close only the P5 sample-1 scalar correction *range* obstruction.

The dependency-preserving repeated tangent core supplies a source-valid upper
bound for the aligned sample-1 accelerometer attitude correction.  The V2
shipping-quaternion primitive independently validates exact homogeneous
quaternion/Cayley composition through 9 rad with radial subdivision above the
old 6-rad monotonicity range.

This bridge proves only that the scalar-core correction family lies inside the
validated deployed-quaternion proof range.  It deliberately does not claim the
complete sample-1 branch: reset, process, tangent-force, sample-1 S, and the
remaining source-family perturbations are still pending.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_repeated_tangent_channel as CORE
import ou3_p5_deployed_quaternion_cayley_cell_v2 as QV2

DEFAULT_DOMAIN = CORE.DEFAULT_DOMAIN
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    core = CORE.build(path, source_pieces=4, source_cell_index=0, p_pieces=32, axial_pieces=32)
    qv2 = QV2.build(path)
    failures = [f"core: {x}" for x in CORE.validate(core)]
    failures += [f"quaternion-v2: {x}" for x in QV2.validate(qv2)]

    dmax = float(core["max_scalar_correction_norm_upper_rad"])
    qmax = float(qv2["maximum_validated_correction_norm_rad"])
    headroom = down(qmax - dmax)
    inside = math.isfinite(dmax) and math.isfinite(qmax) and 0.0 <= dmax < qmax and headroom > 0.0
    if not inside:
        failures.append("scalar sample-1 correction is not inside V2 deployed-quaternion range")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_SCALAR_CORRECTION_RANGE_BRIDGE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "scalar_core_status": core["P5_SAMPLE1_REPEATED_TANGENT_CORE_WITNESS"],
        "scalar_core_is_complete_sample1_certificate": False,
        "scalar_core_correction_norm_upper_rad": dmax,
        "deployed_quaternion_v2_status": qv2["P5_DEPLOYED_QUATERNION_CAYLEY_CELL_V2_PRIMITIVE"],
        "deployed_quaternion_v2_range_upper_rad": qmax,
        "scalar_range_headroom_rad_lower": headroom,
        "scalar_core_inside_validated_quaternion_range": inside,
        "old_six_rad_proof_range_is_active_obstruction": False,
        "reset_process_tangent_force_perturbations_included": False,
        "sample1_S_due_not_due_family_closed_here": False,
        "complete_sample1_branch_closed_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_SCALAR_CORRECTION_RANGE_BRIDGE": "PASS" if not failures else "NOT_ESTABLISHED",
        "next_obligation": "ADD_RESET_PROCESS_TANGENT_FORCE_AND_SAMPLE1_S_PERTURBATIONS_TO_REPEATED_TANGENT_CORE",
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in ("source_generated_not_trajectory_fit", "scalar_core_inside_validated_quaternion_range"):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "scalar_core_is_complete_sample1_certificate",
        "old_six_rad_proof_range_is_active_obstruction",
        "reset_process_tangent_force_perturbations_included",
        "sample1_S_due_not_due_family_closed_here", "complete_sample1_branch_closed_here",
        "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if not float(d.get("scalar_range_headroom_rad_lower", -math.inf)) > 0.0:
        f.append("nonpositive scalar correction range headroom")
    if d.get("deployed_quaternion_v2_status") != "PASS":
        f.append("V2 quaternion prerequisite did not pass")
    if not f and d.get("P5_SAMPLE1_SCALAR_CORRECTION_RANGE_BRIDGE") != "PASS":
        f.append("range bridge did not pass")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    vf = validate(out)
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_SAMPLE1_SCALAR_CORRECTION_RANGE_BRIDGE"],
        "scalar_d_max": out["scalar_core_correction_norm_upper_rad"],
        "validated_range": out["deployed_quaternion_v2_range_upper_rad"],
        "headroom": out["scalar_range_headroom_rad_lower"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
