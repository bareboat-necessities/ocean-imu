#!/usr/bin/env python3
"""Fast attribution of the first coarse H18 correction-range failure.

This companion to ``ou3_p4_h18_correction_subdivision`` deliberately skips the
P3 metric construction and every Jacobian spectral norm.  The correction-range
question does not depend on either object.  It propagates the exact same H=18
shipping prediction/Joseph/reset state and covariance maps on exact P2 source
node zero and reports the first accepted operation whose interval correction
exceeds the already validated deployed-quaternion range.

The parent state box is the same first 0.80-rad Cayley-ball cover cell used by
the complete H18 screen.  This is diagnostic attribution only: no domain is
changed and no P4/P5 theorem can be promoted here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_implementation_word_language as WORDS
import ou3_p4_candidate_full_word as CAND
import ou3_p4_h18_correction_subdivision as SUB
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_source_node_cells as NODES
from ou3_proof_module_state import preserve_module_bindings
import ou3_p5_full_h_prefix_cells as H

DEFAULT_DOMAIN = SUB.DEFAULT_DOMAIN
SCHEMA = 1


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_node_index: int = 0) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("H18 correction attribution must not be trajectory fitted")

    nodes = NODES.build()
    nf = NODES.validate(nodes)
    if nf:
        raise RuntimeError(f"P2 source-node materialization failed: {nf}")
    src = NODES.h18_source_cell(source_node_index, nodes)
    node = NODES.node(source_node_index, nodes)

    sector = SECTOR.build(path)
    sf = SECTOR.validate(sector)
    words = WORDS.build(path)
    wf = WORDS.validate(words)
    q = float(sector["design_cayley_norm_upper"])
    parent = CAND._ball_box_cover(q, max_box_norm_factor=1.5)[0]
    samples = int(words["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"])

    with preserve_module_bindings():
        H._source_cell = lambda: src
        try:
            row = SUB._run_child(path, domain, src, parent, samples)
        except Exception as exc:
            row = {
                "completed": False,
                "entry_cayley_box": parent,
                "first_failure": {
                    "step": None,
                    "operation": "unattributed_exception",
                    "error": f"{type(exc).__name__}: {exc}",
                },
                "maximum_correction": None,
            }

    failures = [f"sector: {x}" for x in sf] + [f"word: {x}" for x in wf]
    first = row.get("first_failure")
    attributed = bool(
        first is not None
        and first.get("operation") in {
            "mandatory_S", "mandatory_accelerometer", "mandatory_magnetometer"
        }
        and first.get("step") is not None
        and float(first.get("correction_norm_upper_rad", 0.0)) > 6.0
    )
    if first is not None and first.get("operation") == "unattributed_exception":
        failures.append("parent correction failure remained unattributed")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_H18_COARSE_PARENT_CORRECTION_FAILURE_ATTRIBUTION",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "P3_metric_constructed_for_this_diagnostic": False,
        "per_prediction_spectral_norm_computed": False,
        "outer_angle_rad": float(sector["design_full_attitude_angle_rad"]),
        "outer_cayley_norm_upper": q,
        "source_node_index": int(source_node_index),
        "source_node": node,
        "full_word_samples": samples,
        "parent_cell": parent,
        "parent_completed_without_correction_range_failure": bool(row.get("completed")),
        "first_failure": first,
        "maximum_correction": row.get("maximum_correction"),
        "first_failure_is_operation_and_step_attributed": attributed,
        "adaptive_state_cell_subdivision_required": not bool(row.get("completed")),
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "bisect the attributed failing state cell, preserving the 0.80-rad outer domain, and rerun only the failing child family until every accepted correction has a validated norm <= 6 rad"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_H18_COARSE_PARENT_CORRECTION_FAILURE_ATTRIBUTION":
        f.append("wrong qualification")
    for key in ("source_generated_not_trajectory_fit",):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "declared_domain_changed",
        "P3_metric_constructed_for_this_diagnostic", "per_prediction_spectral_norm_computed",
        "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("outer_angle_rad", 0.0)) != 0.80:
        f.append("outer angle is not exactly 0.80 rad")
    if d.get("source_node_index") != 0:
        f.append("focused attribution source node is not zero")
    if int(d.get("full_word_samples", 0)) != 202:
        f.append("full word is not 202 samples")
    if d.get("parent_completed_without_correction_range_failure") is False:
        if d.get("first_failure_is_operation_and_step_attributed") is not True:
            f.append("coarse-parent correction-range failure is not attributed")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve(), source_node_index=args.source_node_index)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "validation_pass": not vf,
        "source_node": d["source_node_index"],
        "samples": d["full_word_samples"],
        "parent_completed": d["parent_completed_without_correction_range_failure"],
        "first_failure": d["first_failure"],
        "maximum_correction": d["maximum_correction"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
