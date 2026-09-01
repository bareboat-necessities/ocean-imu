#!/usr/bin/env python3
"""Trial the shared verified innovation inverse on the exact-node H18 failure.

This diagnostic changes no filter operation and no theorem domain.  It runs the
same exact P2-node-zero, first 0.80-rad Cayley cell and 202-sample H word used by
the correction-range diagnostic.  Only the proof enclosure for S^-1 is changed:
try the verified midpoint-Neumann enclosure first, then fail closed to the
existing S>=R spectral-entry fallback when the q<1 criterion cannot be proved.

The result determines whether the million-radian accelerometer correction is an
artifact of the old inverse enclosure or whether additional state/covariance
subdivision is still required.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import ou3_implementation_word_language as WORDS
import ou3_p4_candidate_full_word as CAND
import ou3_p4_h18_correction_subdivision as SUB
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_source_node_cells as NODES
import ou3_verified_spd_inverse as VINV
from ou3_proof_module_state import preserve_module_bindings
import ou3_p5_full_h_prefix_cells as H

DEFAULT_DOMAIN = SUB.DEFAULT_DOMAIN
SCHEMA = 1


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_node_index: int = 0) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    nodes = NODES.build()
    nf = NODES.validate(nodes)
    if nf:
        raise RuntimeError(f"source nodes invalid: {nf}")
    src = NODES.h18_source_cell(source_node_index, nodes)
    sector = SECTOR.build(path)
    sf = SECTOR.validate(sector)
    words = WORDS.build(path)
    wf = WORDS.validate(words)
    q = float(sector["design_cayley_norm_upper"])
    parent = CAND._ball_box_cover(q, max_box_norm_factor=1.5)[0]
    samples = int(words["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"])

    attempts = []
    original_inverse = H._spd_inverse_enclosure

    def improved(S, R):
        try:
            X, meta = VINV.inverse_enclosure(S)
            attempts.append({"accepted": True, **meta})
            return X, "VERIFIED_MIDPOINT_NEUMANN_INVERSE"
        except VINV.VerifiedInverseFailure as exc:
            attempts.append({"accepted": False, "reason": str(exc)})
            return original_inverse(S, R)

    with preserve_module_bindings():
        H._source_cell = lambda: src
        H._spd_inverse_enclosure = improved
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

    accepted = [x for x in attempts if x["accepted"]]
    rejected = [x for x in attempts if not x["accepted"]]
    backend_counts = Counter()
    if row.get("maximum_correction"):
        backend_counts[row["maximum_correction"].get("inverse_backend", "unknown")] += 1
    failures = [f"sector: {x}" for x in sf] + [f"word: {x}" for x in wf]
    if row.get("first_failure", {}).get("operation") == "unattributed_exception":
        failures.append("verified-inverse trial ended in unattributed exception")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_H18_VERIFIED_INNOVATION_INVERSE_TRIAL",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "outer_angle_rad": float(sector["design_full_attitude_angle_rad"]),
        "source_node_index": int(source_node_index),
        "full_word_samples": samples,
        "verified_inverse_attempt_count": len(attempts),
        "verified_inverse_accept_count": len(accepted),
        "verified_inverse_reject_count": len(rejected),
        "verified_inverse_accepted_metadata": accepted,
        "verified_inverse_rejections": rejected,
        "word_completed_without_correction_range_failure": bool(row.get("completed")),
        "first_failure": row.get("first_failure"),
        "maximum_correction": row.get("maximum_correction"),
        "maximum_correction_backend_counts": dict(backend_counts),
        "existing_inverse_fallback_retained_when_neumann_unproved": True,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "if the verified inverse removes the correction-range obstruction, integrate it into the shared P4 measurement enclosure on both routes; otherwise inspect the rejected q bound and subdivide the covariance/state source that prevents q<1 rather than shrinking the 0.80-rad attitude domain"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in ("source_generated_not_trajectory_fit", "existing_inverse_fallback_retained_when_neumann_unproved"):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("source_replay_used", "filter_changed", "declared_domain_changed", "P4_USABLE_CERTIFICATE_PROMOTED"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("outer_angle_rad", 0.0)) != 0.80:
        f.append("outer angle changed")
    if d.get("source_node_index") != 0 or d.get("full_word_samples") != 202:
        f.append("focused exact-node/full-word contract changed")
    if int(d.get("verified_inverse_attempt_count", 0)) <= 0:
        f.append("verified inverse was never attempted")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain.resolve(), source_node_index=a.source_node_index)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "validation_pass": not vf,
        "attempts": d["verified_inverse_attempt_count"],
        "accepted": d["verified_inverse_accept_count"],
        "rejected": d["verified_inverse_reject_count"],
        "completed": d["word_completed_without_correction_range_failure"],
        "first_failure": d["first_failure"],
        "maximum_correction": d["maximum_correction"],
        "last_rejection": d["verified_inverse_rejections"][-1] if d["verified_inverse_rejections"] else None,
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
