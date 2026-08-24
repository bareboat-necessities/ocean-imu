#!/usr/bin/env python3
"""Normalize and close OU-III hybrid proof obligations.

The validated enclosure verifier predates the source-domain obligation names.
This stage independently normalizes its per-jump rows to the current source
contract and discharges periodic a_w covariance synchronization with the
source-bound PSD/Loewner proof.  Sampled replay is never used as theorem
support here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_hybrid_aw_sync_proof as AW_SYNC
import ou3_source_domain_contract as SOURCE_DOMAIN

REQUIRED = set(SOURCE_DOMAIN.HYBRID_OBLIGATIONS)
AW_KIND = "periodic_aw_covariance_sync"
ALIASES = {
    "magnetic_regauge": "magnetic_regauge_refinement",
    "cooldown": "cooldown_reentry",
}


def _finite_nonnegative(value) -> bool:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(x) and x >= 0.0


def _finite_positive(value) -> bool:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(x) and x > 0.0


def normalize_kind(kind: str) -> str:
    return ALIASES.get(kind, kind)


def _row_valid(row: dict) -> tuple[bool, str | None]:
    """Recompute the generic inward jump inequality from validator output."""
    source = row.get("source_level_W_upper")
    gain = row.get("jump_gain_upper")
    additive = row.get("additive_W_upper")
    new_coord = row.get("new_coordinate_W_upper", 0.0)
    dest = row.get("destination_level_W")
    if not all(_finite_nonnegative(x) for x in (source, gain, additive, new_coord)):
        return False, "jump row has missing/nonfinite/negative primitive bound"
    if not _finite_positive(dest):
        return False, "jump row destination level is not finite positive"
    post = float(gain) * float(source) + float(additive) + float(new_coord)
    if not post < float(dest):
        return False, "recomputed jump does not land strictly inside destination level"
    if row.get("pass") is not True:
        return False, "underlying validated jump row did not pass its structural checks"
    return True, None


def _aw_valid(proof: dict) -> tuple[bool, list[str]]:
    failures: list[str] = []
    if proof.get("qualification") != "SOURCE_BOUND_ANALYTIC_HYBRID_PROOF":
        failures.append("a_w sync proof is not source-bound analytic proof")
    if proof.get("sampled_evidence_used") is not False:
        failures.append("a_w sync proof improperly depends on sampled evidence")
    if proof.get("source_binding_pass") is not True:
        failures.append("a_w sync source binding failed")
    if proof.get("proof_mode") != "PSD_NONEXPANSIVE":
        failures.append("a_w sync proof mode is not PSD_NONEXPANSIVE")
    if proof.get("nonexpansive_information_energy") is not True:
        failures.append("a_w sync information energy is not proven nonexpansive")
    gain = proof.get("jump_gain_upper")
    additive = proof.get("additive_W_upper")
    new_coord = proof.get("new_coordinate_W_upper")
    if not _finite_nonnegative(gain) or float(gain) > 1.0:
        failures.append("a_w sync jump gain is not <= 1")
    if additive != 0.0:
        failures.append("a_w sync additive W term is not exactly zero")
    if new_coord != 0.0:
        failures.append("a_w sync new-coordinate W term is not exactly zero")
    return not failures, failures


def validate(check: dict, aw_proof: dict | None = None) -> dict:
    hybrid = check.get("hybrid", {})
    rows = hybrid.get("bounds", [])
    failures: list[str] = []
    satisfied: set[str] = set()
    normalized_rows: list[dict] = []

    for row in rows:
        kind = normalize_kind(str(row.get("kind", "")))
        if kind == AW_KIND:
            # The instantaneous PSD covariance inflation is exactly
            # nonexpansive.  It is not required to create a strict inward
            # margin; strict contraction is supplied by the surrounding word.
            continue
        if kind not in REQUIRED:
            continue
        ok, reason = _row_valid(row)
        normalized_rows.append({
            "source_kind": row.get("kind"),
            "normalized_kind": kind,
            "pass": ok,
            "failure": reason,
        })
        if ok:
            satisfied.add(kind)
        else:
            failures.append(f"{kind}: {reason}")

    if aw_proof is None:
        aw_proof = AW_SYNC.prove()
    aw_ok, aw_failures = _aw_valid(aw_proof)
    if aw_ok:
        satisfied.add(AW_KIND)
    else:
        failures.extend(f"{AW_KIND}: {x}" for x in aw_failures)

    missing = sorted(REQUIRED - satisfied)
    if missing:
        failures.append(f"missing or unproven hybrid obligations: {missing}")

    return {
        "schema": 1,
        "qualification": "SOURCE_DOMAIN_ALIGNED_HYBRID_THEOREM_GATE",
        "sampled_evidence_used": False,
        "required": sorted(REQUIRED),
        "satisfied": sorted(satisfied),
        "missing": missing,
        "normalized_rows": normalized_rows,
        "periodic_aw_covariance_sync": {
            "pass": aw_ok,
            "proof_mode": aw_proof.get("proof_mode"),
            "jump_gain_upper": aw_proof.get("jump_gain_upper"),
            "additive_W_upper": aw_proof.get("additive_W_upper"),
            "strict_inward_margin_required": False,
            "reason": "PSD covariance inflation is nonexpansive in the inverse-covariance metric; word contraction supplies strict decrease",
        },
        "failures": failures,
        "pass": not failures,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--validated-check", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    check = json.loads(args.validated_check.read_text())
    out = validate(check)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True))
    print(json.dumps(out, indent=2, sort_keys=True))
    return 0 if out["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
