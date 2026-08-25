#!/usr/bin/env python3
"""Close OU-III hybrid proof obligations against the current source contract.

There is one hybrid theorem path.  Validated jump rows must use the current
source-domain obligation names directly.  Tilt reset must exclude discarded
pre-reset tilt energy from its multiplicative gain, cooldown reentry must use a
reachable product of certified word factors, and periodic a_w covariance sync
is discharged by the source-bound PSD/Loewner proof.  No replay or legacy-name
fallback participates in promotion.
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


def _row_valid(row: dict) -> tuple[bool, list[str]]:
    failures: list[str] = []
    source = row.get("source_level_W_upper")
    gain = row.get("jump_gain_upper")
    additive = row.get("additive_W_upper")
    new_coord = row.get("new_coordinate_W_upper", 0.0)
    dest = row.get("destination_level_W")
    if not all(_finite_nonnegative(x) for x in (source, gain, additive, new_coord)):
        failures.append("jump row has missing/nonfinite/negative primitive bound")
    if not _finite_positive(dest):
        failures.append("jump row destination level is not finite positive")
    if all(_finite_nonnegative(x) for x in (source, gain, additive, new_coord)) and _finite_positive(dest):
        post = float(gain) * float(source) + float(additive) + float(new_coord)
        if not post < float(dest):
            failures.append("recomputed jump does not land strictly inside destination level")
    if row.get("pass") is not True:
        failures.append("underlying validated jump row did not pass structural checks")

    kind = str(row.get("kind", ""))
    if kind == "tilt_reset":
        if row.get("discarded_pre_reset_tilt_excluded_from_multiplicative_gain") is not True:
            failures.append("tilt reset charges discarded pre-reset tilt in multiplicative gain")
        if row.get("reset_to_funnel_exact_map") is not True:
            failures.append("tilt reset is not an exact reset-to-funnel bound")
    if kind == "cooldown_reentry":
        if row.get("reachable_word_product_used") is not True:
            failures.append("cooldown reentry does not use reachable word products")
        if row.get("global_worst_word_power_used") is not False:
            failures.append("cooldown reentry uses global worst-word power")
    return not failures, failures


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
    if not _finite_nonnegative(gain) or float(gain) > 1.0:
        failures.append("a_w sync jump gain is not <= 1")
    if proof.get("additive_W_upper") != 0.0:
        failures.append("a_w sync additive W term is not exactly zero")
    if proof.get("new_coordinate_W_upper") != 0.0:
        failures.append("a_w sync new-coordinate W term is not exactly zero")
    return not failures, failures


def validate(check: dict, aw_proof: dict | None = None) -> dict:
    rows = check.get("hybrid", {}).get("bounds", [])
    failures: list[str] = []
    satisfied: set[str] = set()
    normalized_rows: list[dict] = []

    for row in rows:
        kind = str(row.get("kind", ""))
        if kind == AW_KIND:
            continue
        if kind not in REQUIRED:
            # Unknown/legacy kinds do not satisfy a current obligation.
            continue
        ok, reasons = _row_valid(row)
        normalized_rows.append({
            "source_kind": kind,
            "normalized_kind": kind,
            "pass": ok,
            "failures": reasons,
        })
        if ok:
            satisfied.add(kind)
        else:
            failures.extend(f"{kind}: {reason}" for reason in reasons)

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
        "schema": 2,
        "qualification": "SOURCE_DOMAIN_ALIGNED_GROUP_FUNNEL_HYBRID_GATE",
        "sampled_evidence_used": False,
        "legacy_name_aliases_used": False,
        "required": sorted(REQUIRED),
        "satisfied": sorted(satisfied),
        "missing": missing,
        "normalized_rows": normalized_rows,
        "tilt_reset_policy": "discarded pre-reset tilt excluded from multiplicative gain",
        "cooldown_policy": "reachable certified word products only",
        "periodic_aw_covariance_sync": {
            "pass": aw_ok,
            "proof_mode": aw_proof.get("proof_mode"),
            "jump_gain_upper": aw_proof.get("jump_gain_upper"),
            "additive_W_upper": aw_proof.get("additive_W_upper"),
            "strict_inward_margin_required": False,
            "reason": "PSD covariance inflation is nonexpansive; surrounding source word supplies strict decrease",
        },
        "failures": failures,
        "pass": not failures,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--validated-check", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = validate(json.loads(args.validated_check.read_text()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True))
    print(json.dumps(out, indent=2, sort_keys=True))
    return 0 if out["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
