#!/usr/bin/env python3
"""Identify the next linear bottleneck after validated full-word translation widening.

This stage does not promote a full-state P4 certificate.  It consumes the
validated one-second translation result and reconstructs the same old limiting
source cell with the direct P3 generalized-matrix backend.  The purpose is to
answer a narrower question before building an 18/21-state Phi/Omega propagator:
after replacing the single-seed translation margin by the complete-word
translation margin, does translation still limit the linear certificate, or do
the attitude/gyro-bias (and active accelerometer-bias) coordinates become the
next quantitative bottleneck?

The comparison is source-only.  The nontranslation value is the already
validated direct generalized lower margin on exactly the same source cell; the
translation value is the outward-validated complete-word result.  Taking their
minimum is only a diagnostic lower bound and is explicitly forbidden from being
used as the final P4 certificate because cross-block complete-word coupling has
not yet been propagated.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_worst_translation_cell as WORST

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"


def build(translation: dict, domain_path: Path = DEFAULT_DOMAIN) -> dict:
    failures: list[str] = []
    modes = {}
    for mode in ("H", "A"):
        c = WORST.build_cell(mode, Path(domain_path).resolve())
        s = WORST.serializable(c)
        t = translation.get("modes", {}).get(mode, {})
        d_trans = float(t.get("complete_word_translation_margin_lower", 0.0))
        d_non = float(s.get("delta_nontranslation_lower", 0.0))
        d_old = float(s.get("delta_full_lower", 0.0))
        if not all(math.isfinite(x) and x > 0.0 for x in (d_trans, d_non, d_old)):
            failures.append(f"{mode}: missing finite positive margin")
            continue
        limiting = "translation_complete_word" if d_trans <= d_non else "nontranslation_existing_direct"
        diag = min(d_trans, d_non)
        modes[mode] = {
            "source_cell": s,
            "validated_complete_word_translation_margin_lower": d_trans,
            "existing_direct_nontranslation_margin_lower": d_non,
            "old_full_single_seed_margin_lower": d_old,
            "diagnostic_blockwise_margin_lower": diag,
            "diagnostic_widening_vs_old_full_margin_lower": diag / d_old,
            "post_translation_limiting_block": limiting,
            "translation_still_limits_after_full_word_widening": d_trans <= d_non,
            "full_state_complete_word_cross_blocks_propagated": False,
            "usable_P4_promoted": False,
        }

    return {
        "qualification": "OU3_P4_POST_TRANSLATION_FULL_STATE_BOTTLENECK_DIAGNOSTIC",
        "source_only": True,
        "trajectory_replay_used": False,
        "blockwise_min_is_final_certificate": False,
        "modes": modes,
        "P4_USABLE_CERTIFICATE_STATUS": "NOT_ESTABLISHED",
        "next_obligation": (
            "propagate full 18/21-state complete-word Phi/Omega, including translation-attitude cross blocks, "
            "on the recurrent worst source cell; use this diagnostic only to prioritize which block requires refinement"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("source_only") is not True or d.get("trajectory_replay_used") is not False:
        f.append("diagnostic is not source-only")
    if d.get("blockwise_min_is_final_certificate") is not False:
        f.append("blockwise diagnostic was promoted as final certificate")
    if d.get("P4_USABLE_CERTIFICATE_STATUS") != "NOT_ESTABLISHED":
        f.append("diagnostic prematurely promoted P4")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        for k in (
            "validated_complete_word_translation_margin_lower",
            "existing_direct_nontranslation_margin_lower",
            "diagnostic_blockwise_margin_lower",
            "diagnostic_widening_vs_old_full_margin_lower",
        ):
            x = m.get(k)
            if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or float(x) <= 0.0:
                f.append(f"{mode}: invalid {k}")
        if m.get("post_translation_limiting_block") not in (
            "translation_complete_word", "nontranslation_existing_direct"
        ):
            f.append(f"{mode}: missing limiting block")
        if m.get("full_state_complete_word_cross_blocks_propagated") is not False:
            f.append(f"{mode}: cross-block propagation incorrectly claimed")
        if m.get("usable_P4_promoted") is not False:
            f.append(f"{mode}: P4 prematurely promoted")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--translation", type=Path, required=True)
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(json.loads(a.translation.read_text(encoding="utf-8")), a.domain)
    f = validate(d)
    d["validation_failures"] = f
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P4_USABLE_CERTIFICATE_STATUS"],
        "modes": {
            mode: {
                "translation": d.get("modes", {}).get(mode, {}).get("validated_complete_word_translation_margin_lower"),
                "nontranslation": d.get("modes", {}).get(mode, {}).get("existing_direct_nontranslation_margin_lower"),
                "diagnostic_margin": d.get("modes", {}).get(mode, {}).get("diagnostic_blockwise_margin_lower"),
                "factor_vs_old": d.get("modes", {}).get(mode, {}).get("diagnostic_widening_vs_old_full_margin_lower"),
                "limiting_block": d.get("modes", {}).get(mode, {}).get("post_translation_limiting_block"),
            }
            for mode in ("H", "A")
        },
        "failures": f,
    }, indent=2, sort_keys=True))
    return 0 if not f else 2


if __name__ == "__main__":
    raise SystemExit(main())
