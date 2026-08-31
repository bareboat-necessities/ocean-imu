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

For a normalized endpoint comparison

    G = [[A, C], [C^T, B]],

with A >= a I and B >= b I, a sufficient Schur/Young condition is

    ||C||_2 < sqrt(a b).

The producer emits a *lower enclosure* of sqrt(a b), not a nearest-rounded
value, so using the number as an open acceptance budget cannot overstate the
true sufficient-condition threshold.  The actual C bound is deliberately not
guessed here and P4 remains fail-closed.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import ou3_p4_worst_translation_cell as WORST

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def _mul_lower(a: float, b: float) -> float:
    if not (a >= 0.0 and b >= 0.0):
        raise ValueError("nonnegative product required")
    return down(float(a) * float(b))


def _sqrt_lower(x: float) -> float:
    """Binary64 lower enclosure of sqrt(x), verified by exact rational square."""
    x = float(x)
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError("finite positive radicand required")
    y = math.sqrt(x)
    qx = Fraction.from_float(x)
    while Fraction.from_float(y) * Fraction.from_float(y) > qx:
        y = math.nextafter(y, -math.inf)
    return down(y)


def _cross_budget_lower(a: float, b: float) -> float:
    return _sqrt_lower(_mul_lower(a, b))


def _ratio_lower(a: float, b: float) -> float:
    if not (a >= 0.0 and b > 0.0):
        raise ValueError("positive denominator required")
    return down(float(a) / float(b))


def build(translation: dict, domain_path: Path = DEFAULT_DOMAIN) -> dict:
    failures: list[str] = []
    modes = {}
    if translation.get("P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS") != "PASS":
        failures.append("validated complete-word translation input is not PASS")

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

        limiting = (
            "translation_complete_word"
            if d_trans <= d_non
            else "nontranslation_existing_direct"
        )
        diag = min(d_trans, d_non)
        widening = _ratio_lower(diag, d_old)
        cross_budget = _cross_budget_lower(d_trans, d_non)
        modes[mode] = {
            "source_cell": s,
            "validated_complete_word_translation_margin_lower": d_trans,
            "existing_direct_nontranslation_margin_lower": d_non,
            "old_full_single_seed_margin_lower": d_old,
            "diagnostic_blockwise_margin_lower": diag,
            "diagnostic_widening_vs_old_full_margin_lower": widening,
            "post_translation_limiting_block": limiting,
            "translation_still_limits_after_full_word_widening": d_trans <= d_non,
            "normalized_full_state_cross_block_spectral_norm_budget_lower_open": cross_budget,
            # compatibility alias; both names are the same conservative lower
            # threshold and validation recomputes them.
            "normalized_full_state_cross_block_spectral_norm_budget_open": cross_budget,
            "cross_block_budget_outward_lower_enclosed": True,
            "full_state_cross_block_sufficient_condition": (
                "validated ||C||_2 < lower_enclosure(sqrt(delta_translation*delta_nontranslation))"
            ),
            "full_state_cross_block_bound_validated": False,
            "full_state_complete_word_cross_blocks_propagated": False,
            "usable_P4_promoted": False,
        }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_POST_TRANSLATION_FULL_STATE_BOTTLENECK_DIAGNOSTIC",
        "source_only": True,
        "trajectory_replay_used": False,
        "blockwise_min_is_final_certificate": False,
        "cross_block_budget_is_final_certificate": False,
        "modes": modes,
        "P4_USABLE_CERTIFICATE_STATUS": "NOT_ESTABLISHED",
        "next_obligation": (
            "propagate full 18/21-state complete-word Phi/Omega on the recurrent source graph and outward-validate "
            "the normalized translation/nontranslation cross-block spectral norm C below the emitted conservative Schur budget; "
            "then validate the exact finite-angle complete return map on the same reachable cells/paths"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("source_only") is not True or d.get("trajectory_replay_used") is not False:
        f.append("diagnostic is not source-only")
    if d.get("blockwise_min_is_final_certificate") is not False:
        f.append("blockwise diagnostic was promoted as final certificate")
    if d.get("cross_block_budget_is_final_certificate") is not False:
        f.append("cross-block budget was promoted as final certificate")
    if d.get("P4_USABLE_CERTIFICATE_STATUS") != "NOT_ESTABLISHED":
        f.append("diagnostic prematurely promoted P4")

    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        for k in (
            "validated_complete_word_translation_margin_lower",
            "existing_direct_nontranslation_margin_lower",
            "old_full_single_seed_margin_lower",
            "diagnostic_blockwise_margin_lower",
            "diagnostic_widening_vs_old_full_margin_lower",
            "normalized_full_state_cross_block_spectral_norm_budget_lower_open",
        ):
            x = m.get(k)
            if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or float(x) <= 0.0:
                f.append(f"{mode}: invalid {k}")

        try:
            d_trans = float(m["validated_complete_word_translation_margin_lower"])
            d_non = float(m["existing_direct_nontranslation_margin_lower"])
            d_old = float(m["old_full_single_seed_margin_lower"])
            expected_limiting = (
                "translation_complete_word"
                if d_trans <= d_non
                else "nontranslation_existing_direct"
            )
            expected_diag = min(d_trans, d_non)
            expected_widening = _ratio_lower(expected_diag, d_old)
            expected_budget = _cross_budget_lower(d_trans, d_non)
            if m.get("post_translation_limiting_block") != expected_limiting:
                f.append(f"{mode}: limiting block does not match margins")
            if m.get("translation_still_limits_after_full_word_widening") is not (d_trans <= d_non):
                f.append(f"{mode}: translation limiting flag does not match margins")
            if float(m.get("diagnostic_blockwise_margin_lower", math.nan)) != expected_diag:
                f.append(f"{mode}: blockwise margin does not equal min diagonal margin")
            if float(m.get("diagnostic_widening_vs_old_full_margin_lower", math.nan)) != expected_widening:
                f.append(f"{mode}: widening factor does not match conservative ratio")
            if float(m.get("normalized_full_state_cross_block_spectral_norm_budget_lower_open", math.nan)) != expected_budget:
                f.append(f"{mode}: cross-block lower budget does not match conservative Schur target")
            if float(m.get("normalized_full_state_cross_block_spectral_norm_budget_open", math.nan)) != expected_budget:
                f.append(f"{mode}: compatibility cross-block budget does not match conservative Schur target")
        except (KeyError, TypeError, ValueError):
            f.append(f"{mode}: cannot recompute diagnostic invariants")

        if m.get("cross_block_budget_outward_lower_enclosed") is not True:
            f.append(f"{mode}: cross-block budget is not outward lower enclosed")
        if m.get("full_state_cross_block_bound_validated") is not False:
            f.append(f"{mode}: cross-block bound incorrectly claimed")
        if m.get("full_state_complete_word_cross_blocks_propagated") is not False:
            f.append(f"{mode}: cross-block propagation incorrectly claimed")
        if m.get("usable_P4_promoted") is not False:
            f.append(f"{mode}: P4 prematurely promoted")
    return list(dict.fromkeys(f))


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
                "cross_block_budget": d.get("modes", {}).get(mode, {}).get("normalized_full_state_cross_block_spectral_norm_budget_lower_open"),
            }
            for mode in ("H", "A")
        },
        "failures": f,
    }, indent=2, sort_keys=True))
    return 0 if not f else 2


if __name__ == "__main__":
    raise SystemExit(main())
