#!/usr/bin/env python3
"""Instantiate the OU-III conditional source-word language for the declared proof domain.

The generic source-word contract intentionally stays blocked without a PE recurrence
window.  This producer supplies that window from the versioned deployment theorem
domain and cross-checks every PE number against the vector-UCO theorem contract.
It does not infer any value from replay.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_source_word_theorem_contract as WORDS
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def close(a, b, tol=1.0e-12):
    return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=tol)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("proof domain must not be trajectory fitted")
    live = domain["normal_live"]
    recurrence = float(live["vector_pe_recurrence_window_s"])
    vector = VECTOR.build()
    pe = vector["operating_envelope"]

    failures: list[str] = []
    checks = {
        "specific_force_norm_lower_mps2": pe["specific_force_norm_lower_mps2"],
        "magnetic_vector_norm_lower_uT": pe["magnetic_vector_norm_lower_uT"],
        "vector_sine_separation_lower": pe["vector_sine_separation_lower"],
        "body_rate_norm_upper_deg_s": pe["body_rate_norm_upper_deg_s"],
    }
    for key, expected in checks.items():
        actual = live.get(key)
        if actual is None or not close(actual, expected):
            failures.append(f"proof domain {key}={actual!r} does not match vector-UCO contract {expected!r}")

    words = WORDS.build(recurrence)
    structural = WORDS.validate(words)
    failures.extend(f"word contract: {x}" for x in structural)
    if words.get("conditional_word_language", {}).get("ready") is not True:
        failures.extend(words.get("failures", []))

    return {
        "schema": SCHEMA,
        "qualification": "DECLARED_DOMAIN_SOURCE_COMPLETE_OU3_NORMAL_LIVE_WORD_LANGUAGE",
        "trajectory_fit": False,
        "operating_domain": live,
        "vector_uco_qualification": vector["qualification"],
        "word_contract": words,
        "H_dimension": 18,
        "A_dimension": 21,
        "source_complete_relative_to_declared_theorem_hypotheses": not failures,
        "continuous_word_enclosed": False,
        "nonlinear_word_enclosed": False,
        "theorem_promotion": "NOT_ESTABLISHED",
        "failures": failures,
        "pass": not failures,
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("trajectory_fit") is not False:
        failures.append("word language is trajectory fitted")
    if d.get("source_complete_relative_to_declared_theorem_hypotheses") is not True:
        failures.append("declared-domain word language is not source complete")
    if d.get("pass") is not True:
        failures.append("word-language producer failed")
    if d.get("continuous_word_enclosed") is not False or d.get("nonlinear_word_enclosed") is not False:
        failures.append("language stage must not masquerade as enclosure")
    if d.get("theorem_promotion") != "NOT_ESTABLISHED":
        failures.append("language stage must not promote theorem")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve())
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "pass": d["pass"],
        "recurrence_window_s": d["operating_domain"]["vector_pe_recurrence_window_s"],
        "word_horizon_lower_s": d["word_contract"]["conditional_word_language"]["word_horizon_lower_s"],
        "word_samples_upper": d["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
