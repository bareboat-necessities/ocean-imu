#!/usr/bin/env python3
"""Instantiate the OU-III conditional source-word language for the declared proof domain.

The generic source-word contract intentionally stays blocked without a PE recurrence
window.  This producer supplies that window from the versioned deployment theorem
domain and cross-checks the declared PE envelope against the generic vector-UCO
theorem contract.

The comparison is monotone, not equality based.  A deployment theorem may assume
stronger observability than the generic weakest-source lemma: larger force/magnetic
norm floors, larger vector-separation floor, and a smaller/equal rate ceiling are
all admissible.  It may not silently weaken any of those generic hypotheses.
No value is inferred from replay.
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
SCHEMA = 2


def finite_positive(value) -> bool:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(x) and x > 0.0


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("proof domain must not be trajectory fitted")
    live = domain["normal_live"]
    recurrence = float(live["vector_pe_recurrence_window_s"])
    vector = VECTOR.build()
    pe = vector["operating_envelope"]

    failures: list[str] = []
    relations = {
        "specific_force_norm_lower_mps2": "ge",
        "magnetic_vector_norm_lower_uT": "ge",
        "vector_sine_separation_lower": "ge",
        "body_rate_norm_upper_deg_s": "le",
    }
    comparison = {}
    for key, relation in relations.items():
        actual = live.get(key)
        generic = pe.get(key)
        ok = finite_positive(actual) and finite_positive(generic)
        if ok:
            a = float(actual)
            g = float(generic)
            ok = a >= g if relation == "ge" else a <= g
        comparison[key] = {
            "declared": actual,
            "generic_contract": generic,
            "required_relation": relation,
            "pass": bool(ok),
        }
        if not ok:
            op = ">=" if relation == "ge" else "<="
            failures.append(
                f"proof domain {key}={actual!r} must be {op} generic vector-UCO bound {generic!r}"
            )

    if not finite_positive(recurrence):
        failures.append("PE recurrence window is not finite positive")

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
        "declared_PE_is_at_least_as_strong_as_generic_contract": not any(
            not row["pass"] for row in comparison.values()
        ),
        "PE_monotone_comparison": comparison,
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
    if d.get("declared_PE_is_at_least_as_strong_as_generic_contract") is not True:
        failures.append("declared PE envelope weakens generic vector-UCO contract")
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
        "PE_monotone_comparison": d["PE_monotone_comparison"],
        "word_horizon_lower_s": d["word_contract"]["conditional_word_language"]["word_horizon_lower_s"],
        "word_samples_upper": d["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
