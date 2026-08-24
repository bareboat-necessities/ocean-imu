#!/usr/bin/env python3
"""Final independent composition gate for the deployed OU-III stability proof.

This sits above the existing deployment-theorem gate.  It adds the two pieces a
Live-only composition cannot establish by itself: source/implementation parity
and the actual pre-Live reset/Mahony/goLive certificate.  It also instantiates
the recurring-PE source-word language from the declared theorem domain.

No upstream PASS bit is sufficient.  The subordinate deployment gate
independently regenerates the source domain and recomputes hybrid, stochastic,
and finite-capture arithmetic; this gate independently rebuilds the source
manifest, startup certificate, and word-language contract.  Only their
conjunction may emit PASS_IMPLEMENTATION_STABLE.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_deployment_gate as DEPLOY
import ou3_implementation_proof_manifest as MANIFEST
import ou3_implementation_word_language as WORDS
import ou3_startup_stability_certificate as STARTUP

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def compose(check: dict, source_domain: dict, primitive: dict,
            domain_path: Path = DEFAULT_DOMAIN) -> dict:
    manifest = MANIFEST.build()
    manifest_failures = MANIFEST.validate(manifest)
    startup = STARTUP.build(domain_path)
    startup_failures = STARTUP.validate(startup)
    words = WORDS.build(domain_path)
    word_failures = WORDS.validate(words)
    deployment = DEPLOY.compose(check, source_domain, primitive)

    downstream_pass = deployment.get("deployment_theorem_certificate") == "PASS"
    final_pass = bool(
        not manifest_failures
        and not startup_failures
        and not word_failures
        and downstream_pass
    )

    failures: list[str] = []
    failures.extend(f"implementation manifest: {x}" for x in manifest_failures)
    failures.extend(f"startup: {x}" for x in startup_failures)
    failures.extend(f"source-word language: {x}" for x in word_failures)
    if not downstream_pass:
        failures.append("continuous nonlinear/hybrid/capture/stochastic deployment certificate did not pass")

    return {
        "schema": SCHEMA,
        "qualification": "INDEPENDENT_ACTUAL_OU3_IMPLEMENTATION_STABILITY_COMPOSITION_GATE",
        "implementation_manifest_pass": not manifest_failures,
        "implementation_bindings": manifest.get("implementation_files"),
        "startup_certificate_pass": not startup_failures,
        "startup": startup,
        "source_complete_word_language_pass": not word_failures,
        "word_language": words,
        "deployment": deployment,
        "downstream_deployment_theorem_pass": downstream_pass,
        "performance_and_replay_role": (
            "retained existing main regression/falsification/evidence layer; not used to derive theorem bounds"
        ),
        "failures": failures,
        "implementation_stability_certificate": (
            "PASS_IMPLEMENTATION_STABLE" if final_pass else "FAIL"
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--validated-check", type=Path, required=True)
    ap.add_argument("--source-domain", type=Path, required=True)
    ap.add_argument("--primitive-bounds", type=Path, required=True)
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    check = json.loads(args.validated_check.read_text(encoding="utf-8"))
    source = json.loads(args.source_domain.read_text(encoding="utf-8"))
    primitive = json.loads(args.primitive_bounds.read_text(encoding="utf-8"))
    out = compose(check, source, primitive, args.domain.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "implementation_stability_certificate": out["implementation_stability_certificate"],
        "implementation_manifest_pass": out["implementation_manifest_pass"],
        "startup_certificate_pass": out["startup_certificate_pass"],
        "source_complete_word_language_pass": out["source_complete_word_language_pass"],
        "downstream_deployment_theorem_pass": out["downstream_deployment_theorem_pass"],
        "failures": out["failures"],
    }, indent=2, sort_keys=True))
    return 0 if out["implementation_stability_certificate"] == "PASS_IMPLEMENTATION_STABLE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
