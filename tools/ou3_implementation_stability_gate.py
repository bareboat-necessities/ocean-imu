#!/usr/bin/env python3
"""Final independent composition gate for the deployed OU-III stability proof.

This sits above the existing deployment-theorem gate.  It adds source and
implementation parity, the actual pre-Live reset/Mahony/goLive certificate, the
recurring-PE source-word language, and the explicit startup-to-inner-H capture
composition required by P5.

The older ``ou3_p5_startup_capture_certificate`` remains an obstruction
identifier: it proves that the tiny local P4 seed cannot be extrapolated to the
P1 handoff family.  It is not the current P5 completion object.  The final gate
also regenerates ``ou3_p5_outer_h_bridge_certificate`` and requires that staged
bridge to close before a finite H-word count can contribute to implementation
stability.  This prevents a stale obstruction diagnostic or the generic
deployment capture arithmetic from standing in for the actual source-faithful
outer bridge.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_deployment_gate as DEPLOY
import ou3_implementation_proof_manifest as MANIFEST
import ou3_implementation_word_language as WORDS
import ou3_p5_outer_h_bridge_certificate as P5BRIDGE
import ou3_p5_startup_capture_certificate as P5ID
import ou3_startup_stability_certificate as STARTUP

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 3


def compose(check: dict, source_domain: dict, primitive: dict,
            domain_path: Path = DEFAULT_DOMAIN) -> dict:
    manifest = MANIFEST.build()
    manifest_failures = MANIFEST.validate(manifest)
    startup = STARTUP.build(domain_path)
    startup_failures = STARTUP.validate(startup)
    words = WORDS.build(domain_path)
    word_failures = WORDS.validate(words)

    p5_id = P5ID.build(domain_path)
    p5_id_failures = P5ID.validate(p5_id)
    p5_bridge = P5BRIDGE.build(domain_path)
    p5_bridge_failures = P5BRIDGE.validate(p5_bridge)
    p5_finite_capture_pass = bool(
        not p5_id_failures
        and not p5_bridge_failures
        and p5_bridge.get("P5_OUTER_H_BRIDGE_CERTIFICATE") == "PASS"
        and isinstance(p5_bridge.get("N_H_words"), int)
        and p5_bridge["N_H_words"] >= 0
    )

    deployment = DEPLOY.compose(check, source_domain, primitive)
    downstream_pass = deployment.get("deployment_theorem_certificate") == "PASS"
    final_pass = bool(
        not manifest_failures
        and not startup_failures
        and not word_failures
        and p5_finite_capture_pass
        and downstream_pass
    )

    failures: list[str] = []
    failures.extend(f"implementation manifest: {x}" for x in manifest_failures)
    failures.extend(f"startup: {x}" for x in startup_failures)
    failures.extend(f"source-word language: {x}" for x in word_failures)
    failures.extend(f"P5 identification validation: {x}" for x in p5_id_failures)
    failures.extend(f"P5 staged bridge validation: {x}" for x in p5_bridge_failures)
    if not p5_finite_capture_pass:
        obstruction = p5_bridge.get(
            "first_failure",
            p5_id.get("first_obstruction", "UNKNOWN_P5_OBSTRUCTION"),
        )
        failures.append(
            "P5 finite startup-to-inner-funnel capture not established: " + str(obstruction)
        )
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
        "P5_identification": p5_id,
        "P5_identification_validation_failures": p5_id_failures,
        "P5_bridge": p5_bridge,
        "P5_bridge_validation_failures": p5_bridge_failures,
        "P5_finite_startup_capture_pass": p5_finite_capture_pass,
        "P5_first_obstruction": p5_bridge.get(
            "first_failure", p5_id.get("first_obstruction")
        ),
        "P5_local_capture_identifier_is_not_completion_certificate": True,
        "deployment": deployment,
        "downstream_deployment_theorem_pass": downstream_pass,
        "generic_deployment_capture_is_not_P5": True,
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
        "P5_finite_startup_capture_pass": out["P5_finite_startup_capture_pass"],
        "P5_first_obstruction": out["P5_first_obstruction"],
        "downstream_deployment_theorem_pass": out["downstream_deployment_theorem_pass"],
        "failures": out["failures"],
    }, indent=2, sort_keys=True))
    return 0 if out["implementation_stability_certificate"] == "PASS_IMPLEMENTATION_STABLE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
