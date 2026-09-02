#!/usr/bin/env python3
"""Single canonical P3 theorem-interface gate for OU-III.

No prerequisite, exact-node probe, premeasurement result, or source-uniform
Cartesian over-approximation is allowed to call itself canonical P3 PASS.
Canonical P3 requires all of the following simultaneously:

* a validated refined P2->P3 correlation certificate;
* source-generated, non-replay bounds;
* time-varying tuner/source coverage over the P3 word;
* interleaved accelerometer and S=0 measurement coverage;
* both H and A fixed-dimensional modes; and
* worst H/A relative Riccati injection margin >= 1e-18.

This module is the sole promotion gate.  Producers may emit useful diagnostics
or candidates, but P4 must consume P3_CANONICAL_PASS rather than their local
PASS-like flags.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p2_p3_correlation_interface as P2I
import ou3_p3_e3_postmeasurement_certificate as CANDIDATE
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P3_CANONICAL_THEOREM_INTERFACE"


def build(domain_path: Path = DEFAULT_DOMAIN,
          p2_correlation_candidate: dict | None = None,
          p3_candidate: dict | None = None) -> dict:
    path = Path(domain_path).resolve()
    p2 = P2I.build(path, p2_correlation_candidate)
    p2f = P2I.validate(p2)
    if p2f:
        raise RuntimeError(f"P2->P3 interface invalid: {p2f}")

    cand = CANDIDATE.build(path) if p3_candidate is None else p3_candidate

    reasons: list[str] = []
    if p2.get("P2_READY_FOR_CANONICAL_P3") is not True:
        reasons.append("refined P2 source-history correlation interface is not ready")
    if cand.get("source_generated_not_trajectory_fit") is not True:
        reasons.append("P3 candidate is not source generated")
    if cand.get("trajectory_replay_used") is not False:
        reasons.append("P3 candidate uses trajectory replay")
    if cand.get("filter_changed") is not False:
        reasons.append("P3 candidate changes the filter")
    if cand.get("zero_lever_arm_branch") is not True:
        reasons.append("P3 candidate is not bound to declared zero-lever-arm branch")
    if cand.get("dormant_transparent_vibration_guard_branch") is not True:
        reasons.append("P3 candidate is not bound to dormant vibration-guard branch")
    if cand.get("same_lifted_measurement_attenuation_as_retained_route") is not True:
        reasons.append("P3 candidate does not cover the retained interleaved measurement lift")

    # The current e3 producer obtains its process lower through the LTV source
    # route, so this is the explicit time-varying-source requirement.  Future
    # candidates must preserve an equivalent positive declaration.
    if cand.get("retained_LTV_determinant_consumed") is not True:
        reasons.append("P3 candidate does not cover time-varying tuner/source process excitation")

    gate = float(BASE.MIN_USEFUL_DELTA)
    mode_margins: dict[str, float | None] = {}
    for mode in ("H", "A"):
        row = cand.get("modes", {}).get(mode, {})
        x = row.get("relative_Riccati_injection_margin_lower")
        if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(float(x)) or float(x) <= 0.0:
            reasons.append(f"{mode}: missing strict post-measurement P3 margin")
            mode_margins[mode] = None
        else:
            mode_margins[mode] = float(x)
            if float(x) < gate:
                reasons.append(f"{mode}: margin below unchanged 1e-18 usefulness gate")

    worst = min((x for x in mode_margins.values() if x is not None), default=0.0)
    passed = not reasons

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_definition_frozen": True,
        "only_this_module_may_promote_P3_for_P4": True,
        "useful_gate": gate,
        "required_properties": {
            "refined_P2_source_history_correlation": True,
            "source_generated_not_replay": True,
            "time_varying_tuner_over_word": True,
            "interleaved_accelerometer_and_S_measurements": True,
            "H_mode": True,
            "A_mode": True,
            "worst_margin_at_least_1e-18": True,
        },
        "P2_interface": {
            "timing_pass": p2["P2_source_timing_certificate_pass"],
            "correlation_ready": p2["P2_CORRELATION_INTERFACE_READY"],
            "ready_for_canonical_P3": p2["P2_READY_FOR_CANONICAL_P3"],
        },
        "candidate_qualification": cand.get("qualification"),
        "mode_margins": mode_margins,
        "worst_H_A_margin": worst,
        "P3_CANONICAL_PASS": passed,
        "P3_CANONICAL_FAIL_REASONS": reasons,
        "P4_MAY_CONSUME_P3": passed,
        "P5_MAY_CONSUME_P4": False,
        "next_obligation": (
            "freeze P3 and feed this exact canonical metric/interface to P4"
            if passed
            else "repair only the listed P2->P3 correlation or fixed P3 numerical obligations; do not redefine the canonical P3 theorem interface"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != QUALIFICATION:
        f.append("wrong qualification")
    if d.get("canonical_definition_frozen") is not True:
        f.append("canonical P3 definition is not frozen")
    if d.get("only_this_module_may_promote_P3_for_P4") is not True:
        f.append("canonical promotion authority is not unique")
    if d.get("useful_gate") != 1.0e-18:
        f.append("canonical P3 useful gate changed")
    reasons = d.get("P3_CANONICAL_FAIL_REASONS", [])
    expected = isinstance(reasons, list) and len(reasons) == 0
    if d.get("P3_CANONICAL_PASS") is not expected:
        f.append("canonical P3 pass flag does not match obligations")
    if d.get("P4_MAY_CONSUME_P3") is not d.get("P3_CANONICAL_PASS"):
        f.append("P4 consumption gate differs from canonical P3 result")
    if d.get("P5_MAY_CONSUME_P4") is not False:
        f.append("P5 was enabled by the P3 gate")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--p2-correlation-candidate", type=Path)
    ap.add_argument("--p3-candidate", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    p2c = json.loads(args.p2_correlation_candidate.read_text(encoding="utf-8")) if args.p2_correlation_candidate else None
    p3c = json.loads(args.p3_candidate.read_text(encoding="utf-8")) if args.p3_candidate else None
    d = build(args.domain, p2c, p3c)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "P3_CANONICAL_PASS": d["P3_CANONICAL_PASS"],
        "worst_H_A_margin": d["worst_H_A_margin"],
        "fail_reasons": d["P3_CANONICAL_FAIL_REASONS"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
