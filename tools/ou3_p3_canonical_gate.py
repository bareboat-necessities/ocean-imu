#!/usr/bin/env python3
"""Single canonical P3 theorem-interface gate for OU-III.

Canonical P3 PASS means exactly:

* the refined, versioned P2 source-history correlation interface is valid;
* the P3 producer explicitly consumes that exact interface version;
* process lower, covariance upper and measurement bounds come from one common
  admissible source history rather than independent Cartesian extrema;
* tuner/source variation over the word is covered;
* interleaved accelerometer and S=0 measurements are covered;
* both H and A modes pass; and
* worst H/A relative Riccati injection margin is >= 1e-18.

Those are theorem obligations, not implementation choices.  A producer may
satisfy the time-varying-source and measurement obligations with the retained
LTV determinant/lift, or with a later equally rigorous source-complete full-
matrix construction.  The gate must not force one proof mechanism after the
semantic theorem interface has been frozen.

No prerequisite or diagnostic may substitute for this verdict.  P4 must consume
P3_CANONICAL_PASS, not producer-local PASS-like flags.
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
SCHEMA = 3
QUALIFICATION = "OU3_P3_CANONICAL_THEOREM_INTERFACE"


def _semantic_coverage(cand: dict, explicit: str, retained_mechanism: str) -> tuple[bool, str | None]:
    """Resolve a theorem obligation without baking one proof device into P3.

    Existing retained producers advertise the historical mechanism flag.  New
    producers may instead advertise the semantic obligation directly.  Either
    path must say True explicitly; absence is fail-closed.
    """
    if cand.get(explicit) is True:
        return True, "semantic"
    if cand.get(retained_mechanism) is True:
        return True, "retained_mechanism"
    return False, None


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

    required_version = p2.get("P3_required_correlation_interface_version")
    if cand.get("P2_correlation_interface_consumed") is not True:
        reasons.append("P3 candidate has not consumed the refined P2 source-history interface")
    if cand.get("P2_correlation_interface_version") != required_version:
        reasons.append("P3 candidate is not bound to the required P2 correlation-interface version")
    if cand.get("process_covariance_measurement_bounds_same_source_history") is not True:
        reasons.append("P3 candidate does not propagate process/covariance/measurement bounds from one source history")

    independent = cand.get("independent_cartesian_tau_sigma_RS_extrema_used")
    if independent is None:
        independent = cand.get("independent_cartesian_tau_sigma_R_S_extrema_used")
    if independent is not False:
        reasons.append("P3 candidate still permits independent Cartesian tau/sigma/R_S extrema")

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

    source_covered, source_method = _semantic_coverage(
        cand, "time_varying_tuner_over_word_covered", "retained_LTV_determinant_consumed"
    )
    if not source_covered:
        reasons.append("P3 candidate does not cover time-varying tuner/source excitation over the word")

    measurements_covered, measurement_method = _semantic_coverage(
        cand, "interleaved_accelerometer_and_S_measurements_covered",
        "same_lifted_measurement_attenuation_as_retained_route",
    )
    if not measurements_covered:
        reasons.append("P3 candidate does not cover interleaved accelerometer and S=0 measurements")

    gate = float(BASE.MIN_USEFUL_DELTA)
    if gate != 1.0e-18:
        raise RuntimeError("canonical P3 useful gate changed")

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
        "proof_mechanism_not_part_of_canonical_definition": True,
        "only_this_module_may_promote_P3_for_P4": True,
        "useful_gate": gate,
        "required_properties": {
            "refined_P2_source_history_correlation": True,
            "exact_P2_correlation_interface_version_consumed": required_version,
            "same_source_history_for_process_covariance_measurement": True,
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
            "correlation_version": p2["correlation_interface_version"],
            "ready_for_canonical_P3": p2["P2_READY_FOR_CANONICAL_P3"],
        },
        "candidate_qualification": cand.get("qualification"),
        "candidate_declared_P2_correlation_version": cand.get("P2_correlation_interface_version"),
        "semantic_coverage": {
            "time_varying_tuner_over_word": source_covered,
            "time_varying_tuner_proof_advertisement": source_method,
            "interleaved_accelerometer_and_S_measurements": measurements_covered,
            "interleaved_measurement_proof_advertisement": measurement_method,
        },
        "mode_margins": mode_margins,
        "worst_H_A_margin": worst,
        "P3_CANONICAL_PASS": passed,
        "P3_CANONICAL_FAIL_REASONS": reasons,
        "P4_MAY_CONSUME_P3": passed,
        "P5_MAY_CONSUME_P4": False,
        "next_obligation": (
            "freeze P3 and feed this exact canonical metric/interface to P4"
            if passed
            else "improve only source-faithful numerical bounds under this frozen theorem interface until the unchanged 1e-18 gate passes; do not redefine P3"
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
    if d.get("proof_mechanism_not_part_of_canonical_definition") is not True:
        f.append("canonical P3 still depends on one proof mechanism")
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
    p2 = d.get("P2_interface", {})
    if p2.get("timing_pass") is not True or p2.get("correlation_ready") is not True:
        f.append("canonical gate did not receive valid P2 timing/correlation interfaces")
    if not isinstance(p2.get("correlation_version"), str) or not p2["correlation_version"]:
        f.append("canonical gate lost P2 correlation interface version")
    required = d.get("required_properties", {})
    for key in (
        "refined_P2_source_history_correlation",
        "same_source_history_for_process_covariance_measurement",
        "source_generated_not_replay",
        "time_varying_tuner_over_word",
        "interleaved_accelerometer_and_S_measurements",
        "H_mode", "A_mode", "worst_margin_at_least_1e-18",
    ):
        if required.get(key) is not True:
            f.append(f"canonical required property changed: {key}")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--p3-candidate", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    p3c = json.loads(args.p3_candidate.read_text(encoding="utf-8")) if args.p3_candidate else None
    d = build(args.domain, None, p3c)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "P2_correlation_version": d["P2_interface"]["correlation_version"],
        "P3_CANONICAL_PASS": d["P3_CANONICAL_PASS"],
        "worst_H_A_margin": d["worst_H_A_margin"],
        "semantic_coverage": d["semantic_coverage"],
        "fail_reasons": d["P3_CANONICAL_FAIL_REASONS"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
