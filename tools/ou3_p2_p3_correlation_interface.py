#!/usr/bin/env python3
"""Frozen P2 -> P3 theorem interface for OU-III.

P2 is deliberately split into two layers:

* source/timing admissibility from the existing EMA/staging/clock model; and
* the frozen correlated stage-transfer interface exported by
  :mod:`ou3_p2_correlation_path_memory`.

The physical 800-cell partition remains useful, but endpoint membership alone is
not sufficient for P3.  Canonical P3 must consume the versioned pair-state
correlation interface and must propagate process, covariance and measurement
bounds from the same admissible source history.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p2_correlation_path_memory as CORR
import ou3_p2_clock_phase_tuner_graph as PHASED
import ou3_p4_sample_clock_source_refinement as CLOCK

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_P2_TO_P3_CORRELATION_INTERFACE"


def build(domain_path: Path = DEFAULT_DOMAIN, candidate: dict | None = None) -> dict:
    # candidate is retained only for API compatibility with the first draft of
    # this interface.  Promotion is now tied to the repository-owned, versioned
    # CORR certificate rather than caller-supplied metadata.
    if candidate is not None:
        raise ValueError("external P2 correlation metadata is not accepted; use the versioned repository certificate")

    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P2->P3 interface must not be trajectory fitted")

    clock = CLOCK.build(path)
    cf = CLOCK.validate(clock)
    if cf:
        raise RuntimeError(f"P2 finite-speed clock refinement failed: {cf}")

    phased = PHASED.build(path)
    pf = PHASED.validate(phased)
    if pf:
        raise RuntimeError(f"P2 staged/committed source graph failed: {pf}")

    corr = CORR.build(path)
    rf = CORR.validate(corr)
    if rf:
        raise RuntimeError(f"P2 correlation path-memory certificate failed: {rf}")

    timing_pass = (
        clock.get("P2_SAMPLE_CLOCK_REFINEMENT_CERTIFICATE") == "PASS"
        and phased.get("P2_PHASED_SOURCE_GRAPH_CERTIFICATE") == "PASS"
    )
    correlation_ready = corr.get("P2_CORRELATION_INTERFACE_CERTIFICATE") == "PASS"

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generated_not_trajectory_fit": True,
        "filter_changed": False,
        "declared_domain_changed": False,
        "P2_source_timing_mathematics_retained": True,
        "P2_source_timing_certificate_pass": timing_pass,
        "physical_800_state_partition_retained": int(corr["physical_source_states"]) == 800,
        "physical_source_states": int(corr["physical_source_states"]),
        "stage_boundary_pair_states": int(corr["stage_boundary_pair_states"]),
        "clock_gap_samples": list(corr["clock_gap_alphabet_samples"]),
        "endpoint_only_800_state_quotient_sufficient_for_P3": False,
        "long_horizon_endpoint_ancestry_may_be_complete": True,
        "independent_cartesian_tau_sigma_RS_extrema_forbidden": True,
        "P3_must_use_one_common_source_history_for_all_bounds": True,
        "correlation_interface_version": corr["interface_version"],
        "correlation_interface_qualification": corr["qualification"],
        "correlation_pair_shift_rule": corr["transition_rule"],
        "correlation_segment_statistics": corr["segment_sufficient_statistics"],
        "correlation_consumer_contract": corr["consumer_contract"],
        "P2_CORRELATION_INTERFACE_READY": correlation_ready,
        "P2_READY_FOR_CANONICAL_P3": timing_pass and correlation_ready,
        "P3_must_declare_correlation_interface_consumed": True,
        "P3_required_correlation_interface_version": CORR.INTERFACE_VERSION,
        "P3_PROMOTED_HERE": False,
        "P4_PROMOTED_HERE": False,
        "next_obligation": (
            "make the canonical P3 covariance/information propagation consume this exact versioned pair-state source-history interface; process lower, covariance upper and measurement R_S bounds must be propagated from one common legal history before any scalar extremization"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != QUALIFICATION:
        f.append("wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit",
        "P2_source_timing_mathematics_retained",
        "P2_source_timing_certificate_pass",
        "physical_800_state_partition_retained",
        "independent_cartesian_tau_sigma_RS_extrema_forbidden",
        "P3_must_use_one_common_source_history_for_all_bounds",
        "P2_CORRELATION_INTERFACE_READY",
        "P2_READY_FOR_CANONICAL_P3",
        "P3_must_declare_correlation_interface_consumed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "filter_changed",
        "declared_domain_changed",
        "endpoint_only_800_state_quotient_sufficient_for_P3",
        "P3_PROMOTED_HERE",
        "P4_PROMOTED_HERE",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if int(d.get("physical_source_states", 0)) != 800:
        f.append("physical P2 source partition changed")
    if d.get("clock_gap_samples") != list(range(13, 27)):
        f.append("finite P2 clock alphabet changed")
    if d.get("correlation_interface_version") != CORR.INTERFACE_VERSION:
        f.append("P2 correlation interface version mismatch")
    if d.get("P3_required_correlation_interface_version") != CORR.INTERFACE_VERSION:
        f.append("canonical P3 correlation version requirement drifted")
    c = d.get("correlation_consumer_contract", {})
    if c.get("correlated_quantities_must_come_from_same_segment_node") is not True:
        f.append("same-segment source correlation requirement lost")
    if c.get("independent_tau_sigma_R_S_extremization_before_propagation") != "FORBIDDEN":
        f.append("independent source extrema became allowed")
    if c.get("global_800_ancestor_hull_as_P3_covariance_information_input") != "FORBIDDEN":
        f.append("flat 800-node ancestor hull became allowed for P3")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "timing_pass": d["P2_source_timing_certificate_pass"],
        "correlation_ready": d["P2_CORRELATION_INTERFACE_READY"],
        "correlation_version": d["correlation_interface_version"],
        "ready_for_canonical_P3": d["P2_READY_FOR_CANONICAL_P3"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
