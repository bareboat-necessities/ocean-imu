#!/usr/bin/env python3
"""Frozen P2 -> P3 correlation interface for the OU-III theorem chain.

P2 has two distinct responsibilities:

1. source/timing admissibility -- shipping EMA, staging/commit clock, physical
   tuner ranges and finite-stage timing; and
2. correlation/path memory -- enough joint history that P3 cannot combine a
   process lower bound, covariance upper bound and measurement R_S bound from
   mutually incompatible source histories.

The existing 800-cell physical quotient is retained as a useful source-state
partition, but endpoint membership alone is explicitly insufficient for P3.
This module freezes the contract that any future P2 correlation certificate
must satisfy before canonical P3 may consume it.

No numerical P3 margin is established here.  In particular, metadata or a
source graph alone cannot make this interface ready.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p2_clock_phase_tuner_graph as PHASED
import ou3_p4_sample_clock_source_refinement as CLOCK
import ou3_p4_source_node_cells as NODES

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P2_TO_P3_CORRELATION_INTERFACE"
CANDIDATE_QUALIFICATION = "OU3_P2_SOURCE_HISTORY_CORRELATION_CERTIFICATE"

_REQUIRED_CORRELATED_QUANTITIES = (
    "tau_applied",
    "sigma_aw_applied",
    "R_S_applied",
    "process_excitation_lower",
    "covariance_upper",
    "measurement_R_S_lower",
)


def _candidate_status(candidate: dict | None, *, physical_states: int,
                      staged_pair_states: int) -> dict:
    if candidate is None:
        return {
            "provided": False,
            "contract_accepted": False,
            "reasons": ["no P2 source-history correlation certificate supplied"],
        }

    reasons: list[str] = []
    if candidate.get("qualification") != CANDIDATE_QUALIFICATION:
        reasons.append("wrong correlation-certificate qualification")
    for key in (
        "source_only",
        "exact_shipping_EMA_semantics",
        "exact_staging_commit_semantics",
        "clock_phase_retained",
        "staged_tuple_retained",
        "committed_tuple_retained",
        "EMA_candidate_memory_retained",
        "all_admissible_source_paths_covered",
        "same_history_used_for_process_covariance_and_measurement_bounds",
        "joint_tau_sigma_RS_correlation_retained",
        "source_conditioned_covariance_upper_available",
        "source_conditioned_process_lower_available",
        "source_conditioned_measurement_R_S_lower_available",
        "full_P3_word_horizon_covered",
        "frozen_clock_branch_covered",
    ):
        if candidate.get(key) is not True:
            reasons.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used",
        "filter_changed",
        "declared_domain_changed",
        "independent_cartesian_extrema_combination_allowed",
        "endpoint_only_800_state_quotient_used_as_complete_history",
    ):
        if candidate.get(key) is not False:
            reasons.append(f"{key} is not false")

    if int(candidate.get("physical_source_states", -1)) != int(physical_states):
        reasons.append("physical source-state count mismatch")
    if int(candidate.get("stage_boundary_pair_states", -1)) != int(staged_pair_states):
        reasons.append("staged/committed pair-state count mismatch")

    quantities = candidate.get("correlated_quantities", [])
    if not isinstance(quantities, list) or any(q not in quantities for q in _REQUIRED_CORRELATED_QUANTITIES):
        reasons.append("required correlated quantity set is incomplete")

    horizon = candidate.get("certified_history_horizon_s")
    if isinstance(horizon, bool) or not isinstance(horizon, (int, float)) or not math.isfinite(float(horizon)) or float(horizon) <= 0.0:
        reasons.append("missing positive certified source-history horizon")

    return {
        "provided": True,
        "contract_accepted": not reasons,
        "reasons": reasons,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, candidate: dict | None = None) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P2->P3 interface must not be trajectory fitted")

    nodes = NODES.build()
    nf = NODES.validate(nodes)
    if nf:
        raise RuntimeError(f"P2 physical source partition failed: {nf}")

    clock = CLOCK.build(path)
    cf = CLOCK.validate(clock)
    if cf:
        raise RuntimeError(f"P2 finite-speed clock refinement failed: {cf}")

    phased = PHASED.build(path)
    pf = PHASED.validate(phased)
    if pf:
        raise RuntimeError(f"P2 staged/committed source graph failed: {pf}")

    physical_states = int(nodes["partition"]["states"])
    pair_states = int(phased["stage_boundary_pair_states"])
    status = _candidate_status(
        candidate,
        physical_states=physical_states,
        staged_pair_states=pair_states,
    )

    timing_pass = (
        clock.get("P2_SAMPLE_CLOCK_REFINEMENT_CERTIFICATE") == "PASS"
        and phased.get("P2_PHASED_SOURCE_GRAPH_CERTIFICATE") == "PASS"
    )
    correlation_ready = bool(status["contract_accepted"])

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generated_not_trajectory_fit": True,
        "filter_changed": False,
        "declared_domain_changed": False,
        "P2_source_timing_mathematics_retained": True,
        "P2_source_timing_certificate_pass": timing_pass,
        "physical_800_state_partition_retained": physical_states == 800,
        "physical_source_states": physical_states,
        "stage_boundary_pair_states": pair_states,
        "clock_gap_samples": phased["clock_phase_gap_alphabet_samples"],
        "staged_committed_memory_available": True,
        "EMA_candidate_memory_available": bool(phased["EMA_candidate_memory_retained_across_commits"]),
        "endpoint_only_800_state_quotient_sufficient_for_P3": False,
        "long_horizon_endpoint_ancestry_may_be_complete": True,
        "independent_cartesian_tau_sigma_RS_extrema_forbidden": True,
        "P3_must_use_one_common_source_history_for_all_bounds": True,
        "required_correlated_quantities": list(_REQUIRED_CORRELATED_QUANTITIES),
        "correlation_candidate": status,
        "P2_CORRELATION_INTERFACE_READY": correlation_ready,
        "P2_READY_FOR_CANONICAL_P3": timing_pass and correlation_ready,
        "P3_PROMOTED_HERE": False,
        "P4_PROMOTED_HERE": False,
        "next_obligation": (
            "construct a source-only history-conditioned enclosure over the retained staged/committed/clock automaton that supplies process excitation, covariance and measurement R_S bounds from the same admissible history; do not replace it by endpoint-only ancestry or independent Cartesian extrema"
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
        "staged_committed_memory_available",
        "EMA_candidate_memory_available",
        "independent_cartesian_tau_sigma_RS_extrema_forbidden",
        "P3_must_use_one_common_source_history_for_all_bounds",
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
    candidate_ok = bool(d.get("correlation_candidate", {}).get("contract_accepted"))
    if d.get("P2_CORRELATION_INTERFACE_READY") is not candidate_ok:
        f.append("P2 correlation-ready flag does not match candidate contract")
    expected_ready = bool(d.get("P2_source_timing_certificate_pass")) and candidate_ok
    if d.get("P2_READY_FOR_CANONICAL_P3") is not expected_ready:
        f.append("P2 canonical-P3 readiness flag is inconsistent")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--candidate", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    candidate = None
    if args.candidate is not None:
        candidate = json.loads(args.candidate.read_text(encoding="utf-8"))
    d = build(args.domain, candidate)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "timing_pass": d["P2_source_timing_certificate_pass"],
        "correlation_ready": d["P2_CORRELATION_INTERFACE_READY"],
        "ready_for_canonical_P3": d["P2_READY_FOR_CANONICAL_P3"],
        "candidate_reasons": d["correlation_candidate"]["reasons"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
