#!/usr/bin/env python3
"""Frozen P2 correlation interface plus the startup-to-Live period handoff.

The stage-transfer mathematics lives in :mod:`ou3_p2_stage_transfer_core` and is
kept byte-for-byte from the previously frozen P2 interface.  This front module
adds only the source-language obligation introduced when WavePeriodEstimator
became startup-usable before Live: the fixed prior may influence the initial
committed/staged pair, but the fixed-prior selector branch is unreachable in
Normal Live after the one-way usable-period latch.

Keeping that handoff outside the transfer core makes the proof change explicit:
no 800-cell partition, 13..26 clock alphabet, pair-shift rule, EMA image, or P3
correlation statistic is changed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p2_stage_transfer_core as CORE
import ou3_p2_live_entry_period_handoff as HANDOFF

REPO = CORE.REPO
DEFAULT_DOMAIN = CORE.DEFAULT_DOMAIN
SCHEMA = CORE.SCHEMA
INTERFACE_VERSION = CORE.INTERFACE_VERSION

down = CORE.down
up = CORE.up
runtime = CORE.runtime
legal_pair = CORE.legal_pair
successors = CORE.successors
segment_kernel = CORE.segment_kernel
transition = CORE.transition


def __getattr__(name):
    return getattr(CORE, name)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    d = CORE.build(path)

    handoff = HANDOFF.build(path)
    hf = HANDOFF.validate(handoff)
    failures = list(d.get("failures", []))
    failures.extend(f"live-entry period handoff: {x}" for x in hf)

    contract = dict(d.get("consumer_contract", {}))
    contract.update({
        "normal_live_fixed_prior_selector_branch": "FORBIDDEN",
        "initial_live_pair_may_retain_pre_live_prior_influence": True,
        "numeric_0p2_hz_value_does_not_identify_prior_provenance": True,
        "future_normal_live_frequency_source":
            "WavePeriodEstimator startup-usable state within the existing source clamp",
    })

    d.update({
        "wave_period_live_entry_handoff_consumed": True,
        "wave_period_live_entry_handoff_certificate_pass": not hf,
        "live_entry_frequency_handoff": handoff,
        "normal_live_fixed_prior_selector_branch_reachable": handoff.get(
            "normal_live_fixed_prior_selector_branch_reachable"
        ),
        "pre_live_prior_influence_retained_in_initial_live_pair": handoff.get(
            "pre_live_prior_influence_retained_in_initial_live_pair"
        ),
        "initial_live_pair_pruned_by_period_handoff": handoff.get(
            "initial_live_committed_staged_pair_pruned_by_this_certificate"
        ),
        "normal_live_future_frequency_source": handoff.get(
            "normal_live_future_frequency_source"
        ),
        "consumer_contract": contract,
        "failures": list(dict.fromkeys(failures)),
    })
    d["P2_CORRELATION_INTERFACE_CERTIFICATE"] = (
        "PASS" if not d["failures"] else "FAIL"
    )
    d["next_obligation"] = (
        "canonical P3 must consume this same pair-state/path-memory interface "
        "together with the Live-entry period handoff: retain possible pre-Live "
        "prior influence in the initial pair, but do not reintroduce a fixed-prior "
        "frequency selector branch inside Normal-Live words"
    )
    return d


def validate(d: dict) -> list[str]:
    f = list(CORE.validate(d))
    for key in (
        "wave_period_live_entry_handoff_consumed",
        "wave_period_live_entry_handoff_certificate_pass",
        "pre_live_prior_influence_retained_in_initial_live_pair",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "normal_live_fixed_prior_selector_branch_reachable",
        "initial_live_pair_pruned_by_period_handoff",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")

    handoff = d.get("live_entry_frequency_handoff", {})
    if handoff.get("P2_LIVE_ENTRY_WAVE_PERIOD_HANDOFF_CERTIFICATE") != "PASS":
        f.append("embedded live-entry period handoff did not pass")
    hf = HANDOFF.validate(handoff) if isinstance(handoff, dict) else [
        "missing live-entry handoff"
    ]
    f.extend(f"live-entry period handoff: {x}" for x in hf)

    c = d.get("consumer_contract", {})
    if c.get("normal_live_fixed_prior_selector_branch") != "FORBIDDEN":
        f.append("consumer contract reintroduced fixed prior in Normal Live")
    if c.get("initial_live_pair_may_retain_pre_live_prior_influence") is not True:
        f.append("consumer contract dropped pre-Live prior influence at handoff")
    if c.get("numeric_0p2_hz_value_does_not_identify_prior_provenance") is not True:
        f.append("consumer contract confuses numeric 0.2 Hz with prior provenance")
    if d.get("P2_CORRELATION_INTERFACE_CERTIFICATE") != "PASS":
        f.append("P2 correlation interface did not pass")
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
        "interface_version": d["interface_version"],
        "status": d["P2_CORRELATION_INTERFACE_CERTIFICATE"],
        "live_entry_handoff_pass": d[
            "wave_period_live_entry_handoff_certificate_pass"
        ],
        "prior_selector_reachable_in_normal_live": d[
            "normal_live_fixed_prior_selector_branch_reachable"
        ],
        "prelive_prior_influence_retained": d[
            "pre_live_prior_influence_retained_in_initial_live_pair"
        ],
        "physical_source_states": d["physical_source_states"],
        "stage_boundary_pair_states": d["stage_boundary_pair_states"],
        "gap_labelled_first_order_edges": d["gap_labelled_first_order_edges"],
        "gap_labelled_pair_state_edges_factorized_count":
            d["gap_labelled_pair_state_edges_factorized_count"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
