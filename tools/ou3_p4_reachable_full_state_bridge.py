#!/usr/bin/env python3
"""Bridge validated translation/full-word results into the next full-state P4 test.

This producer does not promote P4.  It combines three already source-only
objects on the same recurrent worst-cell route:

* the outward-validated one-second complete-word translation margin;
* the direct nontranslation generalized margin diagnostic; and
* the source-dynamic reachability graph.

For a normalized full-state endpoint comparison split into translation and
nontranslation blocks,

    G = [[A, C], [C^T, B]],

with A >= a I and B >= b I, a sufficient Schur/Young condition for G >> 0 is

    ||C||_2 < sqrt(a b).

The returned sqrt(a b) is therefore a concrete rigorous *cross-block budget*
for the next 18/21-state Phi/Omega enclosure.  It replaces the vague
"propagate the full state" blocker with a numerical target.  The actual C
bound is deliberately not guessed here and P4 remains NOT_ESTABLISHED until a
validated complete-word full-state backend proves it on every required
reachable path/cell.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def build(translation: dict, bottleneck: dict, path: dict) -> dict:
    failures: list[str] = []
    if translation.get("P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS") != "PASS":
        failures.append("validated complete-word translation input is not PASS")
    if path.get("path_graph_ready") is not True:
        failures.append("source-path graph is not ready")
    if int(path.get("recurrent_states", 0)) <= 0:
        failures.append("source-path graph has no recurrent states")

    modes: dict[str, dict] = {}
    for mode in ("H", "A"):
        tr = translation.get("modes", {}).get(mode, {})
        bn = bottleneck.get("modes", {}).get(mode, {})
        a = tr.get("complete_word_translation_margin_lower")
        b = bn.get("existing_direct_nontranslation_margin_lower")
        try:
            a = float(a)
            b = float(b)
        except (TypeError, ValueError):
            a = b = float("nan")
        if not (math.isfinite(a) and math.isfinite(b) and a > 0.0 and b > 0.0):
            failures.append(f"{mode}: missing positive block margins")
            continue
        budget = math.sqrt(a * b)
        modes[mode] = {
            "validated_translation_margin_lower": a,
            "direct_nontranslation_margin_lower": b,
            "normalized_cross_block_spectral_norm_budget_upper_open": budget,
            "full_state_sufficient_condition": "validated ||C||_2 < sqrt(delta_translation*delta_nontranslation)",
            "cross_block_bound_validated": False,
            "full_state_linear_certificate_established": False,
        }

    return {
        "qualification": "OU3_P4_REACHABLE_FULL_STATE_CROSS_BLOCK_BRIDGE",
        "source_only": True,
        "trajectory_replay_used": False,
        "reachable_state_count": int(path.get("partition", {}).get("states", 0)),
        "recurrent_state_count": int(path.get("recurrent_states", 0)),
        "source_graph_strongly_connected_components": int(path.get("strongly_connected_components", 0)),
        "old_worst_corner_recurrent": bool(path.get("old_worst_corner_has_internal_recurrent_cycle", False)),
        "modes": modes,
        "P4_USABLE_CERTIFICATE_STATUS": "NOT_ESTABLISHED",
        "next_obligation": (
            "propagate the complete 18/21-state Phi/Omega endpoint on the reachable source graph and "
            "outward-validate the normalized translation/nontranslation cross-block spectral norm C; "
            "prove it is below the emitted sqrt(delta_translation*delta_nontranslation) budget on every required cell/path"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("source_only") is not True or d.get("trajectory_replay_used") is not False:
        f.append("bridge is not source-only")
    if d.get("P4_USABLE_CERTIFICATE_STATUS") != "NOT_ESTABLISHED":
        f.append("bridge prematurely promoted P4")
    if int(d.get("recurrent_state_count", 0)) <= 0:
        f.append("bridge lost recurrent source graph")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        q = m.get("normalized_cross_block_spectral_norm_budget_upper_open")
        if not isinstance(q, (int, float)) or not math.isfinite(float(q)) or float(q) <= 0.0:
            f.append(f"{mode}: invalid cross-block budget")
        if m.get("cross_block_bound_validated") is not False:
            f.append(f"{mode}: bridge falsely claims cross-block validation")
        if m.get("full_state_linear_certificate_established") is not False:
            f.append(f"{mode}: bridge falsely claims full-state certificate")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--translation", type=Path, required=True)
    ap.add_argument("--bottleneck", type=Path, required=True)
    ap.add_argument("--path", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(
        json.loads(a.translation.read_text(encoding="utf-8")),
        json.loads(a.bottleneck.read_text(encoding="utf-8")),
        json.loads(a.path.read_text(encoding="utf-8")),
    )
    f = validate(d)
    d["validation_failures"] = f
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P4_USABLE_CERTIFICATE_STATUS"],
        "reachable_states": d["reachable_state_count"],
        "recurrent_states": d["recurrent_state_count"],
        "old_worst_corner_recurrent": d["old_worst_corner_recurrent"],
        "modes": {
            mode: {
                "translation": d.get("modes", {}).get(mode, {}).get("validated_translation_margin_lower"),
                "nontranslation": d.get("modes", {}).get(mode, {}).get("direct_nontranslation_margin_lower"),
                "cross_block_budget": d.get("modes", {}).get(mode, {}).get("normalized_cross_block_spectral_norm_budget_upper_open"),
            }
            for mode in ("H", "A")
        },
        "failures": f,
    }, indent=2, sort_keys=True))
    return 0 if not f else 2


if __name__ == "__main__":
    raise SystemExit(main())
