#!/usr/bin/env python3
"""Bind the shared H18 nonlinear word screen to one exact P2 source node.

The broadband H18 screen is useful for debugging the nonlinear return map but
cannot represent the eventual path-dependent theorem metric, whose endpoints
are individual P2 source states.  This wrapper selects one of the exact 800
P2 tuner cells from :mod:`ou3_p4_source_node_cells`, injects that cell into the
same shared-operation H18 screen, and restores all process-global bindings on
exit.

This is intentionally still a screen: the P3 computational congruence remains
the whitening diagnostic, source vector orientation cells and optional branch
families are incomplete, and only one source node is evaluated per invocation.
No P4 promotion is possible here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p4_h18_shared_word_screen as SHARED
import ou3_p4_source_node_cells as NODES
from ou3_proof_module_state import preserve_module_bindings

DEFAULT_DOMAIN = SHARED.DEFAULT_DOMAIN
SCHEMA = 1


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_node_index: int = 0,
          samples: int | None = None, cell_limit: int = 1,
          ball_inflation: float = 1.5) -> dict:
    """Run the common H18 word map on one exact P2 source state cell."""
    node_payload = NODES.build()
    nf = NODES.validate(node_payload)
    if nf:
        raise RuntimeError(f"P2 source-node materialization failed: {nf}")
    node = NODES.node(source_node_index, node_payload)
    src = NODES.h18_source_cell(source_node_index, node_payload)

    with preserve_module_bindings():
        # SCREEN._run_cell calls SCREEN.H._source_cell exactly once before
        # constructing transition/covariance objects.  SHARED.build then swaps
        # only the nonlinear operation functions, so both wrappers compose.
        SHARED.SCREEN.H._source_cell = lambda: src
        out = dict(SHARED.build(
            Path(domain_path).resolve(),
            samples=samples,
            cell_limit=cell_limit,
            ball_inflation=ball_inflation,
        ))

    out["schema"] = SCHEMA
    out["qualification"] = "OU3_P4_H18_EXACT_P2_SOURCE_NODE_WORD_SCREEN"
    out["exact_P2_source_node_cell_used"] = True
    out["P2_source_node_count_available"] = node_payload["partition"]["states"]
    out["P2_source_node_index"] = int(source_node_index)
    out["P2_source_node"] = node
    out["all_P2_source_nodes_checked"] = False
    out["actual_per_node_Sigma_KF_whitening_used"] = False
    out["P4_USABLE_CERTIFICATE_PROMOTED"] = False
    out["next_obligation"] = (
        "attach actual source-node covariance/information factors and run the generalized metric on all reachable g->h edges; extend this exact-node screen across all 800 P2 source cells, vector-orientation cells and optional accepted/rejected/not-due branches"
    )
    return out


def validate(d: dict) -> list[str]:
    """Validate exact-node binding while retaining the screen's fail-closed scope."""
    base = dict(d)
    base["qualification"] = "OU3_P4_H18_SHARED_DIFFERENTIAL_WORD_SCREEN"
    failures = SHARED.validate(base)
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_H18_EXACT_P2_SOURCE_NODE_WORD_SCREEN":
        failures.append("wrong source-node H18 qualification")
    if d.get("exact_P2_source_node_cell_used") is not True:
        failures.append("exact P2 source-node cell was not used")
    if d.get("P2_source_node_count_available") != 800:
        failures.append("P2 source-node count is not 800")
    i = d.get("P2_source_node_index")
    if not isinstance(i, int) or not 0 <= i < 800:
        failures.append("invalid P2 source-node index")
    node = d.get("P2_source_node", {})
    if node.get("index") != i:
        failures.append("reported source node does not match selected index")
    if d.get("all_P2_source_nodes_checked") is not False:
        failures.append("single-node screen claimed all P2 nodes")
    if d.get("actual_per_node_Sigma_KF_whitening_used") is not False:
        failures.append("source-cell selection was confused with endpoint metric whitening")
    if d.get("P4_USABLE_CERTIFICATE_PROMOTED") is not False:
        failures.append("source-node screen improperly promoted P4")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-node-index", type=int, default=0)
    ap.add_argument("--samples", type=int, default=None)
    ap.add_argument("--cell-limit", type=int, default=1)
    ap.add_argument("--ball-inflation", type=float, default=1.5)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(
        args.domain.resolve(),
        source_node_index=args.source_node_index,
        samples=args.samples,
        cell_limit=args.cell_limit,
        ball_inflation=args.ball_inflation,
    )
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "validation_pass": not vf,
        "source_node": out["P2_source_node_index"],
        "source_tau_s": out["P2_source_node"]["tau_s"],
        "source_sigma_filter": out["P2_source_node"]["sigma_filter_committed_mps2"],
        "source_R_S": out["P2_source_node"]["R_S_filter_std"],
        "samples": [out["samples_checked"], out["full_word_samples"]],
        "cells": [out["outer_ball_box_cells_completed"], out["outer_ball_box_cells_requested"], out["outer_ball_box_cover_total"]],
        "max_endpoint_conditioned_norm": out["max_endpoint_P3_congruence_conditioned_norm_upper"],
        "max_prefix_conditioned_norm": out["max_prefix_P3_congruence_conditioned_norm_upper"],
        "first_failure": out["first_failure"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
