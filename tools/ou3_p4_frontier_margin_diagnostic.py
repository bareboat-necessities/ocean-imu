#!/usr/bin/env python3
"""Diagnose why the current P4 frontier is unusably small.

This producer does not promote a new nonlinear certificate.  It instruments the
source-cell P3 generalized-matrix calculation, records every source-cell margin
that participates in the source-uniform endpoint theorem, and compares that
linear margin with the current exact-nonlinear P4 frontier.

The purpose is to distinguish two cases before doing more P4 algebra:

1. one or a small family of source cells destroys the uniform margin, in which
   case the next route is a source-correlated/path-dependent metric or source
   subdivision; or
2. every admissible cell genuinely has a machine-epsilon-scale margin, in
   which case a longer recurrent-PE superword/direct nonlinear return map is
   required.

No replay data are used and no usefulness claim is made from a tiny W level.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_source_reachable_matrix_p3_direct as P3D
import ou3_explicit_information_word_certificate as P3
import ou3_p4_frontier_combined_certificate as P4
import ou3_p4_direct_word_contraction_certificate as DIRECT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"


def _ival(x):
    if hasattr(x, "lo") and hasattr(x, "hi"):
        return [float(x.lo), float(x.hi)]
    return x


def _capture_cells(domain_path: Path):
    captured = []
    original = P3D.BASE.mode_cell

    def wrapped(mode, x, rho_trans, sigma, rs, live, vector, process, sched, alpha6):
        row = original(mode, x, rho_trans, sigma, rs, live, vector, process, sched, alpha6)
        captured.append({
            "mode": mode,
            "x_h_over_tau": _ival(x),
            "sigma_aw": _ival(sigma),
            "R_S": _ival(rs),
            "rho_translation_lower": float(rho_trans),
            "alpha6_lower": float(alpha6),
            "delta_full_lower": float(row["relative_Riccati_injection_margin_lower"]),
            "delta_translation_lower": float(row["direct_translation_generalized_margin_lower"]),
            "delta_nontranslation_lower": float(row["direct_nontranslation_margin_lower"]),
            "limiting_block": row["generalized_matrix_inequality"]["limiting_block"],
            "ldlt_downward_shrink_count": int(row["generalized_matrix_inequality"]["rounding_boundary_downward_shrink_count"]),
            "measurement_information_beta_upper": float(row["generalized_matrix_inequality"]["measurement_information_beta_upper"]),
            "translation_posterior_matrix_factor_lower": float(row["generalized_matrix_inequality"]["translation_posterior_matrix_factor_lower"]),
        })
        return row

    P3D.BASE.mode_cell = wrapped
    P3D.BASE._build_cached.cache_clear()
    try:
        p3 = P3.build(domain_path)
    finally:
        P3D.BASE.mode_cell = original
        P3D.BASE._build_cached.cache_clear()
    return p3, captured


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    p = Path(domain_path).resolve()
    p3, cells = _capture_cells(p)
    p4 = P4.build(p)
    failures = [f"P3: {x}" for x in P3.validate(p3)]
    failures += [f"P4: {x}" for x in P4.validate(p4)]

    modes = {}
    for mode in ("H", "A"):
        rows = [r for r in cells if r["mode"] == mode]
        if not rows:
            failures.append(f"{mode}: no P3 source cells captured")
            continue
        rows.sort(key=lambda r: r["delta_full_lower"])
        worst = rows[0]
        best = rows[-1]
        delta = float(p3["modes"][mode]["word_endpoint_relative_Riccati_injection_margin_lower"])
        gap = float(DIRECT.strict_gap(delta))
        fm = p4["modes"][mode]
        W = float(fm["certified_level_W"])
        sqrtW = math.sqrt(W)
        mode_cells = len(rows)
        near_worst = sum(r["delta_full_lower"] <= 10.0 * worst["delta_full_lower"] for r in rows)
        modes[mode] = {
            "p3_word_endpoint_delta_lower": delta,
            "p3_strict_sqrt_endpoint_gap_lower": gap,
            "p3_cell_count": mode_cells,
            "p3_cells_within_10x_worst_margin": near_worst,
            "p3_worst_cell": worst,
            "p3_best_cell": best,
            "p3_best_to_worst_delta_ratio": best["delta_full_lower"] / worst["delta_full_lower"],
            "p4_frontier_W": W,
            "p4_frontier_sqrtW": sqrtW,
            "p4_frontier_prefix_canonical_norm_upper": float(fm["prefix_canonical_error_norm_upper"]),
            "p4_active_cap": fm["frontier_selected_active_cap"],
            "current_scalar_small_gain_route_usable": False,
            "usability_reason": (
                "current P4 is endpoint-limited by a P3 margin at or near binary64 precision; "
                "do not promote this tiny level as a practical basin"
            ),
        }

    out = {
        "qualification": "OU3_P4_USABILITY_AND_P3_MARGIN_DIAGNOSTIC",
        "source_only": True,
        "trajectory_replay_used": False,
        "modes": modes,
        "P4_USABLE_CERTIFICATE_STATUS": "NOT_ESTABLISHED",
        "next_route": (
            "If the worst-cell margin is isolated, replace the source-uniform metric by source-correlated/path-dependent cells. "
            "If margins are uniformly tiny, certify a longer recurrent-PE superword and lift that return map directly with validated subdivision."
        ),
        "failures": failures,
    }
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("P4_USABLE_CERTIFICATE_STATUS") != "NOT_ESTABLISHED":
        f.append("tiny current P4 was incorrectly promoted as usable")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        if not m.get("p3_cell_count", 0) > 0:
            f.append(f"{mode}: missing source-cell diagnostic")
        if m.get("current_scalar_small_gain_route_usable") is not False:
            f.append(f"{mode}: scalar small-gain route incorrectly marked usable")
        if not float(m.get("p3_word_endpoint_delta_lower", 0.0)) > 0.0:
            f.append(f"{mode}: P3 margin not positive")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain)
    f = validate(d)
    d["validation_failures"] = f
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True))
    print(json.dumps({
        "status": d["P4_USABLE_CERTIFICATE_STATUS"],
        "modes": {
            mode: {
                "delta": d.get("modes", {}).get(mode, {}).get("p3_word_endpoint_delta_lower"),
                "gap": d.get("modes", {}).get(mode, {}).get("p3_strict_sqrt_endpoint_gap_lower"),
                "cells": d.get("modes", {}).get(mode, {}).get("p3_cell_count"),
                "near_worst": d.get("modes", {}).get(mode, {}).get("p3_cells_within_10x_worst_margin"),
                "best_worst_ratio": d.get("modes", {}).get(mode, {}).get("p3_best_to_worst_delta_ratio"),
                "worst": d.get("modes", {}).get(mode, {}).get("p3_worst_cell"),
                "W": d.get("modes", {}).get(mode, {}).get("p4_frontier_W"),
                "qprefix": d.get("modes", {}).get(mode, {}).get("p4_frontier_prefix_canonical_norm_upper"),
            }
            for mode in ("H", "A")
        },
        "failures": f,
    }, indent=2, sort_keys=True))
    return 0 if not f else 2


if __name__ == "__main__":
    raise SystemExit(main())
