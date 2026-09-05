#!/usr/bin/env python3
"""Stable source-history-free SEA3 finite-word translation lower.

This facade keeps the v2 theorem but refines instantaneous x cells until the
validated process interval admits a strict deterministic Loewner lower.  The
refinement is arithmetic only: every x-cell identity is discarded after each
sample and no predecessor/source history is created.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, symmetric_positive_definite_ldlt
import ou3_sea3_finite_word_translation_v2_backend as B

SCHEMA = B.SCHEMA
QUALIFICATION = B.QUALIFICATION
USEFUL_GATE = B.USEFUL_GATE
WORD_HORIZON_S = B.WORD_HORIZON_S
DEFAULT_DOMAIN = B.DEFAULT_DOMAIN
MAX_PROCESS_REFINEMENT_DEPTH = 16

_PROFILE = {}


def _refined_process_cells(x: Interval, sigma2: float, depth: int = 0):
    F = B.BASE._transition(x)
    Q = B.BASE._scale(B.TUBE.step_scaled_q(x), sigma2)
    Qlower, eps, route = B.BASE._common_point_lower(Q)
    if symmetric_positive_definite_ldlt(Qlower)[0]:
        return [(x, F, B.BASE._transpose(F), Qlower, eps, route, depth)]
    if depth >= MAX_PROCESS_REFINEMENT_DEPTH:
        raise RuntimeError(
            f"cannot obtain strict deterministic process lower on x cell {x.as_list()}"
        )
    mid = math.sqrt(x.lo * x.hi)
    left = Interval.outward_bounds(x.lo, mid)
    right = Interval.outward_bounds(mid, x.hi)
    return (
        _refined_process_cells(left, sigma2, depth + 1)
        + _refined_process_cells(right, sigma2, depth + 1)
    )


def _cell_data(global_x: Interval, sigma2: float):
    leaves = B.TUBE.split_x_cell(global_x)
    rows = []
    max_eps = 0.0
    max_depth = 0
    relative = 0
    for x, _rho in leaves:
        for cell in _refined_process_cells(x, sigma2):
            xx, F, Ft, Qlower, eps, route, depth = cell
            rows.append((xx, F, Ft, Qlower))
            max_eps = max(max_eps, eps)
            max_depth = max(max_depth, depth)
            relative += int(route == "RELATIVE_DIAGONAL")
    if not rows:
        raise RuntimeError("empty refined instantaneous x cover")
    _PROFILE.clear()
    _PROFILE.update({
        "canonical_x_leaves_before_process_refinement": len(leaves),
        "instantaneous_x_cells_after_process_refinement": len(rows),
        "maximum_process_refinement_depth": max_depth,
        "maximum_process_interval_shave": max_eps,
        "relative_process_interval_shaves": relative,
        "process_refinement_is_instantaneous_arithmetic_only": True,
    })
    return rows


# B.build resolves these helpers in the backend module at runtime.  Replace only
# the instantaneous arithmetic cell constructor; all theorem and validation
# logic remains the audited v2 implementation.
B._cell_data = _cell_data


def build(domain_path: Path = DEFAULT_DOMAIN, tube_path: Path | None = None) -> dict:
    d = B.build(domain_path, tube_path)
    profile = d.setdefault("numerical_profile", {})
    profile.update(_PROFILE)
    d["instantaneous_process_cells_refined_for_validated_SPD"] = True
    d["process_refinement_creates_source_history"] = False
    d["instantaneous_x_cells"] = int(
        _PROFILE.get("instantaneous_x_cells_after_process_refinement", d["instantaneous_x_cells"])
    )
    return d


def validate(d: dict) -> list[str]:
    f = list(B.validate(d))
    if d.get("instantaneous_process_cells_refined_for_validated_SPD") is not True:
        f.append("instantaneous process-cell SPD refinement missing")
    if d.get("process_refinement_creates_source_history") is not False:
        f.append("process-cell refinement acquired source history")
    p = d.get("numerical_profile", {})
    if p.get("process_refinement_is_instantaneous_arithmetic_only") is not True:
        f.append("process refinement is not marked instantaneous-only")
    if int(p.get("instantaneous_x_cells_after_process_refinement", 0)) <= 0:
        f.append("refined instantaneous x cover is empty")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--tube", type=Path, default=None)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.tube)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "x_cells": d["instantaneous_x_cells"],
        "global_x_h_over_tau": d["global_x_h_over_tau"],
        "steps": d["prediction_steps"],
        "translation_delta": d["relative_word_injection_floor_lower"],
        "useful_gate": d["useful_gate"],
        "pass": d["pass"],
        "endpoint_diag": d["endpoint_common_lower_diagonal_z"],
        "profile": d["numerical_profile"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
