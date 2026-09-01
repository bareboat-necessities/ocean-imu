#!/usr/bin/env python3
"""Run the #450 H18 word screen through the shared #449/#450 operations.

The original H18 screen predates extraction of
:mod:`ou3_p4_h18_differential_operations`.  This wrapper makes the numerical
screen consume that common prediction/residual/Joseph-reset implementation
without duplicating the expensive screen harness.  The injected bindings are
scoped and restored after the calculation.

This remains a screening stage: source-orientation coverage, optional branch
coverage and actual endpoint source-node metric whitening are still open, so no
P4 promotion is possible here even if the screening norm is below one.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p4_h18_differential_operations as DOPS
import ou3_p4_h18_interval_ad_word as SCREEN
from ou3_proof_module_state import preserve_module_bindings

DEFAULT_DOMAIN = SCREEN.DEFAULT_DOMAIN
SCHEMA = 1


def build(domain_path: Path = DEFAULT_DOMAIN, *, samples: int | None = None,
          cell_limit: int = 1, ball_inflation: float = 1.5) -> dict:
    with preserve_module_bindings():
        # Install the shared differential layer into the established screen
        # harness.  The harness still owns source scheduling/covariance cells,
        # conditioning diagnostics and fail-closed promotion logic.
        SCREEN._prediction = DOPS.prediction
        SCREEN._rotation_residual_acc = DOPS.accelerometer_residual
        SCREEN._rotation_residual_mag = DOPS.magnetometer_residual
        SCREEN._residual_S = DOPS.S_residual
        SCREEN._H_acc_canonical = DOPS.H_acc_canonical
        SCREEN._H_mag_canonical = DOPS.H_mag_canonical
        SCREEN._accepted_update = DOPS.accepted_update
        out = dict(SCREEN.build(
            Path(domain_path).resolve(),
            samples=samples,
            cell_limit=cell_limit,
            ball_inflation=ball_inflation,
        ))
    out["shared_H18_differential_operations_used"] = True
    out["shared_operation_module"] = "ou3_p4_h18_differential_operations"
    out["qualification"] = "OU3_P4_H18_SHARED_DIFFERENTIAL_WORD_SCREEN"
    return out


def validate(d: dict) -> list[str]:
    base = dict(d)
    base["qualification"] = "OU3_P4_H18_INTERVAL_AD_WORD_SCREEN"
    failures = SCREEN.validate(base)
    if d.get("qualification") != "OU3_P4_H18_SHARED_DIFFERENTIAL_WORD_SCREEN":
        failures.append("wrong shared H18 screen qualification")
    if d.get("shared_H18_differential_operations_used") is not True:
        failures.append("shared H18 differential operations are not active")
    if d.get("shared_operation_module") != "ou3_p4_h18_differential_operations":
        failures.append("wrong shared differential operation module")
    if d.get("P4_USABLE_CERTIFICATE_PROMOTED") is not False:
        failures.append("shared H18 screen improperly promoted P4")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--samples", type=int, default=None)
    ap.add_argument("--cell-limit", type=int, default=1)
    ap.add_argument("--ball-inflation", type=float, default=1.5)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(
        args.domain.resolve(), samples=args.samples,
        cell_limit=args.cell_limit, ball_inflation=args.ball_inflation,
    )
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "validation_pass": not vf,
        "shared_operations": out["shared_H18_differential_operations_used"],
        "samples": [out["samples_checked"], out["full_word_samples"]],
        "cells": [out["outer_ball_box_cells_completed"], out["outer_ball_box_cells_requested"], out["outer_ball_box_cover_total"]],
        "max_endpoint_conditioned_norm": out["max_endpoint_P3_congruence_conditioned_norm_upper"],
        "max_prefix_conditioned_norm": out["max_prefix_P3_congruence_conditioned_norm_upper"],
        "screen_gamma_lt_one": out["H18_SCREEN_GAMMA_LT_ONE"],
        "first_failure": out["first_failure"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
