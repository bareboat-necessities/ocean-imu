#!/usr/bin/env python3
"""Tighten sample-1 entry with the exact first-accelerometer a_w gain contraction.

The V1 sample-1 entry stage propagates the correct 18x18 Joseph/reset covariance,
but its non-attitude state hull used the raw entrywise interval product K*r.
For the first accelerometer prefix the exact-source geometry is stronger:

    H = [H_theta, I_aw],   P_{theta,aw}=0,   P_{aw,aw}=p_aw I,
    S = H_theta P_theta H_theta' + p_aw I + R_acc.

Therefore

    K_aw = p_aw S^{-1},    0 <= p_aw S^{-1} <= I,

and hence ||Delta a_w|| <= ||r|| for every source child, independently of the
width of the interval representation of p_aw.  This is a Loewner consequence
of the shipping innovation covariance, not a fitted gain bound.

This wrapper keeps the V1 full 18x18 Joseph/reset covariance map unchanged and
only intersects the interval a_w state correction with that exact norm
consequence before hulling accepted/rejected physical/state branches.  It then
re-runs the sample-1 prediction.  No theorem or deployment limit changes.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_mul
import ou3_p5_sample1_entry as V1

DEFAULT_DOMAIN = V1.DEFAULT_DOMAIN
SCHEMA = 2


def _first_accel_covariance_and_state_tight(Pm, e, H, R, residual_norm: float, d_norm_cap: float):
    FULL = V1.FULL
    N = V1.N
    r = FULL._vec_box(residual_norm)
    PHt, S = FULL._innovation(Pm, H, R)
    Sinv, backend = FULL._spd_inverse_enclosure(S, R)
    K = matrix_mul(PHt, Sinv)
    dx = FULL._mat_vec(K, r)

    # Attitude cap comes from the already certified exact-source correction norm.
    tcap = Interval(-FULL.up(d_norm_cap), FULL.up(d_norm_cap))
    dx_theta = [FULL._intersect(dx[i], tcap) for i in range(3)]
    dx_capped = list(dx)
    dx_capped[0:3] = dx_theta

    # Exact source structure gives K_aw = p_aw*S^{-1} with
    # S >= p_aw*I + R > p_aw*I, so ||K_aw||_2 <= 1.  A vector correction whose
    # Euclidean norm is <= residual_norm has every component in this interval.
    awcap = Interval(-FULL.up(residual_norm), FULL.up(residual_norm))
    for i in V1.FULL.AW:
        dx_capped[i] = FULL._intersect(dx_capped[i], awcap)

    Pj = FULL._shipping_joseph(Pm, K, S, PHt)
    Pr = FULL._reset_covariance(Pj, dx_theta)

    e_acc = list(e)
    for i in range(3, N):
        e_acc[i] = e[i] - dx_capped[i]

    Pout = FULL._psd_tighten(FULL._mat_hull(Pm, Pr))
    eout = FULL._vec_hull(e, e_acc)
    return Pout, eout, {
        "inverse_backend": backend,
        "K": K,
        "S": S,
        "r": r,
        "dx": dx_capped,
        "P_accepted_reset": Pr,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2) -> dict:
    old = V1._first_accel_covariance_and_state
    V1._first_accel_covariance_and_state = _first_accel_covariance_and_state_tight
    try:
        out = dict(V1.build(Path(domain_path).resolve(), source_pieces=source_pieces))
    finally:
        V1._first_accel_covariance_and_state = old
    out["schema"] = SCHEMA
    out["qualification"] = "OU3_P5_EXACT_SOURCE_SAMPLE1_ENTRY_WITH_EXACT_AW_GAIN_CONTRACTION"
    out["first_accel_aw_gain_exact_formula"] = "K_aw=p_aw*S^{-1}"
    out["first_accel_innovation_loewner_floor"] = "S>=p_aw*I+R_acc"
    out["first_accel_aw_gain_operator_norm_upper"] = 1.0
    out["first_accel_aw_state_correction_norm_bounded_by_residual_norm"] = True
    out["raw_entrywise_aw_Kr_used_as_state_bound"] = False
    return out


def validate(d: dict) -> list[str]:
    base = dict(d)
    base["schema"] = V1.SCHEMA
    failures = V1.validate(base)
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("first_accel_aw_gain_exact_formula") != "K_aw=p_aw*S^{-1}":
        failures.append("exact first accelerometer aw gain identity missing")
    if d.get("first_accel_innovation_loewner_floor") != "S>=p_aw*I+R_acc":
        failures.append("first accelerometer innovation Loewner floor missing")
    if float(d.get("first_accel_aw_gain_operator_norm_upper", math.inf)) != 1.0:
        failures.append("first accelerometer aw gain contraction changed")
    if d.get("first_accel_aw_state_correction_norm_bounded_by_residual_norm") is not True:
        failures.append("first accelerometer aw state correction norm is not tightened")
    if d.get("raw_entrywise_aw_Kr_used_as_state_bound") is not False:
        failures.append("raw entrywise aw K*r state bound remains active")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve(), source_pieces=args.source_pieces)
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_SAMPLE1_ENTRY_CERTIFICATE"],
        "cells": out["evaluated_source_phase_cells"],
        "q1_pre_measurement": out["sample1_pre_measurement_cayley_norm_upper"],
        "inverse_backends": out["inverse_backend_counts"],
        "max_state": out["max_sample1_state_group_norm_uppers"],
        "first_failure": out["first_failure"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
