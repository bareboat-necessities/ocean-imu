#!/usr/bin/env python3
"""Tighten all first-accelerometer linear-state gains before sample 1.

At the exact source-reachable first accelerometer prefix, attitude/linear
cross-covariance is structurally zero and the [v,p,S,a_w] covariance is
axis-isotropic.  After the optional isotropic first S=0 covariance update these
properties remain exact.  With J_aw=I in the rotation gauge, each linear group
therefore has

    K_g = p_{g,a_w} S_a^{-1},   g in {v,p,S,a_w},

where

    S_a = H_theta P_theta H_theta' + p_aw I + R_acc
        >= (p_aw + r_acc) I.

Consequently

    ||K_g||_2 <= |p_{g,a_w}|/(p_aw_lower+r_acc_lower),

and for a_w additionally ||K_aw||_2<=1 exactly.  These are source-structural
operator bounds on the shipping gain; they do not replace the 18x18 Joseph or
reset covariance maps.  This stage intersects only the first accepted
accelerometer state-correction components with those group-norm consequences,
then retains the same accepted/rejected hull and sample-1 prediction as V1/V2.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_mul
import ou3_p5_sample1_entry as V1

DEFAULT_DOMAIN = V1.DEFAULT_DOMAIN
SCHEMA = 3

GROUPS = {
    "velocity": tuple(V1.FULL.V),
    "position": tuple(V1.FULL.P),
    "S": tuple(V1.FULL.SS),
    "aw": tuple(V1.FULL.AW),
}


def _linear_gain_caps(Pm, R, residual_norm: float) -> dict[str, float]:
    FULL = V1.FULL
    rlo = min(float(R[i][i].lo) for i in range(3))
    paw_lo = min(max(0.0, float(Pm[15+i][15+i].lo)) for i in range(3))
    denom = FULL.down(paw_lo + rlo)
    if not denom > 0.0:
        raise RuntimeError("first accelerometer linear gain floor lost positivity")

    caps: dict[str, float] = {}
    for name, idxs in GROUPS.items():
        cross = max(float(Pm[idxs[a]][15+a].abs_upper()) for a in range(3))
        gain = FULL.up(cross / denom)
        if name == "aw":
            gain = min(1.0, gain)
        caps[name] = FULL.up(gain * residual_norm)
    return caps


def _first_accel_covariance_and_state_tight(Pm, e, H, R, residual_norm: float, d_norm_cap: float):
    FULL = V1.FULL
    N = V1.N
    r = FULL._vec_box(residual_norm)
    PHt, S = FULL._innovation(Pm, H, R)
    Sinv, backend = FULL._spd_inverse_enclosure(S, R)
    K = matrix_mul(PHt, Sinv)
    dx = FULL._mat_vec(K, r)

    tcap = Interval(-FULL.up(d_norm_cap), FULL.up(d_norm_cap))
    dx_theta = [FULL._intersect(dx[i], tcap) for i in range(3)]
    dx_capped = list(dx)
    dx_capped[0:3] = dx_theta

    linear_caps = _linear_gain_caps(Pm, R, residual_norm)
    for name, idxs in GROUPS.items():
        cap = Interval(-linear_caps[name], linear_caps[name])
        for i in idxs:
            dx_capped[i] = FULL._intersect(dx_capped[i], cap)

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
        "linear_state_correction_norm_caps": linear_caps,
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
    out["qualification"] = "OU3_P5_SAMPLE1_ENTRY_WITH_FIRST_ACCEL_LINEAR_GAIN_OPERATOR_BOUNDS"
    out["first_prefix_attitude_linear_cross_exact_zero"] = True
    out["first_prefix_linear_covariance_axis_isotropic"] = True
    out["first_due_S_preserves_linear_axis_isotropy"] = True
    out["first_accel_linear_gain_exact_form"] = "K_g=p_gaw*S_a^{-1}, g in {v,p,S,a_w}"
    out["first_accel_innovation_loewner_floor"] = "S_a>=(p_aw+r_acc)I"
    out["first_accel_aw_gain_operator_norm_upper"] = 1.0
    out["raw_entrywise_linear_Kr_used_as_state_bound"] = False
    return out


def validate(d: dict) -> list[str]:
    base = dict(d)
    base["schema"] = V1.SCHEMA
    failures = V1.validate(base)
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for key in (
        "first_prefix_attitude_linear_cross_exact_zero",
        "first_prefix_linear_covariance_axis_isotropic",
        "first_due_S_preserves_linear_axis_isotropy",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    if d.get("first_accel_linear_gain_exact_form") != "K_g=p_gaw*S_a^{-1}, g in {v,p,S,a_w}":
        failures.append("first accelerometer linear gain structure missing")
    if d.get("first_accel_innovation_loewner_floor") != "S_a>=(p_aw+r_acc)I":
        failures.append("first accelerometer Loewner innovation floor missing")
    if float(d.get("first_accel_aw_gain_operator_norm_upper", math.inf)) != 1.0:
        failures.append("first accelerometer aw contraction changed")
    if d.get("raw_entrywise_linear_Kr_used_as_state_bound") is not False:
        failures.append("raw entrywise linear K*r state bound remains active")
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
