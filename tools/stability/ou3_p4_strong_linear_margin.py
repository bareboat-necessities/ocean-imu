#!/usr/bin/env python3
"""Largest source-uniform full-matrix linear margin consumable by P4.

Canonical P3 remains frozen at delta=1e-18.  This module reuses exactly the
same H18 prior-free interval-LDLT sufficient condition, but treats delta as a
query variable and finds a strictly certified lower bound on the largest delta
for which *every* BRMM h/tau cell closes.  The A21 H->A direct-sum release then
caps that H18 margin by its exact active-bias scalar condition.

No P3 matrix, source domain, or acceptance threshold is changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, down, symmetric_positive_definite_ldlt
import ou3_brmm_h18_prior_free_completion as H18
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_full_process_ucc as PROCESS
import ou3_brmm_h18_information_composition as HINFO
import ou3_brmm_live_covariance_seed as LIVE
import ou3_brmm_riccati_metric_p3 as P3

DEFAULT_DOMAIN = H18.DEFAULT_DOMAIN
SCHEMA = 1
QUALIFICATION = "OU3_P4_STRONG_SOURCE_UNIFORM_LINEAR_MARGIN_V1"


def _cell(x: Interval, delta: float, *, process: dict, dynamic: dict, penalty: float):
    if not (0.0 < delta < 1.0):
        return False, -math.inf
    one_minus = H18.I(down(1.0 - delta))
    M = H18._zero(18)
    q_att = float(process["attitude_gyro_bias"]["Q_attitude_gyro_bias_lambda_min_lower"])
    att_diag = down((1.0 - delta) * q_att - penalty)
    if not att_diag > 0.0:
        return False, att_diag
    for i in range(6):
        M[i][i] = H18.I(att_diag)

    inv = dynamic["dynamic_invariant"]
    sigma_floor = float(inv["sigma_aw_filter_mps2"][0])
    h = float(dynamic["validated_rate_and_jump_bounds"]["dt_s"])
    scales = [sigma_floor*h, sigma_floor*h*h, sigma_floor*h*h*h, sigma_floor]
    B = H18.FACTORED.step_scaled_q_over_x(x)
    qscale = H18.I(x.lo)
    Maxis = [[one_minus*qscale*B[i][j] for j in range(4)] for i in range(4)]
    for i in range(4):
        pscaled = H18.TUBE.up(penalty / H18.TUBE.down(scales[i]*scales[i]))
        Maxis[i][i] = Maxis[i][i] - H18.I(pscaled)
    idx = (H18.OFF_V, H18.OFF_P, H18.OFF_S, H18.OFF_AW)
    for axis in range(3):
        for i in range(4):
            for j in range(4):
                M[idx[i]+axis][idx[j]+axis] = Maxis[i][j]
    ok, pivots = symmetric_positive_definite_ldlt(M)
    return bool(ok), min((float(p.lo) for p in pivots), default=-math.inf)


def _h18_test(delta: float, *, leaves, pbar_trace: float, fnorm2: float, process: dict, dynamic: dict):
    penalty = H18.TUBE.up((delta*delta/4.0)*fnorm2*pbar_trace)
    worst = math.inf
    for x in leaves:
        ok, pivot = _cell(x, delta, process=process, dynamic=dynamic, penalty=penalty)
        if not ok:
            return False, pivot, penalty
        worst = min(worst, pivot)
    return True, worst, penalty


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    p3 = P3.build(path)
    if P3.validate(p3) or p3.get("P3_CONDITIONAL_BRMM_PASS") is not True:
        raise RuntimeError("frozen canonical P3 must pass before margin strengthening")
    if float(p3["useful_gate"]) != 1e-18:
        raise RuntimeError("canonical P3 gate changed")

    dynamic = DYNAMIC.build(path)
    process = PROCESS.build()
    hinfo = HINFO.build(path)
    live = LIVE.build(path)
    for name, failures in (
        ("dynamic", DYNAMIC.validate(dynamic)),
        ("process", PROCESS.validate(process)),
        ("H18 information", HINFO.validate(hinfo)),
        ("live seed", LIVE.validate(live)),
    ):
        if failures:
            raise RuntimeError(f"{name} prerequisite failed: {failures}")

    pbar = H18._same_word_covariance_upper(path, dynamic, process, hinfo)
    pbar_trace = float(pbar["Pbar_trace_upper"])
    fnorm2 = float(H18._prediction_norm_sq_upper(process))
    leaves = H18._x_cover(dynamic)

    # Find a failing upper bracket first; 0.5 is certainly below one and the
    # direct matrix condition itself decides whether it remains feasible.
    lo = 1e-18
    ok, pivot, penalty = _h18_test(lo, leaves=leaves, pbar_trace=pbar_trace,
                                  fnorm2=fnorm2, process=process, dynamic=dynamic)
    if not ok:
        raise RuntimeError("frozen P3 delta no longer closes H18")
    hi = 0.5
    if _h18_test(hi, leaves=leaves, pbar_trace=pbar_trace, fnorm2=fnorm2,
                  process=process, dynamic=dynamic)[0]:
        hi = math.nextafter(1.0, 0.0)
    for _ in range(80):
        mid = (lo + hi)/2.0
        ok, _, _ = _h18_test(mid, leaves=leaves, pbar_trace=pbar_trace,
                             fnorm2=fnorm2, process=process, dynamic=dynamic)
        if ok:
            lo = mid
        else:
            hi = mid
    # Reserve two binary64 ulps below the last passing value and verify again.
    h18_delta = math.nextafter(math.nextafter(lo, 0.0), 0.0)
    h18_ok, h18_pivot, h18_penalty = _h18_test(
        h18_delta, leaves=leaves, pbar_trace=pbar_trace, fnorm2=fnorm2,
        process=process, dynamic=dynamic)
    if not h18_ok or not h18_pivot > 0.0:
        raise RuntimeError("strengthened H18 margin did not survive final outward check")

    held = live["held_ba"]
    qba = float(process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"])
    pba0 = float(held["seed_variance"])
    a21_cap = down(qba/(qba+pba0))
    a21_delta = math.nextafter(min(h18_delta, a21_cap), 0.0)
    ba_margin = down((1.0-a21_delta)*qba-a21_delta*pba0)
    if not (a21_delta >= 1e-18 and ba_margin > 0.0):
        raise RuntimeError("strengthened A21 direct-sum margin is not strict")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": p3["canonical_source"],
        "P3_frozen_gate": 1e-18,
        "P3_numerical_matrices_changed": False,
        "source_domain_changed": False,
        "trajectory_replay_used": False,
        "same_H18_prior_free_interval_LDLT_condition_reused": True,
        "same_H18_x_cover_reused": True,
        "x_cell_count": len(leaves),
        "same_word_Pbar_trace_upper": pbar_trace,
        "prediction_F_norm_sq_upper": fnorm2,
        "H18": {
            "strong_relative_margin_lower": h18_delta,
            "worst_interval_LDLT_pivot_lower_at_strong_margin": h18_pivot,
            "delta_squared_penalty_at_strong_margin": h18_penalty,
            "strict": True,
        },
        "A21": {
            "strong_relative_margin_lower": a21_delta,
            "active_ba_direct_sum_delta_cap": a21_cap,
            "first_active_ba_margin_lower": ba_margin,
            "strict": True,
        },
        "both_stronger_than_frozen_gate": h18_delta > 1e-18 and a21_delta > 1e-18,
        "P4_may_consume_strong_margin_without_modifying_P3": True,
    }


def validate(d: dict) -> list[str]:
    f=[]
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if float(d.get("P3_frozen_gate", 0.0)) != 1e-18:
        f.append("frozen P3 gate changed")
    for k in ("same_H18_prior_free_interval_LDLT_condition_reused", "same_H18_x_cover_reused",
              "both_stronger_than_frozen_gate", "P4_may_consume_strong_margin_without_modifying_P3"):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in ("P3_numerical_matrices_changed", "source_domain_changed", "trajectory_replay_used"):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    for mode in ("H18","A21"):
        row=d.get(mode,{})
        if row.get("strict") is not True or not float(row.get("strong_relative_margin_lower",0.0)) > 1e-18:
            f.append(f"{mode} strong margin not strict/useful")
    return f


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args=ap.parse_args()
    d=build(args.domain)
    failures=validate(d)
    d["validation_pass"]=not failures
    d["validation_failures"]=failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "H18_delta":d["H18"]["strong_relative_margin_lower"],
        "A21_delta":d["A21"]["strong_relative_margin_lower"],
        "H18_pivot":d["H18"]["worst_interval_LDLT_pivot_lower_at_strong_margin"],
        "A21_ba_margin":d["A21"]["first_active_ba_margin_lower"],
        "failures":failures,
    },indent=2,sort_keys=True))
    return int(bool(failures))

if __name__=="__main__":
    raise SystemExit(main())
