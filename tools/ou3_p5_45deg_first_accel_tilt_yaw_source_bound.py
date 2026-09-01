#!/usr/bin/env python3
"""Tilt/yaw-separated signed first-accelerometer bound for the 45 deg P5 entrance.

The signed-source V2 stage preserves correction direction but still lets the
Cayley component tangent to gravity range over the entire full-attitude 45 deg
ball.  That is source-infeasible at first Live.  The startup/heading contract
already certifies the gravity-direction error separately, and the exact-source
first-accelerometer producer converts it after the first prediction into

    ||c_tangent|| <= q_tilt,

while the P5 entrance independently supplies ||c||<=q_45.  The intersection is
a two-coordinate attitude chart: large full-angle allowance may be yaw, but the
accelerometer-observable tangent coordinate is limited by q_tilt.

This wrapper feeds that existing q_tilt into the signed/source-correlated
first-accelerometer calculation.  It does not tighten the theorem domain, add a
new assumption, change the filter, or claim yaw contraction from the
accelerometer.  The full q_45 radius is retained in every Cayley composition and
in every PSD-remainder term; only the tangent coordinate uses its already
certified source bound.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_candidate_first_accel_exact_source as FIRST
import ou3_p5_45deg_first_accel_signed_source_bound as V1
import ou3_p5_45deg_first_accel_signed_source_bound_v2 as V2

DEFAULT_DOMAIN = V1.DEFAULT_DOMAIN
SCHEMA = 1
DEFAULT_TANGENT_CELLS = 96


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2,
          tangent_cells: int = DEFAULT_TANGENT_CELLS) -> dict:
    path = Path(domain_path).resolve()
    first = FIRST.build(path, source_pieces=source_pieces)
    ff = FIRST.validate(first)
    if ff:
        raise RuntimeError(f"exact-source first-accelerometer prerequisite failed: {ff}")
    q_tilt = float(first["post_prediction_cayley_tangent_norm_upper"])
    if not (math.isfinite(q_tilt) and q_tilt >= 0.0):
        raise RuntimeError("certified post-prediction tangent Cayley bound is invalid")

    original_tangent = V1._tangent_cells
    def source_tangent_cells(q_full: float, n: int):
        return original_tangent(min(float(q_full), q_tilt), n)
    V1._tangent_cells = source_tangent_cells
    try:
        out = dict(V2.build(path, source_pieces=source_pieces,
                            tangent_cells=tangent_cells))
    finally:
        V1._tangent_cells = original_tangent

    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_45DEG_FIRST_ACCEL_TILT_YAW_SEPARATED_SOURCE_BOUND",
        "full_attitude_q_upper_retained": float(out["pre_update_q_upper"]),
        "certified_gravity_tangent_q_upper": q_tilt,
        "tangent_bound_source": "ou3_p5_first_accel_exact_source.post_prediction_cayley_tangent_norm_upper",
        "P5_full_attitude_45deg_domain_tightened": False,
        "new_tilt_assumption_added": False,
        "accelerometer_claims_yaw_contraction": False,
        "two_coordinate_attitude_chart_used": True,
        "signed_V2_post_update_q_upper": float(out["signed_source_correlated_post_update_q_upper"]),
    })
    # The V2 output already used the patched tangent cells, so compare against
    # the old signed V2 by running it once more after restoring the helper.
    baseline = V2.build(path, source_pieces=source_pieces, tangent_cells=tangent_cells)
    out["signed_full_ball_baseline_post_update_q_upper"] = float(
        baseline["signed_source_correlated_post_update_q_upper"])
    out["tilt_yaw_post_update_q_upper"] = float(out["signed_source_correlated_post_update_q_upper"])
    b = out["signed_full_ball_baseline_post_update_q_upper"]
    qnew = out["tilt_yaw_post_update_q_upper"]
    out["tilt_yaw_vs_signed_full_ball_improvement_factor"] = b / qnew if qnew > 0.0 else math.inf
    out["strictly_improves_signed_full_ball"] = qnew < b
    out["P5_45DEG_FIRST_ACCEL_TILT_YAW_SOURCE_BOUND"] = (
        "PASS" if out.get("P5_45DEG_FIRST_ACCEL_SIGNED_SOURCE_BOUND") == "PASS"
        and out["strictly_improves_signed_full_ball"] else "NOT_ESTABLISHED"
    )
    out["next_obligation"] = (
        "propagate the tilt/yaw-separated accepted and LDLT-fallback identity children to sample1; retain the certified physical H group norms and source-correlated Joseph/reset covariance, then let the first source-reachable magnetometer packet act on the yaw coordinate before testing 30deg recapture"
    )
    return out


def validate(d: dict) -> list[str]:
    f = []
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("P5_45DEG_FIRST_ACCEL_SIGNED_SOURCE_BOUND") != "PASS":
        f.append("underlying signed source bound did not pass")
    for k in ("two_coordinate_attitude_chart_used", "strictly_improves_signed_full_ball"):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in ("P5_full_attitude_45deg_domain_tightened", "new_tilt_assumption_added",
              "accelerometer_claims_yaw_contraction", "source_replay_used", "filter_changed",
              "deployed_correction_limit_increased"):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    qfull = float(d.get("full_attitude_q_upper_retained", math.inf))
    qt = float(d.get("certified_gravity_tangent_q_upper", math.inf))
    qbase = float(d.get("signed_full_ball_baseline_post_update_q_upper", math.inf))
    qnew = float(d.get("tilt_yaw_post_update_q_upper", math.inf))
    if not (0.0 <= qt < qfull < 1.0):
        f.append("tilt/full two-coordinate entrance relation invalid")
    if not (math.isfinite(qnew) and 0.0 < qnew < qbase < 8.0):
        f.append("tilt/yaw q bound is not a strict finite improvement")
    if d.get("P5_45DEG_FIRST_ACCEL_TILT_YAW_SOURCE_BOUND") == "PASS" and f:
        f.append("PASS carries validation failures")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--tangent-cells", type=int, default=DEFAULT_TANGENT_CELLS)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(x.domain.resolve(), source_pieces=x.source_pieces, tangent_cells=x.tangent_cells)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_45DEG_FIRST_ACCEL_TILT_YAW_SOURCE_BOUND"],
        "q_full_pre": d["full_attitude_q_upper_retained"],
        "q_tangent_pre": d["certified_gravity_tangent_q_upper"],
        "q_signed_full_ball": d["signed_full_ball_baseline_post_update_q_upper"],
        "q_tilt_yaw": d["tilt_yaw_post_update_q_upper"],
        "improvement_factor": d["tilt_yaw_vs_signed_full_ball_improvement_factor"],
        "returned30": d["returned_to_30deg_P4_sector_here"],
        "validation_failures": vf,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
