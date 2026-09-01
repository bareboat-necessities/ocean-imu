#!/usr/bin/env python3
"""Deployed-startup tilt/yaw subroute for the 45 deg P5 capture problem.

The generic signed-source P5 stage covers every state in the declared full
SO(3) 45 deg entrance ball.  This producer does something deliberately
narrower and must not be confused with that generic route: it intersects the
same 45 deg full-attitude bound with the *additional* gravity-direction bound
already established by P1 for source-reachable deployed startup handoffs.

The exact-source first-accelerometer certificate converts that startup tilt
information after the first prediction into

    ||c_tangent|| <= q_tilt,

while the full attitude still satisfies ||c||<=q_45.  The resulting two-
coordinate chart is useful when composing P1 into P5, but it does not prove
capture for an arbitrary abstract 45 deg P5 entrance state whose tilt could use
the whole 45 deg allowance.

No new deployment assumption is introduced: q_tilt is imported from P1 through
the source-audited exact first-Live certificate.  Nevertheless it is additional
information relative to the standalone P5 entrance, so this producer is marked
as a startup-source subroute and is forbidden from replacing the generic 45 deg
P5 obligation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_first_accel_exact_source_v2 as EXACT
import ou3_p5_45deg_first_accel_signed_source_bound as V1
import ou3_p5_45deg_first_accel_signed_source_bound_v2 as V2

DEFAULT_DOMAIN = V1.DEFAULT_DOMAIN
SCHEMA = 3
DEFAULT_TANGENT_CELLS = 96


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2,
          tangent_cells: int = DEFAULT_TANGENT_CELLS) -> dict:
    path = Path(domain_path).resolve()
    exact = EXACT.build(path, source_pieces=source_pieces)
    ef = EXACT.validate(exact)
    if ef:
        raise RuntimeError(f"exact-source first-accelerometer prerequisite failed: {ef}")
    q_tilt = float(exact["post_prediction_cayley_tangent_norm_upper"])
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
        "qualification": "OU3_P5_DEPLOYED_STARTUP_TILT_YAW_FIRST_ACCEL_SUBROUTE",
        "full_attitude_q_upper_retained": float(out["pre_update_q_upper"]),
        "certified_gravity_tangent_q_upper": q_tilt,
        "tangent_bound_source": "P1 via ou3_p5_first_accel_exact_source_v2.post_prediction_cayley_tangent_norm_upper",
        "source_tilt_cosine_lower": float(exact["post_prediction_true_gravity_cosine_lower"]),
        "source_reachable_startup_intersection_only": True,
        "uses_additional_P1_tilt_information": True,
        "generic_P5_45deg_entrance_covered_here": False,
        "does_not_replace_generic_P5_45deg_route": True,
        "new_deployment_assumption_added": False,
        "accelerometer_claims_yaw_contraction": False,
        "two_coordinate_attitude_chart_used": True,
    })
    baseline = V2.build(path, source_pieces=source_pieces, tangent_cells=tangent_cells)
    out["signed_full_ball_baseline_post_update_q_upper"] = float(
        baseline["signed_source_correlated_post_update_q_upper"])
    out["startup_intersection_post_update_q_upper"] = float(
        out["signed_source_correlated_post_update_q_upper"])
    b = out["signed_full_ball_baseline_post_update_q_upper"]
    qnew = out["startup_intersection_post_update_q_upper"]
    out["startup_intersection_improvement_factor"] = b / qnew if qnew > 0.0 else math.inf
    out["strictly_improves_source_reachable_startup_subroute"] = qnew < b
    out["P5_DEPLOYED_STARTUP_TILT_YAW_FIRST_ACCEL_SUBROUTE"] = (
        "PASS" if out.get("P5_45DEG_FIRST_ACCEL_SIGNED_SOURCE_BOUND") == "PASS"
        and out["strictly_improves_source_reachable_startup_subroute"] else "NOT_ESTABLISHED"
    )
    out["next_obligation"] = (
        "for deployed startup composition, propagate this P1-intersected child toward magnetic yaw correction; for the standalone generic 45deg P5 entrance, keep the full-ball signed route and preserve the joint attitude-a_w Kalman correction instead of adding the P1 tilt restriction"
    )
    return out


def validate(d: dict) -> list[str]:
    f = []
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("P5_45DEG_FIRST_ACCEL_SIGNED_SOURCE_BOUND") != "PASS":
        f.append("underlying generic signed source bound did not pass")
    for k in (
        "two_coordinate_attitude_chart_used",
        "source_reachable_startup_intersection_only",
        "uses_additional_P1_tilt_information",
        "does_not_replace_generic_P5_45deg_route",
        "strictly_improves_source_reachable_startup_subroute",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "generic_P5_45deg_entrance_covered_here",
        "new_deployment_assumption_added",
        "accelerometer_claims_yaw_contraction",
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    qfull = float(d.get("full_attitude_q_upper_retained", math.inf))
    qt = float(d.get("certified_gravity_tangent_q_upper", math.inf))
    qbase = float(d.get("signed_full_ball_baseline_post_update_q_upper", math.inf))
    qnew = float(d.get("startup_intersection_post_update_q_upper", math.inf))
    ctilt = float(d.get("source_tilt_cosine_lower", -math.inf))
    if not (0.0 <= qt < qfull < 1.0):
        f.append("startup tilt/full two-coordinate relation invalid")
    if not (0.0 < ctilt <= 1.0):
        f.append("source tilt cosine is invalid")
    if not (math.isfinite(qnew) and 0.0 < qnew < qbase < 8.0):
        f.append("startup-intersection q bound is not a strict finite improvement")
    if d.get("P5_DEPLOYED_STARTUP_TILT_YAW_FIRST_ACCEL_SUBROUTE") == "PASS" and f:
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
        "status": d["P5_DEPLOYED_STARTUP_TILT_YAW_FIRST_ACCEL_SUBROUTE"],
        "tilt_cos": d["source_tilt_cosine_lower"],
        "q_full_pre": d["full_attitude_q_upper_retained"],
        "q_tangent_pre": d["certified_gravity_tangent_q_upper"],
        "q_generic_signed": d["signed_full_ball_baseline_post_update_q_upper"],
        "q_startup_intersection": d["startup_intersection_post_update_q_upper"],
        "improvement_factor": d["startup_intersection_improvement_factor"],
        "generic_45deg_covered": d["generic_P5_45deg_entrance_covered_here"],
        "validation_failures": vf,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
