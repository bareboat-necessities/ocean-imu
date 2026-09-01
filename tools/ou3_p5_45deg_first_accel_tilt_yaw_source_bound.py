#!/usr/bin/env python3
"""Tilt/yaw-separated signed first-accelerometer bound for deployed startup P5.

The generic signed 45 deg first-accelerometer stage preserves correction
orientation but still lets the gravity-tangent Cayley coordinate range over the
entire 45 deg ball.  The deployed startup route has additional information from
P1: the gravity-direction error is already certified separately from heading.
Intersecting those two already-proved facts gives a source-reachable two-
coordinate chart:

    ||c|| <= q_45,        ||c_tangent|| <= q_tilt.

This is a deployed startup subroute, not a proof for every abstract state in the
standalone 45 deg P5 entrance set.  No new deployment assumption is introduced
and the generic P5 route remains separate.

The producer also carries the resulting full Cayley bound through the *next*
shipping 5 ms prediction.  Prediction uses only the already-declared gyro-bias
and deterministic transport bounds through the same exact Cayley composition
helper used by the first-prefix proof.  This closes a bookkeeping gap in the
previous continuation: the strong startup intersection is no longer discarded
before sample 1.  Sample-1 S/accelerometer/magnetometer corrections are still
left to the next source-correlated numerical stage.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_first_accel_exact_source_v2 as EXACT
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_45deg_first_accel_signed_source_bound as V1
import ou3_p5_45deg_first_accel_signed_source_bound_v2 as V2

DEFAULT_DOMAIN = V1.DEFAULT_DOMAIN
SCHEMA = 3
DEFAULT_TANGENT_CELLS = 96


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2,
          tangent_cells: int = DEFAULT_TANGENT_CELLS) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
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
        "tangent_bound_source": "ou3_p5_first_accel_exact_source_v2.post_prediction_cayley_tangent_norm_upper",
        "source_tilt_cosine_lower": float(exact["post_prediction_true_gravity_cosine_lower"]),
        "source_reachable_startup_intersection_only": True,
        "uses_additional_P1_tilt_information": True,
        "does_not_replace_generic_P5_45deg_route": True,
        "generic_P5_45deg_entrance_covered_here": False,
        "new_deployment_assumption_added": False,
        "accelerometer_claims_yaw_contraction": False,
        "two_coordinate_attitude_chart_used": True,
        "signed_V2_post_update_q_upper": float(out["signed_source_correlated_post_update_q_upper"]),
    })

    baseline = V2.build(path, source_pieces=source_pieces, tangent_cells=tangent_cells)
    out["signed_full_ball_baseline_post_update_q_upper"] = float(
        baseline["signed_source_correlated_post_update_q_upper"])
    out["startup_intersection_post_update_q_upper"] = float(
        out["signed_source_correlated_post_update_q_upper"])
    b = out["signed_full_ball_baseline_post_update_q_upper"]
    qnew = out["startup_intersection_post_update_q_upper"]
    out["startup_intersection_improvement_factor"] = b / qnew if qnew > 0.0 else math.inf
    out["strictly_improves_generic_signed_route"] = qnew < b

    h = float(FULL._source_cell()["dt_s"])
    q1_startup = RG._q_after_first_prediction(qnew, domain, h)
    q1_generic = RG._q_after_first_prediction(b, domain, h)
    out["next_prediction_dt_s"] = h
    out["startup_subroute_sample1_pre_measurement_q_upper"] = q1_startup
    out["generic_signed_sample1_pre_measurement_q_upper"] = q1_generic
    out["startup_sample1_vs_generic_signed_improvement_factor"] = (
        q1_generic / q1_startup if q1_startup > 0.0 else math.inf
    )
    out["startup_sample1_pre_measurement_inside_q8"] = q1_startup < 8.0
    out["generic_signed_sample1_pre_measurement_inside_q8"] = q1_generic < 8.0
    out["sample1_measurements_evaluated_here"] = False
    out["sample1_source_phase_transition_classified_here"] = False
    out["sample1_to_30deg_recapture_established_here"] = False

    out["P5_DEPLOYED_STARTUP_TILT_YAW_FIRST_ACCEL_SUBROUTE"] = (
        "PASS" if out.get("P5_45DEG_FIRST_ACCEL_SIGNED_SOURCE_BOUND") == "PASS"
        and out["strictly_improves_generic_signed_route"]
        and out["startup_sample1_pre_measurement_inside_q8"] else "NOT_ESTABLISHED"
    )
    out["next_obligation"] = (
        "start the sample1 due/not-due S and accepted/identity vector-prefix enclosure from startup_subroute_sample1_pre_measurement_q_upper rather than the generic q<8 chart; retain the physical H group-norm caps and source-correlated Joseph/reset covariance, then apply the first source-reachable magnetic information before testing 30deg recapture"
    )
    return out


def validate(d: dict) -> list[str]:
    f = []
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("P5_45DEG_FIRST_ACCEL_SIGNED_SOURCE_BOUND") != "PASS":
        f.append("underlying signed source bound did not pass")
    for k in (
        "two_coordinate_attitude_chart_used",
        "strictly_improves_generic_signed_route",
        "source_reachable_startup_intersection_only",
        "uses_additional_P1_tilt_information",
        "does_not_replace_generic_P5_45deg_route",
        "startup_sample1_pre_measurement_inside_q8",
        "generic_signed_sample1_pre_measurement_inside_q8",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "generic_P5_45deg_entrance_covered_here", "new_deployment_assumption_added",
        "accelerometer_claims_yaw_contraction", "source_replay_used", "filter_changed",
        "deployed_correction_limit_increased", "sample1_measurements_evaluated_here",
        "sample1_source_phase_transition_classified_here", "sample1_to_30deg_recapture_established_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    qfull = float(d.get("full_attitude_q_upper_retained", math.inf))
    qt = float(d.get("certified_gravity_tangent_q_upper", math.inf))
    qbase = float(d.get("signed_full_ball_baseline_post_update_q_upper", math.inf))
    qnew = float(d.get("startup_intersection_post_update_q_upper", math.inf))
    q1base = float(d.get("generic_signed_sample1_pre_measurement_q_upper", math.inf))
    q1new = float(d.get("startup_subroute_sample1_pre_measurement_q_upper", math.inf))
    ctilt = float(d.get("source_tilt_cosine_lower", -math.inf))
    if not (0.0 <= qt < qfull < 1.0):
        f.append("tilt/full two-coordinate entrance relation invalid")
    if not (0.0 < ctilt <= 1.0):
        f.append("source tilt cosine is invalid")
    if not (math.isfinite(qnew) and 0.0 < qnew < qbase < 8.0):
        f.append("startup intersection q bound is not a strict finite improvement")
    if not (math.isfinite(q1new) and 0.0 < q1new < q1base < 8.0):
        f.append("startup sample1 prediction did not preserve the strict chart improvement")
    factor = float(d.get("startup_sample1_vs_generic_signed_improvement_factor", 0.0))
    if not (math.isfinite(factor) and factor > 1.0):
        f.append("startup sample1 improvement factor is not strict")
    dt = float(d.get("next_prediction_dt_s", -1.0))
    # FULL._source_cell() deliberately carries outward-rounded source literals.
    # Accept that certified enclosure instead of requiring exactly one binary64
    # ULP around the decimal spelling 0.005.  Four ULPs is still <4e-18 s and
    # cannot hide a different deployed cadence.
    if not math.isclose(dt, 0.005, rel_tol=0.0, abs_tol=4.0 * math.ulp(0.005)):
        f.append("next prediction step is not the deployed 5 ms interval")
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
        "q_generic_signed_post_accel": d["signed_full_ball_baseline_post_update_q_upper"],
        "q_startup_post_accel": d["startup_intersection_post_update_q_upper"],
        "q_generic_sample1_pre": d["generic_signed_sample1_pre_measurement_q_upper"],
        "q_startup_sample1_pre": d["startup_subroute_sample1_pre_measurement_q_upper"],
        "sample1_improvement_factor": d["startup_sample1_vs_generic_signed_improvement_factor"],
        "returned30": d["returned_to_30deg_P4_sector_here"],
        "validation_failures": vf,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
