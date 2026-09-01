#!/usr/bin/env python3
"""Source-audited wrapper for the exact first P5 accelerometer geometry."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p5_first_accel_exact_source as BASE

REPO = Path(__file__).resolve().parents[1]
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = BASE.DEFAULT_DOMAIN
SCHEMA = 2


def _source_semantics() -> tuple[dict, list[str]]:
    k = MEKF.read_text(encoding="utf-8")
    w = WRAPPER.read_text(encoding="utf-8")
    markers = {
        "constructor_zeroes_extended_mean": "xext.setZero();",
        "front_end_does_not_drive_mekf": "updateCore_(dt, gyro, acc, /*tempC=*/35.0f, /*drive_mekf=*/false);",
        "linear_mean_prediction_is_homogeneous": "x_lin_next(i) = sum;",
        "linear_mean_writes_aw_from_homogeneous_result": "xext.template segment<3>(OFF_AW) = x_lin_next.template segment<3>(9);",
        "live_reset_seats_aw_covariance": "reset_aw_covariance_to_stationary();",
        "pseudo_phase_starts_at_zero": "pseudo_update_elapsed_s_ = T(0);",
        "S_pseudo_residual_is_minus_current_S": "const Vector3 r = -xext.template segment<3>(off_S);",
        "S_state_update_is_gain_times_residual": "xext.noalias() += K * r;            // State update",
        "first_accel_reads_aw_mean": "const Vector3 aw = xext.template segment<3>(OFF_AW);",
        "first_accel_force_is_aw_minus_gravity": "f_cog_b = R_wb() * (aw - g_world);",
        "yaw_covariance_axis_is_body_world_down": "Vector3 u_down_body = R_wb() * world_down;",
    }
    joined = k + "\n" + w
    missing = [name for name, marker in markers.items() if marker not in joined]
    ok = lambda name: name not in missing
    return {
        "source_markers": markers,
        "constructor_extended_mean_zero": ok("constructor_zeroes_extended_mean"),
        "mahony_proxy_front_end_leaves_mekf_state_untouched": ok("front_end_does_not_drive_mekf"),
        "first_linear_prediction_is_homogeneous": ok("linear_mean_prediction_is_homogeneous") and ok("linear_mean_writes_aw_from_homogeneous_result"),
        "live_reset_changes_aw_covariance_not_mean": ok("live_reset_seats_aw_covariance") and ok("pseudo_phase_starts_at_zero"),
        "zero_S_mean_makes_first_due_pseudo_mean_correction_zero": ok("S_pseudo_residual_is_minus_current_S") and ok("S_state_update_is_gain_times_residual"),
        "first_accel_force_reads_aw_minus_gravity": ok("first_accel_reads_aw_mean") and ok("first_accel_force_is_aw_minus_gravity"),
        "goLive_yaw_covariance_axis_is_body_gravity_axis": ok("yaw_covariance_axis_is_body_world_down"),
    }, missing


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2) -> dict:
    BASE._source_semantics = _source_semantics
    out = dict(BASE.build(Path(domain_path).resolve(), source_pieces=source_pieces))
    out["schema"] = SCHEMA
    out["qualification"] = "OU3_P5_FIRST_ACCEL_EXACT_STARTUP_SOURCE_GEOMETRY_SOURCE_AUDITED"
    out["first_linear_prediction_homogeneous_zero_mean_certified"] = out["source_semantics"].get("first_linear_prediction_is_homogeneous") is True
    out["first_S_pseudo_zero_residual_mean_identity_certified"] = out["source_semantics"].get("zero_S_mean_makes_first_due_pseudo_mean_correction_zero") is True
    return out


def validate(d: dict) -> list[str]:
    base = dict(d)
    base["schema"] = BASE.SCHEMA
    failures = BASE.validate(base)
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("first_linear_prediction_homogeneous_zero_mean_certified") is not True:
        failures.append("first linear zero-mean propagation is not source certified")
    if d.get("first_S_pseudo_zero_residual_mean_identity_certified") is not True:
        failures.append("first S zero-residual mean identity is not source certified")
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
        "status": out["P5_FIRST_ACCEL_EXACT_SOURCE_CERTIFICATE"],
        "cells": out["evaluated_source_phase_cells"],
        "q_full": out["post_prediction_full_cayley_norm_upper"],
        "q_tangent": out["post_prediction_cayley_tangent_norm_upper"],
        "force": out["first_accel_specific_force_magnitude_mps2"],
        "max_K": out["max_Ktheta_norm_upper"],
        "max_d": out["max_first_accelerometer_correction_norm_upper_rad"],
        "margin": out["minimum_correction_range_margin_rad"],
        "over_limit": out["children_above_validated_correction_limit"],
        "first_unclosed": out["first_unclosed_child"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
