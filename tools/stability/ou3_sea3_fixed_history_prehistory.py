#!/usr/bin/env python3
"""Same-driver startup prehistory for the legal fixed complete-SEA3 member.

The feasibility member in ``ou3_sea3_fixed_history_source_core`` chooses

    a = K* e_0 / sqrt(Q_00),

for one fixed complete-SEA3 lambda/response history.  Therefore its physical
vertical acceleration is the continuum covariance section

    y(t) = Q(t,0) / sqrt(Q(0,0)).

For the fixed real, stationary response used by that member, Q(t,0) depends on
lag through cos(omega*t), so the *negative-time* startup history is determined
by the same exact continuum driver: y(-t)=y(t).  No warm-up source, replay,
finite harmonic realization, or independently chosen startup signal is added.

This module evaluates enough of that same continuum member before t=0 to drive
the deployed measurement-only frontend through WavePeriodEstimator usability
and TunerReady.  Quadrature remains an evaluation of the continuum integral,
not the source definition.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_sea3_fixed_history_source_core as CORE

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_FIXED_HISTORY_SAME_DRIVER_PREHISTORY_V1"
DT = CORE.DT
# 45 s clears the 4/lambda WPE usable floor (~31.83 s) and the 10 s outer
# online-tune warmup while leaving margin for the one-current-period history
# requirement.  Readiness is still checked by the actual C++ frontend next;
# this duration is not itself a proof that TunerReady occurs.
PREHISTORY_S = 45.0
PRE_SAMPLES = int(round(PREHISTORY_S / DT))
PANELS = 8192


def _same_driver_value_at_lag(lag_samples: int, panels: int, scale: float, q00: float) -> float:
    if lag_samples < 0:
        raise ValueError("lag_samples must be nonnegative")
    q = CORE._gram_lag(lag_samples, panels, scale)
    return q / math.sqrt(q00)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    source = CORE.build(Path(domain_path).resolve())
    failures = CORE.validate(source)
    if failures:
        raise RuntimeError(f"fixed-history source core invalid: {failures}")

    scale, _ = CORE._normalized_spectrum_scale(PANELS)
    q00 = CORE._gram_lag(0, PANELS, scale)
    if not (math.isfinite(q00) and q00 > 0.0):
        raise RuntimeError("same-driver Q00 is not positive finite")

    # k=-PRE_SAMPLES,...,-1.  Lag from t=0 is |k|.  The t=0 sample itself is
    # owned by the canonical 601-sample word and is not duplicated here.
    prehistory = []
    max_abs = 0.0
    for j in range(PRE_SAMPLES, 0, -1):
        y = _same_driver_value_at_lag(j, PANELS, scale, q00)
        max_abs = max(max_abs, abs(y))
        prehistory.append({
            "k": -j,
            "t_s": -j * DT,
            "f_cog_vertical_acceleration_mps2": y,
            "gyro_measurement_rad_s": [0.0, 0.0, 0.0],
            "omega_body_corrected_rad_s": [0.0, 0.0, 0.0],
            # Identity attitude, no lever arm: accelerometer specific force is
            # a_non-grav - g e_D in the BODY=NED frame used by the fixed member.
            "specific_force_body_mps2": [0.0, 0.0, y - 9.80665],
            "R_wb": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        })

    first_word = source["source_core"][0]["f_cog_body_mps2"][2]
    y0 = _same_driver_value_at_lag(0, PANELS, scale, q00)
    join_error = abs(float(first_word) - y0)
    live_cap = json.loads(Path(domain_path).read_text(encoding="utf-8"))["normal_live"]
    acc_cap = float(live_cap["non_gravitational_cog_acceleration_norm_upper_mps2"])

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "fixed_history_source_qualification": source["qualification"],
        "driver_choice": "a=K*e0/sqrt(Q00)",
        "same_driver_field_as_word": True,
        "stationary_fixed_history_correlation_extension": True,
        "negative_time_relation": "y(-t)=Q(t,0)/sqrt(Q00)=y(t)",
        "prehistory_duration_s": PREHISTORY_S,
        "dt_s": DT,
        "prehistory_samples": PRE_SAMPLES,
        "prehistory": prehistory,
        "word_join": {
            "same_t0_value": join_error < 2.0e-12,
            "absolute_join_error_mps2": join_error,
            "t0_word_value_mps2": float(first_word),
            "t0_recomputed_same_driver_mps2": y0,
        },
        "normal_live_acceleration_cap_mps2": acc_cap,
        "max_abs_prehistory_non_grav_acceleration_mps2": max_abs,
        "prehistory_inside_declared_acceleration_cap": max_abs <= acc_cap,
        "wpe_usable_time_floor_s": 4.0 / (2.0 * math.pi * 0.02),
        "outer_online_tune_warmup_s": float(
            json.loads(Path(domain_path).read_text(encoding="utf-8"))["startup"]["online_tune_warmup_sec"]
        ),
        "actual_cpp_frontend_TunerReady_verified_here": False,
        "trajectory_replay_used": False,
        "finite_harmonic_source_used": False,
        "independent_startup_signal_used": False,
        "independent_sample_boxes_used": False,
        "P3_changed": False,
        "P4_promoted": False,
        "next_obligation": (
            "feed these exact same-driver negative-time samples to the unchanged C++ updateFrontEnd path, require actual TunerReady, then serialize the source-reachable frontend state and goLive covariance seed before executing the 601-sample word"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("prehistory detached from complete SEA3")
    for key in (
        "same_driver_field_as_word",
        "stationary_fixed_history_correlation_extension",
        "prehistory_inside_declared_acceleration_cap",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    if d.get("word_join", {}).get("same_t0_value") is not True:
        f.append("same-driver prehistory does not join the 601-sample word at t=0")
    if int(d.get("prehistory_samples", 0)) != PRE_SAMPLES:
        f.append("prehistory sample count changed")
    if not math.isclose(float(d.get("dt_s", 0.0)), DT, rel_tol=0.0, abs_tol=1e-15):
        f.append("prehistory sample period changed")
    if float(d.get("prehistory_duration_s", 0.0)) <= float(d.get("wpe_usable_time_floor_s", math.inf)):
        f.append("prehistory does not even cover the WPE usable floor")
    if float(d.get("prehistory_duration_s", 0.0)) <= float(d.get("outer_online_tune_warmup_s", math.inf)):
        f.append("prehistory does not cover outer tuner warmup")
    for key in (
        "actual_cpp_frontend_TunerReady_verified_here",
        "trajectory_replay_used",
        "finite_harmonic_source_used",
        "independent_startup_signal_used",
        "independent_sample_boxes_used",
        "P3_changed",
        "P4_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} must remain false at prehistory stage")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "prehistory_duration_s": d["prehistory_duration_s"],
        "prehistory_samples": d["prehistory_samples"],
        "max_abs_acceleration_mps2": d["max_abs_prehistory_non_grav_acceleration_mps2"],
        "word_join": d["word_join"],
        "next_obligation": d["next_obligation"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
