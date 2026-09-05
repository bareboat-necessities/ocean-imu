#!/usr/bin/env python3
"""Fail-closed executor boundary from complete SEA3 into canonical P3.

The literal H18/A21 Riccati machinery and a connected typed execution kernel
exist downstream of this boundary.  This module is the only consumer path by
which a materialized 3 s SEA3 window may reach them.  It accepts no replay, raw
sample list, independently supplied F/Q/R_S schedule, or precomputed covariance
floor increment.  The input must first pass
``ou3_sea3_hard_finite_window_source.validate_artifact``.

The SEA0 provider is still open, so canonical execution remains fail-closed.
The typed kernel is already implemented and tested; once the provider closes,
the remaining executor work is strict deserialization of the trusted witness
into that kernel followed by full endpoint LDLT checks.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import ou3_sea3_aw_covariance_floor as FLOOR
import ou3_sea3_complete_source as SOURCE
import ou3_sea3_complete_window_execution_kernel as KERNEL
import ou3_sea3_frontend_state_step as FRONTEND
import ou3_sea3_full_normal_live_word as WORD
import ou3_sea3_hard_finite_window_source as SEA0
import ou3_sea3_shipping_prediction_primitives as PRED

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_SEA3_COMPLETE_WINDOW_EXECUTOR_V2"


def validate_window_artifact(d: dict[str, Any]) -> list[str]:
    """Delegate source validity to SEA0 and enforce the exact executor payload."""
    f = list(SEA0.validate_artifact(d))
    transitions = d.get("transitions")
    if isinstance(transitions, list):
        for k, sample in enumerate(transitions):
            if not isinstance(sample, dict):
                continue
            physical = sample.get("joint_physical_output")
            events = sample.get("source_events")
            if isinstance(physical, dict):
                for key in (
                    "gyro_measurement_interval",
                    "omega_body_corrected_interval",
                    "specific_force_body_interval",
                    "f_cog_body_interval",
                    "R_wb_interval",
                ):
                    if key not in physical:
                        f.append(f"sample {k} missing executor physical field {key}")
            if isinstance(events, dict):
                for key in (
                    "magnetometer_events_after_imu",
                    "aw_covariance_floor_requested",
                    "S_zero_due",
                ):
                    if key not in events:
                        f.append(f"sample {k} missing executor event field {key}")
                if "aw_covariance_floor_increment" in events:
                    f.append(
                        f"sample {k} illegally serializes covariance-dependent a_w floor increment"
                    )
    if not isinstance(d.get("front_end_entry"), dict):
        f.append("missing provider-certified front_end_entry")
    if not isinstance(d.get("live_covariance_seed"), dict):
        f.append("missing provider-certified live_covariance_seed")
    return list(dict.fromkeys(f))


def execute_verified_window(d: dict[str, Any]) -> dict[str, Any]:
    """Execute a provider-certified complete window through the trusted kernel."""
    failures = validate_window_artifact(d)
    if failures:
        raise ValueError("SEA3 window rejected before P3 execution: " + "; ".join(failures))
    if not SEA0.PROVIDER_IMPLEMENTATION_CLOSED:
        raise RuntimeError("SEA0 provider gate is open; H18/A21 execution forbidden")

    # The numerical typed kernel is already implemented.  The only remaining
    # path here is a strict parser for provider-owned interval/front-end/P0
    # witnesses.  Keep this hard failure until that parser exists; accepting a
    # structurally plausible dict would recreate the arbitrary-source problem.
    raise NotImplementedError(
        "provider closed but canonical provider-artifact deserialization into the trusted typed kernel is not yet implemented"
    )


def build_status(domain_path: Path = DEFAULT_DOMAIN) -> dict[str, Any]:
    path = Path(domain_path).resolve()
    source = SOURCE.build(path)
    source_failures = SOURCE.validate(source)
    word = WORD.build(path)
    word_failures = WORD.validate(word)
    pred = PRED.build()
    pred_failures = PRED.validate(pred)
    frontend = FRONTEND.build()
    frontend_failures = FRONTEND.validate(frontend)
    floor = FLOOR.build()
    floor_failures = FLOOR.validate(floor)
    kernel = KERNEL.build(path)
    kernel_failures = KERNEL.validate(kernel)
    sea0 = SEA0.build(path)
    sea0_failures = SEA0.validate_status(sea0)
    bad = {
        "source": source_failures,
        "word": word_failures,
        "prediction": pred_failures,
        "frontend": frontend_failures,
        "aw_floor": floor_failures,
        "typed_kernel": kernel_failures,
        "sea0": sea0_failures,
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"complete-window executor prerequisites failed: {bad}")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": SEA0.CANONICAL_SOURCE,
        "source_generator": False,
        "trajectory_replay_used": False,
        "raw_sample_array_accepted": False,
        "independent_F_Q_schedule_accepted": False,
        "independent_RS_schedule_accepted": False,
        "independent_frontend_state_accepted": False,
        "precomputed_aw_floor_increment_accepted": False,
        "only_canonical_SEA0_provider_artifact_accepted": True,
        "provider_implementation_closed": SEA0.PROVIDER_IMPLEMENTATION_CLOSED,
        "provider_window_artifact_accepted": False,
        "typed_execution_kernel_ready": kernel["typed_execution_kernel_ready"],
        "covariance_dependent_aw_floor_enclosure_ready": floor[
            "positive_part_outer_enclosure_closed_in_real_arithmetic"
        ],
        "raw_gyro_and_corrected_rate_remain_distinct": kernel[
            "raw_gyro_and_bias_corrected_rate_are_distinct_same_witness_coordinates"
        ],
        "SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED": False,
        "FULL_H18_WORD_EXECUTED": False,
        "FULL_A21_WORD_EXECUTED": False,
        "FULL_H18_A21_LDLT_CLOSED": False,
        "P3_PROMOTED": False,
        "literal_word_shipping_parity_pass": word["shipping_event_order_parity_pass"],
        "frontend_shipping_parity_pass": frontend["shipping_source_parity_pass"],
        "prediction_primitives_ready": pred["validated_matrix_primitives_ready"],
        "window_horizon_s": SEA0.HORIZON_S,
        "sample_period_s": SEA0.DT_S,
        "complete_window_samples": SEA0.SAMPLES,
        "execution_contract": {
            "same_provider_transition_drives_physical_and_frontend": True,
            "raw_gyro_measurement_feeds_private_Mahony": True,
            "bias_corrected_body_rate_feeds_MEKF_prediction": True,
            "frontend_derives_committed_tau_sigma_RS_TS": True,
            "prediction_derives_F_Q_from_same_committed_schedule_and_body_rate": True,
            "all_due_S_updates_use_actual_applied_per_axis_RS": True,
            "all_valid_accelerometer_updates_execute": True,
            "asynchronous_provider_PE_events_execute": True,
            "provider_carries_floor_request_not_floor_increment": True,
            "floor_increment_computed_from_current_mode_covariance": True,
            "same_window_executes_H18_and_A21": True,
            "endpoint_gate": "Omega_W - delta P_W >= 0 by full validated LDLT",
            "delta_lower_required": 1.0e-18,
        },
        "next_obligation": (
            "close the canonical SEA0 hard finite-window provider, then strictly deserialize its same-history "
            "601-sample witness into the already-tested typed execution kernel and run every endpoint H18/A21 LDLT"
        ),
    }


def validate_status(d: dict[str, Any]) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != SEA0.CANONICAL_SOURCE:
        f.append("canonical source mismatch")
    for key in (
        "only_canonical_SEA0_provider_artifact_accepted",
        "literal_word_shipping_parity_pass",
        "frontend_shipping_parity_pass",
        "prediction_primitives_ready",
        "typed_execution_kernel_ready",
        "covariance_dependent_aw_floor_enclosure_ready",
        "raw_gyro_and_corrected_rate_remain_distinct",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_generator",
        "trajectory_replay_used",
        "raw_sample_array_accepted",
        "independent_F_Q_schedule_accepted",
        "independent_RS_schedule_accepted",
        "independent_frontend_state_accepted",
        "precomputed_aw_floor_increment_accepted",
        "provider_window_artifact_accepted",
        "SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED",
        "FULL_H18_WORD_EXECUTED",
        "FULL_A21_WORD_EXECUTED",
        "FULL_H18_A21_LDLT_CLOSED",
        "P3_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} must remain false before provider/executor closure")
    if d.get("provider_implementation_closed") is not SEA0.PROVIDER_IMPLEMENTATION_CLOSED:
        f.append("provider gate mismatch")
    c = d.get("execution_contract", {})
    for key in (
        "same_provider_transition_drives_physical_and_frontend",
        "raw_gyro_measurement_feeds_private_Mahony",
        "bias_corrected_body_rate_feeds_MEKF_prediction",
        "frontend_derives_committed_tau_sigma_RS_TS",
        "prediction_derives_F_Q_from_same_committed_schedule_and_body_rate",
        "all_due_S_updates_use_actual_applied_per_axis_RS",
        "all_valid_accelerometer_updates_execute",
        "asynchronous_provider_PE_events_execute",
        "provider_carries_floor_request_not_floor_increment",
        "floor_increment_computed_from_current_mode_covariance",
        "same_window_executes_H18_and_A21",
    ):
        if c.get(key) is not True:
            f.append(f"execution contract lost {key}")
    if float(c.get("delta_lower_required", 0.0)) != 1.0e-18:
        f.append("P3 useful gate changed")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--sea0-window", type=Path, default=None)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    if args.sea0_window is None:
        d = build_status(args.domain)
        failures = validate_status(d)
    else:
        candidate = json.loads(args.sea0_window.read_text(encoding="utf-8"))
        artifact_failures = validate_window_artifact(candidate)
        d = build_status(args.domain)
        d["candidate_window_validation_failures"] = artifact_failures
        d["provider_window_artifact_accepted"] = not artifact_failures
        failures = validate_status(d)
        if not artifact_failures:
            d["execution_result"] = execute_verified_window(candidate)
        else:
            failures.extend(artifact_failures)
    d["validation_pass"] = not failures
    d["validation_failures"] = list(dict.fromkeys(failures))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "provider_closed": d["provider_implementation_closed"],
        "typed_kernel_ready": d["typed_execution_kernel_ready"],
        "family_materialized": d["SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"],
        "H18_executed": d["FULL_H18_WORD_EXECUTED"],
        "A21_executed": d["FULL_A21_WORD_EXECUTED"],
        "P3": d["P3_PROMOTED"],
        "failures": d["validation_failures"],
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
