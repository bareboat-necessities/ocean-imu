#!/usr/bin/env python3
"""Canonical fail-closed executor from complete SEA3 into P3.

The numerical H18/A21 word kernel and the strict provider-artifact codec are
complete downstream of this boundary.  This module is the only consumer path by
which a materialized 3 s SEA3 family may reach P3.  It accepts no replay, raw
sample array, independently supplied F/Q/R_S schedule, arbitrary P0, or
precomputed covariance-floor increment.

Canonical execution has exactly this order:

  1. the code-owned SEA0 provider validates the complete same-history witness;
  2. the strict codec deserializes the provider-owned z^t entry, source-reachable
     H18/A21 Live covariance seeds, and all 601 transition payloads;
  3. the connected typed kernel executes every retained front-end branch through
     every shipping prediction/floor/S/accel/mag event in both modes;
  4. every endpoint branch is tested by the full validated matrix inequality

         Omega_W - delta P_W >= 0,  delta = 1e-18,

     using interval LDLT in H18 and A21.

The SEA0 provider is currently still open, so this complete execution path is
intentionally unreachable in canonical CI.  No structural JSON or diagnostic
word can bypass that provider gate.
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
import ou3_sea3_window_artifact_codec as CODEC

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 3
QUALIFICATION = "OU3_SEA3_COMPLETE_WINDOW_EXECUTOR_V3"
DELTA = 1.0e-18


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


def _endpoint_certificate(branch: KERNEL.ExecutionBranch) -> dict[str, Any]:
    H = WORD.certify_literal_endpoint(branch.H, delta=DELTA)
    A = WORD.certify_literal_endpoint(branch.A, delta=DELTA)
    return {
        "source_cell_id": branch.source_cell_id,
        "H18": H,
        "A21": A,
        "H18_complete_601_sample_word": (
            branch.H.imu_samples == SEA0.SAMPLES
            and branch.H.accel_updates == SEA0.SAMPLES
            and branch.H.riccati.predictions == SEA0.SAMPLES
        ),
        "A21_complete_601_sample_word": (
            branch.A.imu_samples == SEA0.SAMPLES
            and branch.A.accel_updates == SEA0.SAMPLES
            and branch.A.riccati.predictions == SEA0.SAMPLES
        ),
        "both_full_matrix_LDLT_closed": bool(H.get("pass") and A.get("pass")),
    }


def execute_verified_window(
    d: dict[str, Any],
    domain_path: Path = DEFAULT_DOMAIN,
) -> dict[str, Any]:
    """Execute one provider-certified complete SEA3 window family.

    There is no fallback path.  Provider acceptance establishes source
    admissibility/reachability; the codec only transports that witness; the
    kernel only executes it; the endpoint certificate alone decides P3.
    """
    failures = validate_window_artifact(d)
    if failures:
        raise ValueError("SEA3 window rejected before P3 execution: " + "; ".join(failures))
    if not SEA0.PROVIDER_IMPLEMENTATION_CLOSED:
        raise RuntimeError("SEA0 provider gate is open; H18/A21 execution forbidden")

    parsed = CODEC.parse_window_artifact(d)
    branches, meta = KERNEL.execute_typed_window(
        frontend_entry=parsed.frontend_entry,
        P0_H=parsed.P0_H,
        P0_A=parsed.P0_A,
        samples=parsed.samples,
        domain_path=Path(domain_path).resolve(),
    )
    if not branches:
        raise RuntimeError("complete SEA3 execution produced no endpoint source cells")

    certs = [_endpoint_certificate(branch) for branch in branches]
    H_complete = all(x["H18_complete_601_sample_word"] for x in certs)
    A_complete = all(x["A21_complete_601_sample_word"] for x in certs)
    ldlt_closed = all(x["both_full_matrix_LDLT_closed"] for x in certs)
    materialized = (
        len(parsed.samples) == SEA0.SAMPLES
        and int(meta.get("samples_executed", -1)) == SEA0.SAMPLES
        and meta.get("same_word_executed_H18_A21") is True
        and meta.get("favorable_frontend_successor_selected") is False
    )
    p3 = bool(materialized and H_complete and A_complete and ldlt_closed)
    return {
        "canonical_source": SEA0.CANONICAL_SOURCE,
        "delta": DELTA,
        "provider_accepted_before_execution": True,
        "strict_codec_used": True,
        "same_history_window_samples": len(parsed.samples),
        "SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED": materialized,
        "FULL_H18_WORD_EXECUTED": H_complete,
        "FULL_A21_WORD_EXECUTED": A_complete,
        "FULL_H18_A21_LDLT_CLOSED": ldlt_closed,
        "endpoint_branch_count": len(branches),
        "kernel_execution": meta,
        "endpoint_certificates": certs,
        "P3_PROMOTED": p3,
    }


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
    codec = CODEC.build()
    codec_failures = CODEC.validate(codec)
    sea0 = SEA0.build(path)
    sea0_failures = SEA0.validate_status(sea0)
    bad = {
        "source": source_failures,
        "word": word_failures,
        "prediction": pred_failures,
        "frontend": frontend_failures,
        "aw_floor": floor_failures,
        "typed_kernel": kernel_failures,
        "artifact_codec": codec_failures,
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
        "arbitrary_P0_accepted": False,
        "precomputed_aw_floor_increment_accepted": False,
        "only_canonical_SEA0_provider_artifact_accepted": True,
        "provider_implementation_closed": SEA0.PROVIDER_IMPLEMENTATION_CLOSED,
        "provider_window_artifact_accepted": False,
        "strict_artifact_codec_ready": codec["strict_codec_ready"],
        "typed_execution_kernel_ready": kernel["typed_execution_kernel_ready"],
        "covariance_dependent_aw_floor_enclosure_ready": floor[
            "positive_part_outer_enclosure_closed_in_real_arithmetic"
        ],
        "raw_gyro_and_corrected_rate_remain_distinct": kernel[
            "raw_gyro_and_bias_corrected_rate_are_distinct_same_witness_coordinates"
        ],
        "frontend_completed_before_async_mag": kernel[
            "frontend_completed_before_async_mag"
        ],
        "gated_601_sample_execution_path_implemented": True,
        "gated_endpoint_full_matrix_LDLT_path_implemented": True,
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
            "frontend_transition_precedes_after_imu_magnetometer_events": True,
            "asynchronous_provider_PE_events_execute": True,
            "provider_carries_floor_request_not_floor_increment": True,
            "floor_increment_computed_from_current_mode_covariance": True,
            "live_P0_must_be_provider_certified_source_reachable": True,
            "same_window_executes_H18_and_A21": True,
            "every_retained_endpoint_branch_must_pass": True,
            "endpoint_gate": "Omega_W - delta P_W >= 0 by full validated LDLT",
            "delta_lower_required": DELTA,
        },
        "next_obligation": (
            "close the canonical SEA0 hard finite-window provider; the downstream strict codec, connected "
            "601-sample H18/A21 execution loop, covariance-dependent floor, and every-branch endpoint LDLT path are implemented"
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
        "strict_artifact_codec_ready",
        "typed_execution_kernel_ready",
        "covariance_dependent_aw_floor_enclosure_ready",
        "raw_gyro_and_corrected_rate_remain_distinct",
        "frontend_completed_before_async_mag",
        "gated_601_sample_execution_path_implemented",
        "gated_endpoint_full_matrix_LDLT_path_implemented",
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
        "arbitrary_P0_accepted",
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
        "frontend_transition_precedes_after_imu_magnetometer_events",
        "asynchronous_provider_PE_events_execute",
        "provider_carries_floor_request_not_floor_increment",
        "floor_increment_computed_from_current_mode_covariance",
        "live_P0_must_be_provider_certified_source_reachable",
        "same_window_executes_H18_and_A21",
        "every_retained_endpoint_branch_must_pass",
    ):
        if c.get(key) is not True:
            f.append(f"execution contract lost {key}")
    if float(c.get("delta_lower_required", 0.0)) != DELTA:
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
            result = execute_verified_window(candidate, args.domain)
            d["execution_result"] = result
            d["SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"] = result[
                "SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"
            ]
            d["FULL_H18_WORD_EXECUTED"] = result["FULL_H18_WORD_EXECUTED"]
            d["FULL_A21_WORD_EXECUTED"] = result["FULL_A21_WORD_EXECUTED"]
            d["FULL_H18_A21_LDLT_CLOSED"] = result["FULL_H18_A21_LDLT_CLOSED"]
            d["P3_PROMOTED"] = result["P3_PROMOTED"]
            # Status validation above intentionally describes the no-artifact
            # fail-closed state.  A provider-accepted execution is validated by
            # the result's own all-branch full-matrix gates instead.
            failures = [] if result["P3_PROMOTED"] else [
                "provider-accepted complete SEA3 execution did not close every H18/A21 endpoint LDLT"
            ]
        else:
            failures.extend(artifact_failures)
    d["validation_pass"] = not failures
    d["validation_failures"] = list(dict.fromkeys(failures))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "provider_closed": d["provider_implementation_closed"],
        "codec_ready": d["strict_artifact_codec_ready"],
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
