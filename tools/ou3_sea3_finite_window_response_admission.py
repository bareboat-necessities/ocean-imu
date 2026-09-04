#!/usr/bin/env python3
"""Fail-closed SEA3 finite-window response admission contract.

This module is the mechanical bridge between the physical SEA3 sea/RAO domain
and the already-declared Normal-Live P1 source bounds.  It intentionally does
*not* turn a JONSWAP spectrum, an RMS value, or an ensemble moment into a
sample-path bound.

A physical finite sea window may enter ``Lhat_SEA3`` only when a separate
validated realization producer supplies outward enclosures that cover every
valid IMU sample in the window after the admitted vessel response operator is
applied.  The currently required sea/response-owned P1 quantities are

* non-gravitational CoG acceleration norm, and
* body-rate norm.

The accelerometer specific-force envelope is then the existing P1 consequence
of gravity plus the non-gravitational acceleration cap.  Magnetic PE, hybrid
reset/regauge obligations, and other non-sea source assumptions remain separate
P1 obligations and are not silently discharged here.

This file therefore closes the *admission predicate/interface*, not the global
left inclusion ``L_actual_sea subset Lhat_SEA3``.  That inclusion remains false
until a replay-free oscillator/IQC (or equivalent validated realization)
producer proves that every physical deployment window being claimed by the
theorem can furnish the required evidence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import ou3_sea3_directional_p2_ha_feasibility as SEA3
import ou3_sea3_p1_compatibility as P1COMPAT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
DEFAULT_RESPONSE_DOMAIN = REPO / "tools" / "ou3_sea3_directional_response_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_FINITE_WINDOW_RESPONSE_ADMISSION_V1"


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _load_domain(path: Path = DEFAULT_DOMAIN) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def build_contract(
    domain_path: Path = DEFAULT_DOMAIN,
    response_domain_path: Path = DEFAULT_RESPONSE_DOMAIN,
    repo: Path = REPO,
) -> dict:
    """Build the non-promoting window-admission interface."""
    domain = _load_domain(domain_path)
    normal = domain["normal_live"]
    runtime = domain["configured_runtime"]

    compatibility = P1COMPAT.build(Path(domain_path), Path(response_domain_path))
    compatibility_failures = P1COMPAT.validate(compatibility)
    if compatibility_failures:
        raise RuntimeError(f"SEA3/P1 compatibility contract invalid: {compatibility_failures}")
    if compatibility.get("coupled_SEA3_domain_required") is not True:
        raise RuntimeError("finite-window admission requires the coupled SEA3 sea/RAO domain")

    response = SEA3.directional_response_enclosure(
        Path(repo).resolve(), Path(response_domain_path).resolve()
    )
    box = response["rao_envelope_parameter_box"]

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_filter_operating_domain_changed": False,
        "coupled_SEA3_domain_required": True,
        "independent_cartesian_sea_x_RAO_domain_rejected": True,
        "response_parameter_box": box,
        "response_parameter_box_sha256": _canonical_sha256(box),
        "normal_live_caps": {
            "non_gravitational_cog_acceleration_norm_upper_mps2": float(
                normal["non_gravitational_cog_acceleration_norm_upper_mps2"]
            ),
            "body_rate_norm_upper_deg_s": float(normal["body_rate_norm_upper_deg_s"]),
        },
        "runtime_branch": {
            "imu_lever_arm_enabled": bool(runtime["imu_lever_arm_enabled"]),
            "accelerometer_vibration_guard_proof_branch": runtime[
                "accelerometer_vibration_guard_proof_branch"
            ],
            "accelerometer_update_required_each_valid_imu_sample_after_live_entry": bool(
                normal["accelerometer_update_required_each_valid_imu_sample_after_live_entry"]
            ),
            "accelerometer_rejection_in_normal_live_scope": bool(
                normal["accelerometer_rejection_in_normal_live_scope"]
            ),
        },
        "required_window_evidence": {
            "validated_arithmetic": True,
            "outward_rounded": True,
            "post_rao_response_enclosed": True,
            "all_valid_imu_samples_covered": True,
            "trajectory_replay_used": False,
            "gaussian_spectrum_only": False,
            "rms_or_psd_only": False,
            "response_parameter_box_sha256": _canonical_sha256(box),
            "window_samples": "positive integer",
            "post_rao_cog_acceleration_norm_upper_mps2": "finite nonnegative scalar",
            "body_rate_norm_upper_deg_s": "finite nonnegative scalar",
        },
        "gaussian_spectrum_alone_can_admit_window": False,
        "RMS_or_PSD_moment_alone_can_admit_window": False,
        "finite_window_admission_predicate_defined": True,
        "finite_window_realization_producer_closed": False,
        "finite_window_realization_certificate_closed": False,
        "L_actual_sea_subset_Lhat_SEA3_closed": False,
        "P2_pruning_promoted": False,
        "P3_promoted": False,
        "P4_promoted": False,
        "P5_promoted": False,
        "next_obligation": (
            "build a replay-free oscillator/IQC or equivalent validated finite-window realization producer that supplies this evidence for the claimed coupled JONSWAP-sea/RAO population; only then may the global left inclusion be promoted"
        ),
    }


def evaluate_window(
    evidence: dict,
    domain_path: Path = DEFAULT_DOMAIN,
    response_domain_path: Path = DEFAULT_RESPONSE_DOMAIN,
    repo: Path = REPO,
) -> dict:
    """Evaluate one finite-window realization against the exact admission contract.

    This predicate accepts only validated pathwise/samplewise evidence.  It does
    not attempt to manufacture such evidence from spectral moments.
    """
    contract = build_contract(domain_path, response_domain_path, repo)
    failures: list[str] = []

    required_true = (
        "validated_arithmetic",
        "outward_rounded",
        "post_rao_response_enclosed",
        "all_valid_imu_samples_covered",
    )
    for key in required_true:
        if evidence.get(key) is not True:
            failures.append(f"{key} must be true")

    if evidence.get("trajectory_replay_used") is not False:
        failures.append("trajectory replay cannot establish finite-window admission")
    if evidence.get("gaussian_spectrum_only") is not False:
        failures.append("Gaussian spectrum alone cannot establish a pathwise P1 bound")
    if evidence.get("rms_or_psd_only") is not False:
        failures.append("RMS/PSD moments alone cannot establish a pathwise P1 bound")

    expected_digest = contract["response_parameter_box_sha256"]
    if evidence.get("response_parameter_box_sha256") != expected_digest:
        failures.append("finite-window evidence is not bound to the certified RAO parameter box")

    samples = evidence.get("window_samples")
    if not isinstance(samples, int) or isinstance(samples, bool) or samples <= 0:
        failures.append("window_samples must be a positive integer")

    caps = contract["normal_live_caps"]
    for key, cap_key in (
        (
            "post_rao_cog_acceleration_norm_upper_mps2",
            "non_gravitational_cog_acceleration_norm_upper_mps2",
        ),
        ("body_rate_norm_upper_deg_s", "body_rate_norm_upper_deg_s"),
    ):
        try:
            value = float(evidence.get(key, math.nan))
        except (TypeError, ValueError):
            value = math.nan
        cap = float(caps[cap_key])
        if not (math.isfinite(value) and value >= 0.0):
            failures.append(f"{key} must be finite and nonnegative")
        elif value > cap:
            failures.append(f"{key} exceeds Normal-Live P1 cap {cap}")

    admitted = not failures
    return {
        "qualification": "OU3_SEA3_FINITE_WINDOW_RESPONSE_ADMISSION_DECISION_V1",
        "window_admitted_to_Lhat_SEA3": admitted,
        "decision": "ADMIT" if admitted else "REJECT",
        "validation_failures": list(dict.fromkeys(failures)),
        "global_left_inclusion_promoted_by_this_decision": False,
        "P1_nonsea_obligations_discharged_by_this_decision": False,
        "contract_response_parameter_box_sha256": expected_digest,
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        failures.append("schema/qualification mismatch")
    for key in (
        "coupled_SEA3_domain_required",
        "independent_cartesian_sea_x_RAO_domain_rejected",
        "finite_window_admission_predicate_defined",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used",
        "filter_changed",
        "declared_filter_operating_domain_changed",
        "gaussian_spectrum_alone_can_admit_window",
        "RMS_or_PSD_moment_alone_can_admit_window",
        "finite_window_realization_producer_closed",
        "finite_window_realization_certificate_closed",
        "L_actual_sea_subset_Lhat_SEA3_closed",
        "P2_pruning_promoted",
        "P3_promoted",
        "P4_promoted",
        "P5_promoted",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")

    box = d.get("response_parameter_box", {})
    if d.get("response_parameter_box_sha256") != _canonical_sha256(box):
        failures.append("RAO parameter-box digest mismatch")

    caps = d.get("normal_live_caps", {})
    if float(caps.get("non_gravitational_cog_acceleration_norm_upper_mps2", math.nan)) != 4.0:
        failures.append("Normal-Live acceleration cap changed from 4 m/s^2")
    if float(caps.get("body_rate_norm_upper_deg_s", math.nan)) != 30.0:
        failures.append("Normal-Live body-rate cap changed from 30 deg/s")

    runtime = d.get("runtime_branch", {})
    if runtime.get("imu_lever_arm_enabled") is not False:
        failures.append("finite-window bridge left zero-lever-arm proof branch")
    if runtime.get("accelerometer_vibration_guard_proof_branch") != "dormant_transparent":
        failures.append("finite-window bridge left dormant vibration-guard proof branch")
    if runtime.get("accelerometer_update_required_each_valid_imu_sample_after_live_entry") is not True:
        failures.append("Normal-Live accelerometer-update requirement disappeared")
    if runtime.get("accelerometer_rejection_in_normal_live_scope") is not False:
        failures.append("accelerometer rejection re-entered Normal-Live proof scope")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--response-domain", type=Path, default=DEFAULT_RESPONSE_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    d = build_contract(args.domain, args.response_domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(
        json.dumps(
            {
                "qualification": d["qualification"],
                "response_parameter_box_sha256": d["response_parameter_box_sha256"],
                "P1_caps": d["normal_live_caps"],
                "finite_window_admission_predicate_defined": d[
                    "finite_window_admission_predicate_defined"
                ],
                "left_inclusion_closed": d["L_actual_sea_subset_Lhat_SEA3_closed"],
                "validation_failures": failures,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
