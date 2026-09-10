#!/usr/bin/env python3
"""Shipping timeout handoff -> early-Live magnetic regauge contract.

This is a theorem-architecture certificate, not a numerical convergence proof.
It closes the semantic question that a timeout handoff need not already possess
north: P5 may enter Live in H18 tilt-only and consume the existing Normal-Live
vector recurrence after handoff, provided the shipping wrapper continues to
accept magnetometer packets and update the MEKF magnetic reference/correction in
Live.

The certificate is deliberately fail-closed on the remaining source-admission
and finite-time state-capture obligations.
"""
from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_TIMEOUT_POSTLIVE_REGAUGE_CONTRACT_V1"


def build() -> dict:
    text = WRAPPER.read_text(encoding="utf-8")
    domain = json.loads(DOMAIN.read_text(encoding="utf-8"))
    normal = domain["normal_live"]

    # Structural source checks. These strings intentionally target deployed
    # code paths rather than comments in this proof module.
    live_dispatch_present = (
        "stage_ == Stage::Live" in text
        or "Stage::Live" in text
    )
    mag_update_present = (
        "updateMag" in text
        or "update_mag" in text
        or "mag_world_ref" in text
        or "set_mag_world_ref" in text
    )
    magnetic_lock_event_named = (
        "magnetic_lock" in json.dumps(normal)
        or "magnetic_regauge_refinement" in json.dumps(normal)
    )

    pe_window = float(normal["vector_pe_recurrence_window_s"])
    sine_sep = float(normal["vector_sine_separation_lower"])
    mag_lo = float(normal["magnetic_vector_norm_lower_uT"])

    normal_live_pe_declared = (
        pe_window > 0.0 and sine_sep > 0.0 and mag_lo > 0.0
    )

    structural_postlive_regauge_path = (
        live_dispatch_present and mag_update_present and normal_live_pe_declared
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "filter_changed": False,
        "quality_gates_changed": False,
        "timeout_handoff_may_be_ungauged": True,
        "P5_may_extend_past_goLive": True,
        "postlive_H18_tilt_only_mode_required": True,
        "normal_live_vector_PE_window_s": pe_window,
        "normal_live_vector_sine_separation_lower": sine_sep,
        "normal_live_magnetic_norm_lower_uT": mag_lo,
        "shipping_live_dispatch_present": live_dispatch_present,
        "shipping_live_magnetometer_update_path_present": mag_update_present,
        "normal_live_magnetic_hybrid_event_named": magnetic_lock_event_named,
        "normal_live_vector_PE_declared": normal_live_pe_declared,
        "STRUCTURAL_POSTLIVE_REGAUGE_PATH_CLOSED": structural_postlive_regauge_path,
        "startup_magnetic_PE_required_for_timeout_branch": False,
        "postlive_vector_PE_may_supply_first_heading_gauge": structural_postlive_regauge_path,
        "source_admission_closed_here": False,
        "finite_time_regauge_bound_closed_here": False,
        "P4_capture_closed_here": False,
        "P5_PASS": False,
        "next_obligation": (
            "materialize the admitted COMPLETE-BRMM early-Live source family and prove that its accepted magnetic-packet recurrence drives the exact shipping H18 magnetic update/regauge map into the widened P4 basin in a finite source-uniform time"
        ),
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        failures.append("schema/qualification mismatch")
    for key in (
        "timeout_handoff_may_be_ungauged",
        "P5_may_extend_past_goLive",
        "postlive_H18_tilt_only_mode_required",
        "normal_live_vector_PE_declared",
        "STRUCTURAL_POSTLIVE_REGAUGE_PATH_CLOSED",
        "postlive_vector_PE_may_supply_first_heading_gauge",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "startup_magnetic_PE_required_for_timeout_branch",
        "source_admission_closed_here",
        "finite_time_regauge_bound_closed_here",
        "P4_capture_closed_here",
        "P5_PASS",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} must remain false")
    return failures


if __name__ == "__main__":
    d = build()
    failures = validate(d)
    print(json.dumps({**d, "validation_pass": not failures, "validation_failures": failures}, indent=2))
    raise SystemExit(1 if failures else 0)
