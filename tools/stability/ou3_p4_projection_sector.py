#!/usr/bin/env python3
"""Canonical adapter exposing the radial projection sector as a P4 prerequisite."""
from __future__ import annotations

import json

import ou3_projection_sector as BASE
import ou3_full_process_ucc as PROCESS


def build() -> dict:
    radius = PROCESS._constants()["accel_bias_tau_s"]  # parity import guard; not the projection radius
    del radius
    report = BASE.build_report(0.4)
    return {
        "qualification": "OU3_P4_PROJECTION_SECTOR_ADAPTER_V1",
        "base_lemma": report["lemma"],
        "global_joint_sector_closed": bool(report["exact_real_operator_sector_closed"]),
        "estimate_ball_invariance_closed": bool(report["estimate_ball_invariant_exact_real"]),
        "saturated_branch_included": bool(report["saturated_branch_included_analytically"]),
        "unsaturated_branch_included": bool(report["unsaturated_branch_included_analytically"]),
        "boundary_clarke_branch_included": bool(report["boundary_clarke_branch_included_analytically"]),
        "projection_radius_mps2": 0.4,
        "finite_precision_closed_here": False,
        "base_report": report,
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("qualification") != "OU3_P4_PROJECTION_SECTOR_ADAPTER_V1":
        f.append("qualification mismatch")
    for key in (
        "global_joint_sector_closed", "estimate_ball_invariance_closed",
        "saturated_branch_included", "unsaturated_branch_included", "boundary_clarke_branch_included",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    if float(d.get("projection_radius_mps2", 0.0)) != 0.4:
        f.append("projection radius changed")
    if d.get("finite_precision_closed_here") is not False:
        f.append("projection adapter falsely claimed finite precision")
    return f


def main() -> int:
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    print(json.dumps(d, sort_keys=True))
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
