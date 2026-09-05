#!/usr/bin/env python3
"""CLI/validation wrapper for the literal 3 s H18/A21 execution engine."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_sea3_full_normal_live_execution as EXEC
import ou3_sea3_full_normal_live_word as WORD
import ou3_sea3_full_word_riccati_backend_tight as TIGHT


def validate_execution(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("canonical_architecture") != "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD":
        failures.append("execution witness detached from canonical architecture")
    for mode in ("H18", "A21"):
        row = d.get(mode, {})
        if row.get("complete_3s_word_executed") is not True:
            failures.append(f"{mode} did not execute a complete 3 s word")
        if int(row.get("predictions", 0)) != int(row.get("samples_executed", -1)):
            failures.append(f"{mode} did not predict every IMU sample")
        if int(row.get("accelerometer_updates", 0)) != int(row.get("samples_executed", -1)):
            failures.append(f"{mode} did not apply every accelerometer update")
        if int(row.get("S_zero_updates", 0)) <= 0:
            failures.append(f"{mode} did not execute pseudo updates")
        if int(row.get("magnetometer_updates", 0)) != 2:
            failures.append(f"{mode} did not execute the two asynchronous PE events")
        if int(row.get("aw_floor_applications", 0)) <= 0:
            failures.append(f"{mode} did not execute state-dependent a_w floor events")
        if row.get("decomposition_identity_preserved_every_sample") is not True:
            failures.append(f"{mode} lost the exact P/Psi/Omega identity")
        if row.get("point_endpoint_may_promote_P3") is not False:
            failures.append(f"{mode} point execution was allowed to promote P3")
    if d.get("both_modes_exact_joint_backend_executed") is not True:
        failures.append("both full-state modes were not executed")
    if d.get("SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED") is not False:
        failures.append("point execution falsely claimed source-family materialization")
    if d.get("FULL_H18_A21_LDLT_CLOSED") is not False:
        failures.append("point execution falsely claimed uniform H/A LDLT closure")
    if d.get("P3_CANONICAL_PASS") is not False:
        failures.append("point execution falsely promoted P3")
    if d.get("no_point_or_reduced_result_may_promote") is not True:
        failures.append("promotion guard is not armed")
    return failures


def _install_tight_backend_and_failure_locator() -> dict[str, int]:
    # Same exact P/Psi/Omega backend and Kalman gain; only redundant validated
    # identities are intersected to stop natural-interval wrapping from erasing
    # the corrective power of recurrent R_S measurements.
    WORD.BACKEND = TIGHT
    EXEC.BACKEND = TIGHT

    counts = {"H": 0, "A": 0}
    orig_imu = WORD.apply_imu_sample
    orig_mag = WORD.apply_magnetometer

    def located_imu(w, *args, **kwargs):
        counts[w.mode] += 1
        try:
            return orig_imu(w, *args, **kwargs)
        except Exception as exc:
            raise RuntimeError(
                f"literal full-word failure mode={w.mode} sample={counts[w.mode]} "
                f"event=imu/prediction-floor-S-accelerometer: {exc}"
            ) from exc

    def located_mag(w, *args, **kwargs):
        try:
            return orig_mag(w, *args, **kwargs)
        except Exception as exc:
            raise RuntimeError(
                f"literal full-word failure mode={w.mode} sample={counts[w.mode]} "
                f"event=magnetometer: {exc}"
            ) from exc

    WORD.apply_imu_sample = located_imu
    WORD.apply_magnetometer = located_mag
    return counts


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=WORD.DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    _install_tight_backend_and_failure_locator()
    manifest = WORD.build(args.domain)
    d = EXEC.build_execution(WORD, float(manifest["word_horizon_s"]))
    d["dependency_tight_joint_backend_used"] = True
    d["R_S_correction_retained_by_Joseph_Schur_intersection"] = True
    failures = validate_execution(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "H18": {k: d["H18"][k] for k in (
            "samples_executed", "S_zero_updates", "magnetometer_updates",
            "aw_floor_applications", "complete_3s_word_executed")},
        "A21": {k: d["A21"][k] for k in (
            "samples_executed", "S_zero_updates", "magnetometer_updates",
            "aw_floor_applications", "complete_3s_word_executed")},
        "source_family_materialized": d["SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"],
        "P3_CANONICAL_PASS": d["P3_CANONICAL_PASS"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
