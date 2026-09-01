#!/usr/bin/env python3
"""Audit the first unclosed 30 deg P4 complete-word candidate cell.

The complete-word candidate backend currently fails closed before signed Cayley
composition when the first H accelerometer correction exceeds the deployed
validated [0,6] rad correction range.  This producer instruments the existing
backend without changing any enclosure or filter operation.  It records the
measurement residual, attitude-gain block, correction norm, innovation box and
the interval operator norm assigned to J_aw=R_wb before the failure.

The exact shipping identity ||R_wb||_2=1 is reported beside the entrywise-box
operator bound only as a diagnostic gap.  This producer does not replace the
box, assume isotropic later-word covariance, or promote P4.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_candidate_full_word as C

DEFAULT_DOMAIN = C.DEFAULT_DOMAIN
SCHEMA = 1


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _op2_upper(A) -> float:
    if not A or not A[0]:
        return 0.0
    rows = []
    for row in A:
        s = 0.0
        for x in row:
            s = up(s + x.abs_upper())
        rows.append(s)
    cols = []
    for j in range(len(A[0])):
        s = 0.0
        for i in range(len(A)):
            s = up(s + A[i][j].abs_upper())
        cols.append(s)
    return up(math.sqrt(up(max(rows) * max(cols))))


def _nonzero(x) -> bool:
    return x.lo != 0.0 or x.hi != 0.0


def _operation(Hm) -> str:
    if any(_nonzero(Hm[i][j]) for i in range(3) for j in C.H.SS):
        return "S"
    if any(_nonzero(Hm[i][j]) for i in range(3) for j in C.H.AW):
        return "accelerometer"
    return "magnetometer_or_other"


def _ser_matrix(A):
    return [[x.as_list() for x in row] for row in A]


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    records = []
    original = C.H._measurement_cell

    def tracked(Pm, Hm, R, r):
        cell = original(Pm, Hm, R, r)
        op = _operation(Hm)
        Ktheta = [list(row) for row in cell["K"][0:3]]
        Jaw = [[Hm[i][j] for j in C.H.AW] for i in range(3)]
        records.append({
            "operation": op,
            "residual_norm_upper": C.H._norm_upper(r),
            "residual_components": [x.as_list() for x in r],
            "attitude_correction_norm_upper_rad": C.H._norm_upper(cell["dx"][0:3]),
            "attitude_gain_operator_norm_upper": _op2_upper(Ktheta),
            "J_aw_entrywise_box_operator_norm_upper": _op2_upper(Jaw) if op == "accelerometer" else None,
            "J_aw_shipping_exact_operator_norm": 1.0 if op == "accelerometer" else None,
            "J_aw_box_inflation_factor_upper": _op2_upper(Jaw) if op == "accelerometer" else None,
            "inverse_backend": cell["inverse_backend"],
            "innovation_covariance": _ser_matrix(cell["S"]),
            "Ktheta": _ser_matrix(Ktheta),
        })
        return cell

    C.H._measurement_cell = tracked
    try:
        parent = C.build(path, max_samples=1, candidate_index=0, cover_factor=1.55)
    finally:
        C.H._measurement_cell = original

    row = parent["candidate_rows"][0]
    ff = row.get("first_unclosed_cell")
    accel = next((x for x in records if x["operation"] == "accelerometer"), None)
    failures = []
    if ff is None:
        failures.append("30deg one-sample audit no longer reproduces an unclosed cell")
    else:
        result = ff.get("result", {})
        fail = result.get("first_failure") or {}
        if int(fail.get("sample", -1)) != 0:
            failures.append("first unclosed 30deg cell moved away from sample 0")
        if "deployed correction norm upper outside validated range" not in str(fail.get("reason", "")):
            failures.append("first unclosed 30deg cell is no longer the correction-range obstruction")
    if accel is None:
        failures.append("accelerometer measurement cell was not observed")
    else:
        if not float(accel["attitude_correction_norm_upper_rad"]) > 6.0:
            failures.append("instrumented accelerometer correction does not reproduce >6 rad obstruction")
        if not float(accel["J_aw_entrywise_box_operator_norm_upper"]) >= 1.0:
            failures.append("invalid J_aw interval operator bound")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_30DEG_FIRST_COMPLETE_WORD_OBSTRUCTION_AUDIT",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "parent_uses_old_q8_chart": bool(parent["old_q8_chart_used"]),
        "candidate_angle_deg": float(row["angle_deg"]),
        "candidate_cover_factor": float(parent["cover_factor"]),
        "candidate_cover_q_upper": float(row["cover_q_upper"]),
        "operation_matched_outer_q_upper": float(parent["operation_matched_outer_q_upper"]),
        "first_unclosed_cell": ff,
        "measurement_records": records,
        "accelerometer_obstruction": accel,
        "shipping_J_aw_is_orthogonal": True,
        "entrywise_J_aw_box_replacement_performed_here": False,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "P4_30DEG_FIRST_OBSTRUCTION_AUDIT": "PASS" if not failures else "FAIL",
        "next_obligation": (
            "replace the sample0 accelerometer entrywise J_aw rotation box by a source-faithful orthogonal/gauged covariance-residual calculation and rerun the same 30deg cell; do not change the candidate angle or deployed correction range"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in ("source_generated_not_trajectory_fit", "shipping_J_aw_is_orthogonal"):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "parent_uses_old_q8_chart",
        "entrywise_J_aw_box_replacement_performed_here",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE", "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("candidate_angle_deg", -1.0)) != 30.0:
        f.append("audit candidate is not 30deg")
    a = d.get("accelerometer_obstruction") or {}
    if not float(a.get("attitude_correction_norm_upper_rad", 0.0)) > 6.0:
        f.append("accelerometer obstruction is not above deployed correction range")
    if d.get("P4_30DEG_FIRST_OBSTRUCTION_AUDIT") == "PASS" and f:
        f.append("PASS carries validation failures")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain.resolve())
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    x = d.get("accelerometer_obstruction") or {}
    print(json.dumps({
        "status": d["P4_30DEG_FIRST_OBSTRUCTION_AUDIT"],
        "cover_q": d["candidate_cover_q_upper"],
        "residual_norm": x.get("residual_norm_upper"),
        "Ktheta_op": x.get("attitude_gain_operator_norm_upper"),
        "correction_norm": x.get("attitude_correction_norm_upper_rad"),
        "Jaw_box_op": x.get("J_aw_entrywise_box_operator_norm_upper"),
        "Jaw_exact_op": x.get("J_aw_shipping_exact_operator_norm"),
        "inverse_backend": x.get("inverse_backend"),
        "validation_failures": vf,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
