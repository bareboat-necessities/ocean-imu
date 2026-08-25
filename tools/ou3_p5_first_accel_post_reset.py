#!/usr/bin/env python3
"""Close the first accepted/rejected accelerometer injection/reset prefix for P5.

The exact-source first-accelerometer stage certifies every deployed correction
at sample zero by d_max < 6 rad.  This producer composes that correction with
the pre-update physical attitude without assuming a favorable correction
sign.  The deployed quaternion map for a finite correction d is an SO(3)
rotation of geodesic angle |d| (the small series branch is evaluated by the
separate deployed primitive and lies far below the present maximum).  Hence
triangle inequality on SO(3) gives

    theta_plus <= theta_minus + |d|.

The pre-update Cayley bound q_minus gives
`theta_minus=2 atan(q_minus/2)`.  If the sum remains below pi, conversion back
to Cayley gives a finite source-complete post-update chart.  A rejected update
is the identity and is automatically contained.

The covariance reset source uses G=I+0.5[d]_x.  Its exact singular-value bound
is sqrt(1+|d|^2/4), so this stage also records a finite reset amplification for
the next full-matrix child propagation.  It does not replace that propagation,
promote the complete word, or set N_H_words.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_deployed_quaternion_cayley_cell as QCOMP

DEFAULT_DOMAIN = FIRST.DEFAULT_DOMAIN
SCHEMA = 1
Q_CHART_TARGET = 8.0


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2) -> dict:
    domain_path = Path(domain_path).resolve()
    first = FIRST.build(domain_path, source_pieces=source_pieces)
    primitive = QCOMP.build(domain_path)
    failures = [f"first-accel: {x}" for x in FIRST.validate(first)]
    failures += [f"deployed-quaternion: {x}" for x in QCOMP.validate(primitive)]
    if first.get("P5_FIRST_ACCEL_EXACT_SOURCE_CERTIFICATE") != "PASS":
        failures.append("exact-source first accelerometer did not pass")

    qminus = float(first["post_prediction_full_cayley_norm_upper"])
    dmax = float(first["max_first_accelerometer_correction_norm_upper_rad"])
    if not (0.0 <= qminus < math.inf and 0.0 <= dmax < math.pi):
        failures.append("invalid pre/reset angle input")
        theta_minus = math.inf
        theta_plus = math.inf
        qplus = math.inf
    else:
        theta_minus = up(2.0 * math.atan(0.5 * qminus))
        theta_plus = up(theta_minus + dmax)
        if theta_plus >= math.pi:
            qplus = math.inf
            failures.append("first accelerometer source family can reach Cayley antipode")
        else:
            qplus = up(2.0 * math.tan(0.5 * theta_plus))

    reset_op = up(math.sqrt(up(1.0 + up(0.25 * dmax * dmax)))) if math.isfinite(dmax) else math.inf
    reset_cov = up(reset_op * reset_op) if math.isfinite(reset_op) else math.inf
    q8_safe = math.isfinite(qplus) and qplus < Q_CHART_TARGET
    if not q8_safe:
        failures.append("first accelerometer post-reset chart is not inside q<8")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_FIRST_ACCEL_DEPLOYED_INJECTION_RESET_PREFIX",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "first_accel_exact_source_status": first["P5_FIRST_ACCEL_EXACT_SOURCE_CERTIFICATE"],
        "pre_accel_cayley_norm_upper": qminus,
        "pre_accel_geodesic_angle_upper_rad": theta_minus,
        "accepted_accel_correction_norm_upper_rad": dmax,
        "rejected_accel_correction_is_identity": True,
        "deployed_correction_limit_rad": float(first["deployed_correction_limit_rad"]),
        "deployed_correction_limit_increased": False,
        "deployed_quaternion_primitive_status": primitive["P5_DEPLOYED_QUATERNION_CAYLEY_CELL_PRIMITIVE"],
        "SO3_triangle_inequality_used_for_sign_complete_composition": True,
        "post_accel_geodesic_angle_upper_rad": theta_plus,
        "post_accel_cayley_norm_upper": qplus,
        "q_chart_target": Q_CHART_TARGET,
        "post_accel_complete_branch_family_inside_q8": q8_safe,
        "reset_operator_norm_upper": reset_op,
        "reset_covariance_quadratic_multiplier_upper": reset_cov,
        "reset_is_nonsingular_for_all_first_accel_children": math.isfinite(reset_op),
        "full_matrix_Joseph_reset_children_propagated_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_FIRST_ACCEL_POST_RESET_PREFIX_CERTIFICATE": "PASS" if not failures else "NOT_ESTABLISHED",
        "next_obligation": (
            "PROPAGATE_EXACT_FIRST_ACCEL_JOSEPH_RESET_COVARIANCE_AND_STATE_CHILDREN_TO_SAMPLE_1"
            if not failures else
            "REFINE_FIRST_ACCEL_SIGNED_CORRECTION_COMPOSITION"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "rejected_accel_correction_is_identity",
        "SO3_triangle_inequality_used_for_sign_complete_composition",
        "post_accel_complete_branch_family_inside_q8",
        "reset_is_nonsingular_for_all_first_accel_children",
    ):
        if d.get(k) is not True:
            failures.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "full_matrix_Joseph_reset_children_propagated_here", "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            failures.append(f"{k} is not false")
    if d.get("first_accel_exact_source_status") != "PASS":
        failures.append("first exact-source correction prerequisite missing")
    if d.get("deployed_quaternion_primitive_status") != "PASS":
        failures.append("deployed quaternion primitive prerequisite missing")
    q = d.get("post_accel_cayley_norm_upper")
    if not isinstance(q, (int, float)) or not math.isfinite(float(q)) or float(q) >= 8.0:
        failures.append("post first-accel q bound is not finite below 8")
    if d.get("P5_FIRST_ACCEL_POST_RESET_PREFIX_CERTIFICATE") != "PASS" and not failures:
        failures.append("first accel post-reset stage did not pass")
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
        "status": out["P5_FIRST_ACCEL_POST_RESET_PREFIX_CERTIFICATE"],
        "q_before": out["pre_accel_cayley_norm_upper"],
        "d_max": out["accepted_accel_correction_norm_upper_rad"],
        "theta_after": out["post_accel_geodesic_angle_upper_rad"],
        "q_after": out["post_accel_cayley_norm_upper"],
        "reset_op": out["reset_operator_norm_upper"],
        "reset_cov_multiplier": out["reset_covariance_quadratic_multiplier_upper"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
