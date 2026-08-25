#!/usr/bin/env python3
"""Numerical/source-semantic audit of the OU-III ungauged H yaw quotient.

The paper's first quotient draft factors the horizontal attitude gauge but leaves
all three body-frame gyro-bias coordinates in the quotient state.  That is not
a detectable gravity-only system under the declared source family.  A source
word with zero body rate, zero a_w, zero translation error and a gyro-bias error
parallel to gravity is admissible.  The bias produces only yaw; gravity is
unchanged by that yaw, the S=0 residual is zero, and no magnetometer correction
is available on the ungauged branch.  Quotienting yaw therefore removes the
attitude motion while the axial gyro-bias error remains exactly unchanged.

This producer turns that structural fact into a numerical word certificate.  It
also checks the covariance side: the shipping H-mode gyro bias is a random walk,
so its unobserved axial covariance receives a positive process increment on
every word.  Retaining that coordinate in a covariance/information quotient
metric cannot yield a uniformly compact recurrent source family by allowing its
information weight to decay instead of contracting the physical error.

The result is intentionally a fail-closed theorem audit.  The correct gravity-
only theorem must remove the neutral axial-bias/yaw zero dynamics from the
strictly contracting quotient (or carry the axial component as a bounded input)
before a numerical nonlinear quotient funnel can be promoted.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_implementation_proof_manifest as MANIFEST
import ou3_p5_heading_handoff_contract as HEADING

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("quotient-word audit domain must not be trajectory fitted")

    manifest = MANIFEST.build()
    heading = HEADING.build(domain_path)
    process = PROCESS.build()
    prereq = [f"manifest: {x}" for x in MANIFEST.validate(manifest)]
    prereq += [f"heading: {x}" for x in HEADING.validate(heading)]
    prereq += [f"process: {x}" for x in PROCESS.validate(process)]
    if prereq:
        return {
            "schema": SCHEMA,
            "qualification": "OU3_P5_UNGAUGED_TIMEOUT_YAW_QUOTIENT_WORD_AUDIT",
            "source_generated_not_trajectory_fit": True,
            "source_replay_used": False,
            "failures": prereq,
            "P5_UNGAUGED_TIMEOUT_QUOTIENT_WORD_CERTIFICATE": "NOT_EVALUATED",
        }

    ungauged = heading["ungauged_timeout_subbranch"]
    if ungauged["required_route"] != "GRAVITY_ONLY_YAW_QUOTIENT_UNTIL_MAGNETIC_GAUGE_HYBRID_EVENT":
        raise RuntimeError("heading contract no longer routes ungauged timeout to yaw quotient")

    dt = float(domain["configured_runtime"]["imu_dt_s"])
    horizon = float(domain["normal_live"]["vector_pe_recurrence_window_s"])
    samples = int(math.ceil(horizon / dt))
    word_time_lo = down(samples * dt)
    word_time_hi = up(samples * dt)

    b_parallel = float(domain["startup"]["physical_handoff_coordinate_bounds"]["gyro_bias_error_norm_upper_rad_s"])
    if not b_parallel > 0.0:
        raise RuntimeError("positive startup gyro-bias witness required")

    # Exact zero-dynamics witness in H mode.  With omega=0 and b_parallel along
    # body/world down, delta_theta_dot=+b creates yaw only.  The paper's quotient
    # removes that yaw representative.  Gravity and S measurements have exactly
    # zero residual on the selected state, and there is no magnetic correction.
    yaw_created_lo = down(b_parallel * word_time_lo)
    yaw_created_hi = up(b_parallel * word_time_hi)
    b_out_lo = down(b_parallel)
    b_out_hi = up(b_parallel)
    physical_bias_ratio_lo = down(b_out_lo / up(b_parallel))
    physical_bias_ratio_hi = up(b_out_hi / down(b_parallel))

    qbg = float(process["source_constants"]["gyro_bias_rw_variance_density"])
    if not qbg > 0.0:
        raise RuntimeError("gyro-bias process density lost strict positivity")
    covariance_growth_lo = down(qbg * word_time_lo)
    covariance_growth_hi = up(qbg * word_time_hi)

    # Repeating the same admissible word leaves the deterministic axial bias
    # unchanged for every N while its covariance grows at least linearly.
    repeated_words_example = 1000
    repeated_time_lo = down(repeated_words_example * word_time_lo)
    repeated_cov_growth_lo = down(qbg * repeated_time_lo)

    zero_dynamics_valid = bool(
        samples >= 1
        and yaw_created_lo > 0.0
        and b_out_lo > 0.0
        and physical_bias_ratio_lo <= 1.0 <= physical_bias_ratio_hi
        and covariance_growth_lo > 0.0
        and ungauged["source_timeout_requires_north_ready"] is False
        and ungauged["full_heading_cayley_bound_available"] is False
    )

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_UNGAUGED_TIMEOUT_YAW_QUOTIENT_WORD_AUDIT",
        "claim": "EXACT_GRAVITY_ONLY_AXIAL_GYRO_BIAS_ZERO_DYNAMICS_WITNESS",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "mode": "H",
        "shipping_dimension": 18,
        "word_horizon_s": horizon,
        "configured_dt_s": dt,
        "word_samples": samples,
        "ungauged_source_route": ungauged["required_route"],
        "paper_yaw_only_quotient_action_leaves_body_bias_coordinates_unchanged": True,
        "witness": {
            "body_rate_rad_s": 0.0,
            "a_w_error_norm_mps2": 0.0,
            "translation_error_norm": 0.0,
            "S_error_norm_m_s": 0.0,
            "initial_attitude_quotient_error": 0.0,
            "gyro_bias_parallel_to_gravity_rad_s": b_parallel,
            "magnetometer_correction_available": False,
            "gravity_residual_from_created_yaw": 0.0,
            "S_zero_residual": 0.0,
            "created_yaw_before_quotient_interval_rad": [yaw_created_lo, yaw_created_hi],
            "post_word_quotient_attitude_error": 0.0,
            "post_word_parallel_gyro_bias_interval_rad_s": [b_out_lo, b_out_hi],
            "parallel_bias_physical_contraction_ratio_interval": [physical_bias_ratio_lo, physical_bias_ratio_hi],
            "zero_dynamics_source_word_valid": zero_dynamics_valid,
        },
        "covariance_obstruction": {
            "gyro_bias_rw_variance_density": qbg,
            "one_word_unobserved_parallel_bias_covariance_growth_interval": [covariance_growth_lo, covariance_growth_hi],
            "repeated_words_example": repeated_words_example,
            "repeated_unobserved_covariance_growth_lower": repeated_cov_growth_lo,
            "uniform_compact_information_metric_possible_if_parallel_bias_retained": False,
            "reason": (
                "the physical parallel-bias error is not corrected while its covariance receives positive process noise; shrinking its information weight is not physical error contraction"
            ),
        },
        "strict_lambda_less_than_one_possible_on_yaw_only_quotient": False,
        "P5_YAW_ONLY_QUOTIENT_OBSTRUCTION_IDENTIFIED": "PASS" if zero_dynamics_valid else "FAIL",
        "P5_UNGAUGED_TIMEOUT_QUOTIENT_WORD_CERTIFICATE": "NOT_ESTABLISHED",
        "required_quotient_correction": {
            "remove_horizontal_attitude_gauge": True,
            "do_not_require_strict_contraction_of_instantaneous_gravity_parallel_gyro_bias": True,
            "parallel_gyro_bias_role": "neutral/unobservable direction in the zero-motion gravity-only word; carry as quotient direction or bounded input until source excitation or magnetic gauge makes it observable",
            "remaining_numerical_obligation": (
                "construct the observable gravity-only source word on the reduced/time-varying tangent subsystem, prove exact nonlinear quotient prefix safety, and compose the magnetic gauge jump into a full-heading P5 node"
            ),
        },
        "paper_theorem_requires_revision": zero_dynamics_valid,
        "next_obligation": (
            "replace the yaw-only quotient theorem by an observable/detectable quotient that does not penalize the axial gyro-bias zero dynamics, then certify its finite nonlinear source-word funnel"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("quotient audit is not source generated")
    if d.get("source_replay_used") is not False:
        failures.append("quotient audit uses replay")
    if d.get("filter_changed") is not False:
        failures.append("quotient audit changes filter")
    if d.get("paper_yaw_only_quotient_action_leaves_body_bias_coordinates_unchanged") is not True:
        failures.append("paper quotient semantics not bound")
    w = d.get("witness", {})
    if w.get("zero_dynamics_source_word_valid") is not True:
        failures.append("gravity-only axial-bias zero-dynamics witness did not validate")
    ratio = w.get("parallel_bias_physical_contraction_ratio_interval")
    if not (isinstance(ratio, list) and len(ratio) == 2 and float(ratio[0]) <= 1.0 <= float(ratio[1])):
        failures.append("parallel gyro-bias error was not certified unchanged")
    c = d.get("covariance_obstruction", {})
    growth = c.get("one_word_unobserved_parallel_bias_covariance_growth_interval")
    if not (isinstance(growth, list) and len(growth) == 2 and float(growth[0]) > 0.0):
        failures.append("unobserved parallel-bias covariance growth is not strict")
    if c.get("uniform_compact_information_metric_possible_if_parallel_bias_retained") is not False:
        failures.append("audit incorrectly permits compact information metric with retained axial bias")
    if d.get("strict_lambda_less_than_one_possible_on_yaw_only_quotient") is not False:
        failures.append("yaw-only quotient incorrectly claims strict contraction")
    if d.get("P5_YAW_ONLY_QUOTIENT_OBSTRUCTION_IDENTIFIED") != "PASS":
        failures.append("yaw-only quotient obstruction was not identified")
    if d.get("P5_UNGAUGED_TIMEOUT_QUOTIENT_WORD_CERTIFICATE") != "NOT_ESTABLISHED":
        failures.append("ungauged quotient was promoted despite exact zero dynamics")
    if d.get("paper_theorem_requires_revision") is not True:
        failures.append("paper theorem revision not required despite validated obstruction")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    failures = validate(out)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out.get("P5_UNGAUGED_TIMEOUT_QUOTIENT_WORD_CERTIFICATE"),
        "obstruction": out.get("P5_YAW_ONLY_QUOTIENT_OBSTRUCTION_IDENTIFIED"),
        "witness": out.get("witness"),
        "covariance_obstruction": out.get("covariance_obstruction"),
        "next_obligation": out.get("next_obligation"),
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
