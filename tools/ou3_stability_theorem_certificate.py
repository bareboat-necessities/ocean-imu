#!/usr/bin/env python3
"""Compose the source-bound OU-III normal-Live stability theorem certificate.

This is an analytical proof-composition stage, not a replay-promotion stage.
It binds the current implementation to the validated scalar OU enclosure,
source-uniform translational UCO/UCC, conditional vector-packet UCO, complete
H/A process UCC, source Gaussian primitive model, and the exact PSD
nonexpansiveness proof for periodic a_w covariance synchronization.

The theorem is deliberately conditional.  Full-heading observability uses the
explicit vector persistent-excitation hypothesis of ``ou3_vector_uco_certificate``:
a proof packet contains accepted accel/mag vectors and the two magnetic packets
are consecutive configured 25 Hz packets.  In addition, the nonlinear local
argument is branch-regular: the nominal source word stays a positive distance
from innovation/gating discontinuities, so a sufficiently small error
neighborhood follows the same finite source branch word.  Arbitrary rejection
runs or nominal points exactly on a gate boundary are not certified here.

Under those hypotheses, uniform detectability/UCC and bounded source
coefficients give the standard bounded stabilizing Riccati family and UES of
each fixed-dimensional H/A normal-Live linearized recursion.  On a sufficiently
small geodesic chart the selected MEKF prediction/correction/reset branches are
C1 with uniformly quadratic remainder; UES plus variation of constants then
gives a nonzero local ISS neighborhood.

This producer does *not* claim an explicit basin radius or promote the stronger
startup/hard-reset/stochastic deployment-funnel certificate.  Those numerical
obligations remain separate and require a declared physical deployment
envelope plus validated word/prefix/hybrid/concentration bounds.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_full_process_ucc as PROCESS
import ou3_hybrid_aw_sync_proof as AW_SYNC
import ou3_scalar_ou_enclosure as SCALAR
import ou3_source_domain_contract as SOURCE
import ou3_source_noise_certificate as NOISE
import ou3_translational_uco_ucc as TRANS
import ou3_validated_transcendentals as VT
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
SCHEMA = 2
QUALIFICATION = "SOURCE_BOUND_CONDITIONAL_BRANCH_REGULAR_NORMAL_LIVE_LOCAL_ISS"


def _finite_positive(value) -> bool:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(x) and x > 0.0


def _finite_positive_interval(bounds) -> bool:
    if not isinstance(bounds, list) or len(bounds) != 2:
        return False
    try:
        lo, hi = map(float, bounds)
    except (TypeError, ValueError):
        return False
    return math.isfinite(lo) and math.isfinite(hi) and 0.0 < lo <= hi


def _source_failures(source: dict) -> tuple[dict, list[str]]:
    box = source.get("validated_parameter_box", {})
    runtime = source.get("configured_runtime_assumption", {})
    checks = {
        "source_complete_parameter_domain": source.get("source_complete_parameter_domain") is True,
        "source_generated_not_trajectory_fit": source.get("source_generated_not_trajectory_fit") is True,
        "validated_parameter_box": box.get("validated_arithmetic") is True,
        "outward_rounded_parameter_box": box.get("outward_rounded") is True,
        "configured_runtime_fixed_source_nominal": runtime.get("sample_period_contract") == "FIXED_SOURCE_NOMINAL",
        "configured_runtime_not_misstated_as_api_guard": runtime.get("api_enforces_this_bound") is False,
        "configured_runtime_positive_dt": _finite_positive(runtime.get("imu_dt_s")),
    }
    failures = [f"source-domain check failed: {name}" for name, ok in checks.items() if not ok]
    cp = box.get("continuous_parameters", {})
    for name in (
        "wave_tune_frequency_hz",
        "tau_aw_s",
        "sigma_aw_mps2",
        "R_S_base",
        "pseudo_update_period_s",
    ):
        if not _finite_positive_interval(cp.get(name)):
            failures.append(f"source continuous parameter is not finite positive: {name}")
    return checks, failures


def _noise_failures(noise: dict) -> tuple[dict, list[str]]:
    z = noise.get("standardized_increment", {})
    checks = {
        "schema": noise.get("schema") == NOISE.SCHEMA,
        "source_generated_not_trajectory_fit": noise.get("source_generated_not_trajectory_fit") is True,
        "standardized_dimension_18": z.get("dimension") == 18,
        "standardized_covariance_le_identity": z.get("covariance_upper_identity") is True,
    }
    failures = [f"source-noise check failed: {name}" for name, ok in checks.items() if not ok]
    for name, value in noise.get("physical_scales", {}).items():
        if not _finite_positive(value):
            failures.append(f"source-noise physical scale is not finite positive: {name}")
    return checks, failures


def _active_bias_stability(source: dict, process: dict) -> dict:
    dt = float(source["configured_runtime_assumption"]["imu_dt_s"])
    tau = float(process["source_constants"]["accel_bias_tau_s"])
    if not (_finite_positive(dt) and _finite_positive(tau)):
        return {"pass": False, "failure": "active accelerometer-bias dt/tau is invalid"}
    x = Interval.outward_bounds(dt / tau, dt / tau)
    if x.hi > VT.MAX_ABS_ARGUMENT:
        return {
            "pass": False,
            "dt_over_tau": x.as_list(),
            "failure": "dt/tau exceeds audited exponential range",
        }
    alpha = VT.exp_interval(-x)
    passed = 0.0 < alpha.lo <= alpha.hi < 1.0
    return {
        "pass": passed,
        "dt_s": dt,
        "tau_s": tau,
        "dt_over_tau": x.as_list(),
        "alpha_interval": alpha.as_list(),
        "failure": None if passed else "active accelerometer-bias decay is not strict",
    }


def _bindings() -> list[dict]:
    paths = (
        REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h",
        REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h",
        REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h",
    )
    return [
        {
            "path": str(path.relative_to(REPO)),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
        for path in paths
    ]


def build(header: Path = SOURCE.DEFAULT_HEADER.resolve()) -> dict:
    header = header.resolve()
    source = SOURCE.build(header)
    scalar = SCALAR.build(header)
    trans = TRANS.build(header)
    vector = VECTOR.build()
    process = PROCESS.build()
    noise = NOISE.build_certificate()
    aw = AW_SYNC.prove()

    source_checks, source_errors = _source_failures(source)
    scalar_errors = SCALAR.validate(scalar)
    trans_errors = TRANS.validate(trans)
    vector_errors = VECTOR.validate(vector)
    process_errors = PROCESS.validate(process)
    noise_checks, noise_errors = _noise_failures(noise)
    ba = _active_bias_stability(source, process)
    aw_errors = (
        [] if aw.get("status") == "PASS"
        else list(aw.get("failures", [])) or ["a_w sync proof failed"]
    )

    pseudo = trans.get("S_observation_uco", {})
    detect = trans.get("integrator_detectability", {})
    vec = vector.get("gyro_bias_two_packet", {})
    modes = process.get("modes", {})
    pe_source = vector.get("operating_envelope", {})

    obligations = {
        "compact_source_domain": not source_errors,
        "validated_exact_scalar_ou_transition": not scalar_errors,
        "translation_uco_ucc_detectability": not trans_errors,
        "full_heading_vector_uco_under_explicit_PE": not vector_errors,
        "full_state_process_ucc_H_A": not process_errors,
        "active_accelerometer_bias_stable_tail": ba.get("pass") is True,
        "finite_source_gaussian_primitives": not noise_errors,
        "periodic_aw_covariance_sync_nonexpansive": not aw_errors,
        "bounded_pseudo_measurement_gap": _finite_positive(pseudo.get("pseudo_gap_max_s")),
        "strict_translation_stable_tail": (
            _finite_positive(detect.get("stable_aw_alpha_upper"))
            and float(detect["stable_aw_alpha_upper"]) < 1.0
        ),
        "strict_vector_information_under_PE": _finite_positive(vec.get("alpha_6_information_lower")),
        "configured_consecutive_vector_packet_gap": _finite_positive_interval(pe_source.get("packet_gap_s")),
        "strict_H_process_excitation": (
            modes.get("H", {}).get("pass") is True
            and _finite_positive(modes.get("H", {}).get("prediction_Q_lambda_min_lower"))
        ),
        "strict_A_process_excitation": (
            modes.get("A", {}).get("pass") is True
            and _finite_positive(modes.get("A", {}).get("prediction_Q_lambda_min_lower"))
        ),
    }

    failures: list[str] = []
    failures += source_errors
    failures += [f"scalar OU: {x}" for x in scalar_errors]
    failures += [f"translation: {x}" for x in trans_errors]
    failures += [f"vector PE: {x}" for x in vector_errors]
    failures += [f"process: {x}" for x in process_errors]
    failures += noise_errors
    failures += [f"a_w sync: {x}" for x in aw_errors]
    if ba.get("pass") is not True:
        failures.append(f"active bias: {ba.get('failure')}")
    failures += [f"theorem obligation failed: {name}" for name, ok in obligations.items() if not ok]
    failures = list(dict.fromkeys(failures))

    theorem_pass = not failures

    pe = dict(pe_source)
    pe.update({
        "persistent_excitation_is_theorem_hypothesis": vector.get("persistent_excitation_is_theorem_hypothesis") is True,
        "accepted_accelerometer_packet_at_vector_times_required": True,
        "accepted_magnetometer_consecutive_pair_required": True,
        "measurement_gate_margin_required": True,
        "qualification": "THEOREM_OPERATING_ENVELOPE_NOT_INFERRED_FROM_EIGHT_REPLAY_TRAJECTORIES",
        "branch_regular_note": (
            "The local nonlinear theorem applies to nominal source words with positive "
            "innovation/gating margin so the same finite branch word persists in a "
            "sufficiently small error neighborhood. Gate-boundary points are excluded."
        ),
    })

    return {
        "schema": SCHEMA,
        "claim": "OU3_SOURCE_BOUND_BRANCH_REGULAR_NORMAL_LIVE_STABILITY_THEOREM_CERTIFICATE",
        "qualification": QUALIFICATION,
        "status": "PASS_CONDITIONAL_LOCAL_ISS" if theorem_pass else "FAIL",
        "sampled_evidence_used": False,
        "trajectory_fit_used": False,
        "implementation_bindings": _bindings(),
        "source_domain_checks": source_checks,
        "source_noise_checks": noise_checks,
        "theorem_obligations": obligations,
        "upstream_certificates": {
            "scalar_ou": scalar["qualification"],
            "translation": trans["qualification"],
            "vector_PE": vector["qualification"],
            "full_process": process["qualification"],
            "source_noise": noise["claim"],
            "periodic_aw_sync": aw["claim"],
        },
        "quantitative_anchors": {
            "pseudo_gap_max_s": pseudo.get("pseudo_gap_max_s"),
            "translation_information_lower": detect.get("information_gramian_lambda_min_lower"),
            "stable_aw_alpha_upper": detect.get("stable_aw_alpha_upper"),
            "vector_alpha_6_information_lower": vec.get("alpha_6_information_lower"),
            "vector_packet_gap_s": pe_source.get("packet_gap_s"),
            "H_prediction_Q_lambda_min_lower": modes.get("H", {}).get("prediction_Q_lambda_min_lower"),
            "A_prediction_Q_lambda_min_lower": modes.get("A", {}).get("prediction_Q_lambda_min_lower"),
            "active_accelerometer_bias": ba,
        },
        "persistent_excitation_operating_envelope": pe,
        "linearized_normal_live": {
            "H_dimension": 18,
            "A_dimension": 21,
            "uniform_detectability_and_stabilizability": theorem_pass,
            "bounded_stabilizing_Riccati": theorem_pass,
            "uniform_exponential_stability": theorem_pass,
            "proof_route": (
                "bounded-gap translational detectability + accepted consecutive vector-packet "
                "UCO + stable tails + full process UCC; standard discrete LTV Kalman/Riccati theorem"
            ),
        },
        "nonlinear_normal_live": {
            "local_iss": theorem_pass,
            "nonzero_neighborhood_exists": theorem_pass,
            "branch_regular_source_word_required": True,
            "explicit_numeric_basin_radius_produced": False,
            "proof_route": (
                "H/A UES + finite branch family + positive gate margin + uniform C1 MEKF "
                "branches on a geodesic chart + quadratic remainder + discrete small gain"
            ),
        },
        "periodic_aw_covariance_sync": {
            "pass": aw.get("status") == "PASS",
            "jump_gain_upper": aw.get("jump_gain_upper"),
            "proof_mode": aw.get("proof_mode"),
        },
        "scope_limits": {
            "configured_runtime_only": True,
            "arbitrary_positive_caller_dt": False,
            "full_heading_requires_persistent_excitation": True,
            "consecutive_accepted_mag_pair_required_for_vector_uco": True,
            "measurement_gate_boundary_points_not_certified": True,
            "permanent_or_unbounded_measurement_rejection_not_certified": True,
            "startup_handoff_and_hard_reset_funnel_numeric_enclosure": "SEPARATE_DEPLOYMENT_CERTIFICATE_OBLIGATION",
            "explicit_nonlinear_basin_radius": "SEPARATE_NUMERICAL_ENCLOSURE_OBLIGATION",
            "stochastic_infinite_horizon_pathwise_bound": "NOT_CLAIMED_FOR_GAUSSIAN_NOISE",
        },
        "claim_separation": {
            "analytical_branch_regular_normal_live_local_ISS": "PASS" if theorem_pass else "FAIL",
            "numerical_source_complete_deployment_funnel": "NOT_ESTABLISHED_BY_THIS_PRODUCER",
            "reason": (
                "the analytical theorem proves existence of a nonzero local ISS neighborhood "
                "under explicit PE/branch-regular hypotheses; the stronger deployment gate "
                "additionally requires explicit word/prefix/hybrid/stochastic/capture constants"
            ),
        },
        "failures": failures,
    }


def validate(payload: dict) -> list[str]:
    failures: list[str] = []
    if payload.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if payload.get("qualification") != QUALIFICATION:
        failures.append("qualification mismatch")
    if payload.get("sampled_evidence_used") is not False:
        failures.append("analytical certificate must not use sampled evidence")
    if payload.get("trajectory_fit_used") is not False:
        failures.append("analytical certificate must not use trajectory fitting")

    obligations = payload.get("theorem_obligations", {})
    if not obligations or not all(v is True for v in obligations.values()):
        failures.append("one or more theorem obligations did not pass")

    linear = payload.get("linearized_normal_live", {})
    if linear.get("uniform_exponential_stability") is not True:
        failures.append("linearized normal-Live UES was not established")

    nonlinear = payload.get("nonlinear_normal_live", {})
    if nonlinear.get("local_iss") is not True or nonlinear.get("nonzero_neighborhood_exists") is not True:
        failures.append("nonlinear normal-Live local ISS was not established")
    if nonlinear.get("branch_regular_source_word_required") is not True:
        failures.append("branch-regular source-word qualification is not explicit")
    if nonlinear.get("explicit_numeric_basin_radius_produced") is not False:
        failures.append("analytical producer must not masquerade as a numerical basin enclosure")

    pe = payload.get("persistent_excitation_operating_envelope", {})
    if pe.get("persistent_excitation_is_theorem_hypothesis") is not True:
        failures.append("full-heading PE qualification is not explicit")
    if pe.get("accepted_accelerometer_packet_at_vector_times_required") is not True:
        failures.append("accepted accelerometer vector-packet hypothesis is not explicit")
    if pe.get("accepted_magnetometer_consecutive_pair_required") is not True:
        failures.append("consecutive accepted magnetic-packet hypothesis is not explicit")
    if pe.get("measurement_gate_margin_required") is not True:
        failures.append("measurement gate-margin hypothesis is not explicit")
    if not _finite_positive_interval(pe.get("packet_gap_s")):
        failures.append("configured vector packet gap is not finite positive")

    if payload.get("status") != "PASS_CONDITIONAL_LOCAL_ISS":
        failures.append("certificate status is not PASS_CONDITIONAL_LOCAL_ISS")
    if payload.get("failures"):
        failures.append("certificate contains upstream failures")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--header", type=Path, default=SOURCE.DEFAULT_HEADER)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    payload = build(args.header.resolve())
    validation_failures = validate(payload)
    payload["validation_pass"] = not validation_failures
    payload["validation_failures"] = validation_failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "validation_pass": payload["validation_pass"],
        "linear_ues": payload["linearized_normal_live"]["uniform_exponential_stability"],
        "nonlinear_local_iss": payload["nonlinear_normal_live"]["local_iss"],
        "failures": payload["failures"],
        "validation_failures": validation_failures,
    }, indent=2, sort_keys=True))
    return 0 if not validation_failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
