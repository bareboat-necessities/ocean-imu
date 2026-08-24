#!/usr/bin/env python3
"""Compose the source-bound OU-III analytical stability certificate.

This producer closes the *analytical* theorem chain from the validated proof
lemmas already present in the repository.  It deliberately does not turn replay
minima, sampled nonlinear radii, or synthetic final-gate fixtures into theorem
evidence.

The composition is:

* the shipping adaptation/timing domain is compact and source-derived;
* the exact scalar OU transition is enclosed with outward-rounded arithmetic;
* the translational chain has source-uniform UCO/detectability and UCC;
* attitude/gyro-bias UCO holds under the explicit vector persistent-excitation
  operating envelope (full heading cannot be unconditional);
* the complete H/A prediction process covariance has a strict source-uniform
  lower bound;
* the active accelerometer-bias Gauss-Markov tail is uniformly exponentially
  stable;
* the source Gaussian primitive-noise model is finite and normalized; and
* periodic a_w covariance synchronization is nonexpansive in the information
  metric by an exact Loewner-order argument.

Uniform detectability plus uniform complete controllability, bounded source
coefficients and positive finite measurement covariances are the standard
Riccati hypotheses used by the manuscript.  They give bounded stabilizing
Riccati solutions and UES of each fixed-dimensional normal-Live linearized
Kalman recursion.  Piecewise-C1 MEKF prediction/correction/reset maps then have
uniform quadratic local remainder on a sufficiently small geodesic chart, so
the usual variation-of-constants/small-gain argument gives a nonzero local ISS
neighborhood.

The resulting PASS is intentionally scoped: it is a conditional normal-Live
local-ISS theorem, not a numerical deployment-funnel certificate.  In
particular it does not claim an explicit basin radius, arbitrary permanent
magnetometer rejection, or the still-separate startup/hard-reset/hybrid funnel
inequalities.  Those stronger quantitative obligations remain owned by
``ou3_validate_enclosure.py``/``ou3_deployment_gate.py``.
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
SCHEMA = 1
QUALIFICATION = "SOURCE_BOUND_CONDITIONAL_NORMAL_LIVE_LOCAL_ISS"


def _finite_positive(x) -> bool:
    try:
        y = float(x)
    except (TypeError, ValueError):
        return False
    return math.isfinite(y) and y > 0.0


def _finite_interval(bounds, *, positive: bool = False) -> bool:
    if not isinstance(bounds, list) or len(bounds) != 2:
        return False
    try:
        lo, hi = map(float, bounds)
    except (TypeError, ValueError):
        return False
    return (
        math.isfinite(lo)
        and math.isfinite(hi)
        and lo <= hi
        and (not positive or lo > 0.0)
    )


def _source_checks(source: dict) -> tuple[dict, list[str]]:
    failures: list[str] = []
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
    for name, ok in checks.items():
        if not ok:
            failures.append(f"source-domain check failed: {name}")

    cp = box.get("continuous_parameters", {})
    for name in (
        "wave_tune_frequency_hz",
        "tau_aw_s",
        "sigma_aw_mps2",
        "R_S_base",
        "pseudo_update_period_s",
    ):
        if not _finite_interval(cp.get(name), positive=True):
            failures.append(f"source continuous parameter is not finite positive: {name}")
    return checks, failures


def _noise_checks(noise: dict) -> tuple[dict, list[str]]:
    failures: list[str] = []
    z = noise.get("standardized_increment", {})
    scales = noise.get("physical_scales", {})
    checks = {
        "schema": noise.get("schema") == NOISE.SCHEMA,
        "source_generated_not_trajectory_fit": noise.get("source_generated_not_trajectory_fit") is True,
        "standardized_dimension_18": z.get("dimension") == 18,
        "standardized_covariance_le_identity": z.get("covariance_upper_identity") is True,
    }
    for name, ok in checks.items():
        if not ok:
            failures.append(f"source-noise check failed: {name}")
    for name, value in scales.items():
        if not _finite_positive(value):
            failures.append(f"source-noise physical scale is not finite positive: {name}")
    return checks, failures


def _active_bias_stability(source: dict, process: dict) -> dict:
    dt = float(source["configured_runtime_assumption"]["imu_dt_s"])
    tau = float(process["source_constants"]["accel_bias_tau_s"])
    if not (dt > 0.0 and tau > 0.0 and math.isfinite(dt) and math.isfinite(tau)):
        return {
            "pass": False,
            "failure": "active accelerometer-bias dt/tau is not finite positive",
        }
    x = Interval.outward_bounds(dt / tau, dt / tau)
    if x.hi > VT.MAX_ABS_ARGUMENT:
        return {
            "pass": False,
            "failure": "active accelerometer-bias dt/tau exceeds validated exponential range",
            "dt_over_tau": x.as_list(),
        }
    alpha = VT.exp_interval(-x)
    passed = 0.0 < alpha.lo <= alpha.hi < 1.0
    return {
        "pass": passed,
        "dt_s": dt,
        "tau_s": tau,
        "dt_over_tau": x.as_list(),
        "alpha_interval": alpha.as_list(),
        "failure": None if passed else "active accelerometer-bias decay is not strictly below one",
    }


def _implementation_binding(paths: tuple[Path, ...]) -> list[dict]:
    rows = []
    for path in paths:
        text = path.read_bytes()
        rows.append({
            "path": str(path.relative_to(REPO)),
            "sha256": hashlib.sha256(text).hexdigest(),
        })
    return rows


def build(header: Path = SOURCE.DEFAULT_HEADER.resolve()) -> dict:
    header = header.resolve()
    source = SOURCE.build(header)
    scalar = SCALAR.build(header)
    trans = TRANS.build(header)
    vector = VECTOR.build()
    process = PROCESS.build()
    noise = NOISE.build_certificate()
    aw = AW_SYNC.prove()

    source_checks, source_failures = _source_checks(source)
    scalar_failures = SCALAR.validate(scalar)
    trans_failures = TRANS.validate(trans)
    vector_failures = VECTOR.validate(vector)
    process_failures = PROCESS.validate(process)
    noise_checks, noise_failures = _noise_checks(noise)
    ba = _active_bias_stability(source, process)

    aw_failures = []
    if aw.get("status") != "PASS":
        aw_failures.extend(str(x) for x in aw.get("failures", []))
        if not aw_failures:
            aw_failures.append("periodic a_w covariance synchronization proof did not pass")

    pseudo = trans.get("S_observation_uco", {})
    detect = trans.get("integrator_detectability", {})
    vec = vector.get("gyro_bias_two_packet", {})
    modes = process.get("modes", {})

    theorem_obligations = {
        "compact_source_domain": not source_failures,
        "validated_exact_scalar_ou_transition": not scalar_failures,
        "translation_uco_ucc_detectability": not trans_failures,
        "full_heading_vector_uco_under_explicit_PE": not vector_failures,
        "full_state_process_ucc_H_A": not process_failures,
        "active_accelerometer_bias_stable_tail": ba.get("pass") is True,
        "finite_source_gaussian_primitives": not noise_failures,
        "periodic_aw_covariance_sync_nonexpansive": not aw_failures,
        "bounded_pseudo_measurement_gap": _finite_positive(pseudo.get("pseudo_gap_max_s")),
        "strict_translation_stable_tail": (
            _finite_positive(detect.get("stable_aw_alpha_upper"))
            and float(detect["stable_aw_alpha_upper"]) < 1.0
        ),
        "strict_vector_information_under_PE": _finite_positive(vec.get("alpha_6_information_lower")),
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
    failures.extend(source_failures)
    failures.extend(f"scalar OU: {x}" for x in scalar_failures)
    failures.extend(f"translation: {x}" for x in trans_failures)
    failures.extend(f"vector PE: {x}" for x in vector_failures)
    failures.extend(f"process: {x}" for x in process_failures)
    failures.extend(noise_failures)
    failures.extend(f"a_w sync: {x}" for x in aw_failures)
    if ba.get("pass") is not True:
        failures.append(f"active bias: {ba.get('failure')}")
    for name, ok in theorem_obligations.items():
        if not ok:
            failures.append(f"theorem obligation failed: {name}")

    # De-duplicate while retaining deterministic order.
    failures = list(dict.fromkeys(failures))
    linear_pass = not failures
    nonlinear_local_pass = linear_pass

    pe = dict(vector.get("operating_envelope", {}))
    pe["accepted_vector_packet_bounded_gap_required"] = True
    pe["qualification"] = (
        "THEOREM_OPERATING_ENVELOPE_NOT_INFERRED_FROM_EIGHT_REPLAY_TRAJECTORIES"
    )

    bindings = _implementation_binding((
        REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h",
        REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h",
        REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h",
    ))

    status = "PASS_CONDITIONAL_LOCAL_ISS" if nonlinear_local_pass else "FAIL"
    return {
        "schema": SCHEMA,
        "claim": "OU3_SOURCE_BOUND_NORMAL_LIVE_STABILITY_THEOREM_CERTIFICATE",
        "qualification": QUALIFICATION,
        "status": status,
        "sampled_evidence_used": False,
        "trajectory_fit_used": False,
        "implementation_bindings": bindings,
        "source_domain_checks": source_checks,
        "source_noise_checks": noise_checks,
        "theorem_obligations": theorem_obligations,
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
            "H_prediction_Q_lambda_min_lower": modes.get("H", {}).get("prediction_Q_lambda_min_lower"),
            "A_prediction_Q_lambda_min_lower": modes.get("A", {}).get("prediction_Q_lambda_min_lower"),
            "active_accelerometer_bias": ba,
        },
        "persistent_excitation_operating_envelope": pe,
        "linearized_normal_live": {
            "H_dimension": 18,
            "A_dimension": 21,
            "uniform_detectability_and_stabilizability": linear_pass,
            "bounded_stabilizing_Riccati": linear_pass,
            "uniform_exponential_stability": linear_pass,
            "proof_route": (
                "source-uniform bounded-gap translational detectability + conditional vector "
                "UCO + stable Gauss-Markov tails + full process UCC; standard discrete LTV "
                "Kalman/Riccati stability theorem"
            ),
        },
        "nonlinear_normal_live": {
            "local_iss": nonlinear_local_pass,
            "nonzero_neighborhood_exists": nonlinear_local_pass,
            "explicit_numeric_basin_radius_produced": False,
            "proof_route": (
                "fixed-mode UES plus uniform piecewise-C1 MEKF map on a geodesic chart; "
                "the nonlinear remainder is quadratic and discrete variation-of-constants "
                "gives a sufficiently small invariant ISS neighborhood"
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
            "permanent_or_unbounded_mag_rejection_not_certified": True,
            "startup_handoff_and_hard_reset_funnel_numeric_enclosure": "SEPARATE_DEPLOYMENT_CERTIFICATE_OBLIGATION",
            "explicit_nonlinear_basin_radius": "SEPARATE_NUMERICAL_ENCLOSURE_OBLIGATION",
            "stochastic_infinite_horizon_pathwise_bound": "NOT_CLAIMED_FOR_GAUSSIAN_NOISE",
        },
        "claim_separation": {
            "analytical_normal_live_local_ISS": "PASS" if nonlinear_local_pass else "FAIL",
            "numerical_source_complete_deployment_funnel": "NOT_ESTABLISHED_BY_THIS_PRODUCER",
            "reason": (
                "the analytical theorem proves existence of a nonzero local ISS neighborhood; "
                "the stronger numerical deployment gate additionally asks for explicit word, "
                "prefix, hybrid, stochastic and finite-capture constants"
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
    nonlinear = payload.get("nonlinear_normal_live", {})
    if linear.get("uniform_exponential_stability") is not True:
        failures.append("linearized normal-Live UES was not established")
    if nonlinear.get("local_iss") is not True or nonlinear.get("nonzero_neighborhood_exists") is not True:
        failures.append("nonlinear normal-Live local ISS was not established")
    if nonlinear.get("explicit_numeric_basin_radius_produced") is not False:
        failures.append("analytical producer must not masquerade as a numerical basin enclosure")
    pe = payload.get("persistent_excitation_operating_envelope", {})
    if pe.get("persistent_excitation_is_theorem_hypothesis") is not True:
        failures.append("full-heading PE qualification is not explicit")
    if pe.get("accepted_vector_packet_bounded_gap_required") is not True:
        failures.append("bounded accepted-vector packet gap hypothesis is not explicit")
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
        "explicit_numeric_basin_radius_produced": payload["nonlinear_normal_live"]["explicit_numeric_basin_radius_produced"],
        "failures": payload["failures"],
        "validation_failures": validation_failures,
    }, indent=2, sort_keys=True))
    return 0 if not validation_failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
