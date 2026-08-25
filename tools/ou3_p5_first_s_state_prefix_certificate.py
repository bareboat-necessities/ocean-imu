#!/usr/bin/env python3
"""Conditional source-staged outer S-state prefix for OU-III P5.

This producer closes the scalar/vector S-error arithmetic needed by the first
due S=0 correction without falling back to the global normal-Live covariance
box.  It is deliberately a *conditional prefix lemma*, not a P5 promotion.
The remaining outer bootstrap must prove that the pre-first-S source prefix
stays in the declared finite-angle/latent-acceleration node.

At goLive the H-mode linear state satisfies the P1 product bounds and the
constructor covariance decomposition used by ``ou3_p5_first_s_gain_certificate``.
Before the first S pseudo, physical accelerometer and magnetometer corrections
have no S measurement column.  Write the state update as prediction plus the
actual additive Kalman correction.  Unrolling to the first S time T gives

  S(T) = S0 + T p0 + T^2 v0/2 + phi_Sa(T) a0
         + sum_j L_j K_j r_j,

with |phi_Sa(T)| <= T^3/6 for every stable OU time constant.  Later corrections
are separate additive terms in the same unrolling, so no hidden one-step
contraction factor is introduced.

The constructor S0/p0/v0 covariance component is invisible to every pre-first-S
physical measurement.  The first-S gain producer already bounds the complete
remaining covariance of the final-S functional by E_T.  For one isotropic
physical measurement R=sigma^2 I and any PSD covariance block, spectral
calculus gives

  ||L K|| <= sqrt(E_T)/(2 sigma).

This is the same sharp scalar maximum sqrt(x)/(x+sigma^2)<=1/(2 sigma), applied
to the source-correlated final-S/measurement block.  It covers accepted
physical corrections without choosing a favorable rejection pattern; rejected
or not-due packets contribute zero.

Residuals are bounded on a declared candidate outer bootstrap q_C<=1 and
||delta a_w||<=the P1 handoff radius.  The accelerometer-bias state is held in H
mode but its P1 error still enters the physical accelerometer residual and is
therefore retained.  Stochastic measurement/model disturbances are not hidden
inside this deterministic capture lemma; they are composed later by the paper's
stochastic localization stage.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_first_s_gain_certificate as FIRST
import ou3_startup_stability_certificate as P1
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
CANDIDATE_CAYLEY_NORM_UPPER = 1.0


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def mul_up(a: float, b: float) -> float:
    return up(float(a) * float(b))


def add_up(a: float, b: float) -> float:
    return up(float(a) + float(b))


def div_up(a: float, b: float) -> float:
    if not b > 0.0:
        raise RuntimeError("positive denominator required")
    return up(float(a) / float(b))


def sqrt_up(x: float) -> float:
    if not (math.isfinite(x) and x >= 0.0):
        raise RuntimeError("finite nonnegative square-root input required")
    return up(math.sqrt(x))


def _rotation_difference_upper(q: float) -> float:
    if not (0.0 <= q < 2.0):
        raise RuntimeError("candidate Cayley norm must lie in [0,2)")
    return div_up(mul_up(2.0, q), down(math.sqrt(4.0 + q*q)))


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("first-S state-prefix domain must not be trajectory fitted")

    first = FIRST.build(domain_path)
    p1 = P1.build(domain_path)
    vector = VECTOR.build()
    failures = [f"first-S-gain: {x}" for x in FIRST.validate(first)]
    failures += [f"P1: {x}" for x in P1.validate(p1)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]
    if failures:
        return {
            "schema": SCHEMA,
            "qualification": "OU3_P5_CONDITIONAL_FIRST_S_STATE_PREFIX",
            "source_generated_not_trajectory_fit": True,
            "source_replay_used": False,
            "P5_FIRST_S_STATE_PREFIX_CERTIFICATE": "NOT_ESTABLISHED",
            "failures": failures,
        }

    b = p1["go_live"]["physical_coordinate_bounds"]
    S0 = float(b["integral_displacement_error_norm_upper_m_s"])
    p0 = float(b["position_error_norm_upper_m"])
    v0 = float(b["velocity_error_norm_upper_mps"])
    aw0 = float(b["latent_acceleration_error_norm_upper_mps2"])
    ba0 = float(b["accelerometer_bias_error_norm_upper_mps2"])
    for label, x in (("S0", S0), ("p0", p0), ("v0", v0), ("aw0", aw0), ("ba0", ba0)):
        if not (math.isfinite(x) and x >= 0.0):
            raise RuntimeError(f"invalid P1 {label} radius")

    timing = first["timing"]
    T = float(timing["first_due_time_upper_s"])
    n_acc = int(timing["first_due_samples_upper"])
    if not (T > 0.0 and n_acc >= 1):
        raise RuntimeError("invalid first-S timing")

    meas = vector["configured_measurement_bounds"]
    acc_sigma = float(meas["acc_measurement_std_mps2"])
    mag_sigma = float(meas["mag_measurement_std_uT"])
    mag_hz = float(meas["mag_odr_hz"])
    if not (acc_sigma > 0.0 and mag_sigma > 0.0 and mag_hz > 0.0):
        raise RuntimeError("physical measurement noise/ODR bounds are not positive")
    n_mag = int(math.ceil(T * mag_hz)) + 1

    # Exact stable-chain open-loop contribution of the declared initial balls.
    t2 = mul_up(T, T)
    t3 = mul_up(t2, T)
    base = add_up(S0, mul_up(T, p0))
    base = add_up(base, mul_up(0.5 * t2, v0))
    base = add_up(base, mul_up(t3 / 6.0, aw0))

    E = float(first["aw_and_process_excess_P_SS_upper"])
    if not (math.isfinite(E) and E >= 0.0):
        raise RuntimeError("invalid first-S correlated final-S covariance excess")
    sqrtE = sqrt_up(E)
    LKa = div_up(sqrtE, mul_up(2.0, acc_sigma))
    LKm = div_up(sqrtE, mul_up(2.0, mag_sigma))

    q = CANDIDATE_CAYLEY_NORM_UPPER
    rot = _rotation_difference_upper(q)
    live = domain["normal_live"]
    fmax = float(live["specific_force_norm_upper_mps2"])
    mmax = float(live["magnetic_vector_norm_upper_uT"])
    if not (fmax > 0.0 and mmax > 0.0):
        raise RuntimeError("normal-Live vector norm ceilings must be positive")

    # Deterministic source residuals on the candidate outer node.  The H-mode
    # accelerometer-bias estimate is held, not nonexistent, so its declared
    # handoff error remains an additive residual source here.
    racc = add_up(mul_up(fmax, rot), add_up(aw0, ba0))
    rmag = mul_up(mmax, rot)

    acc_each = mul_up(LKa, racc)
    mag_each = mul_up(LKm, rmag)
    acc_all = mul_up(float(n_acc), acc_each)
    mag_all = mul_up(float(n_mag), mag_each)
    correction = add_up(acc_all, mag_all)
    Sbound = add_up(base, correction)

    KthetaS = float(first["K_thetaS_operator_norm_upper_first_due"])
    injection = mul_up(KthetaS, Sbound)
    helper_limit = 3.0
    finite_helper_pass = injection < helper_limit

    passed = all(math.isfinite(x) and x >= 0.0 for x in (
        base, E, LKa, LKm, racc, rmag, correction, Sbound, injection
    )) and finite_helper_pass

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_CONDITIONAL_FIRST_S_STATE_PREFIX",
        "claim": "SOURCE_STAGED_FIRST_S_STATE_PREFIX_CONDITIONAL_ON_OUTER_NODE_BOOTSTRAP",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "stochastic_disturbance_used_to_establish_deterministic_capture": False,
        "accepted_rejected_branch_policy": "charge every source-possible accepted correction; rejected/not-due contributes zero",
        "candidate_outer_bootstrap": {
            "cayley_norm_upper": q,
            "rotation_R_minus_I_norm_upper": rot,
            "latent_acceleration_error_norm_upper_mps2": aw0,
            "bootstrap_proved_here": False,
            "required_next": "exact large-angle source correction sector plus complete outer prefix invariance",
        },
        "first_due_timing": {
            "time_upper_s": T,
            "accelerometer_packets_charged_upper": n_acc,
            "magnetometer_packets_charged_upper": n_mag,
            "magnetometer_odr_hz": mag_hz,
        },
        "initial_open_loop_S_bound": {
            "S0_radius": S0,
            "p0_radius": p0,
            "v0_radius": v0,
            "aw0_radius": aw0,
            "phi_Sa_bound": "abs(phi_Sa(T)) <= T^3/6 for stable OU",
            "norm_upper_m_s": base,
        },
        "source_correlated_final_S_gain": {
            "E_T_covariance_excess_upper": E,
            "lemma": "||L K|| <= sqrt(E_T)/(2 sigma) for isotropic R=sigma^2 I",
            "acc_measurement_std_mps2": acc_sigma,
            "mag_measurement_std_uT": mag_sigma,
            "L_K_acc_operator_norm_upper": LKa,
            "L_K_mag_operator_norm_upper": LKm,
            "full_S_cross_gain_retained": True,
        },
        "deterministic_physical_residual_bounds": {
            "specific_force_norm_upper_mps2": fmax,
            "magnetic_vector_norm_upper_uT": mmax,
            "held_accelerometer_bias_error_norm_upper_mps2": ba0,
            "acc_residual_norm_upper_mps2": racc,
            "mag_residual_norm_upper_uT": rmag,
        },
        "physical_correction_final_S_contribution": {
            "per_acc_upper_m_s": acc_each,
            "all_acc_upper_m_s": acc_all,
            "per_mag_upper_m_s": mag_each,
            "all_mag_upper_m_s": mag_all,
            "total_upper_m_s": correction,
        },
        "first_due_S_error_norm_upper_m_s": Sbound,
        "first_due_K_thetaS_operator_norm_upper": KthetaS,
        "first_due_S_induced_attitude_correction_norm_upper_rad": injection,
        "deployed_group_helper_correction_limit_rad": helper_limit,
        "S_induced_correction_inside_group_helper": finite_helper_pass,
        "outer_node_bootstrap_supplied_here": False,
        "exact_large_angle_dissipation_supplied_here": False,
        "P5_FIRST_S_STATE_PREFIX_CERTIFICATE": "PASS_CONDITIONAL" if passed else "NOT_ESTABLISHED",
        "next_obligation": (
            "prove the candidate outer node is prefix invariant with the exact large-angle vector dissipation sector; "
            "then this first-S state/gain pair supplies the complete early S-to-attitude correction bound"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("first-S state prefix is not source generated")
    if d.get("source_replay_used") is not False:
        failures.append("first-S state prefix uses replay")
    if d.get("filter_changed") is not False:
        failures.append("first-S state prefix changes filter")
    if d.get("stochastic_disturbance_used_to_establish_deterministic_capture") is not False:
        failures.append("deterministic capture uses a stochastic fallback")
    if d.get("P5_FIRST_S_STATE_PREFIX_CERTIFICATE") != "PASS_CONDITIONAL":
        failures.append("conditional first-S state prefix did not pass")
    boot = d.get("candidate_outer_bootstrap", {})
    if boot.get("bootstrap_proved_here") is not False:
        failures.append("first-S state prefix silently proves its own outer bootstrap")
    if d.get("outer_node_bootstrap_supplied_here") is not False:
        failures.append("first-S state prefix claims outer-node closure")
    if d.get("exact_large_angle_dissipation_supplied_here") is not False:
        failures.append("first-S state prefix claims large-angle dissipation")
    gain = d.get("source_correlated_final_S_gain", {})
    if gain.get("full_S_cross_gain_retained") is not True:
        failures.append("first-S state prefix drops full S cross gain")
    for key in (
        "first_due_S_error_norm_upper_m_s",
        "first_due_K_thetaS_operator_norm_upper",
        "first_due_S_induced_attitude_correction_norm_upper_rad",
    ):
        x = d.get(key)
        if not (isinstance(x, (int, float)) and math.isfinite(float(x)) and float(x) > 0.0):
            failures.append(f"invalid {key}")
    if d.get("S_induced_correction_inside_group_helper") is not True:
        failures.append("first-S S-induced attitude correction exceeds validated group helper")
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
        "status": out.get("P5_FIRST_S_STATE_PREFIX_CERTIFICATE"),
        "S_first_due_upper": out.get("first_due_S_error_norm_upper_m_s"),
        "S_attitude_correction_upper": out.get("first_due_S_induced_attitude_correction_norm_upper_rad"),
        "bootstrap": out.get("candidate_outer_bootstrap"),
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
