#!/usr/bin/env python3
"""Source-staged first-due S->attitude gain bound for OU-III P5.

The global normal-Live covariance eigenvalue bound is the wrong object at the
startup handoff: using it in a PSD block inequality produced a completely
spurious ~2.8e5 operator bound for K_theta,S.  This producer uses covariance
structure that survives every admissible prefix before the first S=0 update.

At goLive the shipping source gives

    P_theta,S = P_theta,a = P_a,S = 0,
    P_SS = 50^2 I,
    pseudo_elapsed = 0.

Before the first S pseudo, accelerometer and magnetometer measurements have no S
column.  Their Joseph updates may create theta-S correlation indirectly through
a_w, but they cannot remove covariance carried by the independent constructor
S0, p0 and v0 components.  Prediction also never feeds S back into those
measurement-active states.  Thus at time t before the first pseudo the Schur
complement obeys

    P_SS - P_S,theta P_theta,theta^-1 P_theta,S >= D(t) I,
    D(t) = 50^2 + 20^2 t^2 + 1^2 t^4/4.

An open-loop upper for P_SS is D plus the a_w-driven contribution.  Accepted
physical measurements only reduce P_SS; attitude covariance resets leave the S
block unchanged.  Periodic a_w covariance synchronization is covered by adding
one full stationary a_w covariance contribution for every source-possible sync
before the first pseudo.

If P_SS <= U I, the theta/S canonical-correlation norm therefore satisfies

    rho_thetaS <= sqrt((U-D)/D).

For R_S with smallest eigenvalue r and lambda(P_SS)>=D>>r,

    ||K_thetaS||
      <= sqrt(lambda_max(P_theta)) rho_thetaS
         sqrt(D)/(D+r).

sqrt(D)/(D+r) decreases in r, so taking r as the *smallest* eigenvalue of the
deployed S=0 covariance diag(rho_x r_S, rho_y r_S, r_S)^2 keeps the bound valid
whether or not that covariance is isotropic; rho_x = rho_y = 1 is the isotropic
special case this certificate was first written against.

The theta bound is taken from the *directional theta block* of the source-uniform
P3 covariance enclosure, not from the translation-dominated global eigenvalue.
This is still conservative across all accepted/rejected prefixes and covariance
resets, while preserving the full S->attitude gain.

The result certifies the gain coefficient only.  P5 still has to bound the
outer S error at the first due pseudo and then prove the exact finite-angle
correction/funnel prefix.  No favorable rejection pattern and no replay are
used.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import ou3_explicit_information_word_certificate as P3
import ou3_full_process_ucc as PROCESS
import ou3_p5_go_live_covariance_stage as GOLIVE
import ou3_source_domain_contract as SOURCE
import ou3_source_reachable_matrix_p3 as P3BASE

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


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


def _deployed_member_float(text: str, name: str) -> float:
    """Value of a deployed `float <name> = <literal>f;` member of the wrapper."""
    m = re.search(rf"float\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f\s*;", text)
    if m is None:
        raise RuntimeError(f"cannot extract deployed member {name}")
    return float(m.group(1))


def _source_timing() -> dict:
    sched = P3BASE.source_schedule()
    tau_lo, tau_hi = map(float, sched["tau_applied_invariant_s"])
    r = float(sched["pseudo_ratio"])
    pmin = float(sched["pseudo_min_s"])
    pmax = float(sched["pseudo_max_s"])
    cadence_hi = up(min(max(r * tau_hi, pmin), pmax))
    dt = float(sched["dt_s"])
    samples = int(math.ceil(cadence_hi / dt)) + 1
    t_hi = up(samples * dt)

    text = WRAPPER.read_text(encoding="utf-8")
    adapt = float(SOURCE.parse_const(text, "ADAPT_EVERY_SECS"))
    if "last_aw_cov_sync_sec_ = time_;" not in text:
        raise RuntimeError("goLive tune no longer resets the a_w sync clock")
    # periodic_aw_cov_sync_tick_ uses a strict > adapt_every_secs_ test.  Ceil
    # is conservative at boundaries and remains valid if a source sample lands
    # one ulp to either side of the threshold.
    syncs = int(math.ceil(t_hi / adapt))
    return {
        "tau_applied_s": [tau_lo, tau_hi],
        "pseudo_cadence_upper_s": cadence_hi,
        "configured_dt_s": dt,
        "first_due_samples_upper": samples,
        "first_due_time_upper_s": t_hi,
        "aw_sync_period_s": adapt,
        "aw_sync_count_before_first_due_upper": syncs,
        "sigma_aw_std_upper_mps2": float(sched["sigma_aw_applied_safety"][1]),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("first-S gain domain must not be trajectory fitted")

    stage = GOLIVE.build(domain_path)
    sf = GOLIVE.validate(stage)
    p3 = P3.build(domain_path)
    pf = P3.validate(p3)
    proc = PROCESS.build()
    qf = PROCESS.validate(proc)
    failures = [f"goLive-stage: {x}" for x in sf] + [f"P3: {x}" for x in pf] + [f"process: {x}" for x in qf]
    if failures:
        return {
            "schema": SCHEMA,
            "qualification": "OU3_P5_FIRST_DUE_S_TO_ATTITUDE_GAIN",
            "source_generated_not_trajectory_fit": True,
            "source_replay_used": False,
            "P5_FIRST_DUE_S_GAIN_CERTIFICATE": "NOT_ESTABLISHED",
            "failures": failures,
        }

    wrapper = WRAPPER.read_text(encoding="utf-8")
    if "float S_factor_      = 1.0f;" not in wrapper:
        raise RuntimeError("configured first-S proof requires deployed S_factor_=1")
    # The two horizontal axes are independent knobs, so the bound has to take
    # the smaller of them rather than assume one horizontal scale.
    rho_x = _deployed_member_float(wrapper, "R_S_x_factor_")
    rho_y = _deployed_member_float(wrapper, "R_S_y_factor_")
    for name_, value in (("R_S_x_factor_", rho_x), ("R_S_y_factor_", rho_y)):
        if not (0.0 < value <= 4.0):
            raise RuntimeError(f"deployed {name_} out of setter range: {value}")
    rho_h = min(rho_x, rho_y)
    if domain.get("configured_runtime", {}).get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("configured first-S proof requires lever arm disabled")

    timing = _source_timing()
    t = float(timing["first_due_time_upper_s"])
    seed = stage["goLive_H_covariance_seed"]
    Pv = float(seed["P_vv_variance_per_axis"])
    Pp = float(seed["P_pp_variance_per_axis"])
    PS = float(seed["P_SS_variance_per_axis"])
    sigma = float(timing["sigma_aw_std_upper_mps2"])
    tau_lo = float(timing["tau_applied_s"][0])
    syncs = int(timing["aw_sync_count_before_first_due_upper"])

    # Persistent covariance invisible to every pre-first-S physical measurement.
    t2 = mul_up(t, t)
    t4 = mul_up(t2, t2)
    D = add_up(PS, add_up(mul_up(Pp, t2), mul_up(0.25 * Pv, t4)))

    # S contribution from the initial stationary a_w covariance and a grossly
    # conservative full stationary covariance injection at every possible sync.
    t3 = mul_up(t2, t)
    t6 = mul_up(t3, t3)
    a_coeff2 = div_up(t6, 36.0)
    aw_one = mul_up(sigma * sigma, a_coeff2)
    aw_all = mul_up(float(syncs + 1), aw_one)

    # Integrated-OU driving-noise S variance upper used by the established P3
    # translation comparison: q_c t^7 / 252, q_c=2 sigma^2/tau.
    t7 = mul_up(t6, t)
    qc = div_up(mul_up(2.0, sigma * sigma), tau_lo)
    process_S = div_up(mul_up(qc, t7), 252.0)
    excess = add_up(aw_all, process_S)
    PSS_upper = add_up(D, excess)
    rho = sqrt_up(div_up(excess, D))

    # Directional theta upper from P3.  Do not use Sigma_lambda_max_upper,
    # whose limiting coordinate is the translation block.
    mc = p3["modes"]["H"]["matrix_comparison"]
    diag = list(map(float, mc["Sigma_diagonal_upper"]))
    if len(diag) < 6:
        raise RuntimeError("P3 H directional covariance upper is incomplete")
    Ptheta_upper = up(max(diag[0:3]))

    # Smallest per-axis S=0 standard deviation the deployed schedule can reach.
    rs_std = down(min(rho_h, 1.0) * float(SOURCE.parse_const(wrapper, "MIN_R_S")))
    rmin = down(rs_std * rs_std)
    if not D > rmin:
        raise RuntimeError("first-S persistent variance no longer dominates minimum R_S")
    spectral_factor = div_up(sqrt_up(D), down(D + rmin))
    KthetaS = mul_up(sqrt_up(Ptheta_upper), mul_up(rho, spectral_factor))

    old_global = float(p3["modes"]["H"]["Sigma_lambda_max_upper"])
    old_K = div_up(sqrt_up(old_global), mul_up(2.0, rs_std))
    improvement = old_K / KthetaS if KthetaS > 0.0 else math.inf

    passed = all(math.isfinite(x) and x > 0.0 for x in (D, PSS_upper, rho, Ptheta_upper, KthetaS))
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_FIRST_DUE_S_TO_ATTITUDE_GAIN",
        "claim": "SOURCE_STAGED_CANONICAL_CORRELATION_BOUND_FOR_FIRST_S_PSEUDO",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "full_S_to_attitude_gain_retained": True,
        "accepted_and_rejected_physical_prefixes_covered": True,
        "measurement_update_argument": (
            "before first S pseudo H_S=0; Joseph corrections cannot remove the independent S0/p0/v0 covariance component; attitude reset leaves S unchanged"
        ),
        "timing": timing,
        "persistent_S_conditional_covariance_lambda_min_lower": D,
        "P_SS_lambda_max_upper_before_first_due": PSS_upper,
        "aw_and_process_excess_P_SS_upper": excess,
        "theta_S_canonical_correlation_upper": rho,
        "P_theta_theta_directional_lambda_max_upper": Ptheta_upper,
        "R_S_variance_lower": rmin,
        "PSS_spectral_gain_factor_upper": spectral_factor,
        "K_thetaS_operator_norm_upper_first_due": KthetaS,
        "old_global_P3_K_thetaS_operator_norm_upper": old_K,
        "gain_widening_factor_vs_global_P3_bound_lower": improvement,
        "global_translation_dominated_covariance_used_for_theta": False,
        "S_state_error_prefix_bound_supplied_here": False,
        "first_due_attitude_injection_bound_supplied_here": False,
        "P5_FIRST_DUE_S_GAIN_CERTIFICATE": "PASS" if passed else "NOT_ESTABLISHED",
        "next_obligation": (
            "combine this small source-staged K_thetaS coefficient with a validated outer S-state prefix/funnel bound, then evaluate the exact finite-angle correction sector"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("first-S gain is not source generated")
    if d.get("source_replay_used") is not False:
        failures.append("first-S gain uses replay")
    if d.get("filter_changed") is not False:
        failures.append("first-S gain changes filter")
    if d.get("full_S_to_attitude_gain_retained") is not True:
        failures.append("first-S gain drops S-to-attitude")
    if d.get("accepted_and_rejected_physical_prefixes_covered") is not True:
        failures.append("first-S gain selects favorable physical branches")
    if d.get("global_translation_dominated_covariance_used_for_theta") is not False:
        failures.append("first-S gain still uses global translation covariance for theta")
    for key in (
        "persistent_S_conditional_covariance_lambda_min_lower",
        "P_SS_lambda_max_upper_before_first_due",
        "theta_S_canonical_correlation_upper",
        "P_theta_theta_directional_lambda_max_upper",
        "K_thetaS_operator_norm_upper_first_due",
    ):
        x = d.get(key)
        if not (isinstance(x, (int, float)) and math.isfinite(float(x)) and float(x) > 0.0):
            failures.append(f"invalid {key}")
    D = float(d.get("persistent_S_conditional_covariance_lambda_min_lower", 0.0))
    U = float(d.get("P_SS_lambda_max_upper_before_first_due", 0.0))
    if not U >= D > 0.0:
        failures.append("first-S S covariance sandwich invalid")
    rho = float(d.get("theta_S_canonical_correlation_upper", math.inf))
    if not 0.0 < rho < 1.0:
        failures.append("first-S canonical correlation is not strict")
    if d.get("P5_FIRST_DUE_S_GAIN_CERTIFICATE") != "PASS":
        failures.append("first-S gain certificate did not pass")
    if d.get("S_state_error_prefix_bound_supplied_here") is not False:
        failures.append("first-S gain incorrectly claims outer state-prefix closure")
    if d.get("first_due_attitude_injection_bound_supplied_here") is not False:
        failures.append("first-S gain incorrectly promotes finite correction without S-state bound")
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
        "status": out.get("P5_FIRST_DUE_S_GAIN_CERTIFICATE"),
        "timing": out.get("timing"),
        "rho_thetaS": out.get("theta_S_canonical_correlation_upper"),
        "K_thetaS": out.get("K_thetaS_operator_norm_upper_first_due"),
        "widening": out.get("gain_widening_factor_vs_global_P3_bound_lower"),
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
