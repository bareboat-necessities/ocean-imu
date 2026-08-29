#!/usr/bin/env python3
"""Next-generation widened OU-III P4 exact nonlinear word certificate.

This producer is a theorem-preserving refinement of
``ou3_p4_nonlinear_word_certificate``.  The legacy P4 proof transports every
nonlinear state-operation defect with the P3 unit prefix-information bound, but
charges every possible state operation the single largest source-uniform defect
constant C.  With at most prediction, S, accelerometer and magnetometer state
maps per IMU sample this gives the safe but coarse factor

    4 * N_samples * max(C_pred, C_S, C_acc, C_mag).

The source word map already fixes those four operation classes and their order.
Linearity of the residual-defect transport therefore permits the strictly
sharper source-complete sum

    N_samples * (C_pred + C_S + C_acc + C_mag),

without assuming that an update is accepted and without trajectory replay.
Rejected/not-due corrections contribute zero, so charging each correction once
per sample remains an upper bound for every admissible branch.  The S=0
innovation is exactly linear, hence C_S has no vector-residual input defect; it
retains only the exact quaternion/Cayley composition defect.  Accelerometer and
magnetometer corrections receive their own physical vector-residual constants
rather than the maximum of the two.

No P3 margin, information metric, source language, shipping operation, reset,
projection, or chart assumption is weakened.  The result is therefore the same
P4 theorem with a smaller transported nonlinear defect B and a larger certified
W level.  The legacy certificate is retained as a baseline and the refinement
fails closed unless its defect sum is no larger than the legacy four-operation
max budget and its certified level is no smaller in both H and A modes.
"""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path

import ou3_p4_nonlinear_word_certificate as LEGACY

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def _refine_mode(mode: str, base: dict, domain: dict) -> dict:
    m = copy.deepcopy(base)
    live = domain["normal_live"]

    Kmax = float(base["full_gain_norm_upper"])
    Lcorr = float(base["correction_quadratic_bound"]["linear_correction_gain_L"])
    q_design = float(base["correction_quadratic_bound"]["design_error_norm_radius"])
    fmax = float(live["specific_force_norm_upper_mps2"])
    magmax = float(live["magnetic_vector_norm_upper_uT"])

    Cvec_acc = LEGACY.add_up(
        LEGACY.mul_up(LEGACY.ROTATION_REMAINDER_COEFF, fmax), 1.5
    )
    Cvec_mag = LEGACY.mul_up(LEGACY.ROTATION_REMAINDER_COEFF, magmax)
    Cinput_acc = LEGACY.mul_up(Kmax, Cvec_acc)
    Cinput_mag = LEGACY.mul_up(Kmax, Cvec_mag)

    # The S residual is exactly r_S=-delta_S.  It has no nonlinear measurement
    # remainder, but its accepted correction still uses the deployed nonlinear
    # quaternion injection and Cayley composition, so Cinput=0 rather than C_S=0.
    corr_S = LEGACY._composition_quadratic_constant(Lcorr, 0.0, q_design)
    corr_acc = LEGACY._composition_quadratic_constant(Lcorr, Cinput_acc, q_design)
    corr_mag = LEGACY._composition_quadratic_constant(Lcorr, Cinput_mag, q_design)

    Cpred = float(base["prediction_quadratic_bound"]["full_state_quadratic_defect_constant_upper"])
    CS = float(corr_S["full_state_quadratic_defect_constant_upper"])
    Cacc = float(corr_acc["full_state_quadratic_defect_constant_upper"])
    Cmag = float(corr_mag["full_state_quadratic_defect_constant_upper"])
    operation_sum = LEGACY.add_up(LEGACY.add_up(Cpred, CS), LEGACY.add_up(Cacc, Cmag))

    samples = int(base["word_samples_upper"])
    old_C = float(base["uniform_operation_quadratic_defect_constant_upper"])
    old_per_sample = LEGACY.mul_up(4.0, old_C)
    if operation_sum > old_per_sample:
        raise RuntimeError(
            f"{mode}: operation-specific defect sum exceeds legacy four-operation max budget"
        )

    mmin = float(base["metric_lambda_min_lower"])
    mmax = float(base["metric_lambda_max_upper"])
    sqrt_mmax = LEGACY.sqrt_up(mmax)
    delta = float(base["P3_word_endpoint_delta_lower"])

    B = LEGACY.mul_up(LEGACY.PREFIX_BOOTSTRAP_W_FACTOR, float(samples))
    B = LEGACY.mul_up(B, sqrt_mmax)
    B = LEGACY.mul_up(B, operation_sum)
    B = LEGACY.div_up(B, mmin)
    if not (math.isfinite(B) and B > 0.0):
        raise RuntimeError(f"{mode}: next-generation nonlinear word defect gain is invalid")

    legacy_B = float(base["transported_word_defect_B_upper"])
    if B > legacy_B:
        raise RuntimeError(f"{mode}: next-generation B is not a refinement of legacy P4")

    sqrt_W_star = LEGACY.div_down(delta, LEGACY.mul_up(8.0, B))
    W_star = LEGACY.mul_down(sqrt_W_star, sqrt_W_star)
    legacy_W = float(base["certified_level_W"])
    if not (W_star >= legacy_W and W_star > 0.0):
        raise RuntimeError(f"{mode}: next-generation P4 failed to widen legacy W level")

    q_prefix = LEGACY.mul_up(2.0, LEGACY.sqrt_up(LEGACY.div_up(W_star, mmin)))
    if not q_prefix < q_design:
        raise RuntimeError(f"{mode}: widened level does not close design-radius bootstrap")
    if not q_prefix < LEGACY.PROMOTED_CAYLEY_NORM_LIMIT:
        raise RuntimeError(f"{mode}: widened level reaches Cayley chart boundary")

    Cinput_max = max(Cinput_acc, Cinput_mag)
    correction_prefix = LEGACY.add_up(
        LEGACY.mul_up(Lcorr, q_prefix),
        LEGACY.mul_up(Cinput_max, LEGACY.mul_up(q_prefix, q_prefix)),
    )
    if not correction_prefix < 1.0e-2:
        raise RuntimeError(f"{mode}: widened correction can cross source quaternion branch")

    nonlinear_sqrt_fraction = LEGACY.div_up(LEGACY.mul_up(B, sqrt_W_star), delta)
    if nonlinear_sqrt_fraction > 0.125000000000001:
        raise RuntimeError(f"{mode}: widened nonlinear endpoint budget exceeds delta/8")

    projection = copy.deepcopy(base.get("active_bias_projection"))
    if mode == "A" and projection is not None:
        margin = float(projection["interior_margin_lower_mps2"])
        projection["certified_error_norm_prefix_upper"] = q_prefix
        projection["projection_surface_reached_in_certified_funnel"] = not (q_prefix < margin)
        if not q_prefix < margin:
            raise RuntimeError("A: widened P4 funnel reaches accelerometer-bias projection surface")

    m.update({
        "vector_residual_quadratic_constant_acc_upper": Cvec_acc,
        "vector_residual_quadratic_constant_mag_upper": Cvec_mag,
        "operation_specific_quadratic_defect_constants_upper": {
            "prediction": Cpred,
            "S_zero_accepted": CS,
            "accelerometer_accepted": Cacc,
            "magnetometer_accepted": Cmag,
        },
        "operation_specific_defect_sum_per_sample_upper": operation_sum,
        "legacy_four_operation_max_defect_per_sample_upper": old_per_sample,
        "operation_specific_sum_no_larger_than_legacy_max_budget": operation_sum <= old_per_sample,
        "correction_quadratic_bound_S_zero": corr_S,
        "correction_quadratic_bound_accelerometer": corr_acc,
        "correction_quadratic_bound_magnetometer": corr_mag,
        "transported_word_defect_B_upper_legacy": legacy_B,
        "transported_word_defect_B_upper": B,
        "transported_word_defect_B_reduction_factor_lower": LEGACY.div_down(legacy_B, B),
        "certified_level_W_legacy": legacy_W,
        "certified_level_sqrt_W_legacy": float(base["certified_level_sqrt_W"]),
        "certified_level_W": W_star,
        "certified_level_sqrt_W": sqrt_W_star,
        "certified_level_W_widening_factor_lower": LEGACY.div_down(W_star, legacy_W),
        "certified_level_sqrt_W_widening_factor_lower": LEGACY.div_down(
            sqrt_W_star, float(base["certified_level_sqrt_W"])
        ),
        "prefix_canonical_error_norm_upper": q_prefix,
        "accepted_correction_norm_prefix_upper": correction_prefix,
        "nonlinear_sqrt_budget_fraction_of_delta_upper": nonlinear_sqrt_fraction,
        "active_bias_projection": projection,
        "nextgen_operation_specific_defect_transport": True,
        "exact_nonlinear_word_pass": True,
    })
    return m


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("next-generation P4 theorem domain is trajectory fitted")

    legacy = LEGACY.build(domain_path)
    legacy_failures = LEGACY.validate(legacy)
    failures = [f"legacy P4: {x}" for x in legacy_failures]
    modes: dict[str, dict] = {}
    if not failures:
        for mode in ("H", "A"):
            try:
                modes[mode] = _refine_mode(mode, legacy["modes"][mode], domain)
            except Exception as exc:
                failures.append(f"{mode}: {exc}")

    out = copy.deepcopy(legacy)
    out["schema"] = legacy["schema"]
    out["qualification"] = "VALIDATED_NEXTGEN_OPERATION_SPECIFIC_CAYLEY_NONLINEAR_SOURCE_WORD_CERTIFICATE"
    out["claim"] = "P4_NEXTGEN_WIDENED_EXACT_NONLINEAR_H_A_WORD_DISSIPATION_AND_PREFIX_SAFETY"
    out["modes"] = modes
    out["word_branch_coverage"]["method"] = (
        "operation-specific exact nonlinear defect sum plus P3 unit segment information transport"
    )
    out["source_subdivision"] = {
        "kind": "ANALYTIC_OPERATION_CLASS_REFINEMENT",
        "cartesian_source_subdivision_used": False,
        "refinement": "prediction/S/accelerometer/magnetometer nonlinear defects are bounded separately before source-complete word transport",
        "future_widening_allowed": "joint source-node and radius subdivision may further reduce operation-class constants without changing the theorem route",
    }
    passed = not failures and all(
        modes.get(mode, {}).get("exact_nonlinear_word_pass") is True for mode in ("H", "A")
    )
    out["P4_EXACT_NONLINEAR_WORD_CERTIFICATE"] = "PASS" if passed else "FAIL"
    out["P4_NEXTGEN_WIDENED_WORD_CERTIFICATE"] = "PASS" if passed else "FAIL"
    out["theorem_promotion"] = "P4_NEXTGEN_WIDENED_NORMAL_LIVE_EXACT_WORDS" if passed else "NOT_ESTABLISHED"
    out["legacy_P4_retained_as_baseline"] = True
    out["nextgen_refinement_source_only"] = True
    out["failures"] = failures
    out["next_obligation"] = (
        "compare widened P4 inner funnel against P5 startup/outer capture; only continue P5 transport where a gap remains"
    )
    return out


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("nextgen_refinement_source_only") is not True:
        failures.append("next-generation P4 refinement is not source-only")
    if d.get("legacy_P4_retained_as_baseline") is not True:
        failures.append("legacy P4 baseline was not retained")
    if d.get("P4_NEXTGEN_WIDENED_WORD_CERTIFICATE") != "PASS" and not failures:
        failures.append("next-generation P4 status is not PASS")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        if m.get("nextgen_operation_specific_defect_transport") is not True:
            failures.append(f"{mode}: operation-specific defect transport missing")
            continue
        if m.get("operation_specific_sum_no_larger_than_legacy_max_budget") is not True:
            failures.append(f"{mode}: operation-specific defect sum is not a safe refinement")
        if float(m.get("transported_word_defect_B_upper", math.inf)) > float(
            m.get("transported_word_defect_B_upper_legacy", -math.inf)
        ):
            failures.append(f"{mode}: transported defect B did not decrease")
        if float(m.get("certified_level_W", 0.0)) < float(m.get("certified_level_W_legacy", math.inf)):
            failures.append(f"{mode}: certified W level regressed")
        if float(m.get("certified_level_W_widening_factor_lower", 0.0)) < 1.0:
            failures.append(f"{mode}: W widening factor is below one")
        if float(m.get("certified_level_sqrt_W_widening_factor_lower", 0.0)) < 1.0:
            failures.append(f"{mode}: sqrt(W) widening factor is below one")
        if not float(m.get("prefix_canonical_error_norm_upper", math.inf)) < float(
            m.get("cayley_norm_limit", 0.0)
        ):
            failures.append(f"{mode}: widened prefix chart bootstrap failed")
        if not float(m.get("accepted_correction_norm_prefix_upper", math.inf)) < 1.0e-2:
            failures.append(f"{mode}: widened correction leaves deployed quaternion series branch")
        if mode == "A":
            p = m.get("active_bias_projection", {})
            if p.get("projection_surface_reached_in_certified_funnel") is not False:
                failures.append("A: widened funnel reaches bias projection surface")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve())
    failures = validate(d)
    out = dict(d)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    compact = {
        mode: {
            "W_star_legacy": out.get("modes", {}).get(mode, {}).get("certified_level_W_legacy"),
            "W_star_nextgen": out.get("modes", {}).get(mode, {}).get("certified_level_W"),
            "W_widening_factor_lower": out.get("modes", {}).get(mode, {}).get("certified_level_W_widening_factor_lower"),
            "sqrt_W_widening_factor_lower": out.get("modes", {}).get(mode, {}).get("certified_level_sqrt_W_widening_factor_lower"),
            "B_reduction_factor_lower": out.get("modes", {}).get(mode, {}).get("transported_word_defect_B_reduction_factor_lower"),
        }
        for mode in ("H", "A")
    }
    print(json.dumps({
        "P4_NEXTGEN_WIDENED_WORD_CERTIFICATE": out["P4_NEXTGEN_WIDENED_WORD_CERTIFICATE"],
        "numerical": compact,
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
