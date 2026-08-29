#!/usr/bin/env python3
"""Third widening layer for OU-III P4: measurement-specific gain envelopes.

The directional P4 producer already separates the measurement operators, but
it still inherits one global Kalman-gain norm bound based on the smallest
measurement covariance among all correction classes.  The generic covariance
inequality used by P4 is class-local:

    ||K_j|| <= sqrt(lambda_max(Sigma) / lambda_min(R_j)).

Therefore S=0, accelerometer and magnetometer can use their own source-certified
R lower bounds.  Each class-local K bound is no larger than the previous global
K bound because the latter used min_j lambda_min(R_j).  This refinement changes
no source language, branch coverage, metric, P3 margin or deployed operation.
It is source-only and fails closed unless the transported defect decreases and
the certified W level is monotone for both H and A.
"""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path

import ou3_p4_nonlinear_word_certificate as LEGACY
import ou3_p4_nextgen_directional_certificate as P4D

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"


def _refine_mode(mode: str, base: dict) -> dict:
    m = copy.deepcopy(base)
    smax = float(base["Sigma_lambda_max_upper"])
    meas = base["measurement_bounds"]
    r = {
        "S_zero": float(meas["S_zero_variance_lower"]),
        "accelerometer": float(meas["acc_variance_lower"]),
        "magnetometer": float(meas["mag_variance_lower"]),
    }
    Kglobal = float(base["full_gain_norm_upper"])
    K = {name: LEGACY.sqrt_up(LEGACY.div_up(smax, rv)) for name, rv in r.items()}
    for name, kval in K.items():
        if kval > Kglobal:
            raise RuntimeError(f"{mode}: class-local {name} gain exceeds global gain bound")

    h = base["directional_measurement_operator_norm_upper"]
    L = {
        "S_zero": LEGACY.mul_up(K["S_zero"], float(h["S_zero"])),
        "accelerometer": LEGACY.mul_up(K["accelerometer"], float(h["accelerometer"])),
        "magnetometer": LEGACY.mul_up(K["magnetometer"], float(h["magnetometer"])),
    }
    q_design = float(base["correction_quadratic_bound"]["design_error_norm_radius"])
    Cvec_acc = float(base["vector_residual_quadratic_constant_acc_upper"])
    Cvec_mag = float(base["vector_residual_quadratic_constant_mag_upper"])
    Cinput = {
        "accelerometer": LEGACY.mul_up(K["accelerometer"], Cvec_acc),
        "magnetometer": LEGACY.mul_up(K["magnetometer"], Cvec_mag),
    }

    corrS = LEGACY._composition_quadratic_constant(L["S_zero"], 0.0, q_design)
    corrA = LEGACY._composition_quadratic_constant(
        L["accelerometer"], Cinput["accelerometer"], q_design
    )
    corrM = LEGACY._composition_quadratic_constant(
        L["magnetometer"], Cinput["magnetometer"], q_design
    )
    Cpred = float(base["directional_operation_quadratic_defect_constants_upper"]["prediction"])
    C = {
        "prediction": Cpred,
        "S_zero_accepted": float(corrS["full_state_quadratic_defect_constant_upper"]),
        "accelerometer_accepted": float(corrA["full_state_quadratic_defect_constant_upper"]),
        "magnetometer_accepted": float(corrM["full_state_quadratic_defect_constant_upper"]),
    }
    csum = LEGACY.add_up(
        LEGACY.add_up(C["prediction"], C["S_zero_accepted"]),
        LEGACY.add_up(C["accelerometer_accepted"], C["magnetometer_accepted"]),
    )
    prev_sum = float(base["directional_operation_defect_sum_per_sample_upper"])
    if csum > prev_sum:
        raise RuntimeError(f"{mode}: class-local gain defect sum is not monotone")

    samples = int(base["word_samples_upper"])
    mmin = float(base["metric_lambda_min_lower"])
    mmax = float(base["metric_lambda_max_upper"])
    B = LEGACY.mul_up(LEGACY.PREFIX_BOOTSTRAP_W_FACTOR, float(samples))
    B = LEGACY.mul_up(B, LEGACY.sqrt_up(mmax))
    B = LEGACY.mul_up(B, csum)
    B = LEGACY.div_up(B, mmin)
    Bprev = float(base["transported_word_defect_B_upper"])
    if not (0.0 < B <= Bprev and math.isfinite(B)):
        raise RuntimeError(f"{mode}: class-local gain B is not monotone")

    delta = float(base["P3_word_endpoint_delta_lower"])
    sqrt_gap = P4D._endpoint_sqrt_gap_lower(delta)
    sqrtW_endpoint = LEGACY.div_down(sqrt_gap, B)
    sqrtW_bootstrap = LEGACY.div_down(1.0, B)
    # Reuse the previous proven design/projection-safe radius as a floor on
    # available safety, while allowing the new endpoint budget to expand only
    # if it remains inside the same explicit q_design check below.
    sqrtW = min(sqrtW_endpoint, sqrtW_bootstrap)
    Wprev = float(base["certified_level_W"])
    sqrtWprev = float(base["certified_level_sqrt_W"])
    if sqrtW < sqrtWprev:
        raise RuntimeError(f"{mode}: class-local gain refinement regressed sqrt(W)")

    Wstar = LEGACY.mul_down(sqrtW, sqrtW)
    qprefix = LEGACY.mul_up(2.0, LEGACY.sqrt_up(LEGACY.div_up(Wstar, mmin)))
    if not qprefix < q_design:
        # The endpoint improvement is larger than the existing analytic
        # remainder chart supports.  Cap to the previous certified radius;
        # this is fail-safe and leaves later radius subdivision free to widen.
        sqrtW = sqrtWprev
        Wstar = Wprev
        qprefix = float(base["prefix_canonical_error_norm_upper"])

    correction_by_class = {
        "S_zero": LEGACY.mul_up(L["S_zero"], qprefix),
        "accelerometer": LEGACY.add_up(
            LEGACY.mul_up(L["accelerometer"], qprefix),
            LEGACY.mul_up(Cinput["accelerometer"], qprefix * qprefix),
        ),
        "magnetometer": LEGACY.add_up(
            LEGACY.mul_up(L["magnetometer"], qprefix),
            LEGACY.mul_up(Cinput["magnetometer"], qprefix * qprefix),
        ),
    }
    correction_prefix = max(correction_by_class.values())
    if not correction_prefix < 1.0e-2:
        raise RuntimeError(f"{mode}: class-local gain correction leaves quaternion branch")

    projection = copy.deepcopy(base.get("active_bias_projection"))
    if mode == "A" and projection is not None:
        margin = float(projection["interior_margin_lower_mps2"])
        projection["certified_error_norm_prefix_upper"] = qprefix
        projection["projection_surface_reached_in_certified_funnel"] = not (qprefix < margin)
        if not qprefix < margin:
            raise RuntimeError("A: class-local gain refinement reaches bias projection surface")

    if LEGACY.mul_up(B, sqrtW) > sqrt_gap:
        raise RuntimeError(f"{mode}: class-local gain endpoint budget does not close")
    if LEGACY.mul_up(B, sqrtW) > 1.0:
        raise RuntimeError(f"{mode}: class-local gain prefix bootstrap does not close")

    m.update({
        "measurement_specific_R_lambda_min_lower": r,
        "measurement_specific_gain_norm_upper": K,
        "global_gain_norm_upper_previous": Kglobal,
        "measurement_specific_gain_bounds_monotone": all(v <= Kglobal for v in K.values()),
        "gain_refined_linear_correction_gain_upper": L,
        "gain_refined_operation_quadratic_defect_constants_upper": C,
        "gain_refined_operation_defect_sum_per_sample_upper": csum,
        "directional_operation_defect_sum_per_sample_upper_previous": prev_sum,
        "gain_refined_defect_sum_monotone": csum <= prev_sum,
        "transported_word_defect_B_upper_previous_gain_stage": Bprev,
        "transported_word_defect_B_upper": B,
        "gain_stage_B_reduction_factor_lower": LEGACY.div_down(Bprev, B),
        "certified_level_W_previous_gain_stage": Wprev,
        "certified_level_sqrt_W_previous_gain_stage": sqrtWprev,
        "certified_level_W": Wstar,
        "certified_level_sqrt_W": sqrtW,
        "gain_stage_W_widening_factor_lower": LEGACY.div_down(Wstar, Wprev),
        "gain_stage_sqrt_W_widening_factor_lower": LEGACY.div_down(sqrtW, sqrtWprev),
        "total_W_widening_factor_vs_legacy_lower": LEGACY.div_down(
            Wstar, float(base["certified_level_W_legacy"])
        ),
        "prefix_canonical_error_norm_upper": qprefix,
        "accepted_correction_norm_prefix_upper": correction_prefix,
        "accepted_correction_norms_by_class_upper": correction_by_class,
        "active_bias_projection": projection,
        "nextgen_measurement_specific_gain_transport": True,
        "exact_nonlinear_word_pass": True,
    })
    return m


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    previous = P4D.build(Path(domain_path).resolve())
    failures = [f"directional P4: {x}" for x in P4D.validate(previous)]
    modes = {}
    if not failures:
        for mode in ("H", "A"):
            try:
                modes[mode] = _refine_mode(mode, previous["modes"][mode])
            except Exception as exc:
                failures.append(f"{mode}: {exc}")
    out = copy.deepcopy(previous)
    out["qualification"] = "VALIDATED_NEXTGEN_MEASUREMENT_SPECIFIC_GAIN_CAYLEY_SOURCE_WORD_CERTIFICATE"
    out["claim"] = "P4_NEXTGEN_MEASUREMENT_SPECIFIC_GAIN_WIDENED_H_A_WORD_DISSIPATION"
    out["modes"] = modes
    out["nextgen_measurement_specific_gain_refinement"] = True
    passed = not failures and all(modes.get(k, {}).get("exact_nonlinear_word_pass") for k in ("H", "A"))
    out["P4_NEXTGEN_GAIN_WORD_CERTIFICATE"] = "PASS" if passed else "FAIL"
    out["P4_NEXTGEN_DIRECTIONAL_WORD_CERTIFICATE"] = "PASS" if passed else "FAIL"
    out["P4_NEXTGEN_WIDENED_WORD_CERTIFICATE"] = "PASS" if passed else "FAIL"
    out["P4_EXACT_NONLINEAR_WORD_CERTIFICATE"] = "PASS" if passed else "FAIL"
    out["theorem_promotion"] = "P4_NEXTGEN_GAIN_REFINED_NORMAL_LIVE_EXACT_WORDS" if passed else "NOT_ESTABLISHED"
    out["failures"] = failures
    out["next_obligation"] = "support-aware metric transport and joint source-node subdivision remain available P4 widening layers"
    return out


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("nextgen_measurement_specific_gain_refinement") is not True:
        failures.append("measurement-specific gain refinement missing")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        if m.get("nextgen_measurement_specific_gain_transport") is not True:
            failures.append(f"{mode}: class-local gain transport missing")
            continue
        if m.get("measurement_specific_gain_bounds_monotone") is not True:
            failures.append(f"{mode}: class-local gain bounds exceed global bound")
        if m.get("gain_refined_defect_sum_monotone") is not True:
            failures.append(f"{mode}: gain-refined defect sum regressed")
        if float(m.get("transported_word_defect_B_upper", math.inf)) > float(
            m.get("transported_word_defect_B_upper_previous_gain_stage", -math.inf)
        ):
            failures.append(f"{mode}: gain-stage B regressed")
        if float(m.get("certified_level_W", 0.0)) < float(
            m.get("certified_level_W_previous_gain_stage", math.inf)
        ):
            failures.append(f"{mode}: gain-stage W regressed")
        if not float(m.get("accepted_correction_norm_prefix_upper", math.inf)) < 1.0e-2:
            failures.append(f"{mode}: quaternion branch safety failed")
        if mode == "A" and m.get("active_bias_projection", {}).get(
            "projection_surface_reached_in_certified_funnel"
        ) is not False:
            failures.append("A: bias projection interior safety failed")
    if not failures and d.get("P4_NEXTGEN_GAIN_WORD_CERTIFICATE") != "PASS":
        failures.append("gain-refined P4 status is not PASS")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "P4_NEXTGEN_GAIN_WORD_CERTIFICATE": d["P4_NEXTGEN_GAIN_WORD_CERTIFICATE"],
        "numerical": {
            mode: {
                "W_legacy": d.get("modes", {}).get(mode, {}).get("certified_level_W_legacy"),
                "W_before_gain": d.get("modes", {}).get(mode, {}).get("certified_level_W_previous_gain_stage"),
                "W_after_gain": d.get("modes", {}).get(mode, {}).get("certified_level_W"),
                "gain_W_factor": d.get("modes", {}).get(mode, {}).get("gain_stage_W_widening_factor_lower"),
                "total_W_factor": d.get("modes", {}).get(mode, {}).get("total_W_widening_factor_vs_legacy_lower"),
                "gain_B_factor": d.get("modes", {}).get(mode, {}).get("gain_stage_B_reduction_factor_lower"),
            } for mode in ("H", "A")
        },
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
