#!/usr/bin/env python3
"""Second-generation widening of the OU-III P4 nonlinear source-word certificate.

This refinement stacks two theorem-preserving improvements on top of
``ou3_p4_nextgen_widened_certificate``:

1. Measurement-class directional operator bounds.  The first next-generation
   producer already separates prediction/S/accelerometer/magnetometer nonlinear
   defects, but it still gives every accepted correction the one global linear
   correction gain based on the largest measurement operator in the source
   domain.  Here S=0, accelerometer, and magnetometer use separate safe H-norm
   envelopes before Cayley/quaternion composition:

       ||H_S|| <= 1,
       ||H_acc|| <= 2 f_max + 2,
       ||H_mag|| <= 2 m_max + 2.

   Each is no larger than the legacy global envelope 2 max(f_max,m_max)+2.
   No branch is removed and no replay data are used.

2. Exact endpoint budget allocation.  Earlier P4 chose the convenient
   sufficient budget B sqrt(W) <= delta/8.  To retain the same advertised
   endpoint decrease W_1 <= (1-delta/2) W_0, it is sufficient and sharper to
   impose

       sqrt(1-delta) + B sqrt(W) <= sqrt(1-delta/2).

   The positive difference is evaluated without cancellation as

       [delta/2] / [sqrt(1-delta/2)+sqrt(1-delta)].

   The prefix bootstrap remains valid because B sqrt(W) <= 1 implies every
   homogeneous-unit-gain prefix plus accumulated defect stays inside 4 W_0.

The producer fails closed unless both refinements are monotone relative to the
preceding next-generation certificate, preserve chart/quaternion/projection
safety, and keep the same P3 endpoint margin and Cayley information metric.
"""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path

import ou3_p4_nonlinear_word_certificate as LEGACY
import ou3_p4_nextgen_widened_certificate as P4W

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def _sqrt_one_minus_up(x: float) -> float:
    if not (0.0 <= x < 1.0):
        raise ValueError("x must lie in [0,1)")
    # sqrt_up is outward on the supplied binary64 point.  Feeding an upward
    # rounded 1-x makes the resulting value an upper bound.
    return LEGACY.sqrt_up(math.nextafter(1.0 - float(x), math.inf))


def _endpoint_sqrt_gap_lower(delta: float) -> float:
    """Lower bound on sqrt(1-delta/2)-sqrt(1-delta), cancellation free."""
    if not (0.0 < delta < 1.0):
        raise ValueError("P3 delta must lie in (0,1)")
    numerator = LEGACY.mul_down(0.5, delta)
    a = _sqrt_one_minus_up(0.5 * delta)
    b = _sqrt_one_minus_up(delta)
    return LEGACY.div_down(numerator, LEGACY.add_up(a, b))


def _refine_mode(mode: str, base: dict, domain: dict) -> dict:
    m = copy.deepcopy(base)
    live = domain["normal_live"]
    Kmax = float(base["full_gain_norm_upper"])
    q_design = float(base["correction_quadratic_bound"]["design_error_norm_radius"])
    fmax = float(live["specific_force_norm_upper_mps2"])
    magmax = float(live["magnetic_vector_norm_upper_uT"])

    Hglobal = float(base["measurement_linear_operator_norm_upper"])
    HS = 1.0
    Hacc = LEGACY.add_up(LEGACY.mul_up(2.0, fmax), 2.0)
    Hmag = LEGACY.add_up(LEGACY.mul_up(2.0, magmax), 2.0)
    for name, h in (("S_zero", HS), ("accelerometer", Hacc), ("magnetometer", Hmag)):
        if h > Hglobal:
            raise RuntimeError(f"{mode}: {name} directional H bound exceeds global H bound")

    LS = LEGACY.mul_up(Kmax, HS)
    Lacc = LEGACY.mul_up(Kmax, Hacc)
    Lmag = LEGACY.mul_up(Kmax, Hmag)

    Cvec_acc = float(base["vector_residual_quadratic_constant_acc_upper"])
    Cvec_mag = float(base["vector_residual_quadratic_constant_mag_upper"])
    Cinput_acc = LEGACY.mul_up(Kmax, Cvec_acc)
    Cinput_mag = LEGACY.mul_up(Kmax, Cvec_mag)

    corrS = LEGACY._composition_quadratic_constant(LS, 0.0, q_design)
    corrA = LEGACY._composition_quadratic_constant(Lacc, Cinput_acc, q_design)
    corrM = LEGACY._composition_quadratic_constant(Lmag, Cinput_mag, q_design)

    Cpred = float(base["operation_specific_quadratic_defect_constants_upper"]["prediction"])
    CS = float(corrS["full_state_quadratic_defect_constant_upper"])
    Cacc = float(corrA["full_state_quadratic_defect_constant_upper"])
    Cmag = float(corrM["full_state_quadratic_defect_constant_upper"])
    directional_sum = LEGACY.add_up(LEGACY.add_up(Cpred, CS), LEGACY.add_up(Cacc, Cmag))
    previous_sum = float(base["operation_specific_defect_sum_per_sample_upper"])
    if directional_sum > previous_sum:
        raise RuntimeError(f"{mode}: directional operation sum is not monotone")

    samples = int(base["word_samples_upper"])
    mmin = float(base["metric_lambda_min_lower"])
    mmax = float(base["metric_lambda_max_upper"])
    delta = float(base["P3_word_endpoint_delta_lower"])

    B = LEGACY.mul_up(LEGACY.PREFIX_BOOTSTRAP_W_FACTOR, float(samples))
    B = LEGACY.mul_up(B, LEGACY.sqrt_up(mmax))
    B = LEGACY.mul_up(B, directional_sum)
    B = LEGACY.div_up(B, mmin)
    Bprev = float(base["transported_word_defect_B_upper"])
    if not (0.0 < B <= Bprev and math.isfinite(B)):
        raise RuntimeError(f"{mode}: directional B is not a finite monotone refinement")

    # Exact endpoint allocation for the same target relative decrease delta/2.
    sqrt_gap = _endpoint_sqrt_gap_lower(delta)
    sqrtW_endpoint = LEGACY.div_down(sqrt_gap, B)
    if not sqrtW_endpoint > 0.0:
        raise RuntimeError(f"{mode}: exact endpoint budget produced no positive radius")

    # Prefix bootstrap: sqrt(W_prefix) <= sqrt(W0)+B W0 <= 2 sqrt(W0).
    sqrtW_bootstrap = LEGACY.div_down(1.0, B)

    # Keep the same nonlinear-design chart assumption used to derive all C's.
    sqrtW_design = LEGACY.mul_down(
        0.5, LEGACY.sqrt_up(LEGACY.mul_down(mmin, q_design * q_design))
    )
    sqrtW = min(sqrtW_endpoint, sqrtW_bootstrap, sqrtW_design)

    # A-mode projection interior can impose an additional source-faithful cap.
    projection = copy.deepcopy(base.get("active_bias_projection"))
    if mode == "A" and projection is not None:
        margin = float(projection["interior_margin_lower_mps2"])
        sqrtW_projection = LEGACY.mul_down(
            0.5, LEGACY.sqrt_up(LEGACY.mul_down(mmin, margin * margin))
        )
        sqrtW = min(sqrtW, sqrtW_projection)

    Wstar = LEGACY.mul_down(sqrtW, sqrtW)
    Wprev = float(base["certified_level_W"])
    if Wstar < Wprev:
        raise RuntimeError(f"{mode}: second-generation refinement regressed W")

    qprefix = LEGACY.mul_up(2.0, LEGACY.sqrt_up(LEGACY.div_up(Wstar, mmin)))
    if not qprefix < q_design:
        raise RuntimeError(f"{mode}: second-generation level leaves design-radius bootstrap")
    if not qprefix < LEGACY.PROMOTED_CAYLEY_NORM_LIMIT:
        raise RuntimeError(f"{mode}: second-generation level reaches Cayley chart boundary")

    correction_bounds = {
        "S_zero": LEGACY.add_up(LEGACY.mul_up(LS, qprefix), 0.0),
        "accelerometer": LEGACY.add_up(
            LEGACY.mul_up(Lacc, qprefix), LEGACY.mul_up(Cinput_acc, qprefix * qprefix)
        ),
        "magnetometer": LEGACY.add_up(
            LEGACY.mul_up(Lmag, qprefix), LEGACY.mul_up(Cinput_mag, qprefix * qprefix)
        ),
    }
    correction_prefix = max(correction_bounds.values())
    if not correction_prefix < 1.0e-2:
        raise RuntimeError(f"{mode}: second-generation correction crosses quaternion series branch")

    if mode == "A" and projection is not None:
        margin = float(projection["interior_margin_lower_mps2"])
        projection["certified_error_norm_prefix_upper"] = qprefix
        projection["projection_surface_reached_in_certified_funnel"] = not (qprefix < margin)
        if not qprefix < margin:
            raise RuntimeError("A: second-generation P4 reaches bias projection surface")

    # Directly verify the exact target inequality using an upper LHS and lower
    # target gap construction rather than forming a near-one endpoint ratio.
    nonlinear_sqrt = LEGACY.mul_up(B, sqrtW)
    if nonlinear_sqrt > sqrt_gap:
        raise RuntimeError(f"{mode}: exact endpoint sqrt budget does not close")
    if LEGACY.mul_up(B, sqrtW) > 1.0:
        raise RuntimeError(f"{mode}: prefix bootstrap B*sqrt(W)<=1 does not close")

    m.update({
        "directional_measurement_operator_norm_upper": {
            "global_previous": Hglobal,
            "S_zero": HS,
            "accelerometer": Hacc,
            "magnetometer": Hmag,
        },
        "directional_linear_correction_gain_upper": {
            "S_zero": LS,
            "accelerometer": Lacc,
            "magnetometer": Lmag,
        },
        "directional_operation_quadratic_defect_constants_upper": {
            "prediction": Cpred,
            "S_zero_accepted": CS,
            "accelerometer_accepted": Cacc,
            "magnetometer_accepted": Cmag,
        },
        "directional_operation_defect_sum_per_sample_upper": directional_sum,
        "operation_specific_defect_sum_per_sample_upper_previous": previous_sum,
        "directional_defect_sum_monotone": directional_sum <= previous_sum,
        "transported_word_defect_B_upper_previous": Bprev,
        "transported_word_defect_B_upper": B,
        "directional_B_reduction_factor_lower": LEGACY.div_down(Bprev, B),
        "exact_endpoint_sqrt_gap_lower": sqrt_gap,
        "exact_endpoint_budget_target": "sqrt(1-delta)+B*sqrt(W)<=sqrt(1-delta/2)",
        "legacy_delta_over_8_budget_replaced": True,
        "prefix_bootstrap_B_sqrt_W_upper": LEGACY.mul_up(B, sqrtW),
        "certified_level_W_previous_nextgen": Wprev,
        "certified_level_sqrt_W_previous_nextgen": float(base["certified_level_sqrt_W"]),
        "certified_level_W": Wstar,
        "certified_level_sqrt_W": sqrtW,
        "secondgen_W_widening_factor_lower": LEGACY.div_down(Wstar, Wprev),
        "secondgen_sqrt_W_widening_factor_lower": LEGACY.div_down(
            sqrtW, float(base["certified_level_sqrt_W"])
        ),
        "total_W_widening_factor_vs_legacy_lower": LEGACY.div_down(
            Wstar, float(base["certified_level_W_legacy"])
        ),
        "prefix_canonical_error_norm_upper": qprefix,
        "accepted_correction_norm_prefix_upper": correction_prefix,
        "accepted_correction_norms_by_class_upper": correction_bounds,
        "active_bias_projection": projection,
        "nextgen_directional_operator_transport": True,
        "nextgen_exact_endpoint_budget": True,
        "exact_nonlinear_word_pass": True,
    })
    return m


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("second-generation P4 theorem domain is trajectory fitted")

    previous = P4W.build(domain_path)
    failures = [f"previous P4W: {x}" for x in P4W.validate(previous)]
    modes = {}
    if not failures:
        for mode in ("H", "A"):
            try:
                modes[mode] = _refine_mode(mode, previous["modes"][mode], domain)
            except Exception as exc:
                failures.append(f"{mode}: {exc}")

    out = copy.deepcopy(previous)
    out["qualification"] = "VALIDATED_NEXTGEN_DIRECTIONAL_EXACT_BUDGET_CAYLEY_SOURCE_WORD_CERTIFICATE"
    out["claim"] = "P4_NEXTGEN_DIRECTIONAL_AND_EXACT_BUDGET_WIDENED_H_A_WORD_DISSIPATION"
    out["modes"] = modes
    out["nextgen_directional_operator_refinement"] = True
    out["nextgen_exact_endpoint_budget_refinement"] = True
    out["source_subdivision"]["future_widening_allowed"] = (
        "joint source-node subdivision and support-aware metric transport remain independent further refinements"
    )
    passed = not failures and all(modes.get(k, {}).get("exact_nonlinear_word_pass") for k in ("H", "A"))
    out["P4_EXACT_NONLINEAR_WORD_CERTIFICATE"] = "PASS" if passed else "FAIL"
    out["P4_NEXTGEN_WIDENED_WORD_CERTIFICATE"] = "PASS" if passed else "FAIL"
    out["P4_NEXTGEN_DIRECTIONAL_WORD_CERTIFICATE"] = "PASS" if passed else "FAIL"
    out["theorem_promotion"] = "P4_NEXTGEN_DIRECTIONAL_NORMAL_LIVE_EXACT_WORDS" if passed else "NOT_ESTABLISHED"
    out["failures"] = failures
    out["next_obligation"] = (
        "evaluate source-node subdivision and support-aware metric transport, then compare widened P4 to P5 outer capture"
    )
    return out


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("nextgen_directional_operator_refinement") is not True:
        failures.append("directional operator refinement missing")
    if d.get("nextgen_exact_endpoint_budget_refinement") is not True:
        failures.append("exact endpoint budget refinement missing")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        if m.get("nextgen_directional_operator_transport") is not True:
            failures.append(f"{mode}: directional transport missing")
            continue
        if m.get("nextgen_exact_endpoint_budget") is not True:
            failures.append(f"{mode}: exact endpoint budget missing")
        if m.get("directional_defect_sum_monotone") is not True:
            failures.append(f"{mode}: directional defect sum is not monotone")
        if float(m.get("transported_word_defect_B_upper", math.inf)) > float(
            m.get("transported_word_defect_B_upper_previous", -math.inf)
        ):
            failures.append(f"{mode}: directional B regressed")
        if float(m.get("certified_level_W", 0.0)) < float(
            m.get("certified_level_W_previous_nextgen", math.inf)
        ):
            failures.append(f"{mode}: second-generation W regressed")
        if float(m.get("secondgen_W_widening_factor_lower", 0.0)) < 1.0:
            failures.append(f"{mode}: second-generation W widening factor below one")
        if float(m.get("prefix_bootstrap_B_sqrt_W_upper", math.inf)) > 1.0:
            failures.append(f"{mode}: prefix bootstrap does not close")
        if not float(m.get("prefix_canonical_error_norm_upper", math.inf)) < float(
            m.get("cayley_norm_limit", 0.0)
        ):
            failures.append(f"{mode}: Cayley chart safety failed")
        if not float(m.get("accepted_correction_norm_prefix_upper", math.inf)) < 1.0e-2:
            failures.append(f"{mode}: quaternion branch safety failed")
        if mode == "A" and m.get("active_bias_projection", {}).get(
            "projection_surface_reached_in_certified_funnel"
        ) is not False:
            failures.append("A: bias projection interior safety failed")
    if not failures and d.get("P4_NEXTGEN_DIRECTIONAL_WORD_CERTIFICATE") != "PASS":
        failures.append("second-generation P4 status is not PASS")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve())
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "P4_NEXTGEN_DIRECTIONAL_WORD_CERTIFICATE": d["P4_NEXTGEN_DIRECTIONAL_WORD_CERTIFICATE"],
        "numerical": {
            mode: {
                "W_legacy": d.get("modes", {}).get(mode, {}).get("certified_level_W_legacy"),
                "W_stage1": d.get("modes", {}).get(mode, {}).get("certified_level_W_previous_nextgen"),
                "W_stage2": d.get("modes", {}).get(mode, {}).get("certified_level_W"),
                "stage2_factor": d.get("modes", {}).get(mode, {}).get("secondgen_W_widening_factor_lower"),
                "total_factor": d.get("modes", {}).get(mode, {}).get("total_W_widening_factor_vs_legacy_lower"),
                "directional_B_factor": d.get("modes", {}).get(mode, {}).get("directional_B_reduction_factor_lower"),
            } for mode in ("H", "A")
        },
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
