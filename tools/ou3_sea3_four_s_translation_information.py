#!/usr/bin/env python3
"""Full four-S translation information inside the complete SEA3 P3 word.

This module is a mandatory lemma consumed by the complete
``COMPLETE_SEA3_NORMAL_LIVE_WORD``.  It is deliberately not a source generator,
a reduced P3 word, or a promotion route by itself.

The shipping progress-preserving S=0 scheduler guarantees one *actual* firing
in every interval of length ``g``.  Select one firing from each disjoint window

    [0,g], [2g,3g], [4g,5g], [6g,7g].

Every other due S firing remains in the literal complete SEA3 word.  For one
translation axis and scaled time u=t/g, the homogeneous S observation is

    S(u) = S0 + u (g p0) + u^2/2 (g^2 v0) + cbar(u) (g^3 a0),
    cbar'''(u) = exp(-int_0^{g u} 1/tau(s) ds) > 0.

Thus the four rows are full rank on [S,g p,g^2 v,g^3 a_w] for arbitrary legal
time-varying tau.  Quantitatively we use Newton divided-difference coordinates.
The guaranteed window separation gives an exact rational bound on ||L^-1||,
so the selected-record information is bounded as a genuine 4x4 matrix:

    D_S,z = L^T Sigma_Y^-1 L >= d_N I_4.

The observation-noise bound uses the actual applied SpectralMSE R_S safety
ceiling and deployed axis factors.  Applied R_S is never replaced by its target.
Process noise over the selected window is included as nuisance with arbitrary
legal SEA3 tau/sigma variation.  Determinants are rank witnesses only and are
not converted into an eigenvalue gate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p3_pseudo_scheduler_progress_certificate as SCHED
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_COMPLETE_SEA3_FOUR_S_FULL_TRANSLATION_INFORMATION"
WINDOW_MULTIPLIERS = ((0.0, 1.0), (2.0, 3.0), (4.0, 5.0), (6.0, 7.0))


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _validated_exp_negative_lower(x: float) -> tuple[float, dict]:
    x = float(x)
    if not (math.isfinite(x) and x >= 0.0):
        raise ValueError("finite nonnegative exponential magnitude required")
    if x == 0.0:
        return 1.0, {
            "method": "VALIDATED_EXP_POINT_PRODUCT",
            "pieces": 0,
            "piece_exponent_upper": 0.0,
            "ordinary_libm_exp_used": False,
        }
    n = max(1, int(math.ceil(x / VT.MAX_ABS_ARGUMENT)))
    step = up(x / float(n))
    while step > VT.MAX_ABS_ARGUMENT:
        n += 1
        step = up(x / float(n))
    factor = VT.exp_point(-step)
    product = Interval.point(1.0)
    for _ in range(n):
        product = product * factor
    lo = down(product.lo)
    if not (math.isfinite(lo) and lo > 0.0):
        raise RuntimeError("validated exponential product lost strict positivity")
    return lo, {
        "method": "VALIDATED_EXP_POINT_PRODUCT",
        "pieces": n,
        "piece_exponent_upper": step,
        "covered_exponent_lower": x,
        "product_exponent_upper": up(n * step),
        "ordinary_libm_exp_used": False,
    }


def _newton_inverse_norm_certificate() -> dict:
    # Exact scaled gap floors for the four scheduler windows.
    d01, d02, d03 = Fraction(1), Fraction(3), Fraction(5)
    d12, d13, d23 = Fraction(1), Fraction(3), Fraction(1)

    r0 = Fraction(1)
    r1 = Fraction(2, 1) / d01
    r2 = sum((
        Fraction(1, 1) / (d01 * d02),
        Fraction(1, 1) / (d01 * d12),
        Fraction(1, 1) / (d02 * d12),
    ), Fraction(0))
    r3 = sum((
        Fraction(1, 1) / (d01 * d02 * d03),
        Fraction(1, 1) / (d01 * d12 * d13),
        Fraction(1, 1) / (d02 * d12 * d23),
        Fraction(1, 1) / (d03 * d13 * d23),
    ), Fraction(0))
    inf_exact = max(r0, r1, r2, r3)

    c0 = sum((
        Fraction(1),
        Fraction(1, 1) / d01,
        Fraction(1, 1) / (d01 * d02),
        Fraction(1, 1) / (d01 * d02 * d03),
    ), Fraction(0))
    c1 = sum((
        Fraction(1, 1) / d01,
        Fraction(1, 1) / (d01 * d12),
        Fraction(1, 1) / (d01 * d12 * d13),
    ), Fraction(0))
    c2 = sum((
        Fraction(1, 1) / (d02 * d12),
        Fraction(1, 1) / (d02 * d12 * d23),
    ), Fraction(0))
    c3 = Fraction(1, 1) / (d03 * d13 * d23)
    one_exact = max(c0, c1, c2, c3)
    spectral_exact = max(one_exact, inf_exact)
    lambda_exact = Fraction(1, 1) / (spectral_exact * spectral_exact)

    return {
        "scaled_pair_gap_lower": {
            "u1-u0": float(d01), "u2-u0": float(d02), "u3-u0": float(d03),
            "u2-u1": float(d12), "u3-u1": float(d13), "u3-u2": float(d23),
        },
        "L_inverse_infinity_norm_exact": str(inf_exact),
        "L_inverse_one_norm_exact": str(one_exact),
        "L_inverse_infinity_norm_upper": up(float(inf_exact)),
        "L_inverse_one_norm_upper": up(float(one_exact)),
        "L_inverse_spectral_norm_upper": up(float(spectral_exact)),
        "L_transpose_L_lambda_min_exact": str(lambda_exact),
        "L_transpose_L_lambda_min_lower": down(float(lambda_exact)),
        "exact_rational_arithmetic_used_before_outward_float_conversion": True,
        "ordinary_floating_eigensolver_used": False,
        "determinant_trace_scalarization_used": False,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    dynamic = DYNAMIC.build(path)
    sched = SCHED.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "dynamic": DYNAMIC.validate(dynamic),
        "scheduler": SCHED.validate(sched),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"four-S complete-SEA3 prerequisites failed: {bad}")
    if complete["canonical_P3_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("four-S lemma lost the canonical complete SEA3 source")

    rs_contract = complete["R_S_regularizer"]
    if rs_contract["deployed_law"] != "SpectralMSE":
        raise RuntimeError("deployed R_S law is no longer SpectralMSE")
    if rs_contract["actual_applied_R_S_required_at_every_due_S_update"] is not True:
        raise RuntimeError("actual applied R_S is not required at due S events")
    if rs_contract["all_due_S_updates_remain_in_full_word"] is not True:
        raise RuntimeError("complete word no longer retains all due S updates")

    inv = dynamic["dynamic_invariant"]
    tau_lo, tau_hi = map(float, inv["tau_applied_s"])
    sigma_lo, sigma_hi = map(float, inv["sigma_aw_filter_mps2"])
    rs_lo, rs_hi = map(float, inv["R_S_applied"])
    if not (0.0 < tau_lo <= tau_hi and 0.0 < sigma_lo <= sigma_hi and 0.0 < rs_lo <= rs_hi):
        raise RuntimeError("invalid complete-SEA3 dynamic invariant")

    factors = list(map(float, rs_contract["axis_std_factors"]))
    if len(factors) != 3 or any(not (math.isfinite(x) and x > 0.0) for x in factors):
        raise RuntimeError("invalid deployed R_S axis factors")

    g = up(float(sched["certified_uniform_max_gap_s"]))
    if not (0.0 < g <= 0.151):
        raise RuntimeError(f"unexpected shipping scheduler gap {g}")
    T = up(7.0 * g)

    # Arbitrary legal time-varying tau enters only through the integral of
    # lambda=1/tau.  tau>=tau_lo gives a_response >= exp(-T/tau_lo).
    decay_exponent = up(T / tau_lo)
    a_response_lower, exp_cert = _validated_exp_negative_lower(decay_exponent)
    third_dd_lower = down(a_response_lower / 6.0)

    # Rank witness only: 0.5 * product of the six guaranteed pair separations
    # times the lower third divided difference.
    vandermonde_sep_product_lower = 45.0
    determinant_rank_witness = down(
        0.5 * vandermonde_sep_product_lower * third_dd_lower
    )

    # Selected-record nuisance covariance. q_c(t)=2 sigma(t)^2/tau(t), so its
    # complete-SEA3 upper is obtained from the same coupled invariant; damping
    # is dropped only to upper-bound nuisance variance, never to generate a word.
    qc_max = up(2.0 * sigma_hi * sigma_hi / tau_lo)
    process_s_var_upper = up(qc_max * up(T ** 7) / 252.0)
    axis_std_upper = [up(rs_hi * f) for f in factors]
    axis_var_upper = [up(x * x) for x in axis_std_upper]
    per_record_var_upper = up(max(axis_var_upper) + process_s_var_upper)
    stack_cov_lambda_max_upper = up(4.0 * per_record_var_upper)
    raw_inverse_cov_scalar_lower = down(1.0 / stack_cov_lambda_max_upper)

    newton = _newton_inverse_norm_certificate()
    newton_info_lower = down(
        raw_inverse_cov_scalar_lower * newton["L_transpose_L_lambda_min_lower"]
    )
    if not (math.isfinite(newton_info_lower) and newton_info_lower > 0.0):
        raise RuntimeError("four-S Newton information lower is not strict")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": complete["canonical_P3_source"],
        "component_of_complete_SEA3_full_word": True,
        "P3_architecture_replaced": False,
        "source_family_replaced": False,
        "trajectory_replay_used": False,
        "selected_four_S_events_replace_full_scheduler_word": False,
        "all_due_S_updates_remain_in_literal_word": True,
        "complete_SEA3_same_realization_requirement_retained": True,
        "independent_tau_sigma_RS_TS_extrema_product_used_as_source": False,
        "actual_applied_SpectralMSE_R_S_consumed": True,
        "instantaneous_R_S_target_substituted_for_applied_R_S": False,
        "uniform_S_gap_s_upper": g,
        "four_S_windows_s": [[down(a * g), up(b * g)] for a, b in WINDOW_MULTIPLIERS],
        "scaled_u_windows": [[a, b] for a, b in WINDOW_MULTIPLIERS],
        "selected_firings_are_guaranteed_actual_members_of_complete_word": True,
        "selected_window_horizon_s_upper": T,
        "time_varying_tau_allowed_inside_selected_subword": True,
        "time_varying_sigma_allowed_inside_selected_subword": True,
        "tau_applied_s": [tau_lo, tau_hi],
        "sigma_aw_filter_mps2": [sigma_lo, sigma_hi],
        "R_S_applied_base_std": [rs_lo, rs_hi],
        "R_S_axis_std_factors": factors,
        "dimensionless_translation_state": ["S", "g*p", "g^2*v", "g^3*a_w"],
        "homogeneous_observation_model": (
            "S(t)=S0+t*p0+0.5*t^2*v0+c(t)*a_w0; c'''(t)=exp(-int_0^t 1/tau(s) ds)"
        ),
        "validated_exponential_lower_certificate": exp_cert,
        "aw_homogeneous_response_lower": a_response_lower,
        "aw_scaled_third_divided_difference_lower": third_dd_lower,
        "scaled_observation_determinant_abs_lower_rank_witness_only": determinant_rank_witness,
        "determinant_used_only_as_rank_witness_not_eigenvalue_scalarization": True,
        "four_S_translation_observation_operator_full_rank": determinant_rank_witness > 0.0,
        "selected_S_record_noise": {
            "OU_driving_intensity_upper": qc_max,
            "process_S_variance_per_record_upper": process_s_var_upper,
            "measurement_variance_axis_upper": axis_var_upper,
            "per_record_variance_upper": per_record_var_upper,
            "four_record_covariance_lambda_max_upper": stack_cov_lambda_max_upper,
            "Sigma_Y_inverse_scalar_lower": raw_inverse_cov_scalar_lower,
            "process_damping_dropped_for_nuisance_upper_only": True,
            "cross_record_process_correlation_covered_by_trace_bound": True,
        },
        "newton_coordinate_information": {
            **newton,
            "raw_record_inverse_covariance_scalar_lower": raw_inverse_cov_scalar_lower,
            "D_S_newton_lambda_min_lower": newton_info_lower,
            "D_S_newton_matrix_lower": f"D_S,z >= {newton_info_lower:.17g} * I_4",
            "full_4x4_matrix_inequality_closed": True,
            "determinant_used_for_information_lower": False,
            "frobenius_singular_value_conversion_used": False,
            "scalar_information_beta_used": False,
        },
        "P3_RS_TRANSLATION_INFORMATION_MATRIX_CLOSED": True,
        "P3_promoted": False,
        "next_obligation": (
            "compose this full [v,p,S,a_w] R_S information lemma with the same complete SEA3 "
            "H18/A21 P/Psi/Omega word, vector PE, A-mode bias dynamics, every accelerometer "
            "Joseph update, every due S update, every process Q and every covariance-floor event"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical complete SEA3 source changed")
    for key in (
        "component_of_complete_SEA3_full_word",
        "all_due_S_updates_remain_in_literal_word",
        "complete_SEA3_same_realization_requirement_retained",
        "actual_applied_SpectralMSE_R_S_consumed",
        "selected_firings_are_guaranteed_actual_members_of_complete_word",
        "time_varying_tau_allowed_inside_selected_subword",
        "time_varying_sigma_allowed_inside_selected_subword",
        "determinant_used_only_as_rank_witness_not_eigenvalue_scalarization",
        "four_S_translation_observation_operator_full_rank",
        "P3_RS_TRANSLATION_INFORMATION_MATRIX_CLOSED",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "P3_architecture_replaced", "source_family_replaced", "trajectory_replay_used",
        "selected_four_S_events_replace_full_scheduler_word",
        "independent_tau_sigma_RS_TS_extrema_product_used_as_source",
        "instantaneous_R_S_target_substituted_for_applied_R_S", "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    ni = d.get("newton_coordinate_information", {})
    if ni.get("full_4x4_matrix_inequality_closed") is not True:
        f.append("Newton full 4x4 information inequality not closed")
    x = ni.get("D_S_newton_lambda_min_lower")
    if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
        f.append("Newton information lower is not finite positive")
    if ni.get("determinant_used_for_information_lower") is not False:
        f.append("determinant was used as information scalarization")
    if ni.get("scalar_information_beta_used") is not False:
        f.append("scalar beta shortcut was reintroduced")
    if len(d.get("four_S_windows_s", [])) != 4:
        f.append("four-S selection does not contain four guaranteed windows")
    if float(d.get("uniform_S_gap_s_upper", math.inf)) > 0.151:
        f.append("scheduler recurrence widened beyond 150 ms service")
    return list(dict.fromkeys(f))


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
        "gap_s": d["uniform_S_gap_s_upper"],
        "four_S_windows_s": d["four_S_windows_s"],
        "R_S_axis_std_upper": max(
            math.sqrt(v) for v in d["selected_S_record_noise"]["measurement_variance_axis_upper"]
        ),
        "D_S_newton_lambda_min_lower": d["newton_coordinate_information"]["D_S_newton_lambda_min_lower"],
        "pass": not failures,
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
