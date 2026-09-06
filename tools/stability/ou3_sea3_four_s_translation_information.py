#!/usr/bin/env python3
"""Full four-S translation information inside the complete SEA3 P3 word.

This module is a mandatory lemma consumed by the canonical
``COMPLETE_SEA3_NORMAL_LIVE_WORD``.  It is not a source generator, a reduced
P3 word, or a promotion route by itself.

The shipping progress-preserving S=0 scheduler guarantees one *actual* firing
in every interval of length ``g``.  Select one firing from each disjoint window

    [0,g], [2g,3g], [4g,5g], [6g,7g].

Every other due S firing remains in the literal complete SEA3 word.  For one
translation axis and scaled time u=t/g, the homogeneous S observation is

    y(u) = S0 + u (g p0) + u^2/2 (g^2 v0) + cbar(u) (g^3 a0),
    cbar'''(u) = exp(-int_0^{g u} 1/tau(s) ds) > 0.

The quantitative certificate must be in the *physical dimensionless state*

    x = [S, g p, g^2 v, g^3 a_w],

not merely in the four Newton divided differences of y.  Let q0..q3 be the
nested divided differences of the four selected records.  Then

    q0 = S + u0 P + u0^2/2 V + c0 A,
    q1 =     P + (u0+u1)/2 V + c01 A,
    q2 =                         1/2 V + c012 A,
    q3 =                                   c0123 A.

The generalized mean-value theorem and cbar''' in [a_-,1] give

    c0123 >= a_-/6 > 0,
    c0 <= 1/6, c01 <= 9/2, c012 <= 5/2

on the four scheduler windows.  We recover A,V,P,S successively and bound the
row l1 norms of the *full raw-record inverse* M^-1.  Therefore

    ||M^-1||_2^2 <= ||M^-1||_F^2
                  <= sum_i ||row_i(M^-1)||_1^2,

which yields a genuine physical-state matrix bound

    M^T M >= m_phys I4.

This is where the third divided-difference lower enters the quantitative
certificate.  Earlier versions conditioned only the raw-to-Newton evaluation
matrix and used c0123 merely as a rank witness; that omitted the map from
Newton coordinates back to the physical a_w coordinate and overstated the
information lower.

The observation-noise bound uses the actual applied SpectralMSE R_S safety
ceiling and deployed axis factors.  Applied R_S is never replaced by its target.
Process noise over the selected window is included as nuisance with arbitrary
legal SEA3 tau/sigma variation.  The determinant remains a rank witness only;
it is not converted into an eigenvalue bound.
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

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 2
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
    """Exact absolute-row bounds for raw records -> divided differences."""
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
        "divided_difference_raw_row_l1_exact": [
            str(r0), str(r1), str(r2), str(r3)
        ],
        "divided_difference_raw_row_l1_upper": [
            up(float(r0)), up(float(r1)), up(float(r2)), up(float(r3))
        ],
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


def _physical_state_recovery_certificate(newton: dict, third_dd_lower: float) -> dict:
    """Bound the full physical-state inverse, retaining the c''' scale.

    q=T_dd y are nested divided differences.  On the selected windows
    u0<=1,u1<=3,u2<=5 and 0<c'''<=1.  From zero initial c,c',c'' and the
    generalized mean-value theorem:

      c(u0)<=1/6, c[u0,u1]<=9/2, c[u0,u1,u2]<=5/2.

    Recover A,V,P,S from q3,q2,q1,q0 and propagate raw-record row-l1 bounds.
    This is a matrix inverse norm certificate, not determinant/trace conversion.
    """
    if not (math.isfinite(third_dd_lower) and third_dd_lower > 0.0):
        raise ValueError("strict third divided-difference lower required")
    q0, q1, q2, q3 = map(
        float, newton["divided_difference_raw_row_l1_upper"]
    )

    u0_upper = 1.0
    half_u0_sq_upper = 0.5
    half_u0_plus_u1_upper = 2.0
    c0_upper = up(1.0 / 6.0)
    c01_upper = up(9.0 / 2.0)
    c012_upper = up(5.0 / 2.0)

    row_A = up(q3 / third_dd_lower)
    row_V = up(2.0 * up(q2 + up(c012_upper * row_A)))
    row_P = up(q1 + up(half_u0_plus_u1_upper * row_V) + up(c01_upper * row_A))
    row_S = up(
        q0
        + up(u0_upper * row_P)
        + up(half_u0_sq_upper * row_V)
        + up(c0_upper * row_A)
    )
    rows = [row_S, row_P, row_V, row_A]
    inverse_frobenius_sq_upper = 0.0
    for r in rows:
        inverse_frobenius_sq_upper = up(
            inverse_frobenius_sq_upper + up(r * r)
        )
    gram_lambda_min_lower = down(1.0 / inverse_frobenius_sq_upper)
    if not (math.isfinite(gram_lambda_min_lower) and gram_lambda_min_lower > 0.0):
        raise RuntimeError("physical four-S observation matrix lower is not strict")

    return {
        "state_order": ["S", "g*p", "g^2*v", "g^3*a_w"],
        "third_divided_difference_lower": third_dd_lower,
        "third_divided_difference_enters_quantitative_inverse": True,
        "derivative_upper_bounds": {
            "c_u0": c0_upper,
            "c_first_divided_difference": c01_upper,
            "c_second_divided_difference": c012_upper,
            "basis_half_u0_plus_u1": half_u0_plus_u1_upper,
            "basis_u0": u0_upper,
            "basis_half_u0_squared": half_u0_sq_upper,
        },
        "physical_state_inverse_raw_record_row_l1_upper": {
            "S": row_S, "g*p": row_P, "g^2*v": row_V, "g^3*a_w": row_A,
        },
        "physical_state_inverse_frobenius_squared_upper": inverse_frobenius_sq_upper,
        "physical_observation_MtM_lambda_min_lower": gram_lambda_min_lower,
        "inverse_matrix_frobenius_norm_bound_used": True,
        "determinant_used_for_quantitative_bound": False,
        "raw_to_Newton_condition_alone_used_as_physical_bound": False,
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

    # Arbitrary legal time-varying tau enters through the integral of
    # lambda=1/tau.  tau>=tau_lo gives c''' >= exp(-T/tau_lo).
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
    physical = _physical_state_recovery_certificate(newton, third_dd_lower)
    physical_info_lower = down(
        raw_inverse_cov_scalar_lower
        * physical["physical_observation_MtM_lambda_min_lower"]
    )
    if not (math.isfinite(physical_info_lower) and physical_info_lower > 0.0):
        raise RuntimeError("four-S physical-state information lower is not strict")

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
            "physical_state_recovery": physical,
            "raw_record_inverse_covariance_scalar_lower": raw_inverse_cov_scalar_lower,
            "third_divided_difference_lower_enters_quantitative_bound": True,
            "D_S_physical_lambda_min_lower": physical_info_lower,
            # Compatibility name retained for downstream consumers.  It now
            # denotes the corrected physical-state bound, not the old Newton-
            # coordinate-only bound.
            "D_S_newton_lambda_min_lower": physical_info_lower,
            "D_S_newton_matrix_lower": f"D_S,z >= {physical_info_lower:.17g} * I_4",
            "full_4x4_matrix_inequality_closed": True,
            "determinant_used_for_information_lower": False,
            "frobenius_singular_value_conversion_used": False,
            "inverse_matrix_frobenius_norm_bound_used": True,
            "raw_to_Newton_condition_alone_used_as_physical_bound": False,
            "scalar_information_beta_used": False,
        },
        "P3_RS_TRANSLATION_INFORMATION_MATRIX_CLOSED": True,
        "P3_promoted": False,
        "next_obligation": (
            "compose this corrected physical [S,gp,g^2v,g^3a_w] R_S information lemma with "
            "the same complete SEA3 H18/A21 P/Psi/Omega word, vector PE, A-mode bias dynamics, "
            "every accelerometer Joseph update, every due S update, every process Q and every "
            "covariance-floor event"
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
    if ni.get("third_divided_difference_lower_enters_quantitative_bound") is not True:
        f.append("physical a_w divided-difference scale missing from quantitative bound")
    if ni.get("raw_to_Newton_condition_alone_used_as_physical_bound") is not False:
        f.append("raw-to-Newton condition was incorrectly used as physical-state bound")
    physical = ni.get("physical_state_recovery", {})
    if physical.get("third_divided_difference_enters_quantitative_inverse") is not True:
        f.append("third divided difference missing from physical inverse")
    if physical.get("determinant_used_for_quantitative_bound") is not False:
        f.append("determinant entered quantitative physical-state bound")
    for key in ("D_S_physical_lambda_min_lower", "D_S_newton_lambda_min_lower"):
        x = ni.get(key)
        if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
            f.append(f"{key} is not finite positive")
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
    ni = d["newton_coordinate_information"]
    print(json.dumps({
        "gap_s": d["uniform_S_gap_s_upper"],
        "four_S_windows_s": d["four_S_windows_s"],
        "R_S_axis_std_upper": max(
            math.sqrt(v) for v in d["selected_S_record_noise"]["measurement_variance_axis_upper"]
        ),
        "third_divided_difference_lower": d["aw_scaled_third_divided_difference_lower"],
        "physical_MtM_lambda_min_lower": ni["physical_state_recovery"]["physical_observation_MtM_lambda_min_lower"],
        "D_S_physical_lambda_min_lower": ni["D_S_physical_lambda_min_lower"],
        "pass": not failures,
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
