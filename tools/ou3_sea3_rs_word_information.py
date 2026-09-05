#!/usr/bin/env python3
"""SEA3 recurrent R_S four-S translation information certificate.

This is one mandatory component of the complete Normal-Live H18/A21 P3 word.
It does not promote P3 by itself.

Let g be the certified source-independent maximum S=0 recurrence gap. Select
one actual pseudo update from each disjoint window

    [0,g], [2g,3g], [4g,5g], [6g,7g].

For one translation axis and scaled time u=t/g,

    S(u) = S0 + u (g p0) + u^2/2 (g^2 v0) + cbar(u) (g^3 a0),
    cbar'''(u) = exp(-int_0^{g u} 1/tau(s) ds) > 0.

The third divided difference therefore has a source-uniform positive lower for
arbitrary legal time-varying tau.  The determinant is retained only as a rank
witness.

The stronger quantitative step in this version uses Newton divided-difference
coordinates instead of determinant/Frobenius scalarization.  If z is the vector
of nested divided differences of the four homogeneous S values, the raw record
vector is

    y = L(u0,u1,u2,u3) z,

where L is the lower-triangular Newton evaluation matrix.  Its inverse is the
ordinary divided-difference table.  The guaranteed timing gaps imply validated
bounds

    ||L^-1||_inf <= 2,
    ||L^-1||_1   <= 12/5,

hence ||L^-1||_2 <= 12/5 and

    L^T L >= (25/144) I.

The selected raw-record covariance satisfies Sigma_Y <= lambda_Y I from the
already certified applied-R_S/process upper, so in Newton coordinates

    D_S,z = L^T Sigma_Y^-1 L
          >= (1/lambda_Y) L^T L
          >= d_N I,

with d_N>0.  This is a genuine 4x4 matrix inequality in the divided-difference
coordinates; no determinant/trace eigenvalue conversion and no scalar
information-beta attenuation is used.

Applied R_S remains distinct from the instantaneous SpectralMSE target.  The
current safe record-noise upper uses the applied R_S safety ceiling; the final
complete word must carry the actual joint adaptive source path and may tighten
this component, but may not replace applied R_S by its target.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
import re
from pathlib import Path

from ou3_interval import Interval, down, up
import ou3_p3_pseudo_scheduler_progress_certificate as SCHED
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_SEA3_RS_FOUR_S_TRANSLATION_WORD_V2"
WINDOW_MULTIPLIERS = ((0.0, 1.0), (2.0, 3.0), (4.0, 5.0), (6.0, 7.0))


def _axis_factors(text: str) -> list[float]:
    out: list[float] = []
    for name in ("R_S_x_factor_", "R_S_y_factor_"):
        m = re.search(rf"float\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f", text)
        if not m:
            raise RuntimeError(f"cannot extract {name}")
        out.append(float(m.group(1)))
    out.append(1.0)
    if any(not (math.isfinite(x) and x > 0.0) for x in out):
        raise RuntimeError("invalid R_S axis factors")
    return out


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
    n = max(1, int(math.ceil(2.0 * x)))
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
    """Bound L^-1 from the guaranteed scaled timing gaps.

    The separated scaled windows give exact rational gap floors

      d01>=1, d02>=3, d03>=5, d12>=1, d13>=3, d23>=1.

    The absolute row/column sums of the divided-difference inverse are therefore
    evaluated in exact rational arithmetic first.  We outward-round only once
    when converting the exact rational bounds to binary64.  This avoids turning
    the exact 2 and 12/5 bounds into larger numbers through repeated floating
    outward additions.
    """
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
    inf_norm = up(float(inf_exact))
    one_norm = up(float(one_exact))
    spectral_norm_upper = up(float(spectral_exact))
    lambda_exact = Fraction(1, 1) / (spectral_exact * spectral_exact)
    lambda_LtL_lower = down(float(lambda_exact))
    if not lambda_LtL_lower > 0.0:
        raise RuntimeError("Newton evaluation matrix lost a strict conditioning lower")

    return {
        "scaled_pair_gap_lower": {
            "u1-u0": float(d01),
            "u2-u0": float(d02),
            "u3-u0": float(d03),
            "u2-u1": float(d12),
            "u3-u1": float(d13),
            "u3-u2": float(d23),
        },
        "L_inverse_infinity_norm_exact": str(inf_exact),
        "L_inverse_one_norm_exact": str(one_exact),
        "L_inverse_infinity_norm_upper": inf_norm,
        "L_inverse_one_norm_upper": one_norm,
        "L_inverse_spectral_norm_upper": spectral_norm_upper,
        "L_transpose_L_lambda_min_exact": str(lambda_exact),
        "L_transpose_L_lambda_min_lower": lambda_LtL_lower,
        "derivation": "exact rational divided-difference row/column sums; lambda_min(L^T L)=1/||L^-1||_2^2; ||.||_2<=max(||.||_1,||.||_inf)",
        "exact_rational_arithmetic_used_before_outward_float_conversion": True,
        "ordinary_floating_eigensolver_used": False,
        "determinant_trace_scalarization_used": False,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("R_S word certificate may not be trajectory fitted")

    dynamic = DYNAMIC.build(path)
    df = DYNAMIC.validate(dynamic)
    if df:
        raise RuntimeError(f"dynamic SEA3 source invalid: {df}")
    sched = SCHED.build(path)
    sf = SCHED.validate(sched)
    if sf:
        raise RuntimeError(f"pseudo scheduler certificate invalid: {sf}")

    text = WRAPPER.read_text(encoding="utf-8")
    factors = _axis_factors(text)
    source_parity = {
        "spectral_mse_is_deployed_default": (
            "RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;" in text
        ),
        "spectral_mse_uses_realized_pseudo_period": (
            "const float TS = pseudo_update_period_for_(tau);" in text
            and "rs_mse_coeff_ * rs_qeff_pow_" in text
            and "/ std::sqrt(TS)" in text
        ),
        "pseudo_period_is_tau_scaled": (
            "pseudo_update_tau_ratio_ * tau" in text
            and "pseudo_update_period_max_s_" in text
        ),
        "applied_tau_and_RS_have_distinct_emas": (
            "tune_.tau_applied   += alpha" in text
            and "tune_.RS_applied    += alpha_RS" in text
        ),
        "RS_is_applied_as_measurement_standard_deviation": (
            "mekf_->set_RS_noise(Eigen::Vector3f(" in text
        ),
    }
    parity_failures = [k for k, v in source_parity.items() if not v]
    if parity_failures:
        raise RuntimeError(f"R_S shipping-source parity failed: {parity_failures}")

    inv = dynamic["dynamic_invariant"]
    tau_lo, tau_hi = map(float, inv["tau_applied_s"])
    sigma_lo, sigma_hi = map(float, inv["sigma_aw_filter_mps2"])
    rs_lo, rs_hi = map(float, inv["R_S_applied"])
    if not (0.0 < tau_lo <= tau_hi and 0.0 < sigma_lo <= sigma_hi and 0.0 < rs_lo <= rs_hi):
        raise RuntimeError("invalid dynamic source invariant")

    g = float(sched["certified_uniform_max_gap_s"])
    if not (math.isfinite(g) and g > 0.0):
        raise RuntimeError("invalid certified S recurrence gap")
    T = up(7.0 * g)

    decay_exponent = up(T / tau_lo)
    a_response_lower, exp_certificate = _validated_exp_negative_lower(decay_exponent)
    third_dd_lower = down(a_response_lower / 6.0)

    # Rank witness only.  It is no longer used to turn the observation matrix
    # into a smallest-eigenvalue lower.
    vandermonde_sep_product_lower = 45.0
    scaled_observation_det_lower = down(
        0.5 * vandermonde_sep_product_lower * third_dd_lower
    )

    qc_max = up(2.0 * sigma_hi * sigma_hi / tau_lo)
    T7 = up(T ** 7)
    process_S_var_upper = up(qc_max * T7 / 252.0)
    axis_std_upper = [up(rs_hi * f) for f in factors]
    axis_var_upper = [up(x * x) for x in axis_std_upper]
    worst_measurement_var = max(axis_var_upper)
    per_record_var_upper = up(worst_measurement_var + process_S_var_upper)
    stack_lambda_max_upper = up(4.0 * per_record_var_upper)
    stack_information_scalar_lower = down(1.0 / stack_lambda_max_upper)

    newton = _newton_inverse_norm_certificate()
    newton_info_lower = down(
        stack_information_scalar_lower
        * float(newton["L_transpose_L_lambda_min_lower"])
    )
    if not (math.isfinite(newton_info_lower) and newton_info_lower > 0.0):
        raise RuntimeError("four-S Newton-coordinate information lower is not strict")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "old_P2_800_state_graph_consumed": False,
        "source_history_graph_consumed": False,
        "predecessor_path_enumeration_consumed": False,
        "source_parity": source_parity,
        "four_S_windows": [[down(a * g), up(b * g)] for a, b in WINDOW_MULTIPLIERS],
        "scaled_u_windows": [[a, b] for a, b in WINDOW_MULTIPLIERS],
        "uniform_S_gap_s": g,
        "word_horizon_s_upper": T,
        "tau_applied_s": [tau_lo, tau_hi],
        "sigma_aw_filter_mps2": [sigma_lo, sigma_hi],
        "R_S_applied_base_std": [rs_lo, rs_hi],
        "R_S_axis_std_factors": factors,
        "time_varying_tau_allowed_inside_word": True,
        "time_varying_sigma_allowed_inside_word": True,
        "four_S_observation_model": (
            "S(t)=S0+t*p0+0.5*t^2*v0+c(t)*a_w0, c'''(t)=exp(-int_0^t 1/tau(s) ds)"
        ),
        "dimensionless_state": ["S", "g*p", "g^2*v", "g^3*a_w"],
        "newton_divided_difference_state": [
            "S(u0)", "S[u0,u1]", "S[u0,u1,u2]", "S[u0,u1,u2,u3]"
        ],
        "aw_decay_exponent_upper": decay_exponent,
        "validated_exponential_lower_certificate": exp_certificate,
        "aw_homogeneous_response_lower": a_response_lower,
        "aw_scaled_third_divided_difference_lower": third_dd_lower,
        "scaled_time_vandermonde_separation_product_lower": vandermonde_sep_product_lower,
        "scaled_observation_determinant_abs_lower_rank_witness_only": scaled_observation_det_lower,
        "determinant_used_only_as_rank_witness_not_eigenvalue_scalarization": True,
        "four_S_translation_observation_operator_full_rank": scaled_observation_det_lower > 0.0,
        "selected_S_record_noise": {
            "q_c_upper_m2_s5": qc_max,
            "process_S_variance_per_record_upper": process_S_var_upper,
            "measurement_variance_axis_upper": axis_var_upper,
            "per_record_variance_upper": per_record_var_upper,
            "four_record_covariance_lambda_max_trace_upper": stack_lambda_max_upper,
            "Sigma_S_inverse_scalar_lower": stack_information_scalar_lower,
            "process_damping_dropped_for_upper_bound": True,
            "cross_record_process_correlation_covered_by_trace_bound": True,
        },
        "newton_coordinate_information": {
            **newton,
            "raw_record_inverse_covariance_scalar_lower": stack_information_scalar_lower,
            "D_S_newton_matrix_lower": (
                f"D_S,z >= {newton_info_lower:.17g} * I_4"
            ),
            "D_S_newton_lambda_min_lower": newton_info_lower,
            "full_4x4_matrix_inequality_closed": True,
            "determinant_used_for_information_lower": False,
            "frobenius_singular_value_conversion_used": False,
            "scalar_information_beta_used": False,
        },
        "SEA3_target_vs_applied_RS_contract": {
            "SpectralMSE_target_same_tau_sigma_TS_point": True,
            "applied_RS_has_separate_EMA": True,
            "instantaneous_target_formula_substituted_for_applied_RS": False,
            "active_applied_RS_safety_ceiling_used_here": True,
            "complete_word_still_requires_same_joint_source_path": True,
        },
        "P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED": True,
        "P3_RS_BATCH_NOISE_UPPER_CLOSED": True,
        "P3_RS_TRANSLATION_INFORMATION_MATRIX_CLOSED": True,
        "P3_PROMOTED": False,
        "next_obligation": (
            "include this four-S matrix component in the complete source-reachable joint P/Psi/Omega H18/A21 Normal-Live word; do not form a separate information/covariance final ratio"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "time_varying_tau_allowed_inside_word",
        "time_varying_sigma_allowed_inside_word",
        "determinant_used_only_as_rank_witness_not_eigenvalue_scalarization",
        "four_S_translation_observation_operator_full_rank",
        "P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED",
        "P3_RS_BATCH_NOISE_UPPER_CLOSED",
        "P3_RS_TRANSLATION_INFORMATION_MATRIX_CLOSED",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "old_P2_800_state_graph_consumed", "source_history_graph_consumed",
        "predecessor_path_enumeration_consumed", "P3_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if not all(d.get("source_parity", {}).values()):
        f.append("shipping R_S source parity failed")
    if len(d.get("four_S_windows", [])) != 4:
        f.append("four separated S windows missing")
    expc = d.get("validated_exponential_lower_certificate", {})
    if expc.get("method") != "VALIDATED_EXP_POINT_PRODUCT":
        f.append("wide exponential did not use validated product construction")
    if expc.get("ordinary_libm_exp_used") is not False:
        f.append("ordinary libm exponential entered the rank proof")
    if float(expc.get("piece_exponent_upper", math.inf)) > VT.MAX_ABS_ARGUMENT:
        f.append("validated exponential product used an unaudited piece")
    if float(d.get("aw_scaled_third_divided_difference_lower", 0.0)) <= 0.0:
        f.append("a_w divided-difference rank floor is not strict")
    if float(d.get("scaled_observation_determinant_abs_lower_rank_witness_only", 0.0)) <= 0.0:
        f.append("four-S observation rank witness is not strict")
    noise = d.get("selected_S_record_noise", {})
    if float(noise.get("four_record_covariance_lambda_max_trace_upper", 0.0)) <= 0.0:
        f.append("selected S-record covariance upper is invalid")
    if float(noise.get("Sigma_S_inverse_scalar_lower", 0.0)) <= 0.0:
        f.append("selected S-record information scalar is invalid")
    ni = d.get("newton_coordinate_information", {})
    if ni.get("full_4x4_matrix_inequality_closed") is not True:
        f.append("Newton-coordinate information matrix did not close")
    if ni.get("exact_rational_arithmetic_used_before_outward_float_conversion") is not True:
        f.append("Newton norm certificate did not use exact rational arithmetic")
    if ni.get("determinant_used_for_information_lower") is not False:
        f.append("determinant re-entered information conditioning")
    if ni.get("frobenius_singular_value_conversion_used") is not False:
        f.append("Frobenius singular-value conversion re-entered information conditioning")
    if ni.get("scalar_information_beta_used") is not False:
        f.append("scalar information beta re-entered four-S component")
    if ni.get("L_inverse_one_norm_exact") != "12/5":
        f.append("Newton inverse exact one-norm bound changed")
    if ni.get("L_inverse_infinity_norm_exact") != "2":
        f.append("Newton inverse exact infinity-norm bound changed")
    if float(ni.get("L_inverse_one_norm_upper", math.inf)) > up(float(Fraction(12, 5))):
        f.append("Newton inverse one-norm outward bound weakened unexpectedly")
    if float(ni.get("L_inverse_infinity_norm_upper", math.inf)) > up(float(Fraction(2, 1))):
        f.append("Newton inverse infinity-norm outward bound weakened unexpectedly")
    if float(ni.get("D_S_newton_lambda_min_lower", 0.0)) <= 0.0:
        f.append("Newton-coordinate information lower is not strict")
    lag = d.get("SEA3_target_vs_applied_RS_contract", {})
    if lag.get("SpectralMSE_target_same_tau_sigma_TS_point") is not True:
        f.append("SpectralMSE target coupling disappeared")
    if lag.get("applied_RS_has_separate_EMA") is not True:
        f.append("applied R_S lag distinction disappeared")
    if lag.get("instantaneous_target_formula_substituted_for_applied_RS") is not False:
        f.append("instantaneous R_S target was incorrectly substituted for applied R_S")
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
        "qualification": d["qualification"],
        "word_horizon_s": d["word_horizon_s_upper"],
        "aw_dd3_lower": d["aw_scaled_third_divided_difference_lower"],
        "Linv_1_upper": ni["L_inverse_one_norm_upper"],
        "Linv_inf_upper": ni["L_inverse_infinity_norm_upper"],
        "LtL_lambda_lower": ni["L_transpose_L_lambda_min_lower"],
        "D_S_newton_lambda_lower": ni["D_S_newton_lambda_min_lower"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())