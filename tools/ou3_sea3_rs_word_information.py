#!/usr/bin/env python3
"""SEA3 recurrent R_S word certificate for the OU-III translation block.

This producer closes the source-history-free *geometry* and nuisance-covariance
part of the canonical P3 correction word.  It deliberately does not recurse a
Riccati covariance lower and it does not replace the applied R_S state by the
instantaneous SpectralMSE target.

Let g be the certified source-independent maximum S=0 recurrence gap.  Select
one actual pseudo update from each disjoint window

    [0,g], [2g,3g], [4g,5g], [6g,7g].

For one translation axis x=[v,p,S,a_w] at the word start, the homogeneous
S-observation at time t is

    S(t) = S + t p + t^2 v/2 + c(t) a_w,

where for arbitrary time-varying lambda(t)=1/tau(t)

    c'''(t) = exp(-int_0^t lambda(s) ds) > 0.

Scale u=t/g and cbar(u)=c(gu)/g^3.  Then cbar'''(u)=c'''(gu).  The third divided
difference over four distinct u values is therefore bounded below by

    cbar[u1,u2,u3,u4] >= exp(-7g/tau_min)/6.

The polynomial columns [1,u,u^2/2] have zero third divided difference.  Hence
the four S rows are full rank on [S,g p,g^2 v,g^3 a_w] for *every* legal
source path, without freezing tau or enumerating source predecessors.

The same four selected measurements have a source-uniform nuisance covariance
upper.  Measurement variance uses the actually applied R_S safety ceiling; OU
process disturbance is bounded by dropping damping, so

    Var[Delta S_process(t)] <= q_c,max t^7/252,
    q_c,max = 2 sigma_max^2/tau_min.

Trace then gives a safe lambda_max bound on the complete 4x4 selected S-record
noise covariance.  This is intentionally conservative but preserves the
important structure: recurrent R_S correction itself makes translation
observable.  The final H18/A21 P3 gate will combine this observation-space
certificate with a finite-memory covariance lower in the same coordinates.

The deployed default R_S law is SpectralMSE and its target is evaluated at the
same target tau/sigma/T_S point.  Applied R_S, however, has its own EMA.  This
producer records that distinction explicitly: target-law similarity is a valid
source fact, but it is not substituted for the applied R_S state until a
separate lag/reachability theorem proves that tightening.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from ou3_interval import down, up
import ou3_p3_pseudo_scheduler_progress_certificate as SCHED
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_RS_FOUR_S_TRANSLATION_WORD_V1"
WINDOW_MULTIPLIERS = ((0.0, 1.0), (2.0, 3.0), (4.0, 5.0), (6.0, 7.0))


def _axis_factors(text: str) -> list[float]:
    out = []
    for name in ("R_S_x_factor_", "R_S_y_factor_"):
        m = re.search(rf"float\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f", text)
        if not m:
            raise RuntimeError(f"cannot extract {name}")
        out.append(float(m.group(1)))
    out.append(1.0)  # z uses the base standard deviation directly
    if any(not (math.isfinite(x) and x > 0.0) for x in out):
        raise RuntimeError("invalid R_S axis factors")
    return out


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

    # lambda(t)<=1/tau_min.  Use the validated exponential backend so the rank
    # witness does not depend on ordinary libm rounding.
    decay_exponent = up(T / tau_lo)
    decay = VT.exp_point(-decay_exponent)
    a_response_lower = down(decay.lo)
    third_dd_lower = down(a_response_lower / 6.0)

    # In u=t/g, the four windows are [0,1],[2,3],[4,5],[6,7].  Their pairwise
    # minimum separations are 1,3,5,1,3,1, so the Vandermonde product is >=45.
    # The determinant is recorded only as a full-rank witness.  It is NOT used
    # as a determinant/trace eigenvalue scalarization in the P3 quantitative gate.
    vandermonde_sep_product_lower = 45.0
    scaled_observation_det_lower = down(
        0.5 * vandermonde_sep_product_lower * third_dd_lower
    )

    # Selected S-record nuisance covariance.  The process-noise upper drops OU
    # damping, hence remains valid for arbitrary legal time-varying tau/sigma.
    qc_max = up(2.0 * sigma_hi * sigma_hi / tau_lo)
    T7 = T ** 7
    process_S_var_upper = up(qc_max * T7 / 252.0)
    axis_std_upper = [up(rs_hi * f) for f in factors]
    axis_var_upper = [up(x * x) for x in axis_std_upper]
    worst_measurement_var = max(axis_var_upper)
    per_record_var_upper = up(worst_measurement_var + process_S_var_upper)
    stack_lambda_max_upper = up(4.0 * per_record_var_upper)
    stack_information_scalar_lower = down(1.0 / stack_lambda_max_upper)

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
        "SEA3_target_vs_applied_RS_contract": {
            "SpectralMSE_target_same_tau_sigma_TS_point": True,
            "applied_RS_has_separate_EMA": True,
            "instantaneous_target_formula_substituted_for_applied_RS": False,
            "active_applied_RS_safety_ceiling_used_here": True,
            "future_tightening_requires_explicit_lag_or_reachability_theorem": True,
        },
        "P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED": True,
        "P3_RS_BATCH_NOISE_UPPER_CLOSED": True,
        "P3_RS_TRANSLATION_INFORMATION_MATRIX_CLOSED": False,
        "P3_PROMOTED": False,
        "next_obligation": (
            "construct the finite-memory UCC covariance lower in the same four-S divided-difference/observation coordinates and compare it directly with this S-record noise bound; then add the independent vector-PE attitude/gyro-bias block and A-mode bias coupling"
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
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "old_P2_800_state_graph_consumed", "source_history_graph_consumed",
        "predecessor_path_enumeration_consumed",
        "P3_RS_TRANSLATION_INFORMATION_MATRIX_CLOSED", "P3_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if not all(d.get("source_parity", {}).values()):
        f.append("shipping R_S source parity failed")
    if len(d.get("four_S_windows", [])) != 4:
        f.append("four separated S windows missing")
    if float(d.get("aw_scaled_third_divided_difference_lower", 0.0)) <= 0.0:
        f.append("a_w divided-difference rank floor is not strict")
    if float(d.get("scaled_observation_determinant_abs_lower_rank_witness_only", 0.0)) <= 0.0:
        f.append("four-S observation rank witness is not strict")
    noise = d.get("selected_S_record_noise", {})
    if float(noise.get("four_record_covariance_lambda_max_trace_upper", 0.0)) <= 0.0:
        f.append("selected S-record covariance upper is invalid")
    if float(noise.get("Sigma_S_inverse_scalar_lower", 0.0)) <= 0.0:
        f.append("selected S-record information scalar is invalid")
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
    print(json.dumps({
        "four_S_full_rank": d["four_S_translation_observation_operator_full_rank"],
        "word_horizon_s_upper": d["word_horizon_s_upper"],
        "third_divided_difference_lower": d["aw_scaled_third_divided_difference_lower"],
        "rank_witness_det_lower": d["scaled_observation_determinant_abs_lower_rank_witness_only"],
        "S_record_noise_lambda_max_upper": d["selected_S_record_noise"]["four_record_covariance_lambda_max_trace_upper"],
        "S_record_information_scalar_lower": d["selected_S_record_noise"]["Sigma_S_inverse_scalar_lower"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
