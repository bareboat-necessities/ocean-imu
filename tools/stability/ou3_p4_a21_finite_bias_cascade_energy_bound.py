#!/usr/bin/env python3
"""Quantitative complete-SEA3 finite-bias comparison-cascade energy bound.

This closes the *comparison observer* C_Hb bound left open by
``ou3_p4_a21_finite_bias_cascade_bridge`` without using the detectability
helper's deliberately crude ``N*K*L^N`` estimate.

Choose the H18 comparison observer to use the certified complete-SEA3 H18
shipping information geometry and a zero accelerometer-bias correction row.
For h_0=0 the only direct b_a forcing of that H18 observer is the accepted
accelerometer residual.  The exact Joseph signed-energy identity gives

    Delta V_H <= s_H b_j^T R_a^-1 b_j,

while S=0 and vector events have no linear b_a residual input and tangent reset
congruences preserve the H18 information energy.  The comparison bias state is
not corrected and follows the exact shipping homogeneous GM law, so

    ||b_j|| <= exp(-j h_-/tau_b) ||b_0||.

Hence the whole 3 s bias-to-H coupling obeys the closed-form energy bound

    ||C_Hb b_0||_{M_H,+}^2
      <= s_H lambda_max(R_a^-1)
         sum_{j=1}^N exp(-2 j h_-/tau_b) ||b_0||^2.

The geometric sum is the actual finite-tau_b input energy over the required
accepted accelerometer sequence.  It is not an N-times worst nonlinear
remainder and contains no product of suffix norms.

Combining this c^2 with the exact two-block Schur condition produces a finite
full-rank quadratic metric for the *triangular comparison observer*.  This is a
stronger quantitative detectability/UES bridge, but it is deliberately NOT yet
a P4 metric for the actual A21 Kalman map: the actual shipping bias correction
row and H<->b cross terms still need a theorem bridge or a direct source-uniform
full-cross-term metric certificate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_p4_a21_finite_bias_cascade_bridge as BRIDGE
import ou3_sea3_a21_detectability_completion as ADET
import ou3_sea3_full_normal_live_word as WORD
import ou3_sea3_riccati_metric_p3 as P3

DEFAULT_DOMAIN = BRIDGE.DEFAULT_DOMAIN
SCHEMA = 1
QUALIFICATION = "OU3_P4_A21_FINITE_TAUB_COMPARISON_CASCADE_GM_ENERGY_V1"


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def geometric_decay_energy_upper(n: int, h_lower: float, tau_b: float) -> float:
    if n <= 0:
        raise ValueError("sample count must be positive")
    h = float(h_lower)
    tau = float(tau_b)
    if not (math.isfinite(h) and h > 0.0 and math.isfinite(tau) and tau > 0.0):
        raise ValueError("h/tau must be finite positive")
    a = down(2.0 * h / tau)
    # q*(1-q^N)/(1-q), evaluated with expm1 to avoid cancellation near q=1.
    q = up(math.exp(-a))
    numerator = up(-math.expm1(-up(a * n)))
    denominator = down(-math.expm1(-a))
    if not denominator > 0.0:
        raise RuntimeError("GM geometric denominator lost positivity")
    return up(q * up(numerator / denominator))


def sufficient_power10(log10_required: float) -> int:
    x = float(log10_required)
    if not math.isfinite(x):
        raise ValueError("required logarithmic weight must be finite")
    # One full decimal order of slack avoids depending on a rounded equality.
    return int(math.ceil(x)) + 1


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    bridge = BRIDGE.build(path)
    p3 = P3.build(path)
    adet = ADET.build(path)
    word = WORD.build(path)
    process = PROCESS.build()
    bad = {
        "bridge": BRIDGE.validate(bridge),
        "P3": P3.validate(p3),
        "A21_detectability": ADET.validate(adet),
        "word": WORD.validate(word),
        "process": PROCESS.validate(process),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"comparison-cascade prerequisites failed: {bad}")
    if p3["canonical_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("comparison cascade detached from canonical complete SEA3")
    if float(p3["useful_gate"]) != BRIDGE.P3_FROZEN_GATE:
        raise RuntimeError("frozen P3 gate changed")

    runtime = adet["active_bias_process"]
    tau_b = float(runtime["tau_ba_s"])
    # The process certificate owns the configured outward IMU-dt interval; the
    # literal word owns the requirement that every valid sample includes the
    # accelerometer Joseph operation.  Keep those ownerships distinct.
    h_bounds = process["configured_runtime"]["imu_dt_outward_interval_s"]
    h_lower = float(h_bounds[0])
    n_acc = int(word["imu_samples_upper"])
    if n_acc <= 0 or word["every_valid_imu_sample_requires_accelerometer_Joseph"] is not True:
        raise RuntimeError("complete word no longer retains every accelerometer update")

    gsum = geometric_decay_energy_upper(n_acc, h_lower, tau_b)
    mh = bridge["H18_metric_equivalence"]
    s_h = float(mh["normalization_s_H"])
    racc = float(mh["Racc_variance_lower"])
    if not (s_h > 0.0 and racc > 0.0):
        raise RuntimeError("H18 normalization/Racc lower lost positivity")
    c2 = up(up(s_h / racc) * gsum)
    if not (math.isfinite(c2) and c2 > 0.0):
        raise RuntimeError("GM energy C_Hb bound is not finite positive")

    g2_budget = float(bridge["cascade_g2_budget"])
    log_mu = BRIDGE.required_mu_log10(c2, g2_budget)
    p10 = sufficient_power10(log_mu)
    # mu=10^p10 is at least ten times the rounded required value, so
    # c2/mu < 0.1*g2_budget < g2_budget, without materializing huge powers.
    schur_ratio_upper = 0.1
    comparison_closed = bool(
        math.isfinite(log_mu)
        and p10 > 0
        and gsum > 0.0
        and c2 > 0.0
        and g2_budget > 0.0
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "P3_frozen_gate": BRIDGE.P3_FROZEN_GATE,
        "P3_frozen_not_modified": True,
        "paper_active_bias_route": adet["paper_active_bias_route"],
        "comparison_observer_form": "[[E_H,C_Hb],[0,Phi_b]]",
        "comparison_observer_only_not_shipping_estimator": True,
        "comparison_bias_correction_row_is_zero": True,
        "H18_comparison_uses_complete_SEA3_shipping_information_geometry": True,
        "same_complete_SEA3_word_used": True,
        "word_horizon_s": float(word["word_horizon_s"]),
        "required_accelerometer_updates": n_acc,
        "all_valid_accelerometer_updates_retained": True,
        "every_due_S_update_with_actual_RS_retained_inside_H18_word": True,
        "actual_RS_axis_factors": mh["actual_applied_R_S_axis_factors"],
        "actual_RS_regularization_not_replaced": True,
        "tau_ba_s": tau_b,
        "imu_dt_lower_s": h_lower,
        "GM_bias_energy_geometric_sum_upper": gsum,
        "GM_energy_sum_formula": "sum_{j=1}^N exp(-2*j*h_lower/tau_b)",
        "H18_mode_global_normalization_s_H": s_h,
        "Racc_variance_lower": racc,
        "joint_C_Hb_metric_c2_upper": c2,
        "C_Hb_bound_derivation": "exact Joseph bias-input supply plus closed-form GM energy",
        "suffix_norm_product_used": False,
        "N_times_worst_map_coupling_used": False,
        "packet_count_nonlinear_remainder_multiplier_used": False,
        "closed_form_same_word_bias_energy_used": True,
        "cascade_g2_budget": g2_budget,
        "mu_log10_required_upper": log_mu,
        "mu_sufficient_decimal": f"1e{p10}",
        "mu_sufficient_power10": p10,
        "c2_over_mu_to_g2_budget_ratio_upper": schur_ratio_upper,
        "comparison_cascade_Schur_condition_closed": comparison_closed,
        "comparison_observer_full_rank_21_state_metric_closed": comparison_closed,
        "full_ba_state_retained": True,
        "full_H18_state_retained_including_aw": True,
        "state_elimination_used": False,
        "a_w_elimination_used": False,
        "actual_A21_shipping_bias_correction_row_consumed": False,
        "actual_A21_shipping_H_b_cross_terms_closed": False,
        "actual_A21_shipping_full_cross_term_metric_closed": False,
        "nonlinear_residual_sector_attached": False,
        "P4_promoted_here": False,
        "P5_may_start": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "next_obligation": "bridge the quantitative finite-tau_b comparison-cascade storage to the ACTUAL A21 Kalman homogeneous word while retaining its bias correction row and H<->b cross terms, or construct the theorem-equivalent source-indexed full-cross-term metric directly; then attach the signed nonlinear complete-word sector",
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    if float(d.get("P3_frozen_gate", math.nan)) != BRIDGE.P3_FROZEN_GATE:
        f.append("P3 gate changed")
    for key in (
        "P3_frozen_not_modified",
        "comparison_observer_only_not_shipping_estimator",
        "comparison_bias_correction_row_is_zero",
        "H18_comparison_uses_complete_SEA3_shipping_information_geometry",
        "same_complete_SEA3_word_used",
        "all_valid_accelerometer_updates_retained",
        "every_due_S_update_with_actual_RS_retained_inside_H18_word",
        "actual_RS_regularization_not_replaced",
        "closed_form_same_word_bias_energy_used",
        "comparison_cascade_Schur_condition_closed",
        "comparison_observer_full_rank_21_state_metric_closed",
        "full_ba_state_retained",
        "full_H18_state_retained_including_aw",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "suffix_norm_product_used",
        "N_times_worst_map_coupling_used",
        "packet_count_nonlinear_remainder_multiplier_used",
        "state_elimination_used",
        "a_w_elimination_used",
        "actual_A21_shipping_bias_correction_row_consumed",
        "actual_A21_shipping_H_b_cross_terms_closed",
        "actual_A21_shipping_full_cross_term_metric_closed",
        "nonlinear_residual_sector_attached",
        "P4_promoted_here",
        "P5_may_start",
        "filter_changed",
        "declared_domain_changed",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("paper_active_bias_route") != "ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION":
        f.append("active-bias theorem route changed")
    for key in (
        "word_horizon_s",
        "tau_ba_s",
        "imu_dt_lower_s",
        "GM_bias_energy_geometric_sum_upper",
        "H18_mode_global_normalization_s_H",
        "Racc_variance_lower",
        "joint_C_Hb_metric_c2_upper",
        "cascade_g2_budget",
        "mu_log10_required_upper",
    ):
        x = d.get(key)
        if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
            f.append(f"invalid positive field {key}")
    p = d.get("mu_sufficient_power10")
    if not isinstance(p, int) or p <= 0 or d.get("mu_sufficient_decimal") != f"1e{p}":
        f.append("invalid sufficient mu representation")
    if float(d.get("c2_over_mu_to_g2_budget_ratio_upper", math.inf)) >= 1.0:
        f.append("chosen comparison-cascade mu does not satisfy Schur budget")
    if d.get("actual_RS_axis_factors") != [0.72, 0.72, 1.0]:
        f.append("actual R_S anisotropy changed")
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
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"C_Hb_c2_upper": d["joint_C_Hb_metric_c2_upper"], "mu_log10_required": d["mu_log10_required_upper"], "mu_sufficient": d["mu_sufficient_decimal"], "comparison_metric_closed": d["comparison_observer_full_rank_21_state_metric_closed"], "actual_shipping_bridge_closed": d["actual_A21_shipping_full_cross_term_metric_closed"], "P4_promoted": d["P4_promoted_here"], "failures": failures}, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
