#!/usr/bin/env python3
"""Complete-SEA3 four-S information in the shipping Live covariance metric.

This sharpens the existing actual-applied-SpectralMSE-R_S four-S lemma without
creating another source, word, estimator or contraction architecture.  The four
selected S=0 firings are actual members of ``COMPLETE_SEA3_NORMAL_LIVE_WORD``;
every other due S update and every other shipping event remains in the literal
complete word.

For the exact prior-free batch factorization, a sufficient endpoint condition is

    (1-delta) D >= delta P0^-1,

or, equivalently,

    P0^(1/2) D P0^(1/2) >= delta/(1-delta) I.

The existing four-S proof factors the raw-record observation map as

    y --T_dd--> q --R^-1--> z,

where q are Newton divided differences and

    z = [S, g p, g^2 v, g^3 a_w].

Thus

    (M Pz^(1/2))^-1 = Pz^-1/2 R^-1 T_dd.

Keeping that factorization is substantially tighter than first collapsing every
row of M^-1 to a separate l1 bound.  The triangular Newton recovery gives
explicit coefficient bounds for R^-1, while the exact-rational four-S producer
already gives one- and infinity-norm bounds for T_dd.  We therefore certify

    ||Pz^-1/2 R^-1 T_dd||_2^2
      <= ||Pz^-1/2 R^-1||_F^2 ||T_dd||_2^2
      <= F_R^2 ||T_dd||_1 ||T_dd||_inf.

Together with the tight selected-record covariance bound ``Sigma_Y <= rho I``,

    lambda_min(Pz^1/2 D_S Pz^1/2)
      >= 1 / (rho F_R^2 ||T_dd||_1 ||T_dd||_inf).

The Live S,p,v variances are the shipping constructor values.  At Live entry the
a_w covariance is reset to the same committed stationary covariance used by the
word, and the admitted family has sigma_aw >= 0.05 m/s^2; therefore the smallest
Live prior metric uses 0.05.  The upstream four-S nuisance covariance remains
uniform over arbitrary legal time-varying SEA3 tau/sigma and actual applied
SpectralMSE R_S.  No independent parameter box is introduced here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import down, up
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_four_s_translation_information_tight as FOUR
import ou3_sea3_live_covariance_seed as LIVE

DEFAULT_DOMAIN = FOUR.DEFAULT_DOMAIN
SCHEMA = 2
QUALIFICATION = "OU3_COMPLETE_SEA3_FOUR_S_SHIPPING_LIVE_PRIOR_METRIC"
USEFUL_GATE = 1.0e-18


def _sum_squares_upper(values: list[float]) -> float:
    total = 0.0
    for x in values:
        total = up(total + up(x * x))
    return total


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    four = FOUR.build(path)
    live = LIVE.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "four_S": FOUR.validate(four),
        "Live_seed": LIVE.validate(live),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"four-S Live-prior prerequisites failed: {bad}")
    source = complete["canonical_P3_source"]
    if source != "COMPLETE_SEA3_NORMAL_LIVE_WORD" or four["canonical_source"] != source:
        raise RuntimeError("four-S Live-prior metric detached from complete SEA3")

    g = float(four["uniform_S_gap_s_upper"])
    if not (math.isfinite(g) and g > 0.0):
        raise RuntimeError("invalid certified S scheduler gap")

    tr = live["translation_seed"]
    sS = math.sqrt(float(tr["P_S"]))
    sp = g * math.sqrt(float(tr["P_p"]))
    sv = (g * g) * math.sqrt(float(tr["P_v"]))
    aw = live["aw_live_seed"]
    saw_min = float(aw["committed_vertical_std_interval_mps2"][0])
    sa = (g ** 3) * saw_min
    stds = [sS, sp, sv, sa]
    if any(not (math.isfinite(x) and x > 0.0) for x in stds):
        raise RuntimeError("shipping Live prior metric lost positivity")

    ni = four["newton_coordinate_information"]
    recovery = ni["physical_state_recovery"]
    d3 = float(recovery["third_divided_difference_lower"])
    deriv = recovery["derivative_upper_bounds"]
    c0 = float(deriv["c_u0"])
    c01 = float(deriv["c_first_divided_difference"])
    c012 = float(deriv["c_second_divided_difference"])
    windows = four["scaled_u_windows"]
    u0 = float(windows[0][1])
    u1 = float(windows[1][1])
    if not all(math.isfinite(x) and x > 0.0 for x in (d3, c0, c01, c012, u0, u1)):
        raise RuntimeError("Newton recovery bounds lost positivity")
    u01 = up(u0 + u1)
    u0u1 = up(u0 * u1)

    # Exact triangular inverse structure in q coordinates:
    # A = q3/d3
    # V = 2 q2 - 2 c012 q3/d3
    # P = q1 -(u0+u1)q2 + ((u0+u1)c012-c01)q3/d3
    # S = q0-u0 q1+u0*u1 q2 +(u0*c01-u0*u1*c012-c0)q3/d3.
    # Signs/cancellations are not assumed in the interval upper; only the
    # triangular sparsity is retained.
    Rinv_rows = [
        [1.0, u0, u0u1, up((up(u0 * c01) + up(u0u1 * c012) + c0) / d3)],
        [0.0, 1.0, u01, up((up(u01 * c012) + c01) / d3)],
        [0.0, 0.0, 2.0, up(2.0 * c012 / d3)],
        [0.0, 0.0, 0.0, up(1.0 / d3)],
    ]

    weighted_rows = []
    weighted_Rinv_frob_sq = 0.0
    names = ["S", "g*p", "g^2*v", "g^3*a_w"]
    for name, row, std in zip(names, Rinv_rows, stds):
        scaled = [up(abs(x) / std) for x in row]
        row_sq = _sum_squares_upper(scaled)
        weighted_Rinv_frob_sq = up(weighted_Rinv_frob_sq + row_sq)
        weighted_rows.append({
            "state": name,
            "R_inverse_q_coefficient_abs_upper": row,
            "Live_std_lower": std,
            "Pz_inverse_half_R_inverse_row_abs_upper": scaled,
            "row_squared_norm_upper": row_sq,
        })

    Tdd_one = float(ni["L_inverse_one_norm_upper"])
    Tdd_inf = float(ni["L_inverse_infinity_norm_upper"])
    if not (math.isfinite(Tdd_one) and Tdd_one > 0.0 and math.isfinite(Tdd_inf) and Tdd_inf > 0.0):
        raise RuntimeError("raw-record to Newton norm bound invalid")
    Tdd_spectral_sq_upper = up(Tdd_one * Tdd_inf)
    weighted_full_inverse_spectral_sq_upper = up(
        weighted_Rinv_frob_sq * Tdd_spectral_sq_upper
    )
    metric_observation_gram_lower = down(1.0 / weighted_full_inverse_spectral_sq_upper)

    rho = float(four["selected_S_record_noise"]["four_record_covariance_lambda_max_upper"])
    if not (math.isfinite(rho) and rho > 0.0):
        raise RuntimeError("selected four-S covariance upper invalid")
    prior_metric_information_lower = down(metric_observation_gram_lower / rho)
    required = up(USEFUL_GATE / down(1.0 - USEFUL_GATE))
    passed = prior_metric_information_lower >= required

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": source,
        "component_of_complete_SEA3_full_word": True,
        "selected_four_S_events_replace_complete_word": False,
        "all_due_S_updates_remain_in_literal_word": True,
        "actual_applied_SpectralMSE_R_S_consumed": True,
        "same_tight_four_S_record_covariance_consumed": True,
        "shipping_Live_covariance_seed_consumed": True,
        "Live_aw_reset_to_same_committed_stationary_covariance_consumed": True,
        "minimum_Live_aw_std_mps2": saw_min,
        "scaled_translation_state": names,
        "scaled_Live_std_lower": stds,
        "triangular_Newton_recovery_retained": True,
        "weighted_R_inverse_rows": weighted_rows,
        "weighted_R_inverse_frobenius_squared_upper": weighted_Rinv_frob_sq,
        "raw_record_to_Newton_one_norm_upper": Tdd_one,
        "raw_record_to_Newton_infinity_norm_upper": Tdd_inf,
        "raw_record_to_Newton_spectral_norm_squared_upper": Tdd_spectral_sq_upper,
        "weighted_full_inverse_spectral_norm_squared_upper": weighted_full_inverse_spectral_sq_upper,
        "prior_metric_observation_gram_lambda_min_lower": metric_observation_gram_lower,
        "selected_record_covariance_lambda_max_upper": rho,
        "translation_Live_prior_metric_information_lambda_min_lower": prior_metric_information_lower,
        "useful_gate": USEFUL_GATE,
        "batch_required_information_ratio": required,
        "translation_batch_prior_condition_pass": passed,
        "exact_batch_condition": "(1-delta)D >= delta P0^-1",
        "determinant_used_for_metric_bound": False,
        "blockwise_contraction_ratio_used": False,
        "D_W_L_W_product_used": False,
        "trajectory_replay_used": False,
        "source_family_replaced": False,
        "independent_tau_sigma_RS_source_created": False,
        "P3_promoted": False,
        "next_obligation": (
            "retain the accelerometer attitude-a_w row structure and combine this translation Live-prior metric "
            "with eta6 PE in one H18 prior-metric matrix, including process nuisance in the selected vector records; "
            "then extend the same batch construction through the shipping H-to-A accelerometer-bias release"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "component_of_complete_SEA3_full_word",
        "all_due_S_updates_remain_in_literal_word",
        "actual_applied_SpectralMSE_R_S_consumed",
        "same_tight_four_S_record_covariance_consumed",
        "shipping_Live_covariance_seed_consumed",
        "Live_aw_reset_to_same_committed_stationary_covariance_consumed",
        "triangular_Newton_recovery_retained",
        "translation_batch_prior_condition_pass",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "selected_four_S_events_replace_complete_word",
        "determinant_used_for_metric_bound",
        "blockwise_contraction_ratio_used",
        "D_W_L_W_product_used",
        "trajectory_replay_used",
        "source_family_replaced",
        "independent_tau_sigma_RS_source_created",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"forbidden/open flag {key} changed")
    info = d.get("translation_Live_prior_metric_information_lambda_min_lower")
    req = d.get("batch_required_information_ratio")
    if not isinstance(info, (int, float)) or not (math.isfinite(float(info)) and float(info) > 0.0):
        f.append("translation Live-prior information is not strict")
    if not isinstance(req, (int, float)) or not (math.isfinite(float(req)) and float(req) > 0.0):
        f.append("batch required ratio is invalid")
    if isinstance(info, (int, float)) and isinstance(req, (int, float)) and float(info) < float(req):
        f.append("translation Live-prior information does not clear batch condition")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("useful gate changed")
    if float(d.get("weighted_full_inverse_spectral_norm_squared_upper", math.inf)) <= 0.0:
        f.append("weighted inverse spectral upper invalid")
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
        "scaled_Live_std_lower": d["scaled_Live_std_lower"],
        "weighted_R_inverse_frobenius_squared_upper": d["weighted_R_inverse_frobenius_squared_upper"],
        "Tdd_spectral_norm_squared_upper": d["raw_record_to_Newton_spectral_norm_squared_upper"],
        "metric_observation_gram_lower": d["prior_metric_observation_gram_lambda_min_lower"],
        "record_covariance_upper": d["selected_record_covariance_lambda_max_upper"],
        "translation_prior_metric_information_lower": d["translation_Live_prior_metric_information_lambda_min_lower"],
        "batch_required_ratio": d["batch_required_information_ratio"],
        "pass": d["translation_batch_prior_condition_pass"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
