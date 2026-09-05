#!/usr/bin/env python3
"""Complete-SEA3 four-S information in the shipping Live covariance metric.

This is a sharpening of the existing actual-applied-SpectralMSE-R_S four-S
lemma, not a new source or reduced word.  The selected four S=0 firings remain
actual members of ``COMPLETE_SEA3_NORMAL_LIVE_WORD`` and every other shipping
event remains in the complete word.

The prior-free batch identity for one fixed-dimensional word gives the exact
sufficient condition

    (1-delta) D >= delta P0^-1.

Equivalently, with the shipping Live prior metric,

    P0^(1/2) D P0^(1/2) >= delta/(1-delta) I.

The retained four-S proof already bounds the inverse of the raw-record
observation matrix for the scaled translation state

    z = [S, g p, g^2 v, g^3 a_w].

Instead of first discarding the Live covariance geometry and replacing the
observation Gramian by a Euclidean scalar, weight those same inverse rows by
the exact shipping Live standard deviations.  If r_i bounds the l1 norm of row
i of M^-1 and s_i is the corresponding Live standard deviation, then

    ||Pz^-1/2 M^-1||_2^2
      <= ||Pz^-1/2 M^-1||_F^2
      <= sum_i (r_i/s_i)^2.

Together with the tight selected-record covariance bound Sigma_Y <= rho I,
this proves

    lambda_min(Pz^1/2 D_S Pz^1/2)
      >= 1 / (rho * sum_i (r_i/s_i)^2).

The Live S,p,v seed is fixed by shipping initialization.  The Live a_w block is
reset to the same committed stationary covariance before the first prediction;
the theorem family admits sigma_aw >= 0.05 m/s^2, so the smallest prior metric
uses that lower endpoint.  This is a uniform bound over the complete admitted
SEA3 family.  No independent tau/sigma/R_S rectangle is generated here; the
upstream four-S covariance bound remains the same source-uniform nuisance bound
and the metric is taken from the same shipping Live seed.
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
SCHEMA = 1
QUALIFICATION = "OU3_COMPLETE_SEA3_FOUR_S_SHIPPING_LIVE_PRIOR_METRIC"
USEFUL_GATE = 1.0e-18


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

    recovery = four["newton_coordinate_information"]["physical_state_recovery"]
    row_map = recovery["physical_state_inverse_raw_record_row_l1_upper"]
    rows = [float(row_map[name]) for name in ("S", "g*p", "g^2*v", "g^3*a_w")]
    if any(not (math.isfinite(x) and x > 0.0) for x in rows):
        raise RuntimeError("four-S physical inverse row bound lost positivity")

    weighted_inverse_frob_sq = 0.0
    terms = []
    for name, row, std in zip(("S", "g*p", "g^2*v", "g^3*a_w"), rows, stds):
        term = up((row / std) ** 2)
        weighted_inverse_frob_sq = up(weighted_inverse_frob_sq + term)
        terms.append({"state": name, "inverse_row_l1_upper": row, "Live_std_lower": std, "weighted_square_upper": term})
    metric_observation_gram_lower = down(1.0 / weighted_inverse_frob_sq)

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
        "scaled_translation_state": ["S", "g*p", "g^2*v", "g^3*a_w"],
        "scaled_Live_std_lower": stds,
        "weighted_inverse_row_terms": terms,
        "weighted_inverse_frobenius_squared_upper": weighted_inverse_frob_sq,
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
            "retain the accelerometer attitude-a_w row structure and combine this strong translation Live-prior metric "
            "with the eta6 PE information in one H18 prior-metric matrix; then extend the same construction through "
            "the shipping H-to-A accelerometer-bias release"
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
        "weighted_inverse_frobenius_squared_upper": d["weighted_inverse_frobenius_squared_upper"],
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
