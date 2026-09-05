#!/usr/bin/env python3
"""Fail-closed mixed detectable/stable A21 Riccati-completion status.

The active 21-state paper route used by this repository is

    eta6 UCO + finite residual accelerometer-bias correlation time,

not an assumed eta9 point-packet UCO condition.  Consequently the H18
prior-free identity cannot simply be copied to A21: the uniform completion

    Z D Z - delta Z >= -(delta^2/4) D^-1

requires a *full 21-state* strictly positive information matrix D.  Finite
``tau_b`` supplies detectability/stability of a possibly unobserved b_a
direction; it does not manufacture D_A21^-1.  Any certificate that appends the
shipping Q_ba block to the H18 18x18 LDLT while using the full-D inverse identity
has therefore changed the theorem hypothesis and is rejected here.

This module records the valid part of the stable-mode bridge.  At Live/A release
the source gives zero H18<->b_a covariance and a finite bias variance P_b0.  The
first active prediction is block diagonal in b_a and injects the exact shipping
Gauss--Markov process Q_ba.  Hence the homogeneous initial-bias covariance tag Cb
and that first process-noise tag Wb satisfy

    Cb_1 <= c_b Wb_1,
    c_b = phi_b^2 P_b0 / lambda_min(Q_ba).

Every later shipping prediction, Joseph update and left-error reset applies the
same congruence to both tags, while the total Omega receives only additional PSD
noise.  Thus the domination is preserved through the complete word.  This is a
useful quantitative stable-mode lemma, but by itself it does *not* close
Omega-delta P: a correct mixed completion must combine the detectable H18
complete-square term with this finite stable prior without ever invoking a
nonexistent full A21 D^-1.

The canonical P3 gate remains fail-closed until that mixed 21x21 inequality is
proved at delta=1e-18.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_sea3_a21_detectability_completion as ADET
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_full_word_event_algebra as EVENT
import ou3_sea3_h18_prior_free_completion as H18
import ou3_sea3_live_covariance_seed as LIVE

DEFAULT_DOMAIN = COMPLETE.DEFAULT_DOMAIN
SCHEMA = 2
QUALIFICATION = "OU3_COMPLETE_SEA3_A21_MIXED_DETECTABLE_STABLE_COMPLETION_STATUS"
USEFUL_GATE = 1.0e-18
HORIZON_S = 3.0


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    dynamic = DYNAMIC.build(path)
    process = PROCESS.build()
    h18 = H18.build(path)
    adet = ADET.build(path)
    live = LIVE.build(path)
    event = EVENT.build()
    bad = {
        "complete": COMPLETE.validate(complete),
        "dynamic": DYNAMIC.validate(dynamic),
        "process": PROCESS.validate(process),
        "H18": H18.validate(h18),
        "A21_detectability": ADET.validate(adet),
        "Live_seed": LIVE.validate(live),
        "event_algebra": EVENT.validate(event),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"A21 mixed-completion prerequisites failed: {bad}")
    source = complete["canonical_P3_source"]
    if source != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("A21 mixed completion detached from complete SEA3")
    if float(complete["word_horizon_s"]) != HORIZON_S:
        raise RuntimeError("canonical complete-SEA3 word is no longer 3 s")
    if adet["paper_active_bias_route"] != "ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION":
        raise RuntimeError("A21 paper route changed")
    if adet["eta9_point_packet_shortcut_used"] is not False:
        raise RuntimeError("eta9 point-packet shortcut entered A21 route")

    p_b0 = float(adet["H_to_A_release"]["bias_diagonal_floor_variance"])
    q_ba = float(process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"])
    tau_b = float(process["source_constants"]["accel_bias_tau_s"])
    dt = float(dynamic["validated_rate_and_jump_bounds"]["dt_s"])
    if not all(math.isfinite(x) and x > 0.0 for x in (p_b0, q_ba, tau_b, dt)):
        raise RuntimeError("invalid source-bound bias release/process constants")

    # Shipping active-bias prediction is phi_b=exp(-dt/tau_b).  Using 1 as an
    # outward-safe upper for phi_b is enough for the domination ratio and avoids
    # crediting any unvalidated transcendental rounding in this small lemma.
    phi_b_sq_upper = 1.0
    prior_to_first_process_ratio = math.nextafter(
        phi_b_sq_upper * p_b0 / q_ba, math.inf
    )
    delta_times_ratio = math.nextafter(
        USEFUL_GATE * prior_to_first_process_ratio, math.inf
    )
    if not (math.isfinite(prior_to_first_process_ratio) and prior_to_first_process_ratio > 0.0):
        raise RuntimeError("bias prior/process domination ratio invalid")

    preserve = event["full_matrix_margin_preservation"]
    congruence_suffix = all(bool(preserve[k]) for k in (
        "covers_prediction",
        "covers_every_due_S_update",
        "covers_every_Normal_Live_accelerometer_update",
        "covers_asynchronous_magnetometer_update",
        "covers_immediate_left_error_reset",
        "covers_not_due_or_rejected_identity_branches",
    ))

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": source,
        "complete_word_horizon_s": HORIZON_S,
        "component_of_complete_SEA3_full_word": True,
        "paper_active_bias_route": adet["paper_active_bias_route"],
        "eta9_point_packet_shortcut_used": False,
        "H18_prior_free_completion_consumed": True,
        "A21_finite_bias_detectability_consumed": True,
        "finite_tau_detectability_does_not_imply_full_A21_information_inverse": True,
        "full_A21_D_inverse_available": False,
        "full_A21_prior_free_D_inverse_identity_used": False,
        "invalid_append_Qba_to_H18_full_D_completion_rejected": True,
        "source_generated_bias_release": {
            "H18_ba_cross_covariances_zero_at_release": True,
            "bias_prior_variance": p_b0,
            "first_active_prediction_ba_process_lambda_min_lower": q_ba,
            "phi_b_squared_upper": phi_b_sq_upper,
        },
        "stable_bias_prior_process_domination": {
            "form": "Cb_after_first_prediction <= c_b * Wb_first_process",
            "ratio_c_b_upper": prior_to_first_process_ratio,
            "delta_times_ratio_upper": delta_times_ratio,
            "delta_times_ratio_is_small": delta_times_ratio < 1.0,
            "same_congruence_preserves_tag_domination_after_first_prediction": congruence_suffix,
            "later_total_Omega_contains_tag_plus_additional_PSD_noise": True,
            "this_lemma_alone_closes_full_A21_Riccati_inequality": False,
        },
        "actual_applied_SpectralMSE_R_S_retained_through_H18_component": True,
        "full_21x21_Omega_minus_delta_P_LDLT_closed": False,
        "A21_prior_free_completion_closed": False,
        "full_21x21_interval_LDLT_used": False,
        "event_algebra_preserves_margin_after_closure": congruence_suffix,
        "useful_gate": USEFUL_GATE,
        "source_family_replaced": False,
        "trajectory_replay_used": False,
        "independent_tau_sigma_RS_source_created": False,
        "P3_promoted": False,
        "next_obligation": (
            "prove the mixed detectable/stable 21x21 completion: keep the H18 complete-square "
            "term on the same complete SEA3 word, carry the finite b_a prior as a separate "
            "homogeneous covariance tag, and use shipping b_a GM process/noise without ever "
            "requiring full eta9 information or full A21 D^-1"
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
        "H18_prior_free_completion_consumed",
        "A21_finite_bias_detectability_consumed",
        "finite_tau_detectability_does_not_imply_full_A21_information_inverse",
        "invalid_append_Qba_to_H18_full_D_completion_rejected",
        "actual_applied_SpectralMSE_R_S_retained_through_H18_component",
        "event_algebra_preserves_margin_after_closure",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "eta9_point_packet_shortcut_used",
        "full_A21_D_inverse_available",
        "full_A21_prior_free_D_inverse_identity_used",
        "full_21x21_Omega_minus_delta_P_LDLT_closed",
        "A21_prior_free_completion_closed",
        "full_21x21_interval_LDLT_used",
        "source_family_replaced",
        "trajectory_replay_used",
        "independent_tau_sigma_RS_source_created",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"fail-closed/forbidden flag {key} changed")
    if d.get("paper_active_bias_route") != "ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION":
        f.append("paper active-bias route changed")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("useful gate changed")
    release = d.get("source_generated_bias_release", {})
    if release.get("H18_ba_cross_covariances_zero_at_release") is not True:
        f.append("bias release lost zero H18/ba cross covariance")
    dom = d.get("stable_bias_prior_process_domination", {})
    ratio = dom.get("ratio_c_b_upper")
    dr = dom.get("delta_times_ratio_upper")
    if not isinstance(ratio, (int, float)) or not (math.isfinite(float(ratio)) and float(ratio) > 0.0):
        f.append("bias prior/process domination ratio invalid")
    if not isinstance(dr, (int, float)) or not (math.isfinite(float(dr)) and 0.0 < float(dr) < 1.0):
        f.append("delta-times-bias-domination ratio is not small")
    for key in (
        "delta_times_ratio_is_small",
        "same_congruence_preserves_tag_domination_after_first_prediction",
        "later_total_Omega_contains_tag_plus_additional_PSD_noise",
    ):
        if dom.get(key) is not True:
            f.append(f"stable-bias domination lost {key}")
    if dom.get("this_lemma_alone_closes_full_A21_Riccati_inequality") is not False:
        f.append("stable-bias lemma was incorrectly promoted to full A21")
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
    dom = d["stable_bias_prior_process_domination"]
    print(json.dumps({
        "route": d["paper_active_bias_route"],
        "full_A21_D_inverse_available": d["full_A21_D_inverse_available"],
        "invalid_full_D_route_rejected": d["invalid_append_Qba_to_H18_full_D_completion_rejected"],
        "bias_prior_process_ratio_upper": dom["ratio_c_b_upper"],
        "delta_times_ratio_upper": dom["delta_times_ratio_upper"],
        "full_A21_Riccati_matrix_closed": d["full_21x21_Omega_minus_delta_P_LDLT_closed"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())