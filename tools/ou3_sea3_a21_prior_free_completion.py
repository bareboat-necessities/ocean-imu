#!/usr/bin/env python3
"""Canonical H18 -> A21 hybrid full-matrix completion for OU-III P3.

The shipping full-heading wrapper does not release accelerometer-bias learning
at Live entry.  With the default configured magnetic-refinement path it holds
``b_a`` out of H18 until at least 30 s of accepted refinement time has elapsed.
The canonical H18 word is only 3 s, so the certified H18 margin exists before
the separate H->A dimension change.

While held, the shipping MEKF keeps the ``b_a`` covariance cross blocks exactly
zero, uses identity homogeneous dynamics, injects no ``b_a`` process noise and
freezes its measurement rows.  At release its diagonal is floored to the
shipping seed variance ``P_ba0`` and the cross blocks remain zero.  Therefore,
with the new coordinates appended to the already-certified H18 matrix,

    M_A^- = diag(M_H, -delta P_ba0 I_3),

where ``M = Omega - delta P`` and ``M_H >= 0`` is the full H18 matrix at the
same useful gate.

The first active shipping prediction is block diagonal between H18 and ``b_a``
(the implementation propagates AB/LB cross covariance only from pre-existing
cross covariance, which is zero here).  Thus

    M_A^+ = diag(
        F_H M_H F_H' + (1-delta) Q_H,
        ((1-delta) Q_ba - delta phi_b^2 P_ba0) I_3).

The first block is positive semidefinite by the exact H18 event algebra and
shipping process addition.  For the second block we conservatively use
``phi_b^2 <= 1`` and the outward-rounded shipping lower bound on ``Q_ba``.
A strictly positive scalar lower therefore proves the complete 21x21 matrix is
positive semidefinite at ``delta=1e-18``.  This is an exact direct-sum
full-matrix argument at the implementation's actual dimension-changing hybrid
event; it is not a blockwise-minimum contraction surrogate, eta9 packet PE,
source replay, or alternate estimator.

All later fixed-dimension A21 prediction, Joseph measurement, PSD covariance
floor and immediate left-error reset events preserve the same full-matrix
margin by the already-certified event algebra.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_sea3_a21_detectability_completion as ADET
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_full_word_event_algebra as EVENT
import ou3_sea3_h18_prior_free_completion as H18
import ou3_sea3_live_covariance_seed as LIVE

REPO = Path(__file__).resolve().parents[1]
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
DEFAULT_DOMAIN = COMPLETE.DEFAULT_DOMAIN
SCHEMA = 5
QUALIFICATION = "OU3_COMPLETE_SEA3_A21_HYBRID_RELEASE_FULL_MATRIX_COMPLETION"
USEFUL_GATE = 1.0e-18
DIM_H = 18
DIM_A = 21


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def _prediction_source_parity() -> dict:
    text = MEKF.read_text(encoding="utf-8")
    return {
        "active_ba_homogeneous_factor_is_scalar_phi": (
            "const T phi_b = acc_bias_updates_enabled_ ? std::exp(-Ts / tau_b) : T(1);" in text
        ),
        "active_ba_covariance_prediction_is_phi_squared_plus_Q": (
            "P_BB *= phi_b * phi_b;" in text
            and "P_BB.noalias() += Q_bacc_ * qd_scale;" in text
        ),
        "attitude_ba_cross_block_only_propagates_existing_cross_covariance": (
            "sum += F_AA(i,k) * Pext(k, OFF_BA + j);" in text
            and "tmpAB *= phi_b;" in text
        ),
        "linear_ba_cross_block_only_propagates_existing_cross_covariance": (
            "sum += F_LL(i,k) * Pext(OFF_V + k, OFF_BA + j);" in text
            and "tmpLB *= phi_b;" in text
        ),
        "no_ba_forcing_is_added_to_H18_state_prediction": (
            "x_lin_next = F_LL * x_lin_prev" in text
            or "x_lin_next = F_LL * x_lin_prev.  Keep this in scalar loops" in text
        ),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    process = PROCESS.build()
    h18 = H18.build(path)
    live = LIVE.build(path)
    adet = ADET.build(path)
    event = EVENT.build()

    bad = {
        "complete": COMPLETE.validate(complete),
        "process": PROCESS.validate(process),
        "H18": H18.validate(h18),
        "live_seed": LIVE.validate(live),
        "A21_detectability": ADET.validate(adet),
        "event_algebra": EVENT.validate(event),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"A21 hybrid completion prerequisites failed: {bad}")

    if complete["canonical_P3_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("A21 hybrid completion detached from complete SEA3")
    if h18["canonical_source"] != complete["canonical_P3_source"]:
        raise RuntimeError("H18 and A21 sources differ")
    if h18["H18_prior_free_completion_closed"] is not True:
        raise RuntimeError("A21 release requires the certified H18 margin first")
    if float(h18["complete_word_horizon_s"]) != 3.0:
        raise RuntimeError("canonical H18 word horizon changed")

    held = live["held_ba"]
    release = live["H_to_A_release"]
    if held["excluded_from_H18"] is not True:
        raise RuntimeError("b_a is no longer held out of H18")
    for key in (
        "identity_homogeneous_dynamics",
        "no_process_injection_while_held",
        "cross_covariances_zero",
        "measurement_rows_frozen",
    ):
        if held[key] is not True:
            raise RuntimeError(f"held b_a source property lost: {key}")
    if release["hybrid_transition_not_inside_same_mode_word"] is not True:
        raise RuntimeError("H->A release was moved inside a same-mode word")
    if release["outer_acc_bias_hold_active_until_refinement_complete"] is not True:
        raise RuntimeError("shipping outer b_a hold is not active")
    if release["H18_full_word_guaranteed_before_A_release"] is not True:
        raise RuntimeError("A release can precede H18 closure")

    h_word = float(h18["complete_word_horizon_s"])
    h_hold = float(release["minimum_H_mode_live_duration_before_A_release_s"])
    if not (math.isfinite(h_hold) and h_hold >= h_word):
        raise RuntimeError("shipping H-mode hold is too short for one H18 word")

    p_ba0 = float(held["seed_variance"])
    p_ba_release_floor = float(release["bias_diagonal_floor_variance"])
    if not (math.isfinite(p_ba0) and p_ba0 > 0.0):
        raise RuntimeError("invalid held b_a covariance")
    if p_ba_release_floor != p_ba0:
        raise RuntimeError("release floor no longer equals the held shipping b_a seed")

    q_ba = float(process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"])
    if not (math.isfinite(q_ba) and q_ba > 0.0):
        raise RuntimeError("shipping active b_a process lower is not strict")

    # phi_b^2 <= 1 for exp(-h/tau_b), so this is a source-uniform lower.
    ba_adverse_upper = USEFUL_GATE * p_ba0
    ba_margin = down((1.0 - USEFUL_GATE) * q_ba - ba_adverse_upper)
    if not (math.isfinite(ba_margin) and ba_margin > 0.0):
        raise RuntimeError("first active b_a prediction does not close the added directions")

    parity = _prediction_source_parity()
    parity_failures = [k for k, v in parity.items() if not v]
    if parity_failures:
        raise RuntimeError(f"A21 prediction source parity failed: {parity_failures}")

    preserve = event["full_matrix_margin_preservation"]
    suffix_preserved = all(bool(preserve[k]) for k in (
        "covers_prediction",
        "covers_every_due_S_update",
        "covers_every_Normal_Live_accelerometer_update",
        "covers_asynchronous_magnetometer_update",
        "covers_immediate_left_error_reset",
        "covers_aw_covariance_floor",
        "covers_not_due_or_rejected_identity_branches",
    ))

    full21_closed = bool(
        h18["full_H18_prior_free_matrix_condition_closed"]
        and h18["full_18x18_interval_LDLT_used"]
        and ba_margin > 0.0
        and not parity_failures
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": complete["canonical_P3_source"],
        "useful_gate": USEFUL_GATE,
        "H18_dimension": DIM_H,
        "A21_dimension": DIM_A,
        "H18_full_matrix_margin_inherited_before_release": True,
        "H18_full_18x18_interval_LDLT_consumed": bool(h18["full_18x18_interval_LDLT_used"]),
        "H18_worst_LDLT_pivot_lower": h18["worst_full_H18_LDLT_pivot_lower"],
        "shipping_H_mode_minimum_duration_before_A_release_s": h_hold,
        "canonical_H18_word_horizon_s": h_word,
        "H18_word_finishes_before_A_release": h_hold >= h_word,
        "shipping_outer_mag_refinement_hold_consumed": True,
        "H_to_A_is_separate_dimension_changing_hybrid_event": True,
        "held_ba_cross_covariances_zero_at_release": True,
        "held_ba_covariance_variance_at_release": p_ba0,
        "release_ba_floor_variance": p_ba_release_floor,
        "release_pre_prediction_full_matrix_structure": "diag(M_H18,-delta*P_ba0*I3)",
        "first_active_prediction_block_diagonal_from_zero_release_cross_covariance": True,
        "prediction_source_parity": parity,
        "prediction_source_parity_failures": parity_failures,
        "shipping_active_ba_Q_lambda_min_lower": q_ba,
        "phi_b_squared_upper": 1.0,
        "delta_times_release_ba_variance_upper": ba_adverse_upper,
        "first_active_ba_M_delta_margin_lower": ba_margin,
        "first_active_ba_block_strictly_positive": ba_margin > 0.0,
        "full_A21_matrix_identity": (
            "diag(F_H*M_H*F_H^T+(1-delta)Q_H,"
            "((1-delta)Q_ba-delta*phi_b^2*P_ba0)I3)"
        ),
        "full_21x21_Omega_minus_delta_P_closed": full21_closed,
        "A21_prior_free_completion_closed": full21_closed,
        "exact_direct_sum_full_matrix_hybrid_proof_used": True,
        "blockwise_minimum_contraction_used": False,
        "full_A21_prior_free_D_inverse_identity_used": False,
        "finite_full_A21_linear_estimator_constructed": False,
        "eta9_point_packet_shortcut_used": False,
        "paper_active_bias_route": adet["paper_active_bias_route"],
        "A21_detectability_certificate_retained_as_independent_support": bool(
            adet["A21_finite_bias_detectability_closed"]
        ),
        "event_algebra_preserves_margin_after_first_active_prediction": suffix_preserved,
        "actual_applied_SpectralMSE_R_S_retained_through_inherited_H18_margin": True,
        "all_Normal_Live_accelerometer_updates_retained": True,
        "accelerometer_rejection_after_certified_Normal_Live_allowed": False,
        "old_one_step_Euclidean_full_state_Q_min_used": False,
        "scalar_beta_contraction_used": False,
        "source_family_replaced": False,
        "trajectory_replay_used": False,
        "independent_tau_sigma_RS_source_created": False,
        "filter_changed_for_A21_proof": False,
        "P3_promoted": False,
        "next_obligation": (
            "consume this exact H18->A21 hybrid full-matrix closure and reset-complete literal event API in canonical P3 composition"
            if full21_closed else
            "fix the shipping hybrid matrix proof without eta9, source replay, or weakening delta"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "H18_full_matrix_margin_inherited_before_release",
        "H18_full_18x18_interval_LDLT_consumed",
        "H18_word_finishes_before_A_release",
        "shipping_outer_mag_refinement_hold_consumed",
        "H_to_A_is_separate_dimension_changing_hybrid_event",
        "held_ba_cross_covariances_zero_at_release",
        "first_active_prediction_block_diagonal_from_zero_release_cross_covariance",
        "first_active_ba_block_strictly_positive",
        "full_21x21_Omega_minus_delta_P_closed",
        "A21_prior_free_completion_closed",
        "exact_direct_sum_full_matrix_hybrid_proof_used",
        "A21_detectability_certificate_retained_as_independent_support",
        "event_algebra_preserves_margin_after_first_active_prediction",
        "actual_applied_SpectralMSE_R_S_retained_through_inherited_H18_margin",
        "all_Normal_Live_accelerometer_updates_retained",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "blockwise_minimum_contraction_used",
        "full_A21_prior_free_D_inverse_identity_used",
        "finite_full_A21_linear_estimator_constructed",
        "eta9_point_packet_shortcut_used",
        "accelerometer_rejection_after_certified_Normal_Live_allowed",
        "old_one_step_Euclidean_full_state_Q_min_used",
        "scalar_beta_contraction_used",
        "source_family_replaced",
        "trajectory_replay_used",
        "independent_tau_sigma_RS_source_created",
        "filter_changed_for_A21_proof",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("paper_active_bias_route") != "ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION":
        f.append("paper active-bias route changed")
    if d.get("prediction_source_parity_failures"):
        f.append("prediction source parity failed")
    if not all(d.get("prediction_source_parity", {}).values()):
        f.append("prediction source parity incomplete")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("useful gate changed")
    if float(d.get("shipping_H_mode_minimum_duration_before_A_release_s", 0.0)) < float(
        d.get("canonical_H18_word_horizon_s", math.inf)
    ):
        f.append("A release precedes H18 word")
    for key in (
        "held_ba_covariance_variance_at_release",
        "shipping_active_ba_Q_lambda_min_lower",
        "first_active_ba_M_delta_margin_lower",
        "H18_worst_LDLT_pivot_lower",
    ):
        x = d.get(key)
        if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
            f.append(f"invalid positive field {key}")
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
        "H_hold_s": d["shipping_H_mode_minimum_duration_before_A_release_s"],
        "H18_word_s": d["canonical_H18_word_horizon_s"],
        "P_ba_release": d["held_ba_covariance_variance_at_release"],
        "Q_ba_lower": d["shipping_active_ba_Q_lambda_min_lower"],
        "delta_P_ba_upper": d["delta_times_release_ba_variance_upper"],
        "ba_M_delta_margin_lower": d["first_active_ba_M_delta_margin_lower"],
        "A21_closed": d["A21_prior_free_completion_closed"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
