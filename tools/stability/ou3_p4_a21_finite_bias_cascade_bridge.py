#!/usr/bin/env python3
"""Finite-tau_b full-21-state cascade bridge for canonical complete-SEA3 P4.

This module is theorem-facing but deliberately fail-closed.  It formulates the
full-rank active-bias cascade metric that the post-#496 ledger requires without
using the retired longer-point-window route or the detectability helper's crude
packet-count coupling estimate.

Let the H18 complete-word map satisfy, in its source-indexed information metric,

    ||E_H h||_{M_H+}^2 <= (1-delta_H) ||h||_{M_H-}^2,

and let the residual accelerometer-bias homogeneous map satisfy

    ||Phi_b b||^2 <= (1-beta_b) ||b||^2.

For a full 21-state comparison map

    E_A = [[E_H, C_Hb],
           [  0, Phi_b]],

use the full-rank block storage

    M_A(zeta) = diag(M_H(zeta), mu I_3).

If c^2 bounds ||M_H+^(1/2) C_Hb||_2^2 and g^2=c^2/mu, then a target
0 < delta_A < min(delta_H,beta_b) is guaranteed by the exact 2x2 Schur test

    g^2 < ((delta_H-delta_A)(beta_b-delta_A))/(1-delta_A).

Because delta_A << binary64 epsilon, the implementation never evaluates
1-delta_A as a proof of strictness.  It stores positive gaps and determinant
slack directly.

The missing object is intentionally explicit: a source-uniform JOINT bound on
C_Hb over the same complete SEA3 word.  The older detectability helper reports
a finite `N * L^N` style bound only to prove finiteness of a triangular observer;
that packet-count bound is not consumed here and cannot promote P4.

A fresh full-matrix covariance lower bound is also provided for the H18 endpoint
metric.  It starts from the shipping full-process Q lower bound, pessimistically
assimilates every same-sample S/accelerometer/magnetometer Joseph measurement,
and uses the exact reset fact sigma_min(I+0.5[d]_x)=1.  Thus it does not infer
an inverse metric floor from covariance marginals and does not require any
correction-radius assumption.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_sea3_a21_detectability_completion as ADET
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_full_word_event_algebra as EVENT
import ou3_sea3_h18_prior_free_completion as H18
import ou3_sea3_riccati_metric_p3 as P3
import ou3_sea3_windowed_vector_pe as PE

REPO = Path(__file__).resolve().parents[2]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"

SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_A21_FINITE_TAUB_FULL_STATE_CASCADE_BRIDGE_V1"
P3_FROZEN_GATE = 1.0e-18
TARGET_P4_CASCADE_GAP = 5.0e-19


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _positive(x: float, label: str) -> float:
    y = float(x)
    if not (math.isfinite(y) and y > 0.0):
        raise RuntimeError(f"{label} must be finite positive, got {x!r}")
    return y


def cascade_g2_limit_conservative(delta_h: float, beta_b: float, delta_a: float) -> float:
    """Conservative strict-cascade budget for g^2=c^2/mu."""
    dh = _positive(delta_h, "delta_h")
    bb = _positive(beta_b, "beta_b")
    da = _positive(delta_a, "delta_a")
    if not da < min(dh, bb):
        raise ValueError("target cascade gap must be below both diagonal gaps")
    a = down(dh - da)
    b = down(bb - da)
    if not (a > 0.0 and b > 0.0):
        raise RuntimeError("cascade diagonal slack rounded nonpositive")
    return down(a * b)


def cascade_determinant_slack_lower(delta_h: float, beta_b: float, delta_a: float, g2_upper: float) -> float:
    """Lower bound on det((1-delta_A)I - T^T T)."""
    da = _positive(delta_a, "delta_a")
    g2 = float(g2_upper)
    if not (math.isfinite(g2) and g2 >= 0.0):
        raise ValueError("g2 upper must be finite nonnegative")
    a = down(_positive(delta_h, "delta_h") - da)
    b = down(_positive(beta_b, "beta_b") - da)
    return down(down(a * b) - up(g2))


def required_mu_log10(c2_upper: float, g2_budget: float) -> float:
    """Return log10(mu) sufficient for c^2/mu <= g2_budget."""
    c2 = float(c2_upper)
    g2 = _positive(g2_budget, "g2_budget")
    if not (math.isfinite(c2) and c2 >= 0.0):
        raise ValueError("c2 upper must be finite nonnegative")
    if c2 == 0.0:
        return -math.inf
    return up(math.log10(c2) - math.log10(g2))


def _axis_factors() -> list[float]:
    text = WRAPPER.read_text(encoding="utf-8")
    if "RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;" not in text:
        raise RuntimeError("cascade bridge requires deployed SpectralMSE R_S")
    out = []
    for name in ("R_S_x_factor_", "R_S_y_factor_"):
        m = re.search(rf"float\s+{name}\s*=\s*([0-9.eE+-]+)f", text)
        if not m:
            raise RuntimeError(f"cannot extract deployed {name}")
        out.append(float(m.group(1)))
    return out + [1.0]


def _measurement_variance_lower(std_xyz: list[float], label: str) -> float:
    if len(std_xyz) != 3:
        raise RuntimeError(f"{label} std must be length three")
    return min(down(_positive(x, f"{label} std") ** 2) for x in std_xyz)


def _h18_information_metric_upper(domain: dict, h18: dict, pe: dict, process: dict, dynamic: dict, event: dict) -> dict:
    qh = _positive(process["modes"]["H"]["prediction_Q_lambda_min_lower"], "H18 Q lower")
    live = domain["normal_live"]
    fmax = _positive(live["specific_force_norm_upper_mps2"], "force upper")
    mmax = _positive(live["magnetic_vector_norm_upper_uT"], "mag upper")
    meas = pe["measurement_runtime"]
    racc = _measurement_variance_lower(list(map(float, meas["accelerometer_std_mps2"])), "accelerometer")
    rmag = _measurement_variance_lower(list(map(float, meas["magnetometer_std_uT"])), "magnetometer")
    rs_base = _positive(dynamic["dynamic_invariant"]["R_S_applied"][0], "applied R_S lower")
    factors = _axis_factors()
    rs_std_min = down(rs_base * min(factors))
    rs_var_min = down(rs_std_min * rs_std_min)
    if not rs_var_min > 0.0:
        raise RuntimeError("actual-applied R_S lower lost positivity")
    info_s = up(1.0 / rs_var_min)
    info_acc = up(up(fmax * fmax + 1.0) / racc)
    info_mag = up(up(mmax * mmax) / rmag)
    precision_upper = up(up(1.0 / qh) + up(info_s + up(info_acc + info_mag)))
    pmin = down(1.0 / precision_upper)
    if not pmin > 0.0:
        raise RuntimeError("full H18 posterior covariance lower is not strict")
    s_h = _positive(h18["same_word_diffuse_prior_covariance_upper"]["Pbar_trace_upper"], "H18 source-uniform covariance trace upper")
    mplus = up(s_h / pmin)
    if not (math.isfinite(mplus) and mplus > 0.0):
        raise RuntimeError("H18 information metric upper is not finite")
    reset = event["left_error_reset"]
    reset_parity = bool(reset["same_full_matrix_margin_preserved_by_congruence"] and reset["small_angle_needed_for_nonsingularity"] is False and float(reset["determinant_lower"]) == 1.0)
    return {
        "normalization_s_H": s_h,
        "prediction_Q_lambda_min_lower": qh,
        "Racc_variance_lower": racc,
        "Rmag_variance_lower": rmag,
        "actual_applied_R_S_base_lower": rs_base,
        "actual_applied_R_S_axis_factors": factors,
        "actual_applied_R_S_variance_lower": rs_var_min,
        "same_sample_information_spectral_upper": precision_upper - 1.0 / qh,
        "posterior_covariance_lambda_min_lower": pmin,
        "M_H_lambda_max_upper": mplus,
        "M_H_lambda_min_lower_from_Pbar_trace_normalization": 1.0,
        "Joseph_information_identity_used": True,
        "every_possible_same_sample_measurement_assimilated_for_lower_bound": True,
        "reset_source_parity_found": reset_parity,
        "reset_sigma_min_lower": 1.0,
        "reset_lower_bound_requires_correction_radius": False,
        "PSD_covariance_floor_can_only_help_lower_bound": True,
        "marginal_inverse_metric_floor_inference_used": False,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    p3 = P3.build(path)
    h18 = H18.build(path)
    adet = ADET.build(path)
    pe = PE.build(path)
    process = PROCESS.build()
    dynamic = DYNAMIC.build(path)
    event = EVENT.build()
    bad = {
        "P3": P3.validate(p3), "H18": H18.validate(h18),
        "A21_detectability": ADET.validate(adet), "PE": PE.validate(pe),
        "process": PROCESS.validate(process), "dynamic": DYNAMIC.validate(dynamic),
        "event_algebra": EVENT.validate(event),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"finite-tau_b cascade prerequisites failed: {bad}")
    if p3["canonical_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("cascade bridge detached from canonical complete SEA3")
    if p3["P3_CONDITIONAL_SEA3_PASS"] is not True:
        raise RuntimeError("cascade bridge requires frozen conditional P3")
    if float(p3["useful_gate"]) != P3_FROZEN_GATE:
        raise RuntimeError("conditional P3 gate changed")
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("cascade bridge may not be trajectory fitted")
    delta_h = float(p3["modes"]["H18"]["relative_Riccati_injection_margin_lower"])
    if delta_h != P3_FROZEN_GATE:
        raise RuntimeError("H18 P3 gap changed")
    beta_b = _positive(adet["bias_homogeneous_energy_gap_lower"], "finite-tau_b bias energy gap")
    delta_a = TARGET_P4_CASCADE_GAP
    g2_limit = cascade_g2_limit_conservative(delta_h, beta_b, delta_a)
    g2_budget = down(0.25 * g2_limit)
    det_slack = cascade_determinant_slack_lower(delta_h, beta_b, delta_a, g2_budget)
    if not (g2_budget > 0.0 and det_slack > 0.0):
        raise RuntimeError("finite-tau_b cascade Schur budget is not strict")
    metric = _h18_information_metric_upper(domain, h18, pe, process, dynamic, event)
    if metric["reset_source_parity_found"] is not True:
        raise RuntimeError("shipping reset source parity changed")
    rejected_log_c = float(adet["triangular_detectability_observer"]["finite_coupling_log10_upper"])
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "paper_storage_family": "M_A(zeta)=diag(M_H(zeta),mu*I3)",
        "full_A21_dimension": 21, "H18_dimension": 18, "accelerometer_bias_dimension": 3,
        "P3_frozen_gate": P3_FROZEN_GATE,
        "P3_frozen_not_modified": True,
        "P3_CONDITIONAL_SEA3_PASS_consumed": True,
        "finite_tau_b_route": "ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION",
        "bias_energy_gap_lower": beta_b,
        "target_P4_cascade_gap": delta_a,
        "target_gap_stored_directly_not_as_1_minus_delta": True,
        "cascade_master_inequality": "c2/mu < ((delta_H-delta_A)*(beta_b-delta_A))/(1-delta_A)",
        "cascade_g2_limit_conservative_lower": g2_limit,
        "cascade_g2_budget": g2_budget,
        "cascade_determinant_slack_lower_if_budget_met": det_slack,
        "H18_metric_equivalence": metric,
        "actual_applied_R_S_retained": True,
        "all_due_S_updates_remain_in_complete_word": True,
        "all_valid_accelerometer_updates_remain_in_complete_word": True,
        "asynchronous_vector_events_remain_in_complete_word": True,
        "all_process_Q_and_covariance_floor_events_retained": True,
        "immediate_resets_retained": True,
        "full_H18_state_retained_including_aw": True,
        "full_ba_state_retained": True,
        "state_elimination_used": False,
        "a_w_Schur_elimination_used": False,
        "selected_S_word_used": False,
        "independent_tuner_or_RS_box_used": False,
        "trajectory_replay_used": False,
        "finite_harmonic_or_grid_source_used": False,
        "correction_radius_claim_used": False,
        "inverse_metric_floor_claim_used": False,
        "packet_count_coupling_bound_consumed": False,
        "rejected_detectability_packet_count_log10_C_upper": rejected_log_c,
        "joint_same_history_C_Hb_metric_bound_required": True,
        "joint_same_history_C_Hb_metric_bound_closed": False,
        "joint_C_Hb_c2_upper": None,
        "mu_log10_required": None,
        "full_rank_A21_cascade_metric_closed": False,
        "nonlinear_residual_sector_attached": False,
        "P4_promoted_here": False,
        "P5_may_start": False,
        "filter_changed": False,
        "quality_gate_changed": False,
        "next_obligation": "derive one source-uniform bound on ||M_H+^(1/2) C_Hb||^2 from the SAME complete-SEA3 word's joint Joseph/information geometry (including every actual-R_S S event), not from N*L^N; then choose finite mu from the certified Schur budget and attach the exact nonlinear residual sector",
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in ("P3_frozen_not_modified", "P3_CONDITIONAL_SEA3_PASS_consumed", "target_gap_stored_directly_not_as_1_minus_delta", "actual_applied_R_S_retained", "all_due_S_updates_remain_in_complete_word", "all_valid_accelerometer_updates_remain_in_complete_word", "asynchronous_vector_events_remain_in_complete_word", "all_process_Q_and_covariance_floor_events_retained", "immediate_resets_retained", "full_H18_state_retained_including_aw", "full_ba_state_retained", "joint_same_history_C_Hb_metric_bound_required"):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("state_elimination_used", "a_w_Schur_elimination_used", "selected_S_word_used", "independent_tuner_or_RS_box_used", "trajectory_replay_used", "finite_harmonic_or_grid_source_used", "correction_radius_claim_used", "inverse_metric_floor_claim_used", "packet_count_coupling_bound_consumed", "joint_same_history_C_Hb_metric_bound_closed", "full_rank_A21_cascade_metric_closed", "nonlinear_residual_sector_attached", "P4_promoted_here", "P5_may_start", "filter_changed", "quality_gate_changed"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("P3_frozen_gate", math.nan)) != P3_FROZEN_GATE:
        f.append("P3 gate changed")
    da = d.get("target_P4_cascade_gap")
    if not isinstance(da, (int, float)) or not (math.isfinite(float(da)) and 0.0 < float(da) < P3_FROZEN_GATE):
        f.append("target cascade gap is not strict below frozen P3 gap")
    for key in ("bias_energy_gap_lower", "cascade_g2_limit_conservative_lower", "cascade_g2_budget", "cascade_determinant_slack_lower_if_budget_met"):
        x = d.get(key)
        if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
            f.append(f"invalid positive cascade field {key}")
    metric = d.get("H18_metric_equivalence", {})
    for key in ("prediction_Q_lambda_min_lower", "Racc_variance_lower", "Rmag_variance_lower", "actual_applied_R_S_variance_lower", "posterior_covariance_lambda_min_lower", "M_H_lambda_max_upper", "M_H_lambda_min_lower_from_Pbar_trace_normalization", "reset_sigma_min_lower"):
        x = metric.get(key)
        if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
            f.append(f"invalid H18 metric field {key}")
    for key in ("Joseph_information_identity_used", "every_possible_same_sample_measurement_assimilated_for_lower_bound", "reset_source_parity_found", "PSD_covariance_floor_can_only_help_lower_bound"):
        if metric.get(key) is not True:
            f.append(f"H18 metric property lost: {key}")
    for key in ("reset_lower_bound_requires_correction_radius", "marginal_inverse_metric_floor_inference_used"):
        if metric.get(key) is not False:
            f.append(f"forbidden H18 metric shortcut re-entered: {key}")
    if d.get("joint_C_Hb_c2_upper") is not None or d.get("mu_log10_required") is not None:
        f.append("open joint C_Hb bound was populated without closure")
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
    print(json.dumps({"target_gap": d["target_P4_cascade_gap"], "bias_gap": d["bias_energy_gap_lower"], "g2_budget": d["cascade_g2_budget"], "H18_metric_upper": d["H18_metric_equivalence"]["M_H_lambda_max_upper"], "joint_C_Hb_closed": d["joint_same_history_C_Hb_metric_bound_closed"], "P4_promoted": d["P4_promoted_here"], "failures": failures}, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
