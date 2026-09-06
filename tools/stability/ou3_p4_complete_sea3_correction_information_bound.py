#!/usr/bin/env python3
"""Operation-correlated correction/reset bound for complete-SEA3 OU-III P4.

This certificate closes the correction-radius source left open by the exact
reset transport lemma.  It does not invent a correction cap and it does not
consume the obsolete scalar Riccati-tube contraction margin.

For one shipping measurement, in any of the exact P4 congruent coordinates,

    y = H zeta,       J = y^T S^-1 y,
    d = K y,          K = P H^T S^-1.

For every state block selector E_i,

    ||E_i d||^2 <= lambda_max(E_i P E_i^T) J,              (1)

because H P H^T <= S.  The same matrix inequality gives

    J <= zeta^T P^-1 zeta = V.                             (2)

The covariance upper used in (1) is the one that belongs to the *actual P3
closure route*.  The H18 prior-free completion constructs a source-uniform
same-3-second-word diffuse-prior endpoint covariance upper using the directional
vector-PE/four-S information composition.  That composition retains every due
S update and the actual SpectralMSE R_S; selected four-S records are only a PSD
information witness and do not replace the complete word.

Starting from that recurrent endpoint bound, measurements can only reduce
covariance.  Over at most one following 3 s word we therefore bound the
attitude and a_w marginals by prediction-only propagation.  This gives
source-uniform operation ceilings Ptheta_bar and Paw_bar without using the old
full-trace tube scalarization.

The current-source Riccati-tube algebra is retained for one purpose only: its
post-measurement process/Joseph floor is a valid lower bound independently of
whether its old scalar contraction margin passes.  Hence

    P_plus,theta >= mu_theta I,

and an attitude-only exact reset defect rho_theta has metric cost

    rho_theta^T P_plus^-1 rho_theta <= ||rho_theta||^2/mu_theta.  (3)

Combining (1)-(2), every finite-angle candidate q_c induces the *derived*
moving-metric energy radius

    nu_geom = q_c^2 / Ptheta_bar,

inside which both the physical Cayley error and the actual Kalman attitude
correction are <= q_c.  This is a consequence of the same shipping covariance
metric, not a replay fit or a new theorem assumption.

The certificate remains fail-closed: source-indexed transport of the exact
finite-angle e_eta shift still has to be dominated using the same operation
information and actual R_S regularized a_w metric before P4 can promote.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_p4_complete_sea3_invariant_aw_coordinate as INVARIANT
import ou3_p4_exact_reset_transport as RESET
import ou3_p4_moving_metric_rebind as REBIND
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_h18_information_composition as HINFO
import ou3_sea3_h18_prior_free_completion as HPRIOR
import ou3_sea3_riccati_metric_p3 as P3
import ou3_sea3_riccati_tube as TUBE

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_OPERATION_CORRELATED_CORRECTION_INFORMATION_BOUND_V2"
WORD_HORIZON_S = 3.0


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _sum3(values: list[float], offset: int, label: str) -> float:
    xs = [float(x) for x in values[offset:offset + 3]]
    if len(xs) != 3 or any(not (math.isfinite(x) and x > 0.0) for x in xs):
        raise RuntimeError(f"invalid {label} covariance marginal ceiling")
    return up(sum(xs))


def _same_word_operation_ceilings(path: Path, dynamic: dict, process: dict, hinfo: dict) -> dict:
    """Propagate the P3 same-word endpoint upper over at most one next word."""
    pbar = HPRIOR._same_word_covariance_upper(path, dynamic, process, hinfo)
    diag = [float(x) for x in pbar["Pbar_diagonal_variance_upper"]]
    if len(diag) != 18:
        raise RuntimeError("P3 same-word H18 covariance upper changed dimension")

    ptheta0 = _sum3(diag, 0, "endpoint theta")
    pbg0 = _sum3(diag, 3, "endpoint gyro bias")
    paw0 = _sum3(diag, 15, "endpoint aw")

    pc = process["source_constants"]
    sigma_g_axis = pc.get("gyro_noise_density_rad_sqrt_s_per_axis")
    if not isinstance(sigma_g_axis, list) or len(sigma_g_axis) != 3:
        raise RuntimeError("per-axis shipping gyro-noise source constants missing")
    qg = up(max(float(x) for x in sigma_g_axis) ** 2)
    qb = up(float(pc["gyro_bias_rw_variance_density"]))
    if not (qg > 0.0 and qb > 0.0):
        raise RuntimeError("invalid attitude-process upper inputs")

    T = WORD_HORIZON_S
    # theta(t)=R theta0 + B bg0 + process, ||R||=1, ||B||<=T.
    # PSD Cauchy/Frobenius domination gives the factor two on the two initial
    # marginal contributions.  The final term is the direct gyro noise plus
    # integrated gyro-bias random walk, summed over three axes.
    theta_process_trace = up(3.0 * (qg * T + qb * T ** 3 / 3.0))
    ptheta_operation = up(2.0 * ptheta0 + 2.0 * T * T * pbg0 + theta_process_trace)

    sigma_hi = float(dynamic["dynamic_invariant"]["sigma_aw_filter_mps2"][1])
    if not (math.isfinite(sigma_hi) and sigma_hi > 0.0):
        raise RuntimeError("invalid complete-SEA3 sigma_aw ceiling")
    # For a time-varying stable scalar OU recursion with every stationary
    # variance <= sigma_hi^2, accumulated process covariance over any horizon
    # is <= sigma_hi^2 per axis.  Retaining the endpoint marginal as well is a
    # conservative prediction-only upper valid under source motion.
    paw_operation = up(paw0 + 3.0 * sigma_hi * sigma_hi)

    return {
        "same_word_endpoint_covariance_upper": pbar,
        "same_word_endpoint_theta_trace_upper": ptheta0,
        "same_word_endpoint_bg_trace_upper": pbg0,
        "same_word_endpoint_aw_trace_upper": paw0,
        "operation_horizon_s_upper": T,
        "attitude_prediction_process_trace_upper": theta_process_trace,
        "attitude_covariance_lambda_max_upper": ptheta_operation,
        "aw_covariance_lambda_max_upper": paw_operation,
        "measurements_inside_following_word_can_only_reduce_covariance": True,
        "source_motion_inside_word_allowed": True,
        "same_word_directional_information_source": hinfo["qualification"],
        "same_word_actual_RS_consumed": bool(hinfo["actual_applied_SpectralMSE_R_S_consumed"]),
    }


def _mode_row(mode: str, tube: dict, process: dict, ceilings: dict, candidate: dict) -> dict:
    ptheta = float(ceilings["attitude_covariance_lambda_max_upper"])
    paw = float(ceilings["aw_covariance_lambda_max_upper"])
    tm = tube["modes"][mode]
    worst = tm["worst_current_source_cell"]
    rho_post = float(worst["post_measurement_scaled_Omega_lambda_min_lower"])
    qtheta = float(process["attitude_gyro_bias"]["theta_diagonal_lower"])
    if not (math.isfinite(rho_post) and rho_post > 0.0 and math.isfinite(qtheta) and qtheta > 0.0):
        raise RuntimeError(f"invalid {mode} post-measurement attitude floor inputs")
    mu_theta = down(rho_post * qtheta)
    if not mu_theta > 0.0:
        raise RuntimeError(f"{mode} post-measurement attitude covariance floor lost positivity")

    q = float(candidate["cayley_norm_upper"])
    if not (math.isfinite(q) and 0.0 < q < 1.0):
        raise RuntimeError("candidate Cayley radius outside q<1")
    nu = down(q * q / ptheta)
    if not (math.isfinite(nu) and nu > 0.0):
        raise RuntimeError(f"{mode} derived metric-energy radius is not positive")

    rb = RESET.reset_defect_bound(q, q)
    rho = float(rb["reset_attitude_defect_norm_upper"])
    reset_metric_cost = up(rho * rho / mu_theta)

    return {
        "mode": mode,
        "dimension": 18 if mode == "H" else 21,
        "attitude_covariance_lambda_max_upper": ptheta,
        "aw_covariance_lambda_max_upper": paw,
        "attitude_covariance_bound_source": "P3_SAME_WORD_DIRECTIONAL_INFORMATION_ENDPOINT_PLUS_ONE_WORD_PREDICTION",
        "aw_covariance_bound_source": "P3_SAME_WORD_ACTUAL_RS_TRANSLATION_ESTIMATOR_PLUS_OU_PREDICTION",
        "same_operation_attitude_correction_information_inequality": (
            "||d_theta||^2 <= Ptheta_bar * y^T S^-1 y"
        ),
        "same_operation_aw_correction_information_inequality": (
            "||d_aw||^2 <= Paw_bar * y^T S^-1 y"
        ),
        "measurement_information_below_moving_energy": (
            "y^T S^-1 y <= zeta^T P^-1 zeta = V"
        ),
        "state_and_correction_energy_bound": (
            "||c||^2,||d_theta||^2 <= Ptheta_bar*V"
        ),
        "post_measurement_scaled_process_floor": rho_post,
        "attitude_process_scale_qtheta": qtheta,
        "post_measurement_attitude_covariance_floor": mu_theta,
        "reset_defect_metric_inequality": (
            "rho_theta^T Pplus^-1 rho_theta <= ||rho_theta||^2/mu_theta"
        ),
        "candidate_attitude_angle_deg": float(candidate["attitude_angle_deg"]),
        "candidate_cayley_norm_upper": q,
        "derived_metric_energy_radius_upper": nu,
        "derived_not_declared_correction_radius_at_energy_boundary": q,
        "reset_endpoint_diagnostic": {
            "correction_radius": q,
            "reset_attitude_defect_norm_upper": rho,
            "reset_defect_metric_cost_upper": reset_metric_cost,
            "additive_packet_budget_used": False,
        },
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("correction-information certificate must not be trajectory fitted")

    complete = COMPLETE.build(path)
    p3 = P3.build(path)
    dynamic = DYNAMIC.build(path)
    hinfo = HINFO.build(path)
    tube = TUBE.build(path)
    process = PROCESS.build()
    invariant = INVARIANT.build(path)
    rebind = REBIND.build()
    reset = RESET.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "P3": P3.validate(p3),
        "dynamic": DYNAMIC.validate(dynamic),
        "H18_directional_information": HINFO.validate(hinfo),
        "tube_floor_algebra": TUBE.validate(tube),
        "process": PROCESS.validate(process),
        "invariant_aw": INVARIANT.validate(invariant),
        "moving_metric_rebind": REBIND.validate(rebind),
        "reset": RESET.validate(reset),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"correction-information prerequisites failed: {bad}")
    if complete["canonical_P3_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("canonical complete SEA3 source changed")
    if p3["P3_CONDITIONAL_SEA3_PASS"] is not True:
        raise RuntimeError("frozen conditional complete-SEA3 P3 is not closed")
    if hinfo["actual_applied_SpectralMSE_R_S_consumed"] is not True:
        raise RuntimeError("directional same-word covariance route lost actual R_S")

    ceilings = _same_word_operation_ceilings(path, dynamic, process, hinfo)
    candidates = invariant["measurement_linearizing_shift_bounds_reused_without_widening"]
    modes = {
        mode: [_mode_row(mode, tube, process, ceilings, c) for c in candidates]
        for mode in ("H", "A")
    }

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "source_family_replaced": False,
        "P3_frozen_not_modified": True,
        "P3_conditional_complete_SEA3_consumed": True,
        "complete_SEA3_word_retained": True,
        "all_valid_accelerometer_updates_remain_in_complete_word": True,
        "all_due_S_updates_and_actual_RS_remain_in_complete_word": True,
        "invariant_aw_normal_form_consumed": True,
        "same_shipping_P_H_R_K_S_cell_used": True,
        "same_operation_correction_radius_derived_from_information": True,
        "independent_global_correction_radius_assumed": False,
        "endpoint_source_word_enumeration_used": False,
        "packet_count_multiplier_used": False,
        "standalone_eta_Rinv_budget_used": False,
        "measurement_information_operator_inequality": "H^T S^-1 H <= P^-1",
        "attitude_correction_information_inequality": (
            "||E_theta K y||^2 <= lambda_max(P_theta_theta) y^T S^-1 y"
        ),
        "aw_correction_information_inequality": (
            "||E_aw K y||^2 <= lambda_max(P_aw_aw) y^T S^-1 y"
        ),
        "same_word_directional_information_covariance_ceiling_consumed": True,
        "H18_directional_RS_regularizer_consumed": True,
        "actual_RS_enters_same_word_covariance_ceiling": True,
        "old_scalar_Riccati_tube_margin_consumed": False,
        "Riccati_tube_used_only_for_post_measurement_floor": True,
        "complete_SEA3_Riccati_tube_supplies_post_measurement_floor": True,
        "operation_covariance_ceilings": ceilings,
        "modes": modes,
        "candidate_angles_deg": [float(c["attitude_angle_deg"]) for c in candidates],
        "candidate_metric_energy_balls_derived": True,
        "reset_transport_correction_radius_source_closed": True,
        "source_indexed_e_eta_transition_closed_here": False,
        "complete_word_nonlinear_dissipation_closed_here": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "bound e_eta across the same complete-SEA3 correction/reset and source transition using the derived "
            "energy-indexed theta/aw correction bounds and actual-R_S aw metric; scalarize only after the complete word"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "source_generated_not_trajectory_fit",
        "P3_frozen_not_modified",
        "P3_conditional_complete_SEA3_consumed",
        "complete_SEA3_word_retained",
        "all_valid_accelerometer_updates_remain_in_complete_word",
        "all_due_S_updates_and_actual_RS_remain_in_complete_word",
        "invariant_aw_normal_form_consumed",
        "same_shipping_P_H_R_K_S_cell_used",
        "same_operation_correction_radius_derived_from_information",
        "same_word_directional_information_covariance_ceiling_consumed",
        "H18_directional_RS_regularizer_consumed",
        "actual_RS_enters_same_word_covariance_ceiling",
        "Riccati_tube_used_only_for_post_measurement_floor",
        "complete_SEA3_Riccati_tube_supplies_post_measurement_floor",
        "candidate_metric_energy_balls_derived",
        "reset_transport_correction_radius_source_closed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used",
        "filter_changed",
        "declared_domain_changed",
        "source_family_replaced",
        "independent_global_correction_radius_assumed",
        "endpoint_source_word_enumeration_used",
        "packet_count_multiplier_used",
        "standalone_eta_Rinv_budget_used",
        "old_scalar_Riccati_tube_margin_consumed",
        "source_indexed_e_eta_transition_closed_here",
        "complete_word_nonlinear_dissipation_closed_here",
        "P4_promoted_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("candidate_angles_deg") != [30.0, 25.0, 20.0, 15.0]:
        f.append("candidate finite-angle cells changed")
    ceilings = d.get("operation_covariance_ceilings", {})
    for key in ("attitude_covariance_lambda_max_upper", "aw_covariance_lambda_max_upper"):
        x = float(ceilings.get(key, math.nan))
        if not (math.isfinite(x) and x > 0.0):
            f.append(f"invalid operation covariance ceiling {key}")
    if ceilings.get("same_word_actual_RS_consumed") is not True:
        f.append("same-word operation ceiling lost actual R_S")
    modes = d.get("modes", {})
    for mode in ("H", "A"):
        rows = modes.get(mode, [])
        if len(rows) != 4:
            f.append(f"{mode} candidate rows missing")
            continue
        for row in rows:
            for key in (
                "attitude_covariance_lambda_max_upper",
                "aw_covariance_lambda_max_upper",
                "post_measurement_attitude_covariance_floor",
                "derived_metric_energy_radius_upper",
            ):
                x = float(row.get(key, math.nan))
                if not (math.isfinite(x) and x > 0.0):
                    f.append(f"{mode} invalid {key}")
            if row.get("reset_endpoint_diagnostic", {}).get("additive_packet_budget_used") is not False:
                f.append(f"{mode} reset diagnostic became additive packet budget")
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
    print(json.dumps({
        "candidate_angles_deg": d["candidate_angles_deg"],
        "operation_covariance_ceilings": d["operation_covariance_ceilings"],
        "reset_radius_source_closed": d["reset_transport_correction_radius_source_closed"],
        "e_eta_transition_closed": d["source_indexed_e_eta_transition_closed_here"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
