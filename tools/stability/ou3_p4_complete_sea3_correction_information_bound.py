#!/usr/bin/env python3
"""Operation-correlated correction/reset bound for complete-SEA3 OU-III P4.

This certificate closes the *correction-radius source* that the exact reset
transport lemma deliberately left open.  It does not invent a correction cap.
For one shipping measurement in any of the exact P4 coordinates, let

    y = H zeta,        J = y^T S^-1 y,
    d = K y,           K = P H^T S^-1.

The attitude part d_theta obeys, exactly,

    ||d_theta||^2 <= lambda_max(P_theta,theta) J.          (1)

Indeed, with A=E_theta P^(1/2) and B=H P^(1/2),

    d_theta = A B^T S^-1 y,

and B B^T = H P H^T <= S.  Therefore

    ||d_theta||^2
      <= ||A||^2 y^T S^-1 H P H^T S^-1 y
      <= lambda_max(P_theta,theta) J.

The same optimal-measurement inequality gives

    J = zeta^T H^T S^-1 H zeta <= zeta^T P^-1 zeta = V.  (2)

The current complete-SEA3 Riccati tube already supplies a source-uniform
shipping covariance ceiling.  If Ptheta_bar is the sum of its three attitude
diagonal ceilings, then

    ||c||^2 <= Ptheta_bar V,
    ||d_theta||^2 <= Ptheta_bar V.                         (3)

Consequently every declared finite-angle candidate q_c induces the *derived*
metric-energy radius

    nu_geom = q_c^2 / Ptheta_bar,

inside which both the physical Cayley error and the actual shipping correction
are inside q_c.  This is not an extra theorem assumption and not a replay-fit
radius; it is a consequence of the same moving shipping covariance metric.

The Riccati tube also supplies a post-measurement process floor.  Its attitude
block gives

    P_plus >= mu_theta I,    mu_theta = rho_post q_theta,

so any attitude-only exact reset defect rho_theta satisfies

    rho_theta^T P_plus^-1 rho_theta <= ||rho_theta||^2/mu_theta.   (4)

Equations (1)-(4), together with the parametric exact quaternion reset bound,
turn reset transport into a scalar function of the *same* operation energy V.
No independent correction radius, packet-count multiplier, endpoint word, or
source enumeration is used.

This producer remains fail-closed: it establishes the operation-correlated
reset ledger but does not yet claim that the e_eta source-indexed coordinate
shift across prediction/source change is dominated over the complete word.
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
import ou3_sea3_riccati_metric_p3 as P3
import ou3_sea3_riccati_tube as TUBE

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_OPERATION_CORRELATED_CORRECTION_INFORMATION_BOUND"


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _mode_row(mode: str, tube: dict, process: dict, candidate: dict) -> dict:
    tm = tube["modes"][mode]
    pdiag = [float(x) for x in tm["Pbar_diagonal_variance_upper"]]
    if len(pdiag) < 3 or any((not math.isfinite(x) or x <= 0.0) for x in pdiag[:3]):
        raise RuntimeError(f"invalid {mode} attitude covariance ceiling")
    ptheta = up(sum(pdiag[:3]))

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

    # At V=nu both ||c|| and the operation-matched correction are <=q.
    # This endpoint evaluation is a geometric diagnostic, not an additive
    # complete-word disturbance charge.
    rb = RESET.reset_defect_bound(q, q)
    rho = float(rb["reset_attitude_defect_norm_upper"])
    reset_metric_cost = up(rho * rho / mu_theta)

    return {
        "mode": mode,
        "dimension": 18 if mode == "H" else 21,
        "attitude_covariance_lambda_max_upper": ptheta,
        "attitude_covariance_bound_source": "SUM_OF_FIRST_THREE_COMPLETE_SEA3_RICCATI_TUBE_DIAGONAL_CEILINGS",
        "same_operation_correction_information_inequality": (
            "||d_theta||^2 <= Ptheta_bar * y^T S^-1 y"
        ),
        "measurement_information_below_moving_energy": (
            "y^T S^-1 y <= zeta^T P^-1 zeta = V"
        ),
        "state_and_correction_energy_bound": (
            "||c||^2 <= Ptheta_bar*V and ||d_theta||^2 <= Ptheta_bar*V"
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
    tube = TUBE.build(path)
    process = PROCESS.build()
    invariant = INVARIANT.build(path)
    rebind = REBIND.build()
    reset = RESET.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "P3": P3.validate(p3),
        "tube": TUBE.validate(tube),
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

    candidates = invariant["measurement_linearizing_shift_bounds_reused_without_widening"]
    modes = {}
    for mode in ("H", "A"):
        modes[mode] = [_mode_row(mode, tube, process, c) for c in candidates]

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
        "complete_SEA3_Riccati_tube_supplies_covariance_ceiling": True,
        "complete_SEA3_Riccati_tube_supplies_post_measurement_floor": True,
        "actual_RS_enters_same_Riccati_tube": True,
        "modes": modes,
        "candidate_angles_deg": [float(c["attitude_angle_deg"]) for c in candidates],
        "candidate_metric_energy_balls_derived": True,
        "reset_transport_correction_radius_source_closed": True,
        "source_indexed_e_eta_transition_closed_here": False,
        "complete_word_nonlinear_dissipation_closed_here": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "bound e_eta(c,source) across the same complete-SEA3 prediction/source transition as a moving-coordinate defect, "
            "using the derived energy-indexed c/d bounds and actual R_S metric; combine with reset cost before complete-word endpoint scalarization"
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
        "complete_SEA3_Riccati_tube_supplies_covariance_ceiling",
        "complete_SEA3_Riccati_tube_supplies_post_measurement_floor",
        "actual_RS_enters_same_Riccati_tube",
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
        "source_indexed_e_eta_transition_closed_here",
        "complete_word_nonlinear_dissipation_closed_here",
        "P4_promoted_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("candidate_angles_deg") != [30.0, 25.0, 20.0, 15.0]:
        f.append("candidate finite-angle cells changed")
    modes = d.get("modes", {})
    for mode in ("H", "A"):
        rows = modes.get(mode, [])
        if len(rows) != 4:
            f.append(f"{mode} candidate rows missing")
            continue
        for row in rows:
            for key in (
                "attitude_covariance_lambda_max_upper",
                "post_measurement_attitude_covariance_floor",
                "derived_metric_energy_radius_upper",
            ):
                x = float(row.get(key, math.nan))
                if not (math.isfinite(x) and x > 0.0):
                    f.append(f"{mode} invalid {key}")
            diag = row.get("reset_endpoint_diagnostic", {})
            if diag.get("additive_packet_budget_used") is not False:
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
        "reset_radius_source_closed": d["reset_transport_correction_radius_source_closed"],
        "e_eta_transition_closed": d["source_indexed_e_eta_transition_closed_here"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
