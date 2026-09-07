#!/usr/bin/env python3
"""Canonical unpromoted P4: full-state differential contraction on complete SEA3.

Finite endpoint energies in raw or nonlinear-Phi coordinates are diagnostics,
not the theorem architecture.  The canonical nonlinear metric is the full-rank
Riemannian pullback

    M(z,zeta)=D Phi(z,zeta)^T P(zeta)^-1 D Phi(z,zeta),

where Phi is the exact accelerometer-linearizing coordinate using the ORIGINAL
shipping H/P/K/S and the full shift

    epsilon_aw=(Q_aw-I)delta_a_w+e_eta.

D Phi is block triangular with the orthogonal Q_aw block on a_w, so the metric
is nonsingular in all H18/A21 coordinates and reduces exactly to the frozen P3
moving metric at zero error.  No state is eliminated.

For the complete nonlinear word F_W the P4 condition is

    rho M_0 - D F_W(z)^T M_1 D F_W(z) > 0,    0<rho<1,

uniformly over every admitted complete SEA3 word and every state in one
certified finite-angle cell.  The complete Jacobian must retain the same
frontend/tuner/covariance history, full F/Q, every valid accelerometer update,
every due S=0 update with the actual applied anisotropic SpectralMSE R_S,
asynchronous vector events, covariance floors, immediate resets, and the
separate H18->A21 rectangular differential event.

P3 remains frozen at delta=1e-18.  P4/P5 remain open until a source-uniform
outward complete-word Jacobian enclosure passes the full interval-LDLT
inequality.  No replay, packet-count remainder, correction radius, inverse
metric floor, finite-state endpoint optimization or eliminated-state shortcut
can promote this gate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_sea3_riccati_metric_p3 as P3
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_complete_sea3_phi_differential_metric as DIFF

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 8
QUALIFICATION = "OU3_SEA3_FULL_STATE_DIFFERENTIAL_P4_V8"


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    p3 = P3.build(path)
    cayley = CAYLEY.build(path)
    diff = DIFF.build(path)
    failures = (
        [f"P3: {x}" for x in P3.validate(p3)]
        + [f"Cayley: {x}" for x in CAYLEY.validate(cayley)]
        + [f"differential metric: {x}" for x in DIFF.validate(diff)]
    )
    if failures:
        raise RuntimeError(f"canonical differential P4 prerequisites failed: {failures}")

    p3_pass = bool(p3["P3_CONDITIONAL_SEA3_PASS"])
    h_delta = float(p3["modes"]["H18"]["relative_Riccati_injection_margin_lower"])
    a_delta = float(p3["modes"]["A21"]["relative_Riccati_injection_margin_lower"])
    jacobian_closed = bool(diff["source_uniform_complete_word_Jacobian_enclosed"])
    contraction_closed = bool(diff["source_uniform_pullback_differential_contraction_closed"])
    p4_pass = bool(p3_pass and jacobian_closed and contraction_closed)

    fail_reasons: list[str] = []
    if not p3_pass:
        fail_reasons.append("frozen conditional complete-SEA3 P3 prerequisite is not closed")
    if not jacobian_closed:
        fail_reasons.append(
            "source-uniform outward full-state Jacobian of the literal complete SEA3 H18/A21 word is not yet enclosed"
        )
    if not contraction_closed:
        fail_reasons.append(
            "full interval-LDLT pullback differential contraction is not yet strict on a declared [30,25,20,15] degree cell"
        )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_P4_architecture": "FULL_STATE_COMPLETE_SEA3_DIFFERENTIAL_PULLBACK",
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "source_family_replaced": False,
        "P3_CONDITIONAL_SEA3_PASS_consumed": p3_pass,
        "P3_DEPLOYMENT_PASS_consumed_as_if_closed": False,
        "P3_H18_delta_consumed": h_delta,
        "P3_A21_delta_consumed": a_delta,
        "P3_frozen_not_modified": True,
        "outer_angle_rad": float(cayley["outer_angle_rad"]),
        "candidate_angles_deg": [30.0, 25.0, 20.0, 15.0],
        "differential_metric_qualification": diff["qualification"],
        "differential_metric_type": diff["metric_type"],
        "differential_metric_definition": diff["metric_definition"],
        "finite_Phi_storage_used_as_Lyapunov_function": False,
        "finite_raw_endpoint_storage_used_as_P4_certificate": False,
        "all_active_states_retained": bool(diff["all_active_states_retained"]),
        "H18_full_rank": bool(diff["H18_full_rank"]),
        "A21_full_rank": bool(diff["A21_full_rank"]),
        "Phi_jacobian_determinant_exact": diff["Phi_jacobian_determinant_exact"],
        "metric_reduces_exactly_to_P3_at_zero_error": bool(diff["metric_reduces_to_P3_at_zero_error"]),
        "same_complete_SEA3_word_required": True,
        "same_frontend_tuner_covariance_history_required": True,
        "all_due_S_updates_with_actual_applied_RS_required": True,
        "all_valid_accelerometer_updates_required": True,
        "all_process_Q_floor_reset_events_required": True,
        "H_to_A_rectangular_differential_event_required": True,
        "outward_interval_AD_used": True,
        "full_interval_LDLT_required": True,
        "differential_word_inequality": diff["differential_word_inequality"],
        "source_uniform_complete_word_Jacobian_enclosed": jacobian_closed,
        "source_uniform_pullback_differential_contraction_closed": contraction_closed,
        "packet_count_remainder_budget_used": False,
        "packetwise_remainder_norm_sum_used": False,
        "state_elimination_used": False,
        "a_w_Schur_final_certificate_used": False,
        "correction_radius_claim_used": False,
        "inverse_metric_floor_claim_used": False,
        "independent_RS_schedule_used": False,
        "point_word_rho_used_to_promote": False,
        "longer_point_window_optimization_used_to_promote": False,
        "P4_FINITE_WINDOW_CLOSED": p4_pass,
        "P4_CANONICAL_PASS": p4_pass,
        "P5_MAY_START": p4_pass,
        "P4_CANONICAL_FAIL_REASONS": fail_reasons,
        "next_obligation": (
            "construct one source-uniform outward D F_W enclosure for the literal complete SEA3 word using the existing interval-AD primitives; retain actual applied R_S inside every due S event; test rho*M0-J^T*M1*J by full interval LDLT from 30 degrees downward; if no cell closes, perform theorem/metric failure analysis rather than retuning or replacing the source"
            if p3_pass else "repair only frozen P3 prerequisite failure"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_P4_architecture") != "FULL_STATE_COMPLETE_SEA3_DIFFERENTIAL_PULLBACK":
        f.append("canonical P4 architecture changed")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "source_generated_not_trajectory_fit", "P3_CONDITIONAL_SEA3_PASS_consumed",
        "P3_frozen_not_modified", "all_active_states_retained", "H18_full_rank", "A21_full_rank",
        "metric_reduces_exactly_to_P3_at_zero_error", "same_complete_SEA3_word_required",
        "same_frontend_tuner_covariance_history_required",
        "all_due_S_updates_with_actual_applied_RS_required", "all_valid_accelerometer_updates_required",
        "all_process_Q_floor_reset_events_required", "H_to_A_rectangular_differential_event_required",
        "outward_interval_AD_used", "full_interval_LDLT_required",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk", "source_family_replaced",
        "P3_DEPLOYMENT_PASS_consumed_as_if_closed", "finite_Phi_storage_used_as_Lyapunov_function",
        "finite_raw_endpoint_storage_used_as_P4_certificate",
        "source_uniform_complete_word_Jacobian_enclosed",
        "source_uniform_pullback_differential_contraction_closed",
        "packet_count_remainder_budget_used", "packetwise_remainder_norm_sum_used",
        "state_elimination_used", "a_w_Schur_final_certificate_used", "correction_radius_claim_used",
        "inverse_metric_floor_claim_used", "independent_RS_schedule_used", "point_word_rho_used_to_promote",
        "longer_point_window_optimization_used_to_promote", "P4_FINITE_WINDOW_CLOSED",
        "P4_CANONICAL_PASS", "P5_MAY_START",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    for key in ("P3_H18_delta_consumed", "P3_A21_delta_consumed"):
        if float(d.get(key, 0.0)) != 1.0e-18:
            f.append(f"{key} changed from frozen delta")
    if float(d.get("Phi_jacobian_determinant_exact", 0.0)) != 1.0:
        f.append("Phi differential coordinate lost exact full rank")
    if float(d.get("outer_angle_rad", 0.0)) < 0.80:
        f.append("outer Cayley sector shrank")
    if d.get("candidate_angles_deg") != [30.0, 25.0, 20.0, 15.0]:
        f.append("declared P4 candidate sequence changed")
    reasons = d.get("P4_CANONICAL_FAIL_REASONS", [])
    if len(reasons) != 2 or not any("jacobian" in x.lower() for x in reasons) or not any(
        "ldlt" in x.lower() or "contraction" in x.lower() for x in reasons
    ):
        f.append("canonical P4 must expose Jacobian and differential-contraction blockers")
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
        "architecture": d["canonical_P4_architecture"],
        "metric": d["differential_metric_type"],
        "actual_RS_required": d["all_due_S_updates_with_actual_applied_RS_required"],
        "complete_word_Jacobian_closed": d["source_uniform_complete_word_Jacobian_enclosed"],
        "P4_CANONICAL_PASS": d["P4_CANONICAL_PASS"],
        "fail_reasons": d["P4_CANONICAL_FAIL_REASONS"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
