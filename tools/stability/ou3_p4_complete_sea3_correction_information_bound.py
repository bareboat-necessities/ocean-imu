#!/usr/bin/env python3
"""Operation-level information bounds; uniform P4 correction radius remains open.

For the ORIGINAL shipping P,H,R,S,K, y=H Phi(z), and J=y^T S^-1 y,

    ||E_i K y||^2 <= lambda_max(E_i P E_i^T) J,
    J <= Phi(z)^T P^-1 Phi(z).

These identities are valid, but they do not supply a source-uniform P ceiling
at every operation. A prediction-only extension of the H18 endpoint ceiling
omits intervening left-error resets, PSD a_w floor increments and the separate
A21 hybrid extension. It must not produce certified candidate energy balls.

Likewise, selecting a tube cell by its minimum contraction ratio does not
certify the minimum post-measurement process floor over all source cells.
For a defect supported on a state block, the appropriate precision block is
E_i P_J^-1 E_i^T, not (E_i P_J E_i^T)^-1. The primitive below encloses the
FULL posterior inverse and retains cross-covariances.

The bounds here accept actual operation matrices. They create no alternative
source and do not claim coverage of the complete SEA3 family. Uniform numeric
radii remain null until reset/floor-complete source-correlated coverage closes.
P3 is consumed unchanged and cannot be promoted by this module.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, IntervalMatrix, symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, matrix_symmetric_hull
import ou3_p4_complete_sea3_invariant_aw_coordinate as INVARIANT
import ou3_sea3_riccati_metric_p3 as P3

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 3
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_OPERATION_CORRECTION_INFORMATION_BOUND_V3"


def _checked_block(P: IntervalMatrix, indices: tuple[int, ...]) -> IntervalMatrix:
    n = len(P)
    if n not in (18, 21) or any(len(row) != n for row in P):
        raise ValueError("operation covariance must be full H18 or A21")
    if not indices or len(set(indices)) != len(indices) or any(
        not isinstance(i, int) or isinstance(i, bool) or not 0 <= i < n for i in indices
    ):
        raise ValueError("state-block indices must be unique and in range")
    if any(not isinstance(x, Interval) or not (math.isfinite(x.lo) and math.isfinite(x.hi))
           for row in P for x in row):
        raise ValueError("covariance entries must be finite intervals")
    for i in range(n):
        for j in range(i):
            if P[i][j].hi < P[j][i].lo or P[j][i].hi < P[i][j].lo:
                raise ValueError("covariance intervals cannot contain a symmetric matrix")
    symmetric = matrix_symmetric_hull(P)
    ok, _ = symmetric_positive_definite_ldlt(symmetric)
    if not ok:
        raise RuntimeError("full operation covariance did not pass interval LDLT")
    return symmetric


def correction_norm_squared_upper(
    P: IntervalMatrix, information: Interval, indices: tuple[int, ...],
) -> float:
    """Conditional ||E_i K y||^2 upper, with J from this SAME operation.

    J provenance (J=y^T S^-1 y and K=P H^T S^-1) is a caller precondition.
    No source-uniform coverage follows from supplying arbitrary point inputs.
    """
    checked = _checked_block(P, indices)
    if (not isinstance(information, Interval)
            or not (math.isfinite(information.lo) and math.isfinite(information.hi))
            or information.hi < 0.0):
        raise ValueError("information must enclose a finite nonnegative J")
    trace = Interval.point(0.0)
    for i in indices:
        trace = trace + checked[i][i]
    return (trace * Interval(0.0, information.hi)).hi


def defect_precision_trace_upper(P_J: IntervalMatrix, indices: tuple[int, ...]) -> float:
    """Bound lambda_max(E_i P_J^-1 E_i^T), retaining every cross term."""
    inverse = matrix_inverse_gauss_jordan(_checked_block(P_J, indices))
    trace = Interval.point(0.0)
    for i in indices:
        trace = trace + inverse[i][i]
    if not (math.isfinite(trace.hi) and trace.hi > 0.0):
        raise RuntimeError("full posterior precision bound is not finite positive")
    return trace.hi


def _status(candidates: list[dict], p3_pass: bool) -> dict:
    modes = {}
    for mode, dimension in (("H", 18), ("A", 21)):
        modes[mode] = {
            "dimension": dimension,
            "source_uniform_operation_covariance_ceiling_closed": False,
            "source_uniform_post_measurement_precision_bound_closed": False,
            "candidate_cells": [{
                "candidate_attitude_angle_deg": float(c["attitude_angle_deg"]),
                "candidate_cayley_norm_upper": float(c["cayley_norm_upper"]),
                "derived_metric_energy_radius_upper": None,
                "candidate_metric_energy_ball_certified": False,
            } for c in candidates],
        }
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "P3_frozen_not_modified": True,
        "P3_conditional_complete_SEA3_consumed": p3_pass,
        "complete_SEA3_word_retained": True,
        "all_valid_accelerometer_updates_remain_in_complete_word": True,
        "all_due_S_updates_and_actual_RS_remain_in_complete_word": True,
        "same_operation_correction_information_identity_valid": True,
        "same_shipping_P_H_R_K_S_cell_required": True,
        "full_matrix_correction_inequality": "d d^T <= J K S K^T <= J P",
        "measurement_information_operator_inequality": "H^T S^-1 H <= P^-1",
        "attitude_correction_information_inequality": "||E_theta K y||^2 <= lambda_max(P_theta_theta) J",
        "aw_correction_information_inequality": "||E_aw K y||^2 <= lambda_max(P_aw_aw) J",
        "storage_vector": "Phi(z)=z+E_aw*((Q_aw-I)*delta_a_w+e_eta)",
        "storage_is_original_physical_metric_isometry": False,
        "full_posterior_inverse_required_for_defect_cost": True,
        "source_uniform_covariance_reset_and_floor_coverage_closed": False,
        "candidate_metric_energy_balls_derived": False,
        "reset_transport_correction_radius_source_closed": False,
        "source_indexed_e_eta_transition_closed_here": False,
        "complete_word_nonlinear_dissipation_closed_here": False,
        "P4_promoted_here": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "trajectory_replay_used": False,
        "source_family_replaced": False,
        "old_scalar_Riccati_tube_margin_consumed": False,
        "standalone_eta_Rinv_budget_used": False,
        "packet_count_multiplier_used": False,
        "candidate_angles_deg": [float(c["attitude_angle_deg"]) for c in candidates],
        "modes": modes,
        "open_obligations": [
            "reset/floor-complete source-correlated operation covariance ceilings in H18 and A21",
            "source-uniform full posterior precision bound, not a selected-cell marginal floor",
            "full epsilon_aw transport and nonlinear-storage comparison through the same complete word",
        ],
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("correction-information bound must not be trajectory fitted")
    p3 = P3.build(path)
    invariant = INVARIANT.build(path)
    failures = P3.validate(p3) + INVARIANT.validate(invariant)
    if failures or p3["P3_CONDITIONAL_SEA3_PASS"] is not True:
        raise RuntimeError(f"unchanged conditional P3/invariant prerequisites failed: {failures}")
    return _status(invariant["measurement_linearizing_shift_bounds_reused_without_widening"], True)


def validate(d: dict) -> list[str]:
    failures = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        failures.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        failures.append("canonical source changed")
    for key in (
        "P3_frozen_not_modified", "P3_conditional_complete_SEA3_consumed", "complete_SEA3_word_retained",
        "all_valid_accelerometer_updates_remain_in_complete_word",
        "all_due_S_updates_and_actual_RS_remain_in_complete_word",
        "same_operation_correction_information_identity_valid", "full_posterior_inverse_required_for_defect_cost",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "storage_is_original_physical_metric_isometry", "source_uniform_covariance_reset_and_floor_coverage_closed",
        "candidate_metric_energy_balls_derived", "reset_transport_correction_radius_source_closed",
        "source_indexed_e_eta_transition_closed_here", "complete_word_nonlinear_dissipation_closed_here",
        "P4_promoted_here", "filter_changed", "declared_domain_changed", "trajectory_replay_used",
        "source_family_replaced", "old_scalar_Riccati_tube_margin_consumed",
        "standalone_eta_Rinv_budget_used", "packet_count_multiplier_used",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    if d.get("candidate_angles_deg") != [30.0, 25.0, 20.0, 15.0]:
        failures.append("candidate finite-angle cells changed")
    for mode, dimension in (("H", 18), ("A", 21)):
        m = d.get("modes", {}).get(mode, {})
        cells = m.get("candidate_cells", [])
        if m.get("dimension") != dimension or len(cells) != 4:
            failures.append(f"{mode} candidate cells missing")
        for key in ("source_uniform_operation_covariance_ceiling_closed",
                    "source_uniform_post_measurement_precision_bound_closed"):
            if m.get(key) is not False:
                failures.append(f"{mode} unsupported {key}")
        if [row.get("candidate_attitude_angle_deg") for row in cells] != [30.0, 25.0, 20.0, 15.0]:
            failures.append(f"{mode} candidate angles changed")
        for row in cells:
            q = row.get("candidate_cayley_norm_upper")
            if not isinstance(q, (int, float)) or not (math.isfinite(float(q)) and 0.0 < q < 1.0):
                failures.append(f"{mode} candidate Cayley radius invalid")
            if row.get("derived_metric_energy_radius_upper") is not None:
                failures.append(f"{mode} unsupported numeric candidate radius")
            if row.get("candidate_metric_energy_ball_certified") is not False:
                failures.append(f"{mode} candidate radius falsely certified")
    if len(d.get("open_obligations", [])) != 3:
        failures.append("uniform correction-bound obligations lost")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d.update(validation_pass=not failures, validation_failures=failures)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"reset_radius_source_closed": False, "open_obligations": d["open_obligations"],
                      "validation_failures": failures}, indent=2))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
