#!/usr/bin/env python3
"""Exact full a_w shift and correction/reset transport for OU-III P4.

In original physical coordinates z=(c,...,da[,db_a]), let C=[c]x,
E=R_true R_hat^T, Q=R_hat^T E R_hat, and f=R_hat(a_hat-g). Define

    e_eta=R_hat^T ((E-I)-C)f,
    epsilon_aw=(Q-I)da+e_eta,
    phi=Phi(z)=z+E_aw epsilon_aw.

Then y=H phi is the exact accelerometer residual with the ORIGINAL shipping
H,P,R,K,S. The pure finite-angle shift e_eta alone is insufficient: the
mixed wave-error shift (Q-I)da must also be retained.

For A=I-KH, the Joseph identity is exactly

    phi^T P^-1 phi-(A phi)^T P_J^-1(A phi)=y^T S^-1 y.

This is a nonlinear storage choice, not a covariance congruence proving
Phi(z)^T P^-1 Phi(z)=z^T P^-1 z. Uniform metric comparison is still open.

Under the physical correction d=K y, t=z-d, z_plus=G t+rho and P_r=G P_J G^T,
G is identity on the physical a_w block. Thus the exact transport is

    Phi_plus=G A phi+xi,
    xi=rho+E_aw(epsilon_plus-epsilon_minus).

For one measurement da_plus=da-d_aw, its mixed-shift difference is

    (Q_plus-Q_minus)da-(Q_plus-I)d_aw.

That term is required even when the pure e_eta transport is bounded. With
u=A phi and b=G^-1 xi the signed ledger is

    V_plus-V_minus=-J+2 u^T P_J^-1 b+b^T P_J^-1 b.

The same full inverse, including all cross terms, must be used for both
terms; a marginal covariance lower bound cannot replace it.

For a physical prediction z_plus=F z+rho_p, the corresponding identity is

    Phi_plus=F Phi+xi_p,
    xi_p=rho_p+E_aw epsilon_plus-F E_aw epsilon_minus.

The last term includes v,p,S components: F E_aw is not E_aw. Those components
are precisely part of the coupling to the retained S=0/R_S regularizer. For a
source-only step use its actual map; for H18->A21 use the actual rectangular
lift J and separate covariance/initialization operation, not an H18 ceiling.
The evaluator below performs only this algebra. Same-history source admission,
physical defects, covariance/hybrid binding and complete-word dissipation are
not established by it.
Candidate bounds below bound ONLY e_eta, not epsilon_aw. No promotion occurs.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import (
    Interval, IntervalMatrix, matrix_add, matrix_mul, matrix_sub, up,
)
import ou3_p4_complete_sea3_accelerometer_operation_coordinate as ACC
import ou3_p4_complete_sea3_vector_remainder_geometry as GEOM
import ou3_p4_exact_reset_transport as RESET
import ou3_sea3_complete_source as COMPLETE

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_MEASUREMENT_LINEARIZING_AW_COORDINATE_V2"


def evaluate_full_shift_transport(
    linear_map: IntervalMatrix,
    physical_defect: IntervalMatrix,
    epsilon_before: IntervalMatrix,
    epsilon_after: IntervalMatrix,
) -> IntervalMatrix:
    """Transport the full shift through one supplied physical operation.

    For z_after=L z_before+r, this returns
    xi=r+E_aw,out epsilon_after-L E_aw,in epsilon_before.
    For correction/reset, z_before here is t=z-d and Phi_before is A Phi;
    L=G. For prediction L=F, NOT G or I. The v,p,S rows of F E_aw must
    survive. An H18-to-A21 lift is rectangular and must be supplied literally.

    Both shifts include (Q_aw-I) delta_a_w. The supplied map/defect/shifts
    must come from the SAME physical event. Shape checks establish neither
    that provenance nor SEA3 membership, and this routine emits no certificate.
    Covariance transport and actual applied R_S belong to the unchanged word.
    """
    n_out = len(linear_map)
    n_in = len(linear_map[0]) if linear_map else 0
    if (n_out, n_in) not in ((18, 18), (21, 21), (21, 18)):
        raise ValueError("linear_map must be H18, A21, or the 21x18 H-to-A lift")
    for name, value, rows, cols in (
        ("linear_map", linear_map, n_out, n_in),
        ("physical_defect", physical_defect, n_out, 1),
        ("epsilon_before", epsilon_before, 3, 1),
        ("epsilon_after", epsilon_after, 3, 1),
    ):
        if len(value) != rows or any(len(row) != cols for row in value):
            raise ValueError(f"{name} must have shape {rows}x{cols}")
        if any(not isinstance(x, Interval) or not (math.isfinite(x.lo) and math.isfinite(x.hi))
               or x.lo > x.hi for row in value for x in row):
            raise ValueError(f"{name} must contain finite intervals")
    before = [[Interval.point(0.0)] for _ in range(n_in)]
    after = [[Interval.point(0.0)] for _ in range(n_out)]
    before[15:18] = [row[:] for row in epsilon_before]
    after[15:18] = [row[:] for row in epsilon_after]
    return matrix_add(physical_defect, matrix_sub(after, matrix_mul(linear_map, before)))


def _candidate(row: dict, force_upper: float) -> dict:
    """Build outward scalar bounds for one retained finite-angle cell."""
    s = float(row["sin_half_angle_upper"])
    c = float(row["cos_half_angle_lower"])
    if not (0.0 <= s < 1.0 and 0.0 < c <= 1.0):
        raise RuntimeError("invalid validated half-angle geometry")
    q = up(2.0 * s / c)
    # Exact: ||eta||=(q/2)||y_R|| and ||y_R||<=2 sin(theta/2)||f||.
    shift = up(q * s * force_upper)

    # From eta=.5 C y_R:
    # d eta=.5 dC y_R+.5 C d y_R.
    # ||y_R||<=q||f|| and the Cayley map derivative obeys
    # ||dE||<=1+q/2 on q<1, giving
    # ||D eta|| <= q(1+q/4)||f||.
    lip = up(q * up(1.0 + 0.25 * q) * force_upper)
    return {
        "attitude_angle_deg": float(row["attitude_angle_deg"]),
        "cayley_norm_upper": q,
        "e_eta_norm_upper_mps2": shift,
        "e_eta_local_lipschitz_upper_mps2_per_cayley": lip,
        "shift_bound_vanishes_at_zero_angle": True,
        "lipschitz_bound_vanishes_at_zero_angle": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("measurement-linearizing coordinate must not be trajectory fitted")

    complete = COMPLETE.build(path)
    acc = ACC.build(path)
    geom = GEOM.build(path)
    reset = RESET.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "accelerometer_operation_coordinate": ACC.validate(acc),
        "vector_remainder_geometry": GEOM.validate(geom),
        "reset_transport": RESET.validate(reset),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"measurement-linearizing-coordinate prerequisites failed: {bad}")
    if complete["canonical_P3_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("canonical complete SEA3 source changed")

    fmax = float(domain["normal_live"]["specific_force_norm_upper_mps2"])
    if not (math.isfinite(fmax) and fmax > 0.0):
        raise RuntimeError("invalid Normal-Live force upper")
    rows = [_candidate(r, fmax) for r in geom["candidate_cells"]]

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
        "all_valid_accelerometer_updates_remain_in_complete_word": True,
        "all_due_S_updates_and_actual_RS_remain_in_complete_word": True,
        "operation_coordinate_consumed": True,
        "nonlinear_storage_metric_equivalence_to_original_closed": False,
        "full_mixed_aw_shift_retained_in_transport": True,
        "operation_coordinate_transform_back_exact": bool(
            acc["state_coordinate_transform"]["transform_back_after_measurement_exact"]
        ),
        "exact_cayley_remainder_identity": "eta=0.5*[c]x*((E-I)f_hat)",
        "exact_cayley_remainder_identity_closed": True,
        "measurement_linearizing_coordinate": {
            "e_eta": "R_hat^T eta",
            "Phi_aw": "u_aw+e_eta",
            "u_aw": "Q_aw delta_a_w",
            "Q_aw": "R_hat^T E R_hat",
            "triangular_in_aw": True,
            "globally_invertible_in_aw_for_each_finite_c": True,
            "inverse": "delta_a_w=Q_aw^T*(Phi_aw-e_eta(c,zeta))",
            "epsilon_aw": "(Q_aw-I)*delta_a_w+e_eta",
            "full_shift_includes_mixed_aw_error": True,
            "R_hat_orthogonal": True,
            "Q_aw_orthogonal": True,
        },
        "exact_accelerometer_residual_identity": (
            "r_a=[c]x f_hat+R_hat Phi_aw+delta_b_a=H_a Phi(z)"
        ),
        "exact_shipping_tangent_H_used": True,
        "accelerometer_eta_declared_zero_in_original_coordinate": False,
        "accelerometer_eta_dropped_from_physics": False,
        "standalone_eta_Rinv_packet_budget_used": False,
        "packet_count_multiplier_used": False,
        "actual_RS_regularizes_same_aw_coordinate_family": True,
        "joseph_P_H_R_K_S_unchanged": True,
        "phi_storage_exact_Joseph_identity": (
            "phi^T P^-1 phi-((I-KH)phi)^T(P+)^-1((I-KH)phi)"
            "=y^T S^-1 y, y=H phi"
        ),
        "phi_storage_has_no_standalone_eta_penalty": True,
        "combined_correction_reset_transport": {
            "t": "z-Ky",
            "z_plus": "G t+rho",
            "phi_linear_posterior": "(I-KH)phi=t+E_aw epsilon_minus",
            "phi_plus": "G(I-KH)phi+xi",
            "xi": "rho+E_aw(epsilon_plus-epsilon_minus)",
            "mixed_shift_difference": "(Q_plus-Q_minus)*delta_a_w-(Q_plus-I)*d_aw",
            "signed_energy_ledger": "V_plus-V_minus=-J+2*u^T*P_J^-1*b+b^T*P_J^-1*b; b=G^-1*xi",
            "G_identity_on_aw_coordinate": True,
            "reset_covariance_congruence_exact": not RESET.validate(reset) and bool(reset["exact_reset_congruence_identity"]),
            "G_inverse_operator_norm_exact": float(reset["reset_inverse_operator_norm_upper"]),
        },
        "prediction_shift_transport": {
            "physical_map": "z_plus=F*z+rho_p",
            "phi_plus": "F*Phi+xi_p",
            "xi_p": "rho_p+E_aw*epsilon_plus-F*E_aw*epsilon_minus",
            "full_v_p_S_aw_prediction_columns_required": True,
            "reset_only_shift_difference_valid_for_prediction": False,
            "source_uniform_prediction_defect_closed": False,
        },
        "hybrid_shift_transport": {
            "physical_map": "z_A_plus=J_HA*z_H+rho_HA",
            "xi_HA": "rho_HA+E_aw_A*epsilon_A_plus-J_HA*E_aw_H*epsilon_H_minus",
            "actual_rectangular_lift_and_covariance_seed_required": True,
            "H18_ceiling_reused_for_A21": False,
            "source_uniform_hybrid_defect_closed": False,
        },
        "source_indexed_shift_must_persist_across_complete_word": True,
        "packetwise_shift_reset_to_zero_allowed": False,
        "complete_SEA3_source_transition_must_drive_e_plus_minus": True,
        "candidate_cells": rows,
        "candidate_shift_bounds_cover_pure_e_eta_only": True,
        "candidate_shift_bounds_require_nominal_force_bound": True,
        "nominal_force_bound_certified_here": False,
        "shipping_Joseph_binding_closed": False,
        "shipping_Joseph_binding_scope": "SOURCE_UNIFORM_NONLINEAR_TRANSPORT",
        "candidate_shift_bounds_cover_full_epsilon_aw": False,
        "outer_geometry_retained": bool(geom["outer_geometry_cell"]["sector_is_homogeneous_quadratic"]),
        "measurement_remainder_obligation_reduced_to_coordinate_transport": True,
        "complete_source_correlated_transport_defect_closed_here": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "evaluate the signed full-shift xi=rho+E_aw(epsilon_plus-epsilon_minus) on the same complete SEA3 "
            "source-indexed word, retaining the actual R_S directional metric and the "
            "same P,H,R,K,S cell; do not return to standalone eta or packet-count budgets"
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
        "all_valid_accelerometer_updates_remain_in_complete_word",
        "all_due_S_updates_and_actual_RS_remain_in_complete_word",
        "operation_coordinate_consumed",
        "operation_coordinate_transform_back_exact",
        "exact_cayley_remainder_identity_closed",
        "exact_shipping_tangent_H_used",
        "actual_RS_regularizes_same_aw_coordinate_family",
        "joseph_P_H_R_K_S_unchanged",
        "phi_storage_has_no_standalone_eta_penalty",
        "source_indexed_shift_must_persist_across_complete_word",
        "complete_SEA3_source_transition_must_drive_e_plus_minus",
        "full_mixed_aw_shift_retained_in_transport",
        "candidate_shift_bounds_cover_pure_e_eta_only",
        "measurement_remainder_obligation_reduced_to_coordinate_transport",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "nominal_force_bound_certified_here",
        "shipping_Joseph_binding_closed",
        "trajectory_replay_used",
        "filter_changed",
        "declared_domain_changed",
        "source_family_replaced",
        "accelerometer_eta_declared_zero_in_original_coordinate",
        "nonlinear_storage_metric_equivalence_to_original_closed",
        "candidate_shift_bounds_cover_full_epsilon_aw",
        "accelerometer_eta_dropped_from_physics",
        "standalone_eta_Rinv_packet_budget_used",
        "packet_count_multiplier_used",
        "packetwise_shift_reset_to_zero_allowed",
        "complete_source_correlated_transport_defect_closed_here",
        "P4_promoted_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    coord = d.get("measurement_linearizing_coordinate", {})
    for key in (
        "triangular_in_aw",
        "globally_invertible_in_aw_for_each_finite_c",
        "R_hat_orthogonal",
        "Q_aw_orthogonal",
    ):
        if coord.get(key) is not True:
            f.append(f"coordinate {key} is not true")
    tr = d.get("combined_correction_reset_transport", {})
    if tr.get("G_identity_on_aw_coordinate") is not True:
        f.append("reset G does not preserve aw coordinate")
    if tr.get("reset_covariance_congruence_exact") is not True:
        f.append("reset covariance congruence is not exact")
    if float(tr.get("G_inverse_operator_norm_exact", math.inf)) != 1.0:
        f.append("reset inverse norm changed")
    prediction = d.get("prediction_shift_transport", {})
    if prediction.get("xi_p") != "rho_p+E_aw*epsilon_plus-F*E_aw*epsilon_minus":
        f.append("prediction dropped the full F E_aw shift transport")
    if prediction.get("full_v_p_S_aw_prediction_columns_required") is not True:
        f.append("prediction lost its coupling into the S regularizer")
    for key in ("reset_only_shift_difference_valid_for_prediction", "source_uniform_prediction_defect_closed"):
        if prediction.get(key) is not False:
            f.append(f"prediction {key} is not false")
    hybrid = d.get("hybrid_shift_transport", {})
    if hybrid.get("actual_rectangular_lift_and_covariance_seed_required") is not True:
        f.append("hybrid transport detached from its actual lift/covariance seed")
    for key in ("H18_ceiling_reused_for_A21", "source_uniform_hybrid_defect_closed"):
        if hybrid.get(key) is not False:
            f.append(f"hybrid {key} is not false")
    rows = d.get("candidate_cells", [])
    if [r.get("attitude_angle_deg") for r in rows] != [30.0, 25.0, 20.0, 15.0]:
        f.append("candidate finite-angle cells changed")
    for row in rows:
        q = float(row.get("cayley_norm_upper", math.inf))
        s = float(row.get("e_eta_norm_upper_mps2", math.inf))
        L = float(row.get("e_eta_local_lipschitz_upper_mps2_per_cayley", math.inf))
        if not (0.0 < q < 1.0 and math.isfinite(s) and s > 0.0 and math.isfinite(L) and L > 0.0):
            f.append(f"invalid coordinate-shift cell {row.get('attitude_angle_deg')}")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "candidate_cells": d["candidate_cells"],
        "exact_residual": d["exact_accelerometer_residual_identity"],
        "transport_xi": d["combined_correction_reset_transport"]["xi"],
        "actual_RS_retained": d["all_due_S_updates_and_actual_RS_remain_in_complete_word"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
