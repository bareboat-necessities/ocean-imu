#!/usr/bin/env python3
"""Full-rank differential P4 metric induced by the exact nonlinear Phi map.

This is the theorem-facing use of the accelerometer-linearizing coordinate.  It
is deliberately NOT the falsified finite storage ``Phi(z)^T P^-1 Phi(z)``.
Instead let

    T(z,zeta) = D_z Phi(z,zeta),
    M(z,zeta) = T^T P(zeta)^-1 T.

For the retained coordinate

    Phi_aw = Q_aw(c,zeta) delta_a_w + e_eta(c,zeta),

all other state coordinates are unchanged.  Hence T is block triangular; its
``a_w,a_w`` diagonal block is the orthogonal rotation ``Q_aw``.  Therefore

    det T = det Q_aw = +1

for every finite Cayley state and every admitted source rotation.  The metric
is full rank on H18 and A21 and eliminates no active coordinate.  At the zero
error equilibrium T=I, so it reduces exactly to the frozen P3 moving metric.

For a complete nonlinear word F_W with Jacobian J_W=D F_W, the contraction
condition is the full matrix inequality

    rho M_0 - J_W^T M_1 J_W > 0,

with one rho<1 uniformly over the same complete BRMM word and a certified
finite-error region.  Every shipping gain/covariance, and especially every due
S=0 operation with its actual applied per-axis R_S, remains in J_W and P.

The interval routines below construct T using the repository's outward
first-order AD and test the complete differential inequality by full interval
LDL^T.  They do not use finite differences, replay fitting, packet-count
remainders, state elimination, correction radii, or marginal metric floors.
A source-uniform complete-word J_W enclosure is still required before P4 can
promote.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import (
    Interval,
    matrix_identity,
    matrix_mul,
    matrix_sub,
    matrix_transpose,
    symmetric_positive_definite_ldlt,
)
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_interval_ad as AD
import ou3_p4_complete_brmm_accelerometer_operation_coordinate as ACC
import ou3_p4_complete_brmm_measurement_linearizing_aw_coordinate as AW
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_brmm_riccati_metric_p3 as P3

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_BRMM_PHI_PULLBACK_DIFFERENTIAL_METRIC_V1"


def _shape(a) -> tuple[int, int]:
    rows = len(a)
    cols = len(a[0]) if rows else 0
    if any(len(row) != cols for row in a):
        raise ValueError("ragged matrix")
    return rows, cols


def _ad_const_matrix(a: Sequence[Sequence[float]], n: int):
    if len(a) != 3 or any(len(row) != 3 for row in a):
        raise ValueError("R_hat must be 3x3")
    return [[AD.constant(float(a[i][j]), n) for j in range(3)] for i in range(3)]


def _ad_transpose(a):
    return [[a[i][j] for i in range(len(a))] for j in range(len(a[0]))]


def _ad_identity3(n: int):
    return [[AD.constant(1.0 if i == j else 0.0, n) for j in range(3)] for i in range(3)]


def _ad_skew(c):
    z = AD.constant(0.0, c[0].n)
    x, y, w = c
    return [[z, -w, y], [w, z, -x], [-y, x, z]]


def phi_aw_jacobian_interval(
    c: Sequence[Interval],
    delta_aw: Sequence[Interval],
    f_hat: Sequence[Interval],
    R_hat: Sequence[Sequence[float]],
) -> list[list[Interval]]:
    """Outward Jacobian D_(c,delta_aw) Phi_aw on one source/error cell.

    ``R_hat`` is one source-provided orthogonal rotation.  Universal source
    coverage is a caller obligation; this evaluator introduces no independent
    rotation box.
    """
    if len(c) != 3 or len(delta_aw) != 3 or len(f_hat) != 3:
        raise ValueError("c, delta_aw and f_hat must be three-vectors")
    if any(not isinstance(x, Interval) for x in [*c, *delta_aw, *f_hat]):
        raise TypeError("phi differential cells must use outward Interval values")

    n = 6
    cad = AD.independent_vector(c, n=n, offset=0)
    daad = AD.independent_vector(delta_aw, n=n, offset=3)
    fad = [AD.constant(x, n) for x in f_hat]
    R = _ad_const_matrix(R_hat, n)
    Rt = _ad_transpose(R)
    E = AD.rotation_from_cayley(cad)
    Q = AD.matmul(AD.matmul(Rt, E), R)
    C = _ad_skew(cad)
    I3 = _ad_identity3(n)
    EmIminusC = [
        [E[i][j] - I3[i][j] - C[i][j] for j in range(3)]
        for i in range(3)
    ]
    qda = AD.matvec(Q, daad)
    eta = AD.matvec(Rt, AD.matvec(EmIminusC, fad))
    phi_aw = [qda[i] + eta[i] for i in range(3)]
    return AD.jacobian(phi_aw)


def full_phi_jacobian_interval(
    dimension: int,
    c: Sequence[Interval],
    delta_aw: Sequence[Interval],
    f_hat: Sequence[Interval],
    R_hat: Sequence[Sequence[float]],
):
    """Build full H18/A21 T=D Phi; only the a_w row block differs from I."""
    if dimension not in (18, 21):
        raise ValueError("Phi differential metric is only H18 or A21")
    local = phi_aw_jacobian_interval(c, delta_aw, f_hat, R_hat)
    T = matrix_identity(dimension)
    for i in range(3):
        for j in range(3):
            T[15 + i][j] = local[i][j]
            T[15 + i][15 + j] = local[i][3 + j]
    return T


def pullback_metric(P_inverse, T):
    """Return the full interval pullback metric T^T P^-1 T."""
    n, m = _shape(P_inverse)
    if n == 0 or n != m or _shape(T) != (n, n):
        raise ValueError("pullback metric dimension mismatch")
    return matrix_symmetric_hull(matrix_mul(matrix_mul(matrix_transpose(T), P_inverse), T))


def differential_contraction_matrix(
    J_word,
    P0_inverse,
    P1_inverse,
    T0,
    T1,
    rho: float,
):
    """Return rho*M0-J_W^T*M1*J_W for one validated complete-word cell."""
    if not (0.0 < float(rho) < 1.0):
        raise ValueError("differential contraction rho must lie strictly in (0,1)")
    n0, n1 = _shape(P0_inverse)[0], _shape(P1_inverse)[0]
    if _shape(P0_inverse) != (n0, n0) or _shape(P1_inverse) != (n1, n1):
        raise ValueError("endpoint inverse metrics must be square")
    if _shape(T0) != (n0, n0) or _shape(T1) != (n1, n1):
        raise ValueError("Phi pullback Jacobian dimension mismatch")
    if _shape(J_word) != (n1, n0):
        raise ValueError("complete-word Jacobian dimension mismatch")
    M0 = pullback_metric(P0_inverse, T0)
    M1 = pullback_metric(P1_inverse, T1)
    transported = matrix_mul(matrix_mul(matrix_transpose(J_word), M1), J_word)
    ri = Interval.outward_bounds(float(rho), float(rho))
    scaled = [[ri * M0[i][j] for j in range(n0)] for i in range(n0)]
    return matrix_symmetric_hull(matrix_sub(scaled, transported))


def certify_differential_contraction(
    J_word,
    P0_inverse,
    P1_inverse,
    T0,
    T1,
    rho: float,
) -> tuple[bool, list[Interval]]:
    """Full interval-LDLT certificate for the pullback differential inequality."""
    D = differential_contraction_matrix(J_word, P0_inverse, P1_inverse, T0, T1, rho)
    return symmetric_positive_definite_ldlt(D)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    p3 = P3.build(path)
    acc = ACC.build(path)
    aw = AW.build(path)
    cayley = CAYLEY.build(path)
    failures = (
        [f"P3: {x}" for x in P3.validate(p3)]
        + [f"accelerometer: {x}" for x in ACC.validate(acc)]
        + [f"full-shift: {x}" for x in AW.validate(aw)]
        + [f"Cayley: {x}" for x in CAYLEY.validate(cayley)]
    )
    if failures:
        raise RuntimeError(f"Phi differential metric prerequisites failed: {failures}")
    if p3.get("P3_CONDITIONAL_BRMM_PASS") is not True:
        raise RuntimeError("Phi differential metric requires frozen conditional P3")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "P3_frozen_not_modified": True,
        "P3_delta_consumed": 1.0e-18,
        "filter_changed": False,
        "declared_domain_changed": False,
        "source_family_replaced": False,
        "trajectory_replay_used": False,
        "metric_type": "FULL_STATE_SOURCE_DEPENDENT_RIEMANNIAN_PULLBACK",
        "metric_definition": "M(z,zeta)=D_Phi(z,zeta)^T*P(zeta)^-1*D_Phi(z,zeta)",
        "finite_Phi_storage_used_as_Lyapunov_function": False,
        "all_active_states_retained": True,
        "H18_full_rank": True,
        "A21_full_rank": True,
        "Phi_jacobian_block_triangular": True,
        "Phi_aw_diagonal_block": "Q_aw=R_hat^T*E*R_hat",
        "Phi_aw_diagonal_block_orthogonal": True,
        "Phi_jacobian_determinant_exact": 1.0,
        "Phi_jacobian_globally_nonsingular_on_finite_Cayley_chart": True,
        "metric_reduces_to_P3_at_zero_error": True,
        "same_complete_BRMM_word_required": True,
        "same_frontend_tuner_and_covariance_history_required": True,
        "all_due_S_updates_with_actual_applied_RS_required": True,
        "all_valid_accelerometer_updates_required": True,
        "all_process_Q_floor_and_reset_events_required": True,
        "H_to_A_rectangular_differential_event_required": True,
        "packet_count_remainder_budget_used": False,
        "state_elimination_used": False,
        "correction_radius_claim_used": False,
        "inverse_metric_floor_claim_used": False,
        "outward_interval_AD_coordinate_jacobian_available": True,
        "full_interval_LDLT_differential_test_available": True,
        "differential_word_inequality": (
            "rho*M0-D F_W(z)^T*M1*D F_W(z) > 0, 0<rho<1, for every admitted complete BRMM word/state cell"
        ),
        "source_uniform_complete_word_Jacobian_enclosed": False,
        "source_uniform_pullback_differential_contraction_closed": False,
        "P4_promoted_here": False,
        "P5_may_start_here": False,
        "next_obligation": (
            "propagate one outward full-state Jacobian enclosure through the literal complete BRMM H18/A21 word, retaining every actual-R_S S event; test the pullback differential inequality directly by full interval LDLT on [30,25,20,15] degree cells"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    if float(d.get("P3_delta_consumed", 0.0)) != 1.0e-18:
        f.append("frozen P3 delta changed")
    for key in (
        "P3_frozen_not_modified", "all_active_states_retained", "H18_full_rank", "A21_full_rank",
        "Phi_jacobian_block_triangular", "Phi_aw_diagonal_block_orthogonal",
        "Phi_jacobian_globally_nonsingular_on_finite_Cayley_chart",
        "metric_reduces_to_P3_at_zero_error", "same_complete_BRMM_word_required",
        "same_frontend_tuner_and_covariance_history_required",
        "all_due_S_updates_with_actual_applied_RS_required", "all_valid_accelerometer_updates_required",
        "all_process_Q_floor_and_reset_events_required", "H_to_A_rectangular_differential_event_required",
        "outward_interval_AD_coordinate_jacobian_available", "full_interval_LDLT_differential_test_available",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "filter_changed", "declared_domain_changed", "source_family_replaced", "trajectory_replay_used",
        "finite_Phi_storage_used_as_Lyapunov_function", "packet_count_remainder_budget_used",
        "state_elimination_used", "correction_radius_claim_used", "inverse_metric_floor_claim_used",
        "source_uniform_complete_word_Jacobian_enclosed",
        "source_uniform_pullback_differential_contraction_closed", "P4_promoted_here", "P5_may_start_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("Phi_jacobian_determinant_exact", 0.0)) != 1.0:
        f.append("Phi Jacobian determinant is not exactly one")
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
        "metric": d["metric_type"],
        "full_rank_H18_A21": d["H18_full_rank"] and d["A21_full_rank"],
        "actual_RS_required": d["all_due_S_updates_with_actual_applied_RS_required"],
        "complete_word_Jacobian_closed": d["source_uniform_complete_word_Jacobian_enclosed"],
        "P4_promoted_here": d["P4_promoted_here"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
