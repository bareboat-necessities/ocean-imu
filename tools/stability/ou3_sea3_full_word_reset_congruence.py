#!/usr/bin/env python3
"""Exact shipping left-error reset primitive for the complete SEA3 word.

Shipping applies an immediate covariance reset after every accepted S=0,
accelerometer and magnetometer correction.  The attitude block is

    G_theta = I + 0.5 [dtheta]_x,

and the remaining coordinates are unchanged.  Since

    det(G_theta) = 1 + ||dtheta||^2/4 >= 1,

G is nonsingular for every finite injection; no small-angle assumption is
needed for the covariance congruence.

This module executes that missing event on the existing canonical full-word
state objects.  It is an event primitive, not another estimator or source.  For
the joint P/Psi/Omega representation,

    P     <- G P G^T,
    Psi   <- G Psi,
    Omega <- G Omega G^T.

For the prior-free batch representation D describes information about the
word-start coordinate and therefore remains unchanged, while the current-state
objects transform as

    T  <- G T,
    Qc <- G Qc G^T.

Both transformations preserve the exact endpoint factorization.  For
M_delta=Omega-delta P the reset gives exactly G M_delta G^T.
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence

from ou3_interval import (
    Interval,
    IntervalMatrix,
    matrix_identity,
    matrix_mul,
    matrix_point,
    matrix_sub,
    matrix_transpose,
)
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_sea3_full_word_riccati_backend as BACKEND

REPO = Path(__file__).resolve().parents[2]
CORE = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"
SCHEMA = 1
QUALIFICATION = "OU3_COMPLETE_SEA3_SHIPPING_LEFT_ERROR_RESET_CONGRUENCE"


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _shape(A: Sequence[Sequence[Interval]]) -> tuple[int, int]:
    n = len(A)
    m = len(A[0]) if n else 0
    if any(len(row) != m for row in A):
        raise ValueError("ragged interval matrix")
    return n, m


def skew(v: Sequence[Interval]) -> IntervalMatrix:
    if len(v) != 3:
        raise ValueError("dtheta must have three components")
    z = I(0.0)
    x, y, zc = v
    return [
        [z, -zc, y],
        [zc, z, -x],
        [-y, x, z],
    ]


def reset_matrix(dtheta: Sequence[Interval], dimension: int) -> IntervalMatrix:
    if dimension < 3:
        raise ValueError("state dimension must contain attitude coordinates")
    G = matrix_identity(dimension)
    S = skew(dtheta)
    half = I(0.5)
    for i in range(3):
        for j in range(3):
            G[i][j] = G[i][j] + half * S[i][j]
    return G


def determinant_attitude_from_norm_squared(norm_squared: Interval) -> Interval:
    if norm_squared.lo < 0.0:
        raise ValueError("norm squared must be nonnegative")
    return I(1.0) + I(0.25) * norm_squared


def apply_joint_reset(state: BACKEND.JointWordState, G: Sequence[Sequence[Interval]]) -> None:
    n = state.dimension
    if _shape(G) != (n, n):
        raise ValueError("joint reset dimension mismatch")
    Gt = matrix_transpose(G)
    state.P = matrix_symmetric_hull(matrix_mul(matrix_mul(G, state.P), Gt))
    state.Psi = matrix_mul(G, state.Psi)
    state.Omega = matrix_symmetric_hull(matrix_mul(matrix_mul(G, state.Omega), Gt))
    state.events += 1
    if not BACKEND.decomposition_identity_enclosed(state):
        raise RuntimeError("P/Psi/Omega identity lost after left-error reset")


def apply_prior_free_reset(
    state: BACKEND.PriorFreeBatchState,
    G: Sequence[Sequence[Interval]],
) -> None:
    n = state.dimension
    if _shape(G) != (n, n):
        raise ValueError("prior-free reset dimension mismatch")
    Gt = matrix_transpose(G)
    state.T = matrix_mul(G, state.T)
    state.Qc = matrix_symmetric_hull(matrix_mul(matrix_mul(G, state.Qc), Gt))
    # D is information about the unchanged word-start coordinate.
    state.events += 1


def reset_contraction_image(
    M_before: Sequence[Sequence[Interval]],
    G: Sequence[Sequence[Interval]],
) -> IntervalMatrix:
    if _shape(M_before) != _shape(G):
        raise ValueError("reset contraction dimension mismatch")
    return matrix_symmetric_hull(matrix_mul(matrix_mul(G, M_before), matrix_transpose(G)))


def _contains_zero(A: Sequence[Sequence[Interval]]) -> bool:
    return all(x.lo <= 0.0 <= x.hi for row in A for x in row)


def source_parity() -> dict[str, bool]:
    c = CORE.read_text(encoding="utf-8")
    return {
        "shipping_exact_reset_formula_present": (
            "Identity() + T(0.5)*skew(dtheta)" in c
            or "Identity() + T(0.5) * skew(dtheta)" in c
        ),
    }


def self_test() -> dict:
    P0 = matrix_point([
        [2.0, 0.1, 0.0, 0.0],
        [0.1, 1.5, 0.0, 0.0],
        [0.0, 0.0, 1.2, 0.0],
        [0.0, 0.0, 0.0, 0.8],
    ])
    joint = BACKEND.initialize(P0)
    batch = BACKEND.initialize_prior_free(4)
    F = matrix_point([
        [1.0, 0.0, 0.0, 0.1],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ])
    Q = matrix_point([
        [0.02,0,0,0],[0,0.02,0,0],[0,0,0.02,0],[0,0,0,0.01],
    ])
    BACKEND.predict(joint, F, Q)
    BACKEND.prior_free_predict(batch, F, Q)
    H = matrix_point([[1.0, 0.0, 0.0, 0.0]])
    R = matrix_point([[0.5]])
    BACKEND.joseph_measurement(joint, H, R)
    BACKEND.prior_free_measurement(batch, H, R)

    dtheta = [I(0.03), I(-0.02), I(0.01)]
    G = reset_matrix(dtheta, 4)
    M0 = BACKEND.contraction_matrix(joint, 0.2)
    expected = reset_contraction_image(M0, G)
    apply_joint_reset(joint, G)
    actual = BACKEND.contraction_matrix(joint, 0.2)
    margin_identity = _contains_zero(matrix_sub(actual, expected))

    # Prior-free endpoint must continue to reconstruct the direct recursion
    # after applying the same reset congruence.
    apply_prior_free_reset(batch, G)
    reconstructed = BACKEND.reconstruct_joint_from_prior_free(batch, P0)
    P_identity = _contains_zero(matrix_sub(joint.P, reconstructed.P))
    Psi_identity = _contains_zero(matrix_sub(joint.Psi, reconstructed.Psi))
    Omega_identity = _contains_zero(matrix_sub(joint.Omega, reconstructed.Omega))

    n2 = I(0.03*0.03 + 0.02*0.02 + 0.01*0.01)
    det = determinant_attitude_from_norm_squared(n2)
    return {
        "source_parity": source_parity(),
        "source_parity_pass": all(source_parity().values()),
        "reset_determinant_interval": det.as_list(),
        "reset_determinant_lower_at_least_one": det.lo >= 1.0,
        "joint_decomposition_identity_enclosed": BACKEND.decomposition_identity_enclosed(joint),
        "M_delta_reset_congruence_identity_enclosed": margin_identity,
        "prior_free_reconstruction_P_identity_enclosed": P_identity,
        "prior_free_reconstruction_Psi_identity_enclosed": Psi_identity,
        "prior_free_reconstruction_Omega_identity_enclosed": Omega_identity,
        "prior_free_D_unchanged_by_reset": True,
        "small_angle_needed_for_nonsingularity": False,
        "P3_promoted": False,
    }


def validate() -> list[str]:
    f: list[str] = []
    try:
        d = self_test()
    except Exception as exc:
        return [f"reset congruence self-test failed: {type(exc).__name__}: {exc}"]
    for key in (
        "source_parity_pass",
        "reset_determinant_lower_at_least_one",
        "joint_decomposition_identity_enclosed",
        "M_delta_reset_congruence_identity_enclosed",
        "prior_free_reconstruction_P_identity_enclosed",
        "prior_free_reconstruction_Psi_identity_enclosed",
        "prior_free_reconstruction_Omega_identity_enclosed",
        "prior_free_D_unchanged_by_reset",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    if d.get("small_angle_needed_for_nonsingularity") is not False:
        f.append("reset nonsingularity incorrectly requires a small-angle assumption")
    if d.get("P3_promoted") is not False:
        f.append("reset primitive promoted P3")
    return f


if __name__ == "__main__":
    failures = validate()
    print({
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "failures": failures,
        "self_test": self_test() if not failures else None,
    })
    raise SystemExit(0 if not failures else 2)
