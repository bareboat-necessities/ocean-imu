#!/usr/bin/env python3
"""Dependency-tight facade for the canonical full OU-III joint Riccati backend.

This is not a second certificate.  It executes the same P/Psi/Omega recursion
and the same Kalman gain as ``ou3_sea3_full_word_riccati_backend``.  The only
change is enclosure tightening by intersecting redundant *exact* identities.

For each measurement the shipping Joseph covariance enclosure

    P+ = A P A^T + K R K^T

is intersected with the algebraically identical Schur enclosure

    P+ = P - P H^T S^-1 H P.

The full cross-covariance action is retained.  This matters especially for the
S=0 pseudo update: the actual applied R_S appears in S=P_SS+R_S and therefore
regularizes the complete P[:,S] correction instead of being lost to repeated
natural-interval dependency wrapping.

Omega is likewise intersected with the exact decomposition identity

    Omega = P - Psi P0 Psi^T.

Both operands are independently validated interval enclosures of the same exact
matrix.  Their intersection is therefore a strictly tighter validated enclosure;
empty intersection is treated as an implementation error, never widened away.
"""
from __future__ import annotations

from typing import Sequence

import ou3_sea3_full_word_riccati_backend as BASE
from ou3_interval import (
    Interval,
    IntervalMatrix,
    matrix_add,
    matrix_identity,
    matrix_mul,
    matrix_sub,
    matrix_transpose,
)
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, matrix_symmetric_hull

JointWordState = BASE.JointWordState
PriorFreeBatchState = BASE.PriorFreeBatchState
USEFUL_GATE = BASE.USEFUL_GATE


def _shape(A: Sequence[Sequence[Interval]]) -> tuple[int, int]:
    return BASE._shape(A)


def _intersection(A: Sequence[Sequence[Interval]], B: Sequence[Sequence[Interval]]) -> IntervalMatrix:
    if _shape(A) != _shape(B):
        raise ValueError("interval matrix intersection shape mismatch")
    r, c = _shape(A)
    out: IntervalMatrix = []
    for i in range(r):
        row: list[Interval] = []
        for j in range(c):
            lo = max(A[i][j].lo, B[i][j].lo)
            hi = min(A[i][j].hi, B[i][j].hi)
            if lo > hi:
                raise RuntimeError(
                    f"redundant exact enclosures disagree at ({i},{j}): "
                    f"[{A[i][j].lo},{A[i][j].hi}] vs [{B[i][j].lo},{B[i][j].hi}]"
                )
            row.append(Interval(lo, hi))
        out.append(row)
    return out


def joseph_measurement(
    state: JointWordState,
    H: Sequence[Sequence[Interval]],
    R: Sequence[Sequence[Interval]],
) -> dict:
    n = state.dimension
    rows, cols = _shape(H)
    if cols != n or rows == 0 or _shape(R) != (rows, rows):
        raise ValueError("measurement H/R dimension mismatch")

    Pm = [[x for x in row] for row in state.P]
    Om = [[x for x in row] for row in state.Omega]
    Ht = matrix_transpose(H)
    PHt = matrix_mul(Pm, Ht)
    S = matrix_symmetric_hull(matrix_add(matrix_mul(H, PHt), R))
    Sinv = matrix_inverse_gauss_jordan(S)
    K = matrix_mul(PHt, Sinv)
    A = matrix_sub(matrix_identity(n), matrix_mul(K, H))
    At = matrix_transpose(A)
    KRKt = matrix_mul(matrix_mul(K, R), matrix_transpose(K))

    Pj = matrix_symmetric_hull(
        matrix_add(matrix_mul(matrix_mul(A, Pm), At), KRKt)
    )
    # Exact covariance-form Kalman update.  H P is PHt^T because P is
    # symmetric; retaining PHt from the same prior also avoids an extra natural
    # interval multiplication before the strongest R_S correction.
    Ps = matrix_symmetric_hull(
        matrix_sub(Pm, matrix_mul(K, matrix_transpose(PHt)))
    )
    state.P = _intersection(Pj, Ps)

    state.Psi = matrix_mul(A, state.Psi)
    Oj = matrix_symmetric_hull(
        matrix_add(matrix_mul(matrix_mul(A, Om), At), KRKt)
    )
    Odecomp = matrix_symmetric_hull(
        matrix_sub(
            state.P,
            matrix_mul(
                matrix_mul(state.Psi, state.P0),
                matrix_transpose(state.Psi),
            ),
        )
    )
    state.Omega = _intersection(Oj, Odecomp)

    state.events += 1
    state.measurements += 1
    if not BASE.decomposition_identity_enclosed(state):
        raise RuntimeError("P/Psi/Omega identity lost after tightened Joseph measurement")
    return {
        "S": S,
        "K": K,
        "A": A,
        "shipping_Joseph_enclosure_intersected": True,
        "exact_Schur_covariance_enclosure_intersected": True,
        "exact_decomposition_Omega_enclosure_intersected": True,
    }


# Everything else is the canonical backend unchanged.
initialize = BASE.initialize
initialize_prior_free = BASE.initialize_prior_free
predict = BASE.predict
prior_free_predict = BASE.prior_free_predict
prior_free_measurement = BASE.prior_free_measurement
add_psd_floor = BASE.add_psd_floor
prior_free_add_psd_floor = BASE.prior_free_add_psd_floor
certify_contraction = BASE.certify_contraction
decomposition_identity_enclosed = BASE.decomposition_identity_enclosed
decomposition_residual = BASE.decomposition_residual
contraction_matrix = BASE.contraction_matrix
prediction_contraction_image = BASE.prediction_contraction_image
joseph_contraction_image = BASE.joseph_contraction_image
floor_contraction_image = BASE.floor_contraction_image
contraction_identity_enclosed = BASE.contraction_identity_enclosed
shipping_source_parity = BASE.shipping_source_parity
contraction_preservation_identities = BASE.contraction_preservation_identities
validate_backend = BASE.validate_backend
_self_test = BASE._self_test
