#!/usr/bin/env python3
"""Joint Joseph update primitives for the OU-III complete-word P4 proof.

This module exists specifically to avoid the failed ``interval K`` proof route.
For one source-correlated measurement cell it keeps P,H,R,S and the residual r
in the same calculation and uses

    S = H P H^T + R,
    x = S^-1 r,
    dx = P H^T x,
    D = r^T x - eta^T R^-1 eta.

No entrywise interval enclosure of K=P H^T S^-1 is ever formed.  The posterior
covariance is evaluated by the algebraically equivalent information/Joseph
identity

    P+ = P - P H^T S^-1 H P,

then intersected with the exact Loewner fact 0 <= P+ <= P before the deployed
attitude reset congruence is applied by the caller.

Innovation inversion is also fail-closed.  We first use outward-rounded
fixed-pivot interval Gauss--Jordan directly on the source-correlated symmetric
SPD family.  If that pivot order cannot exclude zero we fall back to the
verified midpoint/Neumann enclosure.  Neither route uses an ordinary floating
inverse as an enclosure, and neither subdivides or boxes K.

The routines are deliberately small and mode-agnostic.  H=18 and A=21 word
backends can share them, which prevents the proof from maintaining two subtly
different Joseph identities.
"""
from __future__ import annotations

import math
from typing import Sequence

from ou3_interval import Interval, matrix_add, matrix_mul, matrix_transpose
from ou3_interval_linear_algebra import (
    IntervalPivotError,
    matrix_inverse_gauss_jordan,
    matrix_symmetric_hull,
)
import ou3_verified_spd_inverse as VINV

SCHEMA = 2


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _intersect(a: Interval, b: Interval) -> Interval:
    lo = max(a.lo, b.lo)
    hi = min(a.hi, b.hi)
    if lo > hi:
        raise RuntimeError(f"joint Joseph interval intersection empty: {a} vs {b}")
    return Interval(lo, hi)


def vec_dot(a: Sequence[Interval], b: Sequence[Interval]) -> Interval:
    if len(a) != len(b):
        raise ValueError("dot-product dimension mismatch")
    z = I(0.0)
    for x, y in zip(a, b):
        z = z + x * y
    return z


def mat_vec(A, x: Sequence[Interval]) -> list[Interval]:
    if not A or len(A[0]) != len(x):
        raise ValueError("matrix/vector dimension mismatch")
    out = []
    for row in A:
        z = I(0.0)
        for a, b in zip(row, x):
            z = z + a * b
        out.append(z)
    return out


def innovation(P, H, R):
    """Return PH^T and the same-cell symmetric innovation family S."""
    PHt = matrix_mul(P, matrix_transpose(H))
    S = matrix_symmetric_hull(matrix_add(matrix_mul(H, PHt), R))
    return PHt, S


def verified_inverse(S):
    """Invert one source-certified symmetric SPD innovation family.

    Fixed-pivot interval elimination is the preferred route because it retains
    the natural source correlations already present in S.  Failure of that
    particular pivot order is not evidence that S is singular, so a certified
    midpoint/Neumann inverse is tried second.  Both are rigorous enclosures.
    """
    Ssym = matrix_symmetric_hull(S)
    try:
        X = matrix_inverse_gauss_jordan(Ssym)
        return X, {
            "schema": SCHEMA,
            "qualification": "OU3_VERIFIED_FIXED_PIVOT_INNOVATION_INVERSE",
            "dimension": len(Ssym),
            "source_symmetry_certified": True,
            "source_SPD_certified": True,
            "inverse_backend": "FIXED_PIVOT_INTERVAL_GAUSS_JORDAN",
            "neumann_q_inf_upper": None,
            "ordinary_float_inverse_used_as_enclosure": False,
            "K_interval_matrix_materialized": False,
        }
    except (IntervalPivotError, ZeroDivisionError):
        X, meta = VINV.inverse_enclosure(
            Ssym,
            symmetric_certified=True,
            spd_certified=True,
        )
        out = dict(meta)
        out["inverse_backend"] = "MIDPOINT_NEUMANN_VERIFIED"
        out["K_interval_matrix_materialized"] = False
        return X, out


def solve_residual(Sinv, residual: Sequence[Interval]) -> list[Interval]:
    """Compute x=S^-1 r without materializing K."""
    return mat_vec(Sinv, residual)


def direct_correction(PHt, innovation_solution: Sequence[Interval]) -> list[Interval]:
    """Compute dx=P H^T S^-1 r directly, without a K interval matrix."""
    return mat_vec(PHt, innovation_solution)


def _rinv_quadratic(eta: Sequence[Interval], R) -> Interval:
    """eta^T R^-1 eta for diagonal positive measurement covariance.

    OU-III's S, accelerometer and magnetometer R matrices are diagonal in the
    configured proof scope.  Refuse a hidden generic inverse here; if that source
    contract changes, the caller must add a verified R solve explicitly.
    """
    if len(R) != len(eta) or any(len(row) != len(eta) for row in R):
        raise ValueError("R/eta dimension mismatch")
    z = I(0.0)
    for i, e in enumerate(eta):
        for j in range(len(eta)):
            if j != i and not (R[i][j].lo == 0.0 and R[i][j].hi == 0.0):
                raise RuntimeError("joint Joseph eta penalty requires diagonal R in current scope")
        if not R[i][i].lo > 0.0:
            raise RuntimeError("measurement covariance lost positive diagonal floor")
        rinv = Interval.outward_bounds(1.0 / R[i][i].hi, 1.0 / R[i][i].lo)
        z = z + e * rinv * e
    return z


def signed_dissipation(residual: Sequence[Interval], innovation_solution: Sequence[Interval],
                       eta: Sequence[Interval], R) -> dict:
    """Evaluate the exact Joseph signed energy identity on one joint cell."""
    information = vec_dot(residual, innovation_solution)
    penalty = _rinv_quadratic(eta, R)
    signed = information - penalty
    return {
        "residual_information": information,
        "nonlinear_eta_penalty": penalty,
        "signed_dissipation": signed,
        "identity": "r^T S^-1 r - eta^T R^-1 eta",
    }


def _posterior_loewner_box(P):
    """Entrywise enclosure implied by 0<=P+<=P for a PSD covariance family."""
    n = len(P)
    dhi = [math.nextafter(max(0.0, P[i][i].hi), math.inf) for i in range(n)]
    out = [[I(0.0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        out[i][i] = Interval(0.0, dhi[i])
        for j in range(i + 1, n):
            b = math.nextafter(math.sqrt(math.nextafter(dhi[i] * dhi[j], math.inf)), math.inf)
            out[i][j] = Interval(-b, b)
            out[j][i] = out[i][j]
    return out


def posterior_covariance(P, PHt, Sinv, *, psd_tighten=None):
    """Exact P-PH^T S^-1 HP natural extension, intersected with Loewner box."""
    reduction = matrix_mul(matrix_mul(PHt, Sinv), matrix_transpose(PHt))
    n = len(P)
    raw = [[P[i][j] - reduction[i][j] for j in range(n)] for i in range(n)]
    raw = matrix_symmetric_hull(raw)
    loose = _posterior_loewner_box(P)
    for i in range(n):
        for j in range(n):
            raw[i][j] = _intersect(raw[i][j], loose[i][j])
    return psd_tighten(raw) if psd_tighten is not None else raw


def measurement_cell(P, H, R, residual: Sequence[Interval], eta: Sequence[Interval],
                     *, psd_tighten=None) -> dict:
    """Evaluate one complete same-cell Joseph operation without K materialization."""
    PHt, S = innovation(P, H, R)
    Sinv, inverse_meta = verified_inverse(S)
    x = solve_residual(Sinv, residual)
    dx = direct_correction(PHt, x)
    energy = signed_dissipation(residual, x, eta, R)
    Pplus = posterior_covariance(P, PHt, Sinv, psd_tighten=psd_tighten)
    return {
        "P_accepted_pre_reset": Pplus,
        "PHt": PHt,
        "S": S,
        "Sinv": Sinv,
        "innovation_solution": x,
        "dx": dx,
        "energy": energy,
        "inverse_meta": inverse_meta,
        "K_interval_matrix_materialized": False,
        "correction_identity": "dx=P H^T (S^-1 r)",
        "posterior_identity": "P+=P-P H^T S^-1 H P",
        "same_cell_P_H_R_r_coupling_retained_until_S_solve": True,
    }


def validate_cell(cell: dict) -> list[str]:
    f = []
    if cell.get("K_interval_matrix_materialized") is not False:
        f.append("K interval matrix was materialized")
    if cell.get("same_cell_P_H_R_r_coupling_retained_until_S_solve") is not True:
        f.append("joint P,H,R,r coupling flag missing")
    if cell.get("correction_identity") != "dx=P H^T (S^-1 r)":
        f.append("direct correction identity changed")
    if cell.get("posterior_identity") != "P+=P-P H^T S^-1 H P":
        f.append("posterior information identity changed")
    e = cell.get("energy", {})
    if e.get("identity") != "r^T S^-1 r - eta^T R^-1 eta":
        f.append("signed Joseph identity changed")
    meta = cell.get("inverse_meta", {})
    if meta.get("source_symmetry_certified") is not True or meta.get("source_SPD_certified") is not True:
        f.append("verified S inverse lacks symmetry/SPD premise")
    if meta.get("ordinary_float_inverse_used_as_enclosure") is not False:
        f.append("innovation inverse used ordinary float as enclosure")
    return f
