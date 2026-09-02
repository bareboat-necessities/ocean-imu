#!/usr/bin/env python3
"""Augmented error/covariance interval-AD primitives for complete-word OU-III P4.

The shipping covariance is not an exogenous constant after an accepted attitude
correction: the Joseph posterior is followed by the deployed left-error reset
congruence G(dtheta) P G(dtheta)^T, and dtheta depends on the innovation.  A
complete-word Jacobian with respect to the word-entry error must therefore carry
that covariance dependence into every later S^-1 and correction.

This module lifts the covariance matrix itself to the existing outward interval
AD type.  For one fixed P2/source/acceptance branch it propagates

    P^- = F P F^T + Q,
    S   = H P^- H^T + R,
    X   = S^-1,
    dx  = P^- H^T X r,
    P+  = P^- - P^- H^T X H P^-,
    P_r = G(dx_theta) P+ G(dx_theta)^T.

The inverse derivative is enclosed analytically by

    dX = -X (dS) X,

using the same verified interval enclosure X of the actual SPD innovation
inverse.  No interval Kalman-gain matrix K=P H^T S^-1 is materialized.

Source/tuner/discrete branch variables are parameters, not AD coordinates.  A
final theorem must enumerate them through the P2/source-word automaton.  Within
each fixed branch this module includes every continuous state-dependent
covariance/gain/reset path required by the finite-cell mean-value lemma.
"""
from __future__ import annotations

import argparse
import json
from typing import Sequence

from ou3_interval import Interval, matrix_mul
import ou3_interval_ad as AD
import ou3_p4_joint_joseph as JJ

SCHEMA = 1


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _zero_ad(n: int) -> AD.AD:
    return AD.constant(0.0, n)


def constant_matrix(A, nder: int):
    return [[AD.constant(x, nder) for x in row] for row in A]


def values(A):
    return [[x.val for x in row] for row in A]


def transpose(A):
    return [list(row) for row in zip(*A)]


def add(A, B):
    if len(A) != len(B) or (A and len(A[0]) != len(B[0])):
        raise ValueError("AD matrix add shape mismatch")
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def sub(A, B):
    if len(A) != len(B) or (A and len(A[0]) != len(B[0])):
        raise ValueError("AD matrix subtract shape mismatch")
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def matmul(A, B):
    return AD.matmul(A, B)


def matvec(A, x):
    return AD.matvec(A, x)


def symmetric_hull(A):
    if len(A) != len(A[0]):
        raise ValueError("symmetric hull requires square matrix")
    n = len(A)
    out = [[A[i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            z = AD.hull_ad(A[i][j], A[j][i])
            out[i][j] = z
            out[j][i] = z
    return out


def tighten_values(Pad, value_tighten=None):
    """Intersect value boxes with a proven covariance enclosure, keep derivatives.

    A PSD/Loewner tightening callback acts only on interval values.  Since the
    derivative intervals were produced by the exact algebra before the value
    intersection, retaining them unchanged remains an enclosure of the actual
    derivative while benefiting later value-only innovation calculations.
    """
    if value_tighten is None:
        return symmetric_hull(Pad)
    old = symmetric_hull(Pad)
    tight = value_tighten(values(old))
    n = len(old)
    out = [[old[i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            out[i][j] = AD.AD(tight[i][j], old[i][j].der)
    return symmetric_hull(out)


def predict(Pad, F, Q, *, value_tighten=None):
    nder = Pad[0][0].n
    Fad = constant_matrix(F, nder)
    Qad = constant_matrix(Q, nder)
    raw = add(matmul(matmul(Fad, Pad), transpose(Fad)), Qad)
    return tighten_values(raw, value_tighten)


def innovation(Pad, H, R):
    """Return AD PH^T and S=HPH^T+R for a fixed source H/R cell."""
    nder = Pad[0][0].n
    Had = constant_matrix(H, nder)
    Rad = constant_matrix(R, nder)
    PHt = matmul(Pad, transpose(Had))
    S = symmetric_hull(add(matmul(Had, PHt), Rad))
    return PHt, S


def verified_inverse_ad(Sad):
    """Verified S^-1 value and full entry-error derivative enclosure."""
    Sval = values(Sad)
    Xval, meta = JJ.verified_inverse(Sval)
    n = len(Sad)
    nder = Sad[0][0].n
    dX_by_k = []
    for k in range(nder):
        dS = [[Sad[i][j].der[k] for j in range(n)] for i in range(n)]
        # dX = -X dS X.  Xval encloses the same actual inverse used by the
        # source cell; natural interval multiplication safely retains all
        # correlations as an over-approximation without an ordinary float inverse.
        tmp = matrix_mul(matrix_mul(Xval, dS), Xval)
        dX_by_k.append([[-tmp[i][j] for j in range(n)] for i in range(n)])
    Xad = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(AD.AD(
                Xval[i][j],
                tuple(dX_by_k[k][i][j] for k in range(nder)),
            ))
        Xad.append(row)
    out_meta = dict(meta)
    out_meta.update({
        "inverse_derivative_identity": "dS^-1=-S^-1(dS)S^-1",
        "state_dependent_innovation_inverse_derivative_enclosed": True,
        "K_interval_matrix_materialized": False,
    })
    return symmetric_hull(Xad), out_meta


def direct_measurement(Pad, H, R, residual: Sequence[AD.AD], *, value_tighten=None):
    """AD Joseph measurement update without materializing K."""
    PHt, Sad = innovation(Pad, H, R)
    Xad, meta = verified_inverse_ad(Sad)
    sol = matvec(Xad, residual)
    dx = matvec(PHt, sol)
    reduction = matmul(matmul(PHt, Xad), transpose(PHt))
    Pj = tighten_values(sub(Pad, reduction), value_tighten)
    return {
        "PHt": PHt,
        "S": Sad,
        "Sinv": Xad,
        "innovation_solution": sol,
        "dx": dx,
        "P_accepted_pre_reset": Pj,
        "inverse_meta": meta,
        "K_interval_matrix_materialized": False,
        "correction_identity": "dx=P H^T (S^-1 r)",
        "posterior_identity": "P+=P-P H^T S^-1 H P",
        "covariance_path_derivative_included": True,
    }


def reset_matrix(dx_theta: Sequence[AD.AD], dimension: int):
    if len(dx_theta) != 3:
        raise ValueError("attitude reset requires three correction coordinates")
    nder = dx_theta[0].n
    G = [[_zero_ad(nder) for _ in range(dimension)] for _ in range(dimension)]
    for i in range(dimension):
        G[i][i] = AD.constant(1.0, nder)
    h = AD.constant(0.5, nder)
    x, y, z = dx_theta
    G[0][1] = -h * z; G[0][2] = h * y
    G[1][0] = h * z;  G[1][2] = -h * x
    G[2][0] = -h * y; G[2][1] = h * x
    return G


def reset_covariance(Pj, dx_theta: Sequence[AD.AD], *, value_tighten=None):
    G = reset_matrix(dx_theta, len(Pj))
    raw = matmul(matmul(G, Pj), transpose(G))
    return tighten_values(raw, value_tighten)


def measurement_and_reset(Pad, H, R, residual: Sequence[AD.AD], *, value_tighten=None):
    cell = direct_measurement(Pad, H, R, residual, value_tighten=value_tighten)
    Pout = reset_covariance(
        cell["P_accepted_pre_reset"], cell["dx"][:3], value_tighten=value_tighten
    )
    cell["P_accepted_post_reset"] = Pout
    cell["reset_covariance_depends_on_correction_AD"] = True
    return cell


def _self_test() -> dict:
    failures = []

    # Scalar inverse family S(x)=2+x, x in [0,0.1].  Exact derivative is
    # -1/(2+x)^2 and must lie inside the analytic interval inverse derivative.
    x = AD.independent(Interval.outward_bounds(0.0, 0.1), 0, 1)
    S = [[AD.constant(2.0, 1) + x]]
    X, meta = verified_inverse_ad(S)
    d = X[0][0].der[0]
    exact0 = -1.0 / (2.0 ** 2)
    exact1 = -1.0 / (2.1 ** 2)
    for v in (exact0, exact1):
        if not d.lo <= v <= d.hi:
            failures.append(f"inverse derivative misses exact scalar endpoint {v}")
    if meta.get("state_dependent_innovation_inverse_derivative_enclosed") is not True:
        failures.append("inverse derivative metadata missing")
    if meta.get("K_interval_matrix_materialized") is not False:
        failures.append("inverse primitive materialized K")

    # One-state Joseph algebra with P(x)=1+0.1x, H=R=1 and residual=x.
    # This checks that covariance derivatives survive innovation, inverse and
    # posterior propagation rather than being silently frozen.
    p = AD.constant(1.0, 1) + AD.constant(0.1, 1) * x
    cell = direct_measurement(
        [[p]], [[I(1.0)]], [[I(1.0)]], [x], value_tighten=None
    )
    dp = cell["P_accepted_pre_reset"][0][0].der[0]
    if dp.lo == 0.0 and dp.hi == 0.0:
        failures.append("Joseph posterior covariance derivative was frozen to zero")
    if cell.get("covariance_path_derivative_included") is not True:
        failures.append("measurement covariance-path derivative flag missing")
    if cell.get("K_interval_matrix_materialized") is not False:
        failures.append("measurement primitive materialized K")

    return {
        "schema": SCHEMA,
        "pass": not failures,
        "failures": failures,
        "scalar_inverse_value": [X[0][0].val.lo, X[0][0].val.hi],
        "scalar_inverse_derivative": [d.lo, d.hi],
        "joseph_posterior_derivative": [dp.lo, dp.hi],
        "K_interval_matrix_materialized": False,
        "state_dependent_covariance_gain_reset_path_targeted": True,
        "P4_PROMOTED_HERE": False,
        "P5_PROMOTED_HERE": False,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if not a.self_test:
        ap.error("standalone mode currently supports only --self-test")
    d = _self_test()
    print(json.dumps(d, indent=2, sort_keys=True))
    return 0 if d["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
