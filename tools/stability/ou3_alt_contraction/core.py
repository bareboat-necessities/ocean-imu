"""Exploratory inverse-free word algebra, NOT a shipping stability certificate.

All arrays are host floating-point candidate objects. The exact algebra is
specified in docs/ou3-alt-contraction.md. No routine here sets a proof PASS.
In particular, H is the residual sensitivity, N is the ACTUAL gain numerator,
and S is the ACTUAL innovation; a held-row mask must not be reconstructed away.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import math

import numpy as np


def matrix(value, name: str, shape=None) -> np.ndarray:
    out = np.array(value, dtype=float, copy=True)
    if out.ndim != 2 or not np.isfinite(out).all():
        raise ValueError(f"{name} must be a finite matrix")
    if shape is not None and out.shape != shape:
        raise ValueError(f"{name}: expected {shape}, got {out.shape}")
    return out


def symmetric(value, name: str, size=None) -> np.ndarray:
    out = matrix(value, name)
    if out.shape[0] != out.shape[1] or (size is not None and out.shape != (size, size)):
        raise ValueError(f"{name} must be square with the requested dimension")
    if not np.allclose(out, out.T, rtol=0, atol=1e-12):
        raise ValueError(f"{name} must be symmetric")
    return (out + out.T) / 2


@dataclass(frozen=True)
class MeasurementLift:
    """Variables chi=(e, w, q), with S q = H e + w and e+ = e - N q."""
    equality: np.ndarray
    before: np.ndarray
    after: np.ndarray
    disturbance: np.ndarray


def measurement_lift(H, N, S) -> MeasurementLift:
    """Assemble the actual solve graph without constructing K or S**-1.

    This is a fixed-coefficient algebra primitive. For an endogenous gain the
    differential graph ALSO needs dS*q and dN*q; it is not the full Jacobian.
    SPD is checked numerically only, never claimed as an outward certificate.
    """
    H = matrix(H, "H")
    m, n = H.shape
    if min(m, n) == 0:
        raise ValueError("measurement dimensions must be nonzero")
    N = matrix(N, "N", (n, m))
    S = symmetric(S, "S", m)
    np.linalg.cholesky(S)
    return MeasurementLift(
        np.hstack((-H, -np.eye(m), S)),
        np.hstack((np.eye(n), np.zeros((n, 2*m)))),
        np.hstack((np.eye(n), np.zeros((n, m)), -N)),
        np.hstack((np.zeros((m, n)), np.eye(m), np.zeros((m, m)))),
    )


def dissipation_form(before, after, disturbance, M0, M1, rho, Gamma,
                     bounded=None, Beta=None) -> np.ndarray:
    """V+ - rho*V - w'Gamma*w - c'Beta*c; joint storage is not marginalized."""
    if not math.isfinite(rho) or not 0 < rho < 1:
        raise ValueError("rho must lie strictly between zero and one")
    before = matrix(before, "before")
    n, k = before.shape
    after = matrix(after, "after", (n, k))
    disturbance = matrix(disturbance, "disturbance")
    if disturbance.shape[1] != k:
        raise ValueError("disturbance must use the same lifted variables")
    M0, M1 = symmetric(M0, "M0", n), symmetric(M1, "M1", n)
    np.linalg.cholesky(M0)
    np.linalg.cholesky(M1)
    Gamma = symmetric(Gamma, "Gamma", disturbance.shape[0])
    if np.linalg.eigvalsh(Gamma).min(initial=0) < -1e-12:
        raise ValueError("Gamma must be positive semidefinite")
    Q = after.T @ M1 @ after - rho * before.T @ M0 @ before
    Q -= disturbance.T @ Gamma @ disturbance
    if (bounded is None) != (Beta is None):
        raise ValueError("bounded and Beta must be supplied together")
    if bounded is not None:
        bounded = matrix(bounded, "bounded")
        if bounded.shape[1] != k:
            raise ValueError("bounded ports must use the same lifted variables")
        Beta = symmetric(Beta, "Beta", bounded.shape[0])
        if np.linalg.eigvalsh(Beta).min(initial=0) < -1e-12:
            raise ValueError("Beta must be positive semidefinite")
        Q -= bounded.T @ Beta @ bounded
    return (Q + Q.T) / 2


def sector_iqc(dimension: int, lower: float = 0, upper: float = 1) -> np.ndarray:
    """q=(dy-lower*dx)'(upper*dx-dy)>=0, variables (dx,dy).

    lower=0, upper=1 is the firmly nonexpansive Euclidean projection IQC.
    Only nonnegative SCALAR multipliers are admitted by the generic assembler;
    arbitrary matrix multipliers need a separate valid nonlinearity theorem.
    """
    if not isinstance(dimension, int) or dimension <= 0:
        raise ValueError("dimension must be a positive integer")
    if not (math.isfinite(lower) and math.isfinite(upper) and lower <= upper):
        raise ValueError("invalid sector endpoints")
    I = np.eye(dimension)
    return np.block([[-lower*upper*I, (lower+upper)*I/2],
                     [(lower+upper)*I/2, -I]])


def iqc_master(Q, equality, multiplier, iqcs=()) -> np.ndarray:
    """Q + He(Y E) + sum(lambda_i Q_i), with chi'Q_i chi >= 0.

    Nonpositivity implies the desired dissipation on E chi=0 and the IQC
    graphs. The PLUS sign is essential. Parameter-dependent coefficients make
    a family of LMIs, not automatically one convex LMI or a vertex theorem.
    """
    Q = symmetric(Q, "Q")
    k = Q.shape[0]
    E = matrix(equality, "equality")
    if E.shape[1] != k:
        raise ValueError("equality and Q dimensions disagree")
    Y = matrix(multiplier, "multiplier", (k, E.shape[0]))
    out = Q + Y @ E + E.T @ Y.T
    for weight, term in iqcs:
        if not math.isfinite(weight) or weight < 0:
            raise ValueError("IQC multipliers must be nonnegative")
        out += weight * symmetric(term, "IQC", k)
    return (out + out.T) / 2


def structural_findings() -> dict:
    """Exact rational obstructions to two over-strong architectural shortcuts.

    The 2-state example is a masking/algebra regression, not a claim that it
    constitutes a full admissible physical BRMM trajectory of the vessel.
    """
    A = ((F(2, 3), F(-1, 3)), (F(0), F(1)))
    v = (F(-1), F(1))
    assert tuple(sum(A[i][j]*v[j] for j in range(2)) for i in range(2)) == v
    # P=I, H_residual=[1,1], actual N=[1,0]', S=3.
    # The shipping-style algebraic Joseph expression gives diag(2/3,1).
    actual_information = ((F(3, 2), F(0)), (F(0), F(1)))
    unmasked_information = ((F(2), F(1)), (F(1), F(2)))
    gap = [[str(actual_information[i][j]-unmasked_information[i][j])
            for j in range(2)] for i in range(2)]
    return {
        "held_mode_neutral_witness": {
            "arithmetic": "exact rational",
            "A": [[str(x) for x in row] for row in A],
            "v": [str(x) for x in v], "A_v_equals_v": True,
            "any_SPD_common_metric_gap": "(1-rho) * v^T M v > 0 for rho < 1",
            "classification": "over-strong strict common-metric target, not filter instability",
        },
        "masked_information_identity": {
            "P": [[1, 0], [0, 1]], "H_residual": [[1, 1]],
            "actual_N": [[1], [0]], "actual_S": [[3]],
            "actual_minus_unmasked_information": gap,
            "unqualified_identity_valid": False,
            "effective_model_bridge_still_possible": True,
        },
    }


def open_obligations() -> dict:
    """Fixed fail-closed status, not user-supplied booleans posing as evidence."""
    return {
        "ALT_LIVE_PASS": False, "ALT_STARTUP_PASS": False,
        "ALT_END_TO_END_PASS": False, "P4_promoted": False, "P5_promoted": False,
        "missing": [
            "source-uniform actual joint24 word graph including true-bias recurrence",
            "BIAS0, BIAS1 and BIAS2 admission with retained same-history coupling",
            "endogenous gain/tuner derivatives or a justified common-storage alternative",
            "guard-complete H18/A21 and H18-to-A21 transport",
            "uniform coercive storage and hard finite-horizon nonlinear IQCs",
            "endpoint dissipation, compatible edges and every-prefix chart retention",
            "actual Mahony/proxy finite-time basin capture at the padded physical envelope",
            "deployment finite-precision enclosure and canonical P3 prerequisite",
        ],
    }


def bias_prediction_lift(F21, phi_hat: float, phi_true: float):
    """Exact joint [true-minus-estimate errors; true bias] prediction graph.

    Keep (phi_true-phi_hat)*beta INSIDE A, not as an independent disturbance.
    One shared driver w enters both e_b and beta. This function is a point
    algebra primitive; interval/source-family qualification is separate.
    """
    F21 = matrix(F21, "F21", (21, 21))
    if not (0 < phi_hat <= 1 and 0 < phi_true <= 1):
        raise ValueError("bias factors must be in (0,1]")
    if not np.allclose(F21[18:21, 18:21], phi_hat*np.eye(3), rtol=0, atol=1e-14):
        raise ValueError("prediction bias block does not match phi_hat")
    if not np.allclose(F21[18:21, :18], 0, rtol=0, atol=1e-14):
        raise ValueError("unexpected prediction coupling into bias")
    A = np.zeros((24, 24))
    A[:21, :21] = F21
    A[18:21, 21:24] = (phi_true-phi_hat)*np.eye(3)
    A[21:24, 21:24] = phi_true*np.eye(3)
    B = np.zeros((24, 3))
    B[18:21, :] = B[21:24, :] = np.eye(3)
    return A, B


def joint_bias_projection_iqc() -> np.ndarray:
    """For chi=(de_before, dbeta, de_after), retain e+=beta-Pi(beta-e).

    du=dbeta-de_before; dy=dbeta-de_after, so dy'(du-dy)>=0.
    The same beta occurs on both sides; it is not an independent reset input.
    """
    I, Z = np.eye(3), np.zeros((3, 3))
    U, Y = np.hstack((-I, I, Z)), np.hstack((Z, I, -I))
    return (Y.T @ U + U.T @ Y)/2 - Y.T @ Y


@dataclass(frozen=True)
class IncrementLift:
    """chi=(delta_r, row_vec(delta_N), row_vec(delta_S), delta_q)."""
    equality: np.ndarray
    correction_difference: np.ndarray


def measurement_increment_lift(N_next, S_next, q_base) -> IncrementLift:
    """Exact FINITE increment of the endogenous innovation solve, not a Jacobian.

    S0*q0=r0, S1*q1=r1 imply
        S1*dq + dS*q0 - dr = 0,
        N1*q1-N0*q0 = N1*dq + dN*q0.
    Thus gain and covariance dependence are retained without intervalizing an
    inverse or assuming a zero residual. dN/dS must still be bound to the same
    shipping state/source, and floating-point solve defects need supply ports.
    Arrays are numerical candidate coefficients, not certified enclosures.
    """
    N = matrix(N_next, "N_next")
    n, m = N.shape
    if min(n, m) == 0:
        raise ValueError("innovation dimensions must be nonzero")
    S = symmetric(S_next, "S_next", m)
    np.linalg.cholesky(S)
    q = np.asarray(q_base, dtype=float)
    if q.shape != (m,) or not np.isfinite(q).all():
        raise ValueError("q_base must be a finite innovation vector")
    Nq = np.kron(np.eye(n), q.reshape(1, m))
    Sq = np.kron(np.eye(m), q.reshape(1, m))
    E = np.hstack((-np.eye(m), np.zeros((m, n*m)), Sq, S))
    C = np.hstack((np.zeros((n, m)), Nq, np.zeros((n, m*m)), N))
    return IncrementLift(E, C)


def masked_supply_certificate() -> dict:
    """Exact positive result for the 2-state masking regression, NOT ocean-IMU.

    A neutral held coordinate precludes strict homogeneous contraction, but a
    joint SPD storage with a bounded held-coordinate supply can close exactly.
    Sylvester's 2x2 criterion is checked with rational arithmetic.
    """
    A = [[F(2, 3), F(-1, 3)], [F(0), F(1)]]
    M = [[F(1), F(1, 4)], [F(1, 4), F(1)]]
    rho = F(9, 10)
    Q = [[sum(A[k][i]*M[k][l]*A[l][j] for k in range(2) for l in range(2))
          -rho*M[i][j]-F(i == 1 and j == 1) for j in range(2)] for i in range(2)]
    det_m = M[0][0]*M[1][1]-M[0][1]*M[1][0]
    det_q = Q[0][0]*Q[1][1]-Q[0][1]*Q[1][0]
    assert M[0][0] > 0 and det_m > 0
    assert Q[0][0] < 0 and det_q > 0 and Q[0][1] == Q[1][0]
    return {'scope': 'exact two-state masked-update analogue only',
            'M': [[str(x) for x in row] for row in M], 'rho': str(rho),
            'bounded_held_supply_weight': '1',
            'Q': [[str(x) for x in row] for row in Q],
            'M_determinant': str(det_m), 'Q_determinant': str(det_q),
            'negative_definite_exact': True, 'joint_cross_terms_retained': True,
            'shipping_certificate': False}
