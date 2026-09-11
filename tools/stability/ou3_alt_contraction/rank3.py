"""Exact thin-factor primitives for ALT's three-component measurements.

These are arithmetic/representation optimizations only. They preserve the full
H18/A21/joint24 state and all cross terms; they do not freeze endogenous gains.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


def _matrix(x, name, shape=None):
    a = np.array(x, dtype=float, copy=True)
    if a.ndim != 2 or not np.isfinite(a).all():
        raise ValueError(f"{name} must be finite matrix")
    if shape is not None and a.shape != shape:
        raise ValueError(f"{name}: expected {shape}, got {a.shape}")
    return a


def _symmetric(x, name, size=None):
    a = _matrix(x, name)
    if a.shape[0] != a.shape[1] or (size is not None and a.shape != (size, size)):
        raise ValueError(f"{name} must be square")
    if not np.allclose(a, a.T, rtol=0, atol=1e-12):
        raise ValueError(f"{name} must be symmetric")
    return (a + a.T) / 2


def state_apply(Psi, K, H):
    """Exact (I-KH)Psi as Psi-K(H Psi), retaining the full state."""
    Psi = _matrix(Psi, "Psi")
    n, _ = Psi.shape
    K = _matrix(K, "K", (n, 3))
    H = _matrix(H, "H", (3, n))
    return Psi - K @ (H @ Psi)


def joseph(P, K, S, PCt):
    """Actual-numerator Joseph polynomial with n x 3 thin factors."""
    P = _symmetric(P, "P")
    n = P.shape[0]
    K = _matrix(K, "K", (n, 3))
    PCt = _matrix(PCt, "PCt", (n, 3))
    S = _symmetric(S, "S", 3)
    out = P - K @ PCt.T - PCt @ K.T + K @ S @ K.T
    return (out + out.T) / 2


def storage_delta(M, K, H):
    """Exact (I-KH)'M(I-KH)-M; each event has rank at most six."""
    M = _symmetric(M, "M")
    n = M.shape[0]
    K = _matrix(K, "K", (n, 3))
    H = _matrix(H, "H", (3, n))
    MK = M @ K
    core = K.T @ MK
    out = -H.T @ MK.T - MK @ H + H.T @ core @ H
    return (out + out.T) / 2


@dataclass(frozen=True)
class ProductPortLift:
    equality: np.ndarray
    correction_difference: np.ndarray


def increment_product_ports(N1, S1):
    """Thin exact descriptor conditional on hard product graphs.

    chi=(dr,dq,uS,uN), uS=dS*q0, uN=dN*q0. Then
      S1*dq+uS-dr=0,
      dcorrection=N1*dq+uN.
    uS/uN are NOT free disturbances; the source-uniform master must enforce
    their same-history product identities.
    """
    N1 = _matrix(N1, "N1")
    n, m = N1.shape
    S1 = _symmetric(S1, "S1", m)
    np.linalg.cholesky(S1)
    k = 3 * m + n
    Dr = np.hstack((np.eye(m), np.zeros((m, k - m))))
    Dq = np.hstack((np.zeros((m, m)), np.eye(m), np.zeros((m, m + n))))
    Us = np.hstack((np.zeros((m, 2 * m)), np.eye(m), np.zeros((m, n))))
    Un = np.hstack((np.zeros((n, 3 * m)), np.eye(n)))
    return ProductPortLift(-Dr + S1 @ Dq + Us, N1 @ Dq + Un)
