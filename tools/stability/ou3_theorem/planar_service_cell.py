"""Analytical parity-cell transport identities for the literal planar audit.

No invariant cell or all-time admission is asserted here. Numerical evaluations
of the formulas are diagnostics; exact rational regression evidence is separate.

Two issues matter before propagating scheduler cells:
* the linked S-shift has rank <= 2m (4 in EVEN, 2 in ODD), and is O(step),
  whereas separately bounded removal terms need not be small;
* shipping additive AW-floor synchronization is NOT Loewner monotone. It is
  nonexpansive in a Frobenius norm (also after block scaling scalar on AW).
"""
from __future__ import annotations
from fractions import Fraction as F
import numpy as np
from .matrix_certificates import add, matmul, transpose, ldlt, encoded


def _sym(a):
    return (a + a.T) / 2


def linked_shift(P, transition, process, H, R, weight=None):
    """Evaluate the exact linked rank-2m difference, without subtracting two Ps.

    U=F P H', Z=[F P (F'-I)+Q]H', A=H P H'+R,
    E=H[(F-I)P F'+P(F'-I)+Q]H', B=A+E.
    D=U A^-1 E B^-1 U' - U B^-1 Z' - Z B^-1 U' - Z B^-1 Z'.
    All P appearances denote the SAME P; all F/Q/R come from one operation.
    """
    P, transition, process, H, R = map(np.asarray, (P, transition, process, H, R))
    n = len(P)
    dF = transition - np.eye(n)
    U = transition @ P @ H.T
    Z = (transition @ P @ dF.T + process) @ H.T
    A = _sym(H @ P @ H.T + R)
    E = _sym(H @ (dF @ P @ transition.T + P @ dF.T + process) @ H.T)
    B = _sym(A + E)
    ai = np.linalg.inv(A)
    bi = np.linalg.inv(B)
    middle = np.block([[ai @ E @ bi, -bi], [-bi, -bi]])
    L = np.concatenate((U, Z), axis=1)
    if weight is not None:
        L = np.asarray(weight) @ L
    out = _sym(L @ middle @ L.T)
    # Thin QR reduces the norm computation to <=2m dimensions exactly.
    _, qr = np.linalg.qr(L, mode="reduced")
    core = _sym(qr @ middle @ qr.T)
    norm = float(np.max(np.abs(np.linalg.eigvalsh(core))))
    return out, {"rank_upper": 2 * H.shape[0], "norm_diagnostic": norm,
                 "linked_measurement_change_norm": float(np.linalg.norm(E, 2)),
                 "linked_cross_change_norm": float(np.linalg.norm(Z, 2))}


def linked_relative_cell_bound(center, root, radius, transition, process, H, R,
                               weight=None):
    """Analytically valid norm formula on P=center+root E root', ||E||_2<=r.

    The returned floating evaluation is NOT directed rounding. The certificate
    records this distinction. This is a declared covariance-cell relaxation,
    not a claim that every covariance in the cell is shipping-reachable.
    """
    if not np.isfinite(radius) or radius < 0:
        raise ValueError("nonnegative finite cell radius required")
    center, root, transition, process, H, R = map(np.asarray,
        (center, root, transition, process, H, R))
    n = len(center)
    W = np.eye(n) if weight is None else np.asarray(weight)
    dF = transition - np.eye(n)
    lower = _sym(center - radius * root @ root.T)
    if np.linalg.eigvalsh(lower).min() <= 0:
        raise ValueError("cell lower covariance must be SPD")
    A_lower = _sym(H @ lower @ H.T + R)
    B_lower = _sym(H @ (transition @ lower @ transition.T + process) @ H.T + R)
    a, b = np.linalg.eigvalsh(A_lower).min(), np.linalg.eigvalsh(B_lower).min()
    if min(a, b) <= 0:
        raise ValueError("positive innovation lower faces required")
    norm = lambda x: float(np.linalg.norm(x, 2))
    U = transition @ center @ H.T
    Z = (transition @ center @ dF.T + process) @ H.T
    E = H @ (dF @ center @ transition.T + center @ dF.T + process) @ H.T
    u = norm(W @ U) + radius * norm(W @ transition @ root) * norm(H @ root)
    z = norm(W @ Z) + radius * norm(W @ transition @ root) * norm(H @ dF @ root)
    e = norm(E) + radius * norm(H @ dF @ root) * (norm(H @ transition @ root) + norm(H @ root))
    return float(u*u*e/(a*b) + (2*u*z + z*z)/b)


def positive_part(a):
    values, vectors = np.linalg.eigh(_sym(np.asarray(a)))
    return (vectors * np.maximum(values, 0)) @ vectors.T


def aw_floor(P, target, aw_indices):
    """Literal real-arithmetic additive PSD-floor map, not marginal congruence."""
    out = np.array(P, dtype=float, copy=True)
    ix = np.ix_(aw_indices, aw_indices)
    out[ix] += positive_part(np.asarray(target) - out[ix])
    return _sym(out)


def _inverse2(a):
    det = a[0][0]*a[1][1] - a[0][1]*a[1][0]
    if det == 0:
        raise ValueError("singular exact matrix")
    return [[a[1][1]/det, -a[0][1]/det], [-a[1][0]/det, a[0][0]/det]]


def certificate():
    """Exact rational regressions, with the analytical theorem stated in prose.

    Sampled identities are regression tests, not a substitute for the general
    algebraic proof in the appendix. The counterexamples concern proof bounds,
    NOT an admitted physical OU-III execution.
    """
    # Both P and Phi-P are positive definite. P <= Phi does not order PH'H P.
    P = [[F(5000), F(5,2)], [F(5,2), F(1,200)]]
    upper = [[F(10000), F(0)], [F(0), F(1,100)]]
    _, pivots_p = ldlt(P)
    _, pivots_gap = ldlt(add(upper, P, F(-1)))
    # F=I, Q=diag(0,.01), H=[0,1], R=.01.
    # D11 = (5/2)^2*(1/.015 - 1/.025) = 500/3.
    d11 = F(25,4) * (1/F(3,200) - 1/F(1,40))
    invalid_old_bound = F(1,100)**2/F(1,100) + F(1,50)**2/F(1,100)
    assert d11 > invalid_old_bound
    # Full-state AW floor is not Loewner monotone even for one AW coordinate.
    # P0=I, P1=I+ones(2,2); target=3 on coordinate 1. Both AW faces are active.
    floor_gap = [[F(1), F(1)], [F(1), F(0)]]
    floor_gap_det = floor_gap[0][0]*floor_gap[1][1] - floor_gap[0][1]**2
    assert floor_gap_det < 0
    # Noncommuting innovation matrices: exact resolvent linked difference.
    A = [[F(2),F(1,3)],[F(1,3),F(3)]]
    B = [[F(3),F(1,5)],[F(1,5),F(4)]]
    U = [[F(1),F(2)],[F(3),F(4)],[F(2),F(-1)]]
    Z = [[F(1,20),F(-1,30)],[F(1,10),F(1,40)],[F(-1,20),F(1,25)]]
    ai, bi = _inverse2(A), _inverse2(B)
    mm = lambda x,y,z: matmul(matmul(x,y),z)
    V = add(U,Z)
    direct = add(mm(U,ai,transpose(U)),mm(V,bi,transpose(V)),F(-1))
    linked = mm(U,matmul(matmul(ai,add(B,A,F(-1))),bi),transpose(U))
    for x,y in ((U,Z),(Z,U),(Z,Z)):
        linked = add(linked,mm(x,bi,transpose(y)),F(-1))
    assert linked == direct
    return {
        "qualification": "OU3_PLANAR_SERVICE_CELL_ALGEBRA_V1",
        "result_type": "PROVED analytical theorem; exact rational algebra regressions",
        "linked_shift_identity": "D=U A^-1 (B-A) B^-1 U' - U B^-1 Z' - Z B^-1 U' - Z B^-1 Z'",
        "linked_shift_rank_upper": {"even":4,"odd":2},
        "linked_identity_exact_rational_residual": "0",
        "noncommuting_innovation_identity_matrix": encoded(direct),
        "old_bound_counterexample": {"P":encoded(P),"Phi":encoded(upper),
            "P_ldlt_pivots":[str(x) for x in pivots_p],
            "Phi_minus_P_ldlt_pivots":[str(x) for x in pivots_gap],
            "defect_D11":str(d11),"invalid_bound":str(invalid_old_bound),
            "classification":"E_IMPLEMENTATION_DIAGNOSTIC_FAILURE",
            "shipping_counterexample":False},
        "aw_floor": {
            "loewner_monotone":False,
            "ordered_input_difference":[["1","1"],["1","1"]],
            "output_difference":encoded(floor_gap),
            "output_difference_determinant":str(floor_gap_det),
            "frobenius_nonexpansive_for_fixed_target":True,
            "scaled_frobenius_condition":"block-diagonal coordinate scaling, scalar on each invariant AW block",
            "proof":"A+(T-A)_+=T+(A-T)_+; PSD-cone projection is Frobenius nonexpansive; other matrix blocks are unchanged",
            "varying_target_bound":"||F_T(P)-F_T0(P0)||_F <= ||P-P0||_F+||T-T0||_F",
        },
        "numerical_norms_are_directed_rounding":False,
        "scheduler_phase_cell_forward_invariant":False,
        "all_time_magnetic_service_verified":False,
        "full_shipping_counterexample_admitted":False,
        "theorem_closed":False,
    }


if __name__ == "__main__":
    import json
    print(json.dumps(certificate(),indent=2,sort_keys=True))
