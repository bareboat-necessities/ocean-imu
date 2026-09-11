"""Exact finite measurement descriptor for the physical joint24 ALT state.

This is a real-arithmetic graph, NOT a Jacobian or a shipping stability gate.
The proof is in docs/ou3-alt-finite-measurement-proof.md.  Rational operations
below allow exact regression; polynomial identities also hold for real inputs.
The scalar quaternion and projection parameters remain hard graph variables.
No sample/parameter table here qualifies their reachable source family.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from typing import Sequence

NSTATE = 24
NEST = 21


def rational(x):
    if isinstance(x, (float, bool)):
        raise TypeError('exact graph accepts integers, Fraction or rational strings, not float/bool')
    return F(x)


def vec(x: Sequence, n: int):
    if len(x) != n:
        raise ValueError(f'expected vector of length {n}')
    return [rational(v) for v in x]


def mat(x: Sequence[Sequence], rows: int, cols: int):
    if len(x) != rows or any(len(r) != cols for r in x):
        raise ValueError(f'expected {rows} by {cols} matrix')
    return [vec(r, cols) for r in x]


def zeros(rows, cols):
    return [[F(0) for _ in range(cols)] for _ in range(rows)]


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return list(map(list, zip(*a)))


def mm(a, b):
    if not a or not b or len(a[0]) != len(b):
        raise ValueError('matrix product shape mismatch')
    return [[sum((x*y for x, y in zip(row, col)), 0)
             for col in zip(*b)] for row in a]


def mv(a, x):
    return [r[0] for r in mm(a, [[v] for v in x])]


def plus(a, b, scale=1):
    if len(a) != len(b) or any(len(ar) != len(br) for ar, br in zip(a, b)):
        raise ValueError('matrix addition shape mismatch')
    return [[x+scale*y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scaled(a, s):
    return [[s*x for x in row] for row in a]


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), 0)


def skew(c):
    x, y, z = c
    return [[0, -z, y], [z, 0, -x], [-y, x, 0]]


def cayley_polynomials(c):
    """Polynomial U,D,E_num; accepts formal indeterminates as well as rationals.

    T=(I-[c]x/2)^-1=U/D, E=E_num/D, E-I=T[c]x.
    c=2*tan(theta/2)*axis, not a rotation vector.
    """
    C = skew(c)
    D = 4+dot(c, c)
    U = plus(plus(scaled(eye(3), 4), scaled(C, 2)),
             [[x*y for y in c] for x in c])
    E_num = plus(plus(scaled(eye(3), D), scaled(C, 4)), scaled(mm(C, C), 2))
    return U, D, E_num


def right_reset_polynomials(c, d, w, k):
    """For E+=E(c)Q(d)^T: W(c+-c)+k U d=0, no differentiation.

    Q is the normalization of (w,k*d). Normalization cancels in this identity.
    Both the polynomial and axis-angle shipping branches satisfy it.
    """
    U, _, _ = cayley_polynomials(c)
    v = [k*x for x in d]
    W = 2*w+dot(c, v)
    cv = mv(skew(c), v)
    twice_V = [2*(w*c[i]-2*v[i]-cv[i]) for i in range(3)]
    return W, twice_V, scaled(U, k)


def small_quaternion_parameters(d):
    """Exact real branch only; do not replace the axis branch by this helper."""
    d = vec(d, 3)
    u = dot(d, d)
    if u >= F(1, 10000):
        raise ValueError('not in the strict theta<1/100 polynomial branch')
    return 1-u/F(8)+u*u/F(384), F(1, 2)-u/F(48)+u*u/F(3840)


def residual_factors(kind, c, *, f_hat=None, R_hat=None, m_hat=None):
    """Return actual linearization H and exact finite secant Hbar.

    r=Hbar*z+nu for sensors. For S, nu=-S_phys, with its single Live
    primitive origin; nu is NOT a new independently admissible source input.
    Hbar need not be the Jacobian. Its coefficients retain their c/source graph.
    Both H18 and A21 accel residuals include held/active e_ba.
    """
    c = vec(c, 3)
    U, D, E_num = cayley_polynomials(c)
    H = zeros(3, NEST)
    Hbar = zeros(3, NSTATE)
    if kind in ('accelerometer', 'magnetometer'):
        f = vec(f_hat if kind == 'accelerometer' else m_hat, 3)
        Htheta = scaled(skew(f), -1)
        Hfinite = scaled(mm(U, Htheta), 1/D)
        for i in range(3):
            H[i][:3], Hbar[i][:3] = Htheta[i], Hfinite[i]
        if kind == 'accelerometer':
            R = mat(R_hat, 3, 3)
            ER = scaled(mm(E_num, R), 1/D)
            for i in range(3):
                H[i][15:18], Hbar[i][15:18] = R[i], ER[i]
                H[i][18+i] = Hbar[i][18+i] = F(1)
    elif kind == 'S_zero':
        for i in range(3):
            H[i][12+i] = Hbar[i][12+i] = F(1)
    else:
        raise ValueError('unsupported measurement kind')
    return H, Hbar


def measurement_operands(P21, R, H, *, active_bias):
    """Actual zero-lever symmetric-real branch, including held BA uncertainty.

    All 21 covariance coordinates are required also in H18. S uses the full
    residual model; N uses the held column AND row masks. This function does
    not prove PSD/acceptance or reachable covariance membership. The parent
    source relation must do so, and attach any safe-LDLT repair separately.
    """
    if type(active_bias) is not bool:
        raise TypeError('active_bias must be a literal branch boolean')
    P, R, H = mat(P21, 21, 21), mat(R, 3, 3), mat(H, 3, 21)
    if P != transpose(P) or R != transpose(R):
        raise ValueError('real symmetric branch required; numerical asymmetry is a separate defect')
    Hg = [row[:] for row in H]
    if not active_bias:
        for row in Hg:
            row[18:21] = [F(0)]*3
    N = mm(P, transpose(Hg))
    if not active_bias:
        N[18:21] = zeros(3, 3)
    S = plus(mm(mm(H, P), transpose(H)), R)
    return N, S


def projection_matrix(alpha):
    """Exact e_ba+=alpha*e_ba_pre+(1-alpha)*beta; keep its hard radial graph."""
    alpha = rational(alpha)
    if not 0 < alpha <= 1:
        raise ValueError('positive radial projection factor in (0,1] required')
    P = eye(24)
    for i in range(3):
        P[18+i][18+i] = alpha
        P[18+i][21+i] = 1-alpha
    return P


def projection_graph_holds(pre_bias_error, beta, alpha, radius):
    """Exact branch predicates; no arbitrary alpha can satisfy this admission."""
    e, b = vec(pre_bias_error, 3), vec(beta, 3)
    a, R = rational(alpha), rational(radius)
    if R <= 0 or not 0 < a <= 1:
        return False
    n2 = sum((b[i]-e[i])**2 for i in range(3))
    return (a == 1 and n2 <= R*R) or (a < 1 and n2 > R*R and a*a*n2 == R*R)


@dataclass(frozen=True)
class FiniteDescriptor:
    # chi=(z24,q3,nu3,h); h=1. E*chi=0 retains c=z[:3], d=N_theta*q.
    equality: list
    before: list
    after: list
    B: list
    denominator: F


def descriptor(c, d, w, k, N, S, Hbar, *, alpha):
    """Exact finite descriptor, conditional on d=N_theta*q and branch graphs.

    Returns no closed/pass flag. Supplying a numerical coefficient matrix is
    not proof of these hard graphs or of a complete physical word.
    """
    c, d = vec(c, 3), vec(d, 3)
    w, k = rational(w), rational(k)
    N, S, Hbar = mat(N, 21, 3), mat(S, 3, 3), mat(Hbar, 3, 24)
    W, _, Lnum = right_reset_polynomials(c, d, w, k)
    if W == 0:
        raise ValueError('finite Cayley chart denominator is zero')
    B = zeros(24, 3)
    B[:3] = scaled(mm(Lnum, N[:3]), 1/W)
    B[3:21] = [row[:] for row in N[3:21]]
    P = projection_matrix(alpha)
    before = [eye(24)[i]+[F(0)]*7 for i in range(24)]
    pre = [eye(24)[i]+[-x for x in B[i]]+[F(0)]*4 for i in range(24)]
    equality = [[-x for x in Hbar[i]]+S[i]+[-F(i == j) for j in range(3)]+[F(0)] for i in range(3)]
    for i in range(3):
        equality.append([F(0)]*24+N[i]+[F(0)]*3+[-d[i]])
        equality.append([F(i == j) for j in range(24)]+[F(0)]*6+[-c[i]])
    return FiniteDescriptor(equality, before, mm(P, pre), B, W)


def projected_metric(M, alpha):
    """P_alpha' M P_alpha via the three bias columns/rows, no dense congruence."""
    M = mat(M, 24, 24)
    a = rational(alpha)
    projection_matrix(a)  # Validate the branch factor, not its radial graph.
    MP = [row[:] for row in M]
    for i in range(24):
        for j in range(3):
            MP[i][18+j] = a*M[i][18+j]
            MP[i][21+j] = M[i][21+j]+(1-a)*M[i][18+j]
    out = [row[:] for row in MP]
    for i in range(3):
        for j in range(24):
            out[18+i][j] = a*MP[18+i][j]
            out[21+i][j] = MP[21+i][j]+(1-a)*MP[18+i][j]
    return out


def finite_storage_change(M, z, q, B, *, alpha):
    """Exact joint-storage ledger using a 24x3 factor and a 3x3 core.

    Measurement and radial projection, not covariance transport, are accounted
    for here. M may contain ALL motion/bias/truth cross terms.
    """
    M, B, z, q = mat(M, 24, 24), mat(B, 24, 3), vec(z, 24), vec(q, 3)
    if M != transpose(M):
        raise ValueError('symmetric storage required')
    Mp = projected_metric(M, alpha)
    MB = mm(Mp, B)
    return dot(z, mv(plus(Mp, M, -1), z))-2*dot(z, mv(MB, q))+dot(q, mv(mm(transpose(B), MB), q))
