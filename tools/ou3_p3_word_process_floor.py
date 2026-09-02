#!/usr/bin/env python3
"""Whole-word process-noise floor for the retained OU-III P3 route.

The P3 covariance ceiling is a multi-second word quantity.  Comparing it with
one 5-ms process covariance is rigorous but useless, especially in the
integrated displacement coordinate.  This module supplies only the missing
word-scale process algebra; it does not scan source paths or promote P3.

For the integrated-OU chain [v,p,S,a_w], N consecutive predictions with fixed
source parameters inject exactly Q(Nh), i.e. the same analytic family read at
X=Nh/tau.  In word scaling diag(sigma*T,sigma*T^2,sigma*T^3,sigma), the family
is dimensionless.  We evaluate its exact correlated exponential series through
X=2.5, apply the exact rational conditioning congruence C=R L^-1, subdivide X
until interval LDLT/inversion certifies each family, and pull the information
bound back to the original word coordinates.

The exact-series coefficients are immutable algebra, so they are memoized.
The expensive conditioned subcell inverse is memoized by its binary64 interval
endpoints.  These caches change no enclosure; they only make a complete source
scan practical.

For each scaled (attitude, gyro-bias) axis pair, the word process floor uses the
exact doubling recursion Omega_2N = Phi_N Omega_N Phi_N' + Omega_N followed by
the conservative endpoint measurement-information update.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import math

from ou3_interval import (
    Interval,
    hull,
    matrix_add,
    matrix_identity,
    matrix_mul,
    matrix_transpose,
    symmetric_gershgorin_upper,
    symmetric_positive_definite_ldlt,
)
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, matrix_symmetric_hull
import ou3_p3_scaled_process as SCALED
import ou3_validated_transcendentals as VT

F = Fraction
WORD_EXACT_SERIES_MAX_X = 2.5
DEFAULT_INFORMATION_SPLIT_DEPTH = 14


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def IF(q: Fraction) -> Interval:
    return I(float(q))


# Exact x->0 LDL conditioning transform of the integrated-OU process shape.
_L_INV = (
    (F(1), F(0), F(0), F(0)),
    (-F(3, 8), F(1), F(0), F(0)),
    (F(1, 15), -F(4, 9), F(1), F(0)),
    (-F(15, 2), F(30), -F(105, 2), F(1)),
)
_R = (F(1), F(10), F(100), F(2))
_CQ = tuple(tuple(_R[i] * _L_INV[i][j] for j in range(4)) for i in range(4))
C = [[IF(v) for v in row] for row in _CQ]
CT = matrix_transpose(C)


@lru_cache(maxsize=None)
def _series_terms(name: str):
    """Exact coefficients of B_name(X) after symbolic leading-power cancel."""
    if name not in SCALED._EXACT:
        raise KeyError(name)
    p = SCALED._DEN_POWER[name]
    N = SCALED.NEAR_EXACT_SERIES_ORDER
    for n in range(p + 1):
        if SCALED._series_coefficient(name, n) != 0:
            raise RuntimeError(
                f"exact OU series {name} lost leading-power cancellation at n={n}"
            )
    # B removes the D_h power p plus the common positive X factor.
    return tuple(
        SCALED._series_coefficient(name, n)
        for n in range(p + 1, N + 1)
    )


def _horner_fraction_interval(x: Interval, coeffs) -> Interval:
    y = IF(coeffs[-1])
    for c in reversed(coeffs[:-1]):
        y = IF(c) + x * y
    return y


def _validated_exp_upper_fraction(q: Fraction) -> Fraction:
    """Rigorous rational upper enclosure of exp(q), q>=0, by range reduction."""
    if q < 0:
        raise ValueError("nonnegative exponential majorant required")
    scale = 1
    while True:
        z = q / F(scale)
        zf = float(z)
        if F.from_float(zf) < z:
            zf = math.nextafter(zf, math.inf)
        if math.isfinite(zf) and zf <= VT.MAX_ABS_ARGUMENT:
            break
        scale *= 2
    upper = F.from_float(VT.exp_point(zf).hi)
    s = scale
    while s > 1:
        upper *= upper
        s //= 2
    return upper


def _exp_tail_bound(rate: int, xmax: Fraction, order: int) -> Fraction:
    if rate < 0 or xmax < 0 or order < 0:
        raise ValueError("invalid exponential remainder arguments")
    q = F(rate) * xmax
    return _validated_exp_upper_fraction(q) * q ** (order + 1) / F(math.factorial(order + 1))


def _large_word_entry(name: str, x: Interval) -> Interval:
    """Exact exponential-branch B entry on 0.01 <= X <= 2.5."""
    if not (SCALED.BRANCH_X <= x.lo <= x.hi <= WORD_EXACT_SERIES_MAX_X):
        raise ValueError("word exact-series interval outside [0.01,2.5]")
    p = SCALED._DEN_POWER[name]
    N = SCALED.NEAR_EXACT_SERIES_ORDER
    y = _horner_fraction_interval(x, _series_terms(name))

    xmax = F.from_float(float(x.hi))
    xmin = F.from_float(float(x.lo))
    tail_num = F(0)
    exp_terms, _pure = SCALED._EXACT[name]
    for rate, pcoeffs in exp_terms:
        for j, pj in pcoeffs:
            M = N - j
            if M < 0:
                tail_num += abs(pj) * xmax ** j * _validated_exp_upper_fraction(F(rate) * xmax)
            else:
                tail_num += abs(pj) * xmax ** j * _exp_tail_bound(rate, xmax, M)
    tail_B = tail_num / (xmin ** (p + 1))
    t = up(float(tail_B))
    return Interval(math.nextafter(y.lo - t, -math.inf), math.nextafter(y.hi + t, math.inf))


def _large_word_normalized_matrix(x: Interval):
    names = (
        ("vv", "vp", "vS", "va"),
        ("vp", "pp", "pS", "pa"),
        ("vS", "pS", "SS", "Sa"),
        ("va", "pa", "Sa", "aa"),
    )
    return [[_large_word_entry(name, x) for name in row] for row in names]


def word_normalized_matrix(x: Interval):
    """Enclose B(X)=Q_scaled(X)/X for positive X<=2.5."""
    if not (0.0 < x.lo <= x.hi <= WORD_EXACT_SERIES_MAX_X):
        raise ValueError("word X interval outside audited positive range")
    if x.hi < SCALED.BRANCH_X:
        return SCALED.small_normalized_matrix(x)
    if x.lo >= SCALED.BRANCH_X:
        return _large_word_normalized_matrix(x)
    left = SCALED.small_normalized_matrix(
        Interval(x.lo, math.nextafter(SCALED.BRANCH_X, -math.inf))
    )
    right = _large_word_normalized_matrix(Interval(SCALED.BRANCH_X, x.hi))
    return [[hull(left[i][j], right[i][j]) for j in range(4)] for i in range(4)]


@lru_cache(maxsize=200000)
def _translation_information_cached(lo: float, hi: float):
    x = Interval(float(lo), float(hi))
    B = word_normalized_matrix(x)
    xlo = I(x.lo)
    Qlo = matrix_symmetric_hull(
        [[xlo * B[i][j] for j in range(4)] for i in range(4)]
    )
    conditioned = matrix_symmetric_hull(matrix_mul(matrix_mul(C, Qlo), CT))
    if not symmetric_positive_definite_ldlt(conditioned)[0]:
        return None
    try:
        inv = matrix_inverse_gauss_jordan(conditioned)
    except Exception:
        return None
    info = matrix_symmetric_hull(matrix_mul(matrix_mul(CT, inv), C))
    return tuple(tuple(info[i][j] for j in range(4)) for i in range(4))


def translation_information_upper(x: Interval):
    """Dominating Q_scaled(X)^-1 family on one certifiable X subcell."""
    block = _translation_information_cached(float(x.lo), float(x.hi))
    if block is None:
        return None
    return [list(row) for row in block]


def translation_information_cover(
    x: Interval, *, max_depth: int = DEFAULT_INFORMATION_SPLIT_DEPTH
):
    """Cover X by subcells whose conditioned information inverse validates."""
    if not (0.0 < x.lo <= x.hi <= WORD_EXACT_SERIES_MAX_X):
        raise ValueError("translation information cover outside audited X range")
    out = []

    def rec(cell: Interval, depth: int) -> None:
        info = translation_information_upper(cell)
        if info is not None:
            out.append((cell, info))
            return
        if depth >= max_depth:
            raise RuntimeError(
                f"cannot certify word translation information cell {cell.as_list()} "
                f"within split depth {max_depth}"
            )
        mid = math.sqrt(cell.lo * cell.hi)
        if not (cell.lo < mid < cell.hi):
            mid = 0.5 * (cell.lo + cell.hi)
        if not (cell.lo < mid < cell.hi):
            raise RuntimeError(f"word translation information cell cannot split: {cell.as_list()}")
        # Closed children overlap at mid and therefore cover the real parent.
        rec(Interval(cell.lo, mid), depth + 1)
        rec(Interval(mid, cell.hi), depth + 1)

    rec(x, 0)
    if not out:
        raise RuntimeError("translation information cover is empty")
    return out


def translation_margin_from_information(
    information, sigma_root: list[float], measurement_information_diag: list[float]
) -> float:
    """Lower delta for (Omega^-1+D)^-1 >= delta*Sigma."""
    if len(sigma_root) != 4 or len(measurement_information_diag) != 4:
        raise ValueError("translation comparison requires four coordinates")
    if any(not (math.isfinite(v) and v > 0.0) for v in sigma_root):
        raise ValueError("Sigma roots must be finite positive")
    if any(not (math.isfinite(v) and v >= 0.0) for v in measurement_information_diag):
        raise ValueError("measurement information must be finite nonnegative")
    g = [I(v) for v in sigma_root]
    M = [[g[i] * information[i][j] * g[j] for j in range(4)] for i in range(4)]
    for i in range(4):
        M[i][i] = M[i][i] + I(measurement_information_diag[i])
    top = symmetric_gershgorin_upper(matrix_symmetric_hull(M))
    if not (math.isfinite(top) and top > 0.0):
        return 0.0
    return down(1.0 / top)


def word_step_doublings(word_horizon_lower_s: float, dt_s: float) -> int:
    if not (word_horizon_lower_s > 0.0 and dt_s > 0.0):
        raise ValueError("positive word horizon and dt required")
    steps = math.floor(down(word_horizon_lower_s / up(dt_s)))
    if steps < 1:
        raise RuntimeError("word does not certainly contain one prediction")
    return int(math.floor(math.log2(steps)))


def _diag(values: list[float]):
    n = len(values)
    return [[I(values[i]) if i == j else I(0.0) for j in range(n)] for i in range(n)]


def _measurement_posterior(Omega, information: list[float]):
    n = len(Omega)
    if len(information) != n:
        raise ValueError("information dimension mismatch")
    D = _diag(information)
    M = matrix_add(matrix_identity(n), matrix_mul(D, Omega))
    Minv = matrix_inverse_gauss_jordan(M)
    return matrix_symmetric_hull(matrix_mul(Omega, Minv))


def attitude_bias_word_noise(
    rho: float, coupling_per_step: float, doublings: int, attitude_information: float
):
    """Posterior lower floor for one scaled (theta,b_g) axis pair."""
    if not (0.0 < rho <= 1.0):
        raise ValueError("rho must lie in (0,1]")
    if doublings < 0 or not (math.isfinite(attitude_information) and attitude_information >= 0.0):
        raise ValueError("invalid doubling count or information")
    Phi = [[I(1.0), I(-abs(coupling_per_step))], [I(0.0), I(1.0)]]
    Omega = [[I(rho), I(0.0)], [I(0.0), I(rho)]]
    A = Phi
    for _ in range(doublings):
        Omega = matrix_symmetric_hull(
            matrix_add(matrix_mul(matrix_mul(A, Omega), matrix_transpose(A)), Omega)
        )
        A = matrix_mul(A, A)
    return _measurement_posterior(Omega, [attitude_information, 0.0])


def generalized_delta(Omega, Sigma, gate: float = 1.0e-18) -> float:
    """Validated generalized lower margin Omega >= delta*Sigma."""
    if gate <= 0.0:
        raise ValueError("positive search gate required")

    def holds(delta: float) -> bool:
        q = I(delta)
        A = matrix_symmetric_hull([
            [Omega[i][j] - q * Sigma[i][j] for j in range(len(Omega))]
            for i in range(len(Omega))
        ])
        return symmetric_positive_definite_ldlt(A)[0]

    if holds(gate):
        lo = hi = gate
        while hi < 1.0:
            trial = min(1.0, hi * 10.0)
            if not holds(trial):
                hi = trial
                break
            lo = hi = trial
            if trial == 1.0:
                return down(lo)
    else:
        hi = gate
        lo = gate
        for _ in range(80):
            lo /= 10.0
            if holds(lo):
                break
        else:
            return 0.0

    for _ in range(48):
        mid = math.sqrt(lo * hi)
        if holds(mid):
            lo = mid
        else:
            hi = mid
    return down(lo)
