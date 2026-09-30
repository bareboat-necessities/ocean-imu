"""Linked soft returns with the full covariance baseline retained.

The exact rational identities are proof algebra, not a shipping certificate.
A rank-one *information* prior does not mean a rank-one covariance prior.
See docs/ou3-linked-soft-return.md. No next-root precision is applied as a
new observation when checking invariance at that same root.
"""
from __future__ import annotations

from fractions import Fraction as F
import math

from .lin_path_certificate import inverse
from .matrix_certificates import add, identity, is_psd, ldlt, matmul, transpose


def _q(value):
    try:
        return F(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError("finite rational-compatible scalar required") from exc


def _matrix(value):
    result = [[_q(x) for x in row] for row in value]
    if not result or not result[0] or any(len(row) != len(result[0]) for row in result):
        raise ValueError("nonempty rectangular matrix required")
    return result


def _scale(matrix, scalar):
    return [[scalar*x for x in row] for row in matrix]


def _column(value, size):
    result = [[_q(x)] for x in value]
    if len(result) != size:
        raise ValueError("vector dimension mismatch")
    return result


def soft_scalar_return(d_perp, coupling_sq, information, ceiling, prior_norm_sq=1):
    """Exact Schur return d + c*ell^2/(|nu|^2+c*j); no source promotion.

    d_perp includes BOTH the known-root covariance and the eliminated quotient
    contribution. prior_norm_sq carries the declared physical covector scale.
    The fixed-point polynomial has the sign of c - return, not the reverse.
    """
    d, b, j, c, s = map(_q, (d_perp, coupling_sq, information, ceiling, prior_norm_sq))
    if min(d, b, j) < 0 or min(c, s) <= 0:
        raise ValueError("d,b,j must be nonnegative; c,s must be positive")
    value = d + c*b/(s+c*j)
    return {
        "variance": value,
        "ratio": value/c,
        "invariance_margin": c-value,
        "derivative": b*s/(s+c*j)**2,
        "fixed_point_polynomial": j*c*c+(s-b-j*d)*c-s*d,
        "source_uniform_O2_closed": False,
    }


def _conditional_baseline(baseline, information):
    b, j = _matrix(baseline), _matrix(information)
    ldlt(b)
    if len(j) != len(b) or not is_psd(j):
        raise ValueError("information must be PSD with the baseline dimension")
    c = inverse(add(inverse(b), j))
    # L=(I+B J)^-1; H=J-J C J is symmetric PSD even for singular J.
    transport = add(identity(len(b)), matmul(c, j), F(-1))
    action = add(j, matmul(matmul(j, c), j), F(-1))
    return b, c, transport, action


def conditional_rank_one(baseline, information, direction, variance):
    """Condition B+a*y*y' on REAL information J, retaining B and its cross terms.

    Returns C+a/(1+a*h)*z*z', where C=(B^-1+J)^-1,
    z=(I+B J)^-1 y and h=y'(J-J C J)y. No fictitious c_next is inserted.
    """
    b, c, transport, action = _conditional_baseline(baseline, information)
    a = _q(variance)
    if a < 0:
        raise ValueError("nonnegative variance required")
    y = _column(direction, len(b))
    z = matmul(transport, y)
    h = matmul(transpose(y), matmul(action, y))[0][0]
    coefficient = a/(1+a*h)
    return {
        "baseline": c,
        "transport": transport,
        "action": action,
        "filtered_direction": z,
        "effective_information": h,
        "excess_coefficient": coefficient,
        "posterior": add(c, _scale(matmul(z, transpose(z)), coefficient)),
        "source_uniform_O2_closed": False,
    }


def linked_directional_max(baseline, information, functional, variance, basis=None):
    """Sharp directional maximum for B+p*u*u', |u|=1, u in Range(E).

    E can be a nonorthonormal full-column-rank homogeneous chart. The exact
    generalized-eigenvalue metric is E'E+p E'H E, not I unless E'E=I.
    The returned extremizer is a covariance, so no irrational normalization
    or numerical eigensolver is required for this rational identity.
    """
    b, c, transport, action = _conditional_baseline(baseline, information)
    p = _q(variance)
    if p < 0:
        raise ValueError("nonnegative variance required")
    n = _column(functional, len(b))
    e = identity(len(b)) if basis is None else _matrix(basis)
    if len(e) != len(b):
        raise ValueError("basis dimension mismatch")
    et = transpose(e)
    gram = matmul(et, e)
    ldlt(gram)
    metric = add(gram, _scale(matmul(et, matmul(action, e)), p))
    ldlt(metric)
    v = matmul(et, matmul(transpose(transport), n))
    coefficients = matmul(inverse(metric), v)
    raw = matmul(e, coefficients)
    norm_sq = matmul(transpose(raw), raw)[0][0]
    extra = p*matmul(transpose(v), coefficients)[0][0]
    base = matmul(transpose(n), matmul(c, n))[0][0]
    if norm_sq == 0:
        # Every allowed direction attains the same functional; retain trace p.
        raw = [[row[0]] for row in e]
        norm_sq = matmul(transpose(raw), raw)[0][0]
    extremizer = add(b, _scale(matmul(raw, transpose(raw)), p/norm_sq))
    return {
        "baseline_variance": base,
        "max_excess": extra,
        "max_variance": base+extra,
        "generalized_metric": metric,
        "extremizing_prior": extremizer,
        "source_uniform_O2_closed": False,
    }


def linked_cauchy(d: float, h: float, direct: float):
    if not all(math.isfinite(v) for v in (d, h, direct)) or d < 0 or h < 0:
        raise ValueError("finite nonnegative d,h required")
    return {"linked_product": d*h, "direct_sq": direct*direct,
            "cauchy_holds": direct*direct <= d*h*(1+1e-14)}


def rank_one_soft_return(a: float, g2: float, rho: float, lambda2: float, c_next: float):
    """Legacy SYNTHETIC observation update, not a shipping soft-return proof.

    This algebra applies an additional precision 1/c_next to an otherwise
    rank-one covariance. Neither operation follows from Corollary K. Retained
    to distinguish its valid scalar algebra from the invalid O2 inference.
    """
    if (not all(math.isfinite(v) for v in (a, g2, rho, lambda2, c_next))
            or a < 0 or g2 < 0 or lambda2 <= 0 or c_next <= 0 or abs(rho) > 1+1e-15):
        raise ValueError("invalid covariance/information parameters")
    x = min(1.0, max(0.0, rho*rho))
    prior = a*g2
    info = x/c_next+lambda2*(1.0-x)
    post = prior/(1.0+prior*info)
    return {"propagated_prior_variance": prior, "effective_information": info,
            "posterior_rank_one_variance": post,
            "next_kernel_variance": post*x,
            "next_quotient_variance": post*(1.0-x),
            "synthetic_observation_update_only": True,
            "proves_shipping_invariance": False}


def linked_product_status():
    return {"independent_dbar_Hbar_required": False,
            "same_word_Pi_cancels_generally": False,
            "full_covariance_baseline_required": True,
            "next_ceiling_is_an_applied_measurement": False,
            "source_uniform_O2_closed": False,
            "replacement": "full-baseline rank-one excess and linked generalized eigenvalue"}
