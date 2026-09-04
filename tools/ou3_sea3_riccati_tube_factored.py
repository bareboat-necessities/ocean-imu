#!/usr/bin/env python3
"""Numerically stable backend for the canonical SEA3 moving-Riccati tube.

This is not a second theorem route.  It delegates the complete tube construction
and validation to ``ou3_sea3_riccati_tube`` and replaces only the numerical
evaluation order of the scaled integrated-OU process covariance.

For the implementation's small-x branch, powers of x are cancelled
symbolically before interval evaluation.  For the exact branch used at
x >= 0.01, the same closed-form expressions are expanded about x=0 using exact
rational coefficients through order 24 and a validated exponential Taylor
remainder.  On the deployed proof domain x=h/tau <= about 0.012 this removes
catastrophic interval cancellation while rigorously enclosing the same exact
formula.

The scaled process matrix has one additional common positive factor x.  The
adaptive x-cell certifier removes that factor symbolically, proves the slowly
varying normalized matrix B(x)=Qscaled(x)/x positive definite, and then returns
rho = x_lo * lambda_min(B).  This preserves the exact matrix inequality while
preventing interval LDLT from independently widening a scale factor shared by
every entry.

The source contracts, covariance ceiling, Joseph comparison, H/A aggregation,
1e-18 gate and zero-history-depth current-source cover all remain the canonical
BASE implementation.
"""
from __future__ import annotations

from fractions import Fraction
from math import factorial

import ou3_sea3_riccati_tube as BASE


N_EXACT = 24


def _p(x, terms):
    return BASE.poly(x, tuple(terms))


def _small_scaled_q_factored(x):
    """Exact symbolic x-power cancellation for BASE's small-x Qbar series."""
    qvv = _p(x, ((1,2/3),(2,-1/2),(3,7/30),(4,-1/12),(5,31/1260),(6,-1/160),(7,127/90720)))
    qvp = _p(x, ((1,1/4),(2,-1/6),(3,5/72),(4,-1/45),(5,17/2880),(6,-41/30240)))
    qvS = _p(x, ((1,1/15),(2,-1/24),(3,41/2520),(4,-7/1440),(5,109/90720)))
    qva = _p(x, ((1,1),(2,-1),(3,7/12),(4,-1/4),(5,31/360),(6,-1/40),(7,127/20160),(8,-17/12096)))

    qpp = _p(x, ((1,1/10),(2,-1/18),(3,5/252),(4,-1/180),(5,17/12960)))
    qpS = _p(x, ((1,1/36),(2,-1/72),(3,13/2880),(4,-1/864)))
    qpa = _p(x, ((1,1/3),(2,-1/3),(3,11/60),(4,-13/180),(5,19/840),(6,-1/168),(7,247/181440)))

    qSS = _p(x, ((1,1/126),(2,-1/288),(3,13/12960)))
    qSa = _p(x, ((1,1/12),(2,-1/12),(3,2/45),(4,-1/60),(5,11/2240),(6,-73/60480)))

    qaa = _p(x, ((1,2),(2,-2),(3,4/3),(4,-2/3),(5,4/15),(6,-4/45),(7,8/315),(8,-2/315),(9,4/2835)))

    return [
        [qvv, qvp, qvS, qva],
        [qvp, qpp, qpS, qpa],
        [qvS, qpS, qSS, qSa],
        [qva, qpa, qSa, qaa],
    ]


def _small_scaled_q_over_x(x):
    """Same scaled covariance with the common positive factor x removed."""
    qvv = _p(x, ((0,2/3),(1,-1/2),(2,7/30),(3,-1/12),(4,31/1260),(5,-1/160),(6,127/90720)))
    qvp = _p(x, ((0,1/4),(1,-1/6),(2,5/72),(3,-1/45),(4,17/2880),(5,-41/30240)))
    qvS = _p(x, ((0,1/15),(1,-1/24),(2,41/2520),(3,-7/1440),(4,109/90720)))
    qva = _p(x, ((0,1),(1,-1),(2,7/12),(3,-1/4),(4,31/360),(5,-1/40),(6,127/20160),(7,-17/12096)))

    qpp = _p(x, ((0,1/10),(1,-1/18),(2,5/252),(3,-1/180),(4,17/12960)))
    qpS = _p(x, ((0,1/36),(1,-1/72),(2,13/2880),(3,-1/864)))
    qpa = _p(x, ((0,1/3),(1,-1/3),(2,11/60),(3,-13/180),(4,19/840),(5,-1/168),(6,247/181440)))

    qSS = _p(x, ((0,1/126),(1,-1/288),(2,13/12960)))
    qSa = _p(x, ((0,1/12),(1,-1/12),(2,2/45),(3,-1/60),(4,11/2240),(5,-73/60480)))

    qaa = _p(x, ((0,2),(1,-2),(2,4/3),(3,-2/3),(4,4/15),(5,-4/45),(6,8/315),(7,-2/315),(8,4/2835)))

    return [
        [qvv, qvp, qvS, qva],
        [qvp, qpp, qpS, qpa],
        [qvS, qpS, qSS, qSa],
        [qva, qpa, qSa, qaa],
    ]


# Each exact Qbar entry has form
#   c2 exp(-2x) + A(x) exp(-x) + P(x).
# The final integer is the power x^k removed by BASE's scaled coordinate.
# Coefficients are exact rationals, so the cancellations through degree k are
# checked symbolically rather than inferred from floating point.
_EXACT = {
    "vv": (Fraction(-1), [Fraction(4)], {0:Fraction(-3),1:Fraction(2)}, 2),
    "vp": (Fraction(1), [Fraction(-2),Fraction(2)], {0:Fraction(1),1:Fraction(-2),2:Fraction(1)}, 3),
    "vS": (Fraction(-1), [Fraction(4),Fraction(0),Fraction(1)], {0:Fraction(-3),1:Fraction(2),2:Fraction(-1),3:Fraction(1,3)}, 4),
    "va": (Fraction(1), [Fraction(-2)], {0:Fraction(1)}, 1),
    "pp": (Fraction(-1), [Fraction(0),Fraction(-4)], {0:Fraction(1),1:Fraction(2),2:Fraction(-2),3:Fraction(2,3)}, 4),
    "pS": (Fraction(1), [Fraction(-2),Fraction(2),Fraction(-1)], {0:Fraction(1),1:Fraction(-2),2:Fraction(2),3:Fraction(-1),4:Fraction(1,4)}, 5),
    "pa": (Fraction(-1), [Fraction(0),Fraction(-2)], {0:Fraction(1)}, 2),
    "SS": (Fraction(-1), [Fraction(4),Fraction(0),Fraction(2)], {0:Fraction(-3),1:Fraction(2),2:Fraction(-2),3:Fraction(4,3),4:Fraction(-1,2),5:Fraction(1,10)}, 6),
    "Sa": (Fraction(1), [Fraction(-2),Fraction(0),Fraction(-1)], {0:Fraction(1)}, 3),
    "aa": (Fraction(-1), [], {0:Fraction(1)}, 0),
}


def _exact_coefficients(c2, a, poly, nmax=N_EXACT):
    coeff = []
    for n in range(nmax + 1):
        c = c2 * Fraction((-2) ** n, factorial(n))
        for j, aj in enumerate(a):
            if n >= j:
                c += aj * Fraction((-1) ** (n - j), factorial(n - j))
        c += poly.get(n, Fraction(0))
        coeff.append(c)
    return coeff


_EXACT_COEFFS = {
    name: _exact_coefficients(spec[0], spec[1], spec[2])
    for name, spec in _EXACT.items()
}


def _exact_context(x):
    """Shared validated exponential/remainder ingredients for one x cell."""
    X = BASE.Interval.outward_bounds(0.0, x.hi)
    X2 = BASE.Interval.outward_bounds(0.0, 2.0 * x.hi)
    return {
        "X": X,
        "X2": X2,
        "exp1": BASE.VT.exp_interval(X),
        "exp2": BASE.VT.exp_interval(X2),
        "xN1": BASE.ipow(X, N_EXACT + 1),
        "x2N1": BASE.ipow(X2, N_EXACT + 1),
        "factorialN1": BASE.I(float(factorial(N_EXACT + 1))),
    }


def _exact_scaled_entry(x, name, ctx, remove_common_x=False):
    c2, a, _poly, k = _EXACT[name]
    coeff = _EXACT_COEFFS[name]
    if any(coeff[n] != 0 for n in range(k + 1)):
        raise RuntimeError("exact OU scaled-series cancellation identity failed")

    divide_power = k + (1 if remove_common_x else 0)
    terms = [
        (n - divide_power, float(coeff[n]))
        for n in range(k + 1, N_EXACT + 1)
        if coeff[n] != 0
    ]
    y = _p(x, terms)

    # Rigorous remainder of the exact exponential representation.  Reuse the
    # validated exp/power intervals for all ten matrix entries in this x cell.
    rem = BASE.I(0.0)
    if c2:
        rem = rem + (
            BASE.I(abs(float(c2))) * ctx["exp2"] * ctx["x2N1"]
            / ctx["factorialN1"]
        )
    for j, aj in enumerate(a):
        if not aj:
            continue
        rem = rem + (
            BASE.I(abs(float(aj))) * ctx["exp1"] * ctx["xN1"]
            / BASE.I(float(factorial(N_EXACT - j + 1)))
        )

    denom = BASE.ipow(BASE.I(x.lo), divide_power).lo if divide_power else 1.0
    if not (denom > 0.0):
        raise RuntimeError("exact OU scaled-series denominator lost positivity")
    r = BASE.up(rem.hi / denom)
    return y + BASE.Interval.outward_bounds(-r, r)


def _exact_scaled_q_stable(x, remove_common_x=False):
    ctx = _exact_context(x)
    vv = _exact_scaled_entry(x, "vv", ctx, remove_common_x)
    vp = _exact_scaled_entry(x, "vp", ctx, remove_common_x)
    vS = _exact_scaled_entry(x, "vS", ctx, remove_common_x)
    va = _exact_scaled_entry(x, "va", ctx, remove_common_x)
    pp = _exact_scaled_entry(x, "pp", ctx, remove_common_x)
    pS = _exact_scaled_entry(x, "pS", ctx, remove_common_x)
    pa = _exact_scaled_entry(x, "pa", ctx, remove_common_x)
    SS = _exact_scaled_entry(x, "SS", ctx, remove_common_x)
    Sa = _exact_scaled_entry(x, "Sa", ctx, remove_common_x)
    aa = _exact_scaled_entry(x, "aa", ctx, remove_common_x)
    return [
        [vv, vp, vS, va],
        [vp, pp, pS, pa],
        [vS, pS, SS, Sa],
        [va, pa, Sa, aa],
    ]


def _branch_matrix(x, remove_common_x):
    if x.hi < BASE.BRANCH_X:
        return _small_scaled_q_over_x(x) if remove_common_x else _small_scaled_q_factored(x)
    if x.lo >= BASE.BRANCH_X:
        return _exact_scaled_q_stable(x, remove_common_x)

    # Outward cells touching the source branch threshold must cover both actual
    # C++ branches.  Evaluate them separately with the stable algebra and hull.
    left_hi = BASE.math.nextafter(BASE.BRANCH_X, -BASE.math.inf)
    families = []
    if x.lo <= left_hi:
        left = BASE.Interval.outward_bounds(x.lo, left_hi)
        families.append(
            _small_scaled_q_over_x(left)
            if remove_common_x
            else _small_scaled_q_factored(left)
        )
    if x.hi >= BASE.BRANCH_X:
        right = BASE.Interval.outward_bounds(BASE.BRANCH_X, x.hi)
        families.append(_exact_scaled_q_stable(right, remove_common_x))
    return [
        [BASE.hull(*(A[i][j] for A in families)) for j in range(4)]
        for i in range(4)
    ]


def step_scaled_q(x):
    return _branch_matrix(x, False)


def step_scaled_q_over_x(x):
    return _branch_matrix(x, True)


_PROFILE = {
    "x_certifications": 0,
    "x_certified_leaves": 0,
    "x_splits": 0,
    "max_split_depth": 0,
}


def _reset_profile():
    for key in _PROFILE:
        _PROFILE[key] = 0


def split_x_cell(x, depth=0):
    """Certify Qscaled(x) >= rho I after removing its common x factor."""
    _PROFILE["x_certifications"] += 1
    _PROFILE["max_split_depth"] = max(_PROFILE["max_split_depth"], depth)

    # Pointwise Qscaled(x) = x B(x), x>0.  Therefore a validated
    # B(x) >= rho_B I on the whole interval implies
    # Qscaled(x) >= x.lo * rho_B I.  Factoring x before interval LDLT avoids
    # treating the same scale uncertainty independently in all matrix entries.
    rho_b = BASE.certified_rho(step_scaled_q_over_x(x))
    rho = BASE.down(x.lo * rho_b) if rho_b > 0.0 else 0.0
    if rho > 0.0:
        _PROFILE["x_certified_leaves"] += 1
        return [(x, rho)]
    if depth >= 14:
        raise RuntimeError(
            f"cannot certify normalized scaled OU process cell {x.as_list()}"
        )

    _PROFILE["x_splits"] += 1
    mid = BASE.math.sqrt(x.lo * x.hi)
    return split_x_cell(
        BASE.Interval.outward_bounds(x.lo, mid), depth + 1
    ) + split_x_cell(
        BASE.Interval.outward_bounds(mid, x.hi), depth + 1
    )


# Patch only numerical primitives used by BASE.build.  All theorem/source
# semantics remain BASE's implementation.
_BASE_BUILD = BASE.build
BASE.step_scaled_q = step_scaled_q
BASE.split_x_cell = split_x_cell

SCHEMA = BASE.SCHEMA
QUALIFICATION = BASE.QUALIFICATION
USEFUL_GATE = BASE.USEFUL_GATE
DEFAULT_DOMAIN = BASE.DEFAULT_DOMAIN


def build(domain_path=DEFAULT_DOMAIN):
    _reset_profile()
    d = _BASE_BUILD(domain_path)
    d["numerical_profile"] = {
        **_PROFILE,
        "common_positive_x_factor_removed_before_ldlt": True,
        "exact_exponential_context_shared_per_x_cell": True,
        "exact_series_coefficients_precomputed": True,
    }
    return d


# Make BASE.main use the profiled/stabilized build as well.
BASE.build = build


def validate(payload):
    return BASE.validate(payload)


def main():
    return BASE.main()


if __name__ == "__main__":
    raise SystemExit(main())
