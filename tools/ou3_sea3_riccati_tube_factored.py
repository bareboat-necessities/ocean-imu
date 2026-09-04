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


def _exact_scaled_entry(x, spec):
    c2, a, poly, k = spec
    coeff = _exact_coefficients(c2, a, poly)
    if any(coeff[n] != 0 for n in range(k + 1)):
        raise RuntimeError("exact OU scaled-series cancellation identity failed")

    terms = [
        (n - k, float(coeff[n]))
        for n in range(k + 1, N_EXACT + 1)
        if coeff[n] != 0
    ]
    y = _p(x, terms)

    # Rigorous remainder of the exact exponential representation.  All exact-
    # branch cells in this proof satisfy 0 < x <= 0.012..., well inside the
    # validated exponential backend's audited argument range.
    X = BASE.Interval.outward_bounds(0.0, x.hi)
    X2 = BASE.Interval.outward_bounds(0.0, 2.0 * x.hi)
    exp1 = BASE.VT.exp_interval(X)
    exp2 = BASE.VT.exp_interval(X2)

    rem = BASE.I(0.0)
    if c2:
        rem = rem + (
            BASE.I(abs(float(c2))) * exp2 * BASE.ipow(X2, N_EXACT + 1)
            / BASE.I(float(factorial(N_EXACT + 1)))
        )
    xN1 = BASE.ipow(X, N_EXACT + 1)
    for j, aj in enumerate(a):
        if not aj:
            continue
        rem = rem + (
            BASE.I(abs(float(aj))) * exp1 * xN1
            / BASE.I(float(factorial(N_EXACT - j + 1)))
        )

    denom = BASE.ipow(BASE.I(x.lo), k).lo if k else 1.0
    if not (denom > 0.0):
        raise RuntimeError("exact OU scaled-series denominator lost positivity")
    r = BASE.up(rem.hi / denom)
    return y + BASE.Interval.outward_bounds(-r, r)


def _exact_scaled_q_stable(x):
    vv = _exact_scaled_entry(x, _EXACT["vv"])
    vp = _exact_scaled_entry(x, _EXACT["vp"])
    vS = _exact_scaled_entry(x, _EXACT["vS"])
    va = _exact_scaled_entry(x, _EXACT["va"])
    pp = _exact_scaled_entry(x, _EXACT["pp"])
    pS = _exact_scaled_entry(x, _EXACT["pS"])
    pa = _exact_scaled_entry(x, _EXACT["pa"])
    SS = _exact_scaled_entry(x, _EXACT["SS"])
    Sa = _exact_scaled_entry(x, _EXACT["Sa"])
    aa = _exact_scaled_entry(x, _EXACT["aa"])
    return [
        [vv, vp, vS, va],
        [vp, pp, pS, pa],
        [vS, pS, SS, Sa],
        [va, pa, Sa, aa],
    ]


def step_scaled_q(x):
    if x.hi < BASE.BRANCH_X:
        return _small_scaled_q_factored(x)
    if x.lo >= BASE.BRANCH_X:
        return _exact_scaled_q_stable(x)

    # Outward cells touching the source branch threshold must cover both actual
    # C++ branches.  Evaluate them separately with the stable algebra and hull.
    left_hi = BASE.math.nextafter(BASE.BRANCH_X, -BASE.math.inf)
    families = []
    if x.lo <= left_hi:
        families.append(
            _small_scaled_q_factored(BASE.Interval.outward_bounds(x.lo, left_hi))
        )
    if x.hi >= BASE.BRANCH_X:
        families.append(
            _exact_scaled_q_stable(BASE.Interval.outward_bounds(BASE.BRANCH_X, x.hi))
        )
    return [
        [BASE.hull(*(A[i][j] for A in families)) for j in range(4)]
        for i in range(4)
    ]


# Patch only the numerical primitive used by BASE.build/split_x_cell.  All
# theorem/source semantics remain BASE's implementation.
BASE.step_scaled_q = step_scaled_q

SCHEMA = BASE.SCHEMA
QUALIFICATION = BASE.QUALIFICATION
USEFUL_GATE = BASE.USEFUL_GATE
DEFAULT_DOMAIN = BASE.DEFAULT_DOMAIN


def build(domain_path=DEFAULT_DOMAIN):
    return BASE.build(domain_path)


def validate(payload):
    return BASE.validate(payload)


def main():
    return BASE.main()


if __name__ == "__main__":
    raise SystemExit(main())
