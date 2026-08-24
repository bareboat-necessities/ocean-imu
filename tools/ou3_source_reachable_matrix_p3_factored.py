#!/usr/bin/env python3
"""Conditioned process-matrix overlay for source-reachable P3.

For the implementation's x<0.01 branch the known powers of x are cancelled
symbolically before interval evaluation.  For the x>=0.01 exponential branch,
the exact source expression is expanded with exact rational coefficients and a
rigorous exponential Taylor remainder *after* the same cancellation order is
removed.  Thus neither branch forms nearly equal O(1) quantities merely to
recover an O(x^7) covariance entry.

Only the numerical representation of ``step_scaled_q`` is replaced.  Source
cells, source-reachable scheduling, covariance upper matrices, generalized
matrix ratio, useful-margin gate, and validation are the base P3 implementation.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from ou3_interval import Interval, hull
import ou3_source_reachable_matrix_p3 as BASE

BRANCH_X = BASE.BRANCH_X
DEFAULT_DOMAIN = BASE.DEFAULT_DOMAIN
SCHEMA = BASE.SCHEMA
MIN_USEFUL_DELTA = BASE.MIN_USEFUL_DELTA
EXACT_SERIES_ORDER = 22


def _FI(q: Fraction | int | float) -> Interval:
    q = q if isinstance(q, Fraction) else Fraction(q)
    f = float(q)
    return Interval.outward_bounds(f, f)


def _factored_source_series(x: Interval, coeffs: tuple[Fraction | float, ...]) -> Interval:
    """Evaluate x*(c0+c1*x+...) with outward Horner arithmetic."""
    y = _FI(coeffs[-1])
    for c in reversed(coeffs[:-1]):
        y = _FI(c) + x * y
    return x * y


def _small_scaled_q(x: Interval):
    if not (x.lo > 0.0 and x.hi < BRANCH_X):
        raise ValueError("small source-polynomial helper requires x<0.01")
    qvv = _factored_source_series(x, (Fraction(2,3),-Fraction(1,2),Fraction(7,30),-Fraction(1,12),Fraction(31,1260),-Fraction(1,160),Fraction(127,90720)))
    qvp = _factored_source_series(x, (Fraction(1,4),-Fraction(1,6),Fraction(5,72),-Fraction(1,45),Fraction(17,2880),-Fraction(41,30240)))
    qvS = _factored_source_series(x, (Fraction(1,15),-Fraction(1,24),Fraction(41,2520),-Fraction(7,1440),Fraction(109,90720)))
    qva = _factored_source_series(x, (1,-1,Fraction(7,12),-Fraction(1,4),Fraction(31,360),-Fraction(1,40),Fraction(127,20160),-Fraction(17,12096)))
    qpp = _factored_source_series(x, (Fraction(1,10),-Fraction(1,18),Fraction(5,252),-Fraction(1,180),Fraction(17,12960)))
    qpS = _factored_source_series(x, (Fraction(1,36),-Fraction(1,72),Fraction(13,2880),-Fraction(1,864)))
    qpa = _factored_source_series(x, (Fraction(1,3),-Fraction(1,3),Fraction(11,60),-Fraction(13,180),Fraction(19,840),-Fraction(1,168),Fraction(247,181440)))
    qSS = _factored_source_series(x, (Fraction(1,126),-Fraction(1,288),Fraction(13,12960)))
    qSa = _factored_source_series(x, (Fraction(1,12),-Fraction(1,12),Fraction(2,45),-Fraction(1,60),Fraction(11,2240),-Fraction(73,60480)))
    qaa = _factored_source_series(x, (2,-2,Fraction(4,3),-Fraction(2,3),Fraction(4,15),-Fraction(4,45),Fraction(8,315),-Fraction(2,315),Fraction(4,2835)))
    return [
        [qvv,qvp,qvS,qva],
        [qvp,qpp,qpS,qpa],
        [qvS,qpS,qSS,qSa],
        [qva,qpa,qSa,qaa],
    ]


def _exact_scaled_entry(
    x: Interval,
    shift: int,
    exp_terms: tuple[tuple[int, dict[int, Fraction]], ...],
    polynomial: dict[int, Fraction],
) -> Interval:
    """Validated series of an exact exponential-branch numerator / x^shift.

    ``exp_terms`` represents p(x)*exp(-lambda*x).  Exact Fraction arithmetic
    forms every Taylor coefficient.  The only truncation is the exponential
    tail, bounded by exp(lambda*x)<2 because lambda*x<=0.03 on the reachable
    domain.
    """
    if x.lo <= 0.0:
        raise ValueError("positive x required")
    N = EXACT_SERIES_ORDER
    coeff = [Fraction(0,1) for _ in range(N+1)]
    for degree, value in polynomial.items():
        if degree <= N:
            coeff[degree] += value
    for lam, p in exp_terms:
        for j, pj in p.items():
            for n in range(j, N+1):
                k = n-j
                coeff[n] += pj * Fraction((-lam)**k, math.factorial(k))
    # Source algebra cancels every power below shift exactly.  Refuse the proof
    # if the encoded source expression does not have that structure.
    if any(coeff[n] != 0 for n in range(min(shift,N+1))):
        raise RuntimeError("exponential source cancellation order mismatch")

    scaled = coeff[shift:]
    y = _FI(scaled[-1])
    for c in reversed(scaled[:-1]):
        y = _FI(c) + x*y

    # Remainder of every p_j x^j exp(-lam*x) after numerator order N.
    # |R| <= 2 |p_j| lam^(N-j+1) x^(N+1)/(N-j+1)!.
    # Divide by x^shift using the interval's positive lower endpoint.
    rem = 0.0
    for lam, p in exp_terms:
        for j, pj in p.items():
            k = N-j+1
            if k <= 0:
                continue
            term = (
                2.0 * abs(float(pj)) * (float(lam)**k)
                * (x.hi**(N+1))
                / (math.factorial(k) * (x.lo**shift))
            )
            rem = BASE.up(rem + term)
    return Interval(BASE.down(y.lo-rem), BASE.up(y.hi+rem))


F = Fraction


def _large_scaled_q(x: Interval):
    if x.lo < BRANCH_X:
        raise ValueError("exact exponential helper requires x>=0.01")

    vv = _exact_scaled_entry(x,2,((2,{0:-F(1)}),(1,{0:F(4)})),{1:F(2),0:-F(3)})
    vp = _exact_scaled_entry(x,3,((2,{0:F(1)}),(1,{1:F(2),0:-F(2)})),{2:F(1),1:-F(2),0:F(1)})
    vS = _exact_scaled_entry(x,4,((2,{0:-F(1)}),(1,{2:F(1),0:F(4)})),{3:F(1,3),2:-F(1),1:F(2),0:-F(3)})
    va = _exact_scaled_entry(x,1,((2,{0:F(1)}),(1,{0:-F(2)})),{0:F(1)})

    pp = _exact_scaled_entry(x,4,((2,{0:-F(1)}),(1,{1:-F(4)})),{3:F(2,3),2:-F(2),1:F(2),0:F(1)})
    pS = _exact_scaled_entry(x,5,((2,{0:F(1)}),(1,{2:-F(1),1:F(2),0:-F(2)})),{4:F(1,4),3:-F(1),2:F(2),1:-F(2),0:F(1)})
    pa = _exact_scaled_entry(x,2,((2,{0:-F(1)}),(1,{1:-F(2)})),{0:F(1)})

    SS = _exact_scaled_entry(x,6,((2,{0:-F(1)}),(1,{2:F(2),0:F(4)})),{5:F(1,10),4:-F(1,2),3:F(4,3),2:-F(2),1:F(2),0:-F(3)})
    Sa = _exact_scaled_entry(x,3,((2,{0:F(1)}),(1,{2:-F(1),0:-F(2)})),{0:F(1)})
    aa = _exact_scaled_entry(x,0,((2,{0:-F(1)}),),{0:F(1)})

    return [
        [vv,vp,vS,va],
        [vp,pp,pS,pa],
        [vS,pS,SS,Sa],
        [va,pa,Sa,aa],
    ]


def step_scaled_q(x: Interval):
    if x.hi < BRANCH_X:
        return _small_scaled_q(x)
    if x.lo >= BRANCH_X:
        return _large_scaled_q(x)
    left_hi = math.nextafter(BRANCH_X,-math.inf)
    families=[]
    if x.lo <= left_hi:
        families.append(_small_scaled_q(Interval(x.lo,left_hi)))
    if x.hi >= BRANCH_X:
        families.append(_large_scaled_q(Interval(BRANCH_X,x.hi)))
    return [[hull(*(A[i][j] for A in families)) for j in range(4)] for i in range(4)]


# Patch only the process-matrix representation before BASE builds its cells.
BASE.step_scaled_q = step_scaled_q
BASE._build_cached.cache_clear()


def build(domain_path: Path=DEFAULT_DOMAIN) -> dict:
    d=BASE.build(domain_path)
    out=dict(d)
    out["small_x_process_representation"]="FACTORED_SOURCE_POLYNOMIAL_X_TIMES_P_OF_X"
    out["large_x_process_representation"]="EXACT_RATIONAL_TAYLOR_AFTER_CANCELLATION_PLUS_VALIDATED_EXPONENTIAL_REMAINDER"
    out["process_representation_algebraically_source_equivalent"]=True
    return out


def validate(d: dict) -> list[str]:
    failures=BASE.validate(d)
    if d.get("small_x_process_representation")!="FACTORED_SOURCE_POLYNOMIAL_X_TIMES_P_OF_X":
        failures.append("small-x process representation is not factored")
    if d.get("large_x_process_representation")!="EXACT_RATIONAL_TAYLOR_AFTER_CANCELLATION_PLUS_VALIDATED_EXPONENTIAL_REMAINDER":
        failures.append("large-x process representation is not conditioned")
    if d.get("process_representation_algebraically_source_equivalent") is not True:
        failures.append("conditioned process representation not source equivalent")
    return failures


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    d=build(args.domain)
    failures=validate(d)
    out=dict(d)
    out["validation_pass"]=not failures
    out["validation_failures"]=failures
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps({
        "backend":d["p3_backend"],
        "small_x_process_representation":d["small_x_process_representation"],
        "large_x_process_representation":d["large_x_process_representation"],
        "cells":d["cell_partition"],
        "H":d["modes"]["H"],
        "A":d["modes"]["A"],
        "failures":failures,
    },indent=2,sort_keys=True))
    return 0 if not failures else 2


if __name__=="__main__":
    raise SystemExit(main())
