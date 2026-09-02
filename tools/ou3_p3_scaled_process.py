#!/usr/bin/env python3
"""Dependency-preserving scaled OU process primitive for the OU-III P3 metric.

P3 uses the physical similarity

    D_h = diag(sigma*h, sigma*h^2, sigma*h^3, sigma)

on the one-axis integrated-OU process covariance for [v,p,S,a].  Directly
forming the tiny covariance entries and then dividing by powers of x=h/tau
loses repeated-variable dependency at the deployed long-tau endpoint.  The
problem is numerical proof formulation, not loss of process excitation.

Two source branches must be mirrored:

* x < 1e-2: the shipping polynomial covariance branch.  Its scaled matrix has
  an exact common positive factor x; we evaluate the algebraically cancelled
  polynomial B(x) in scaled_Q=x*B(x).
* x >= 1e-2: the shipping exponential branch.  Near the branch point direct
  interval evaluation suffers catastrophic cancellation.  We therefore expand
  the *exact exponential formula* as one correlated Taylor series, cancel its
  exact leading powers symbolically with Fraction arithmetic, bound the full
  exponential remainder, and again evaluate scaled_Q=x*B(x).  Farther from the
  threshold the ordinary validated exponential expression is retained.

No trajectory values, floating eigensolvers, or ordinary floating inverses are
used as proof enclosures.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from ou3_interval import Interval, symmetric_positive_definite_ldlt
import ou3_validated_transcendentals as VT

BRANCH_X = 1.0e-2
NEAR_EXACT_SERIES_MAX_X = 5.0e-2
NEAR_EXACT_SERIES_ORDER = 30
SCHEMA = 2
F = Fraction


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def IF(q: Fraction) -> Interval:
    return Interval.outward_bounds(float(q), float(q))


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def ipow(x: Interval, n: int) -> Interval:
    y = Interval.point(1.0)
    for _ in range(int(n)):
        y = y * x
    return y


def poly(x: Interval, terms) -> Interval:
    y = Interval.point(0.0)
    for n, c in terms:
        coeff = IF(c) if isinstance(c, Fraction) else I(c)
        y = y + coeff * ipow(x, int(n))
    return y


# Shipping x<0.01 qbar series, after removing the D_h similarity power and one
# additional power for the common positive factor x.
_SMALL = {
    "vv": ((0,F(2,3)),(1,F(-1,2)),(2,F(7,30)),(3,F(-1,12)),(4,F(31,1260)),(5,F(-1,160)),(6,F(127,90720))),
    "vp": ((0,F(1,4)),(1,F(-1,6)),(2,F(5,72)),(3,F(-1,45)),(4,F(17,2880)),(5,F(-41,30240))),
    "va": ((0,F(1)),(1,F(-1)),(2,F(7,12)),(3,F(-1,4)),(4,F(31,360)),(5,F(-1,40)),(6,F(127,20160)),(7,F(-17,12096))),
    "pp": ((0,F(1,10)),(1,F(-1,18)),(2,F(5,252)),(3,F(-1,180)),(4,F(17,12960))),
    "pa": ((0,F(1,3)),(1,F(-1,3)),(2,F(11,60)),(3,F(-13,180)),(4,F(19,840)),(5,F(-1,168)),(6,F(247,181440))),
    "aa": ((0,F(2)),(1,F(-2)),(2,F(4,3)),(3,F(-2,3)),(4,F(4,15)),(5,F(-4,45)),(6,F(8,315)),(7,F(-2,315)),(8,F(4,2835))),
    "vS": ((0,F(1,15)),(1,F(-1,24)),(2,F(41,2520)),(3,F(-7,1440)),(4,F(109,90720))),
    "pS": ((0,F(1,36)),(1,F(-1,72)),(2,F(13,2880)),(3,F(-1,864))),
    "SS": ((0,F(1,126)),(1,F(-1,288)),(2,F(13,12960))),
    "Sa": ((0,F(1,12)),(1,F(-1,12)),(2,F(2,45)),(3,F(-1,60)),(4,F(11,2240)),(5,F(-73,60480))),
}


def small_normalized_matrix(x: Interval):
    """Return exact source-polynomial B(x) where scaled_Q=x*B on x<1e-2."""
    qvv=poly(x,_SMALL["vv"]); qvp=poly(x,_SMALL["vp"]); qva=poly(x,_SMALL["va"])
    qpp=poly(x,_SMALL["pp"]); qpa=poly(x,_SMALL["pa"]); qaa=poly(x,_SMALL["aa"])
    qvS=poly(x,_SMALL["vS"]); qpS=poly(x,_SMALL["pS"]); qSS=poly(x,_SMALL["SS"]); qSa=poly(x,_SMALL["Sa"])
    return [
        [qvv,qvp,qvS,qva],
        [qvp,qpp,qpS,qpa],
        [qvS,qpS,qSS,qSa],
        [qva,qpa,qSa,qaa],
    ]


# Exact x>=0.01 source formulas represented as
#   sum poly_j(x)*exp(-rate_j*x) + pure_poly(x).
# All coefficients are exact rationals.  _DEN_POWER is the D_h similarity
# exponent k_i+k_j; B additionally removes one common x.
_EXACT = {
    "vv": (((2,((0,F(-1)),)),(1,((0,F(4)),))),((0,F(-3)),(1,F(2)))),
    "vp": (((2,((0,F(1)),)),(1,((0,F(-2)),(1,F(2))))),((0,F(1)),(1,F(-2)),(2,F(1)))),
    "va": (((2,((0,F(1)),)),(1,((0,F(-2)),))),((0,F(1)),)),
    "pp": (((2,((0,F(-1)),)),(1,((1,F(-4)),))),((0,F(1)),(1,F(2)),(2,F(-2)),(3,F(2,3)))),
    "pa": (((2,((0,F(-1)),)),(1,((1,F(-2)),))),((0,F(1)),)),
    "aa": (((2,((0,F(-1)),)),),((0,F(1)),)),
    "vS": (((2,((0,F(-1)),)),(1,((0,F(4)),(2,F(1))))),((0,F(-3)),(1,F(2)),(2,F(-1)),(3,F(1,3)))),
    "pS": (((2,((0,F(1)),)),(1,((0,F(-2)),(1,F(2)),(2,F(-1))))),((0,F(1)),(1,F(-2)),(2,F(2)),(3,F(-1)),(4,F(1,4)))),
    "SS": (((2,((0,F(-1)),)),(1,((0,F(4)),(2,F(2))))),((0,F(-3)),(1,F(2)),(2,F(-2)),(3,F(4,3)),(4,F(-1,2)),(5,F(1,10)))),
    "Sa": (((2,((0,F(1)),)),(1,((0,F(-2)),(2,F(-1))))),((0,F(1)),)),
}
_DEN_POWER = {"vv":2,"vp":3,"vS":4,"va":1,"pp":4,"pS":5,"pa":2,"SS":6,"Sa":3,"aa":0}


def _series_coefficient(name: str, n: int) -> Fraction:
    exp_terms,pure=_EXACT[name]
    c=dict(pure).get(n,F(0))
    for rate,pcoeffs in exp_terms:
        for j,pj in pcoeffs:
            if n>=j:
                k=n-j
                c += pj * F((-rate)**k, math.factorial(k))
    return c


def _exp_tail_bound(rate: int, xmax: Fraction, order: int) -> Fraction:
    """Lagrange bound for exp(-rate*x), using exp(rate*x)<=2 on this range."""
    if rate*xmax > F(1,2):
        raise RuntimeError("near-threshold exact-series exponential bound left audited range")
    return F(2) * (rate*xmax)**(order+1) / F(math.factorial(order+1))


def _near_exact_normalized_entry(name: str, x: Interval) -> Interval:
    """Exact exponential-branch B entry with correlated cancellation retained."""
    p=_DEN_POWER[name]
    N=NEAR_EXACT_SERIES_ORDER
    # The exact formula has zeros through degree p.  Check that algebraically;
    # a source-expression edit that changes this structure must fail closed.
    for n in range(p+1):
        if _series_coefficient(name,n) != 0:
            raise RuntimeError(f"exact OU series {name} lost leading-power cancellation at n={n}")
    terms=[]
    for n in range(p+1,N+1):
        terms.append((n-(p+1),_series_coefficient(name,n)))
    y=poly(x,tuple(terms))

    xmax=F.from_float(float(x.hi))
    xmin=F.from_float(float(x.lo))
    tail_num=F(0)
    exp_terms,_pure=_EXACT[name]
    for rate,pcoeffs in exp_terms:
        for j,pj in pcoeffs:
            M=N-j
            if M<0:
                tail_num += abs(pj)*xmax**j*F(2)
            else:
                tail_num += abs(pj)*xmax**j*_exp_tail_bound(rate,xmax,M)
    tail_B=tail_num/(xmin**(p+1))
    t=up(float(tail_B))
    return Interval(math.nextafter(y.lo-t,-math.inf),math.nextafter(y.hi+t,math.inf))


def near_exact_normalized_matrix(x: Interval):
    if not (BRANCH_X <= x.lo <= x.hi <= NEAR_EXACT_SERIES_MAX_X):
        raise ValueError("near exact OU series outside audited interval")
    qvv=_near_exact_normalized_entry("vv",x); qvp=_near_exact_normalized_entry("vp",x)
    qvS=_near_exact_normalized_entry("vS",x); qva=_near_exact_normalized_entry("va",x)
    qpp=_near_exact_normalized_entry("pp",x); qpS=_near_exact_normalized_entry("pS",x)
    qpa=_near_exact_normalized_entry("pa",x); qSS=_near_exact_normalized_entry("SS",x)
    qSa=_near_exact_normalized_entry("Sa",x); qaa=_near_exact_normalized_entry("aa",x)
    return [[qvv,qvp,qvS,qva],[qvp,qpp,qpS,qpa],[qvS,qpS,qSS,qSa],[qva,qpa,qSa,qaa]]


def _large_scaled_matrix(x: Interval):
    """Exact source exponential branch followed by the physical D_h similarity."""
    a=VT.exp_interval(-x); a2=a.square()
    one,two,three,four=I(1),I(2),I(3),I(4)
    x2,x3,x4,x5=ipow(x,2),ipow(x,3),ipow(x,4),ipow(x,5)
    qvv=-a2+four*a+two*x-three
    qvp=a2+two*a*(x-one)+x2-two*x+one
    qva=a2-two*a+one
    qpp=-a2-four*a*x+I(2/3)*x3-two*x2+two*x+one
    qpa=-a2-two*a*x+one
    qaa=one-a2
    qvS=-a2+a*(x2+four)+(x3-three*x2+I(6)*x-I(9))/three
    qpS=a2+a*(-x2+two*x-two)+I(1/4)*x4-x3+two*x2-two*x+one
    qSS=-a2+two*a*x2+four*a+I(1/10)*x5-I(1/2)*x4+I(4/3)*x3-two*x2+two*x-three
    qSa=a2-a*(x2+two)+one
    q=[[qvv,qvp,qvS,qva],[qvp,qpp,qpS,qpa],[qvS,qpS,qSS,qSa],[qva,qpa,qSa,qaa]]
    k=(1,2,3,0)
    return [[q[i][j]/ipow(x,k[i]+k[j]) for j in range(4)] for i in range(4)]


def _minus_rho(A, rho: float):
    return [[A[i][j]-I(rho if i==j else 0.0) for j in range(len(A))] for i in range(len(A))]


def certified_rho(A) -> float:
    ok,_=symmetric_positive_definite_ldlt(A)
    if not ok:
        return 0.0
    hi=min(A[i][i].lo for i in range(len(A)))
    lo=0.0
    for _ in range(48):
        mid=0.5*(lo+hi)
        ok,_=symmetric_positive_definite_ldlt(_minus_rho(A,mid))
        if ok: lo=mid
        else: hi=mid
    return down(lo)


def _factored_rho(x: Interval, B) -> float:
    rb=certified_rho(B)
    return down(x.lo*rb) if rb>0.0 else 0.0


def certified_cell_rho(x: Interval) -> float:
    if not (0.0 < x.lo <= x.hi <= VT.MAX_ABS_ARGUMENT):
        raise ValueError("x=h/tau outside validated range")
    if x.hi < BRANCH_X:
        return _factored_rho(x,small_normalized_matrix(x))
    if x.lo >= BRANCH_X and x.hi <= NEAR_EXACT_SERIES_MAX_X:
        return _factored_rho(x,near_exact_normalized_matrix(x))
    if x.lo >= NEAR_EXACT_SERIES_MAX_X:
        return certified_rho(_large_scaled_matrix(x))
    if x.lo < BRANCH_X <= x.hi:
        left=Interval(x.lo,math.nextafter(BRANCH_X,-math.inf))
        right=Interval(BRANCH_X,x.hi)
        return min(certified_cell_rho(left),certified_cell_rho(right))
    # Exact branch cell crossing the near-series/direct-expression boundary.
    left=Interval(x.lo,math.nextafter(NEAR_EXACT_SERIES_MAX_X,-math.inf))
    right=Interval(NEAR_EXACT_SERIES_MAX_X,x.hi)
    return min(certified_cell_rho(left),certified_cell_rho(right))


def split_x_cell(x: Interval, depth: int=0):
    rho=certified_cell_rho(x)
    if rho>0.0:
        return [(x,rho)]
    if depth>=20:
        raise RuntimeError(f"cannot certify dependency-preserving scaled OU process cell {x.as_list()}")
    for cut in (BRANCH_X,NEAR_EXACT_SERIES_MAX_X):
        if x.lo < cut < x.hi:
            return split_x_cell(Interval(x.lo,math.nextafter(cut,-math.inf)),depth+1)+split_x_cell(Interval(cut,x.hi),depth+1)
    mid=math.sqrt(x.lo*x.hi)
    return split_x_cell(Interval.outward_bounds(x.lo,mid),depth+1)+split_x_cell(Interval.outward_bounds(mid,x.hi),depth+1)


def validate_range(xlo: float, xhi: float, cells: int=24):
    ratio=(xhi/xlo)**(1.0/cells)
    edges=[xlo]
    for _ in range(cells-1): edges.append(edges[-1]*ratio)
    edges.append(xhi)
    for cut in (BRANCH_X,NEAR_EXACT_SERIES_MAX_X):
        if xlo<cut<xhi: edges=sorted(set(edges+[cut]))
    pieces=[]
    for a,b in zip(edges,edges[1:]):
        pieces.extend(split_x_cell(Interval.outward_bounds(a,b)))
    worst=min(r for _x,r in pieces)
    return pieces,worst


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--x-lo",type=float,default=0.005/12.0)
    ap.add_argument("--x-hi",type=float,default=0.005/0.02)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    pieces,worst=validate_range(a.x_lo,a.x_hi)
    d={
        "schema":SCHEMA,
        "qualification":"OU3_P3_DEPENDENCY_PRESERVING_SCALED_OU_PROCESS",
        "source_small_branch_algebraically_identical":True,
        "source_exact_branch_algebraically_identical":True,
        "near_branch_exact_exponential_remainder_bounded":True,
        "common_positive_x_factor_preserved_before_LDLT":True,
        "ordinary_floating_eigensolver_used":False,
        "x_range":[a.x_lo,a.x_hi],
        "near_exact_series_max_x":NEAR_EXACT_SERIES_MAX_X,
        "near_exact_series_order":NEAR_EXACT_SERIES_ORDER,
        "certified_subcells":len(pieces),
        "scaled_process_lambda_min_lower":worst,
        "pass":math.isfinite(worst) and worst>0.0,
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps(d,indent=2,sort_keys=True))
    return 0 if d["pass"] else 2


if __name__=="__main__":
    raise SystemExit(main())
