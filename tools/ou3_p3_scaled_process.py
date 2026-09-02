#!/usr/bin/env python3
"""Dependency-preserving scaled OU process primitive for the OU-III P3 metric.

The P3 similarity uses D_h=diag(sigma*h, sigma*h^2, sigma*h^3, sigma).
On the shipping small-x branch every entry of D_h^-1 Q D_h^-T has a common
positive factor x=h/tau.  Evaluating Q first and then interval-dividing by
powers of x destroys that dependency near the deployed long-tau endpoint.

This module evaluates the algebraically cancelled series directly and factors
out the common positive x before the interval LDLT test.  It is mathematically
identical to the source polynomial branch; it changes only the validated
arithmetic expression.  The >=1e-2 source branch retains the exact exponential
formula.  No replay values, floating eigensolvers, or non-rigorous inverses are
used.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, hull, symmetric_positive_definite_ldlt
import ou3_validated_transcendentals as VT

BRANCH_X = 1.0e-2
SCHEMA = 1


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def ipow(x: Interval, n: int) -> Interval:
    y = Interval.point(1.0)
    for _ in range(int(n)):
        y = y * x
    return y


def poly(x: Interval, terms) -> Interval:
    y = Interval.point(0.0)
    for n, c in terms:
        y = y + I(c) * ipow(x, int(n))
    return y


# Source small-x qbar series, with the D_h similarity power removed and then
# one additional power removed for the common positive x factor.
_SMALL = {
    "vv": ((0,2/3),(1,-1/2),(2,7/30),(3,-1/12),(4,31/1260),(5,-1/160),(6,127/90720)),
    "vp": ((0,1/4),(1,-1/6),(2,5/72),(3,-1/45),(4,17/2880),(5,-41/30240)),
    "va": ((0,1),(1,-1),(2,7/12),(3,-1/4),(4,31/360),(5,-1/40),(6,127/20160),(7,-17/12096)),
    "pp": ((0,1/10),(1,-1/18),(2,5/252),(3,-1/180),(4,17/12960)),
    "pa": ((0,1/3),(1,-1/3),(2,11/60),(3,-13/180),(4,19/840),(5,-1/168),(6,247/181440)),
    "aa": ((0,2),(1,-2),(2,4/3),(3,-2/3),(4,4/15),(5,-4/45),(6,8/315),(7,-2/315),(8,4/2835)),
    "vS": ((0,1/15),(1,-1/24),(2,41/2520),(3,-7/1440),(4,109/90720)),
    "pS": ((0,1/36),(1,-1/72),(2,13/2880),(3,-1/864)),
    "SS": ((0,1/126),(1,-1/288),(2,13/12960)),
    "Sa": ((0,1/12),(1,-1/12),(2,2/45),(3,-1/60),(4,11/2240),(5,-73/60480)),
}


def small_normalized_matrix(x: Interval):
    """Return B(x) where scaled_Q(x) = x * B(x) on x<1e-2."""
    qvv=poly(x,_SMALL["vv"]); qvp=poly(x,_SMALL["vp"]); qva=poly(x,_SMALL["va"])
    qpp=poly(x,_SMALL["pp"]); qpa=poly(x,_SMALL["pa"]); qaa=poly(x,_SMALL["aa"])
    qvS=poly(x,_SMALL["vS"]); qpS=poly(x,_SMALL["pS"]); qSS=poly(x,_SMALL["SS"]); qSa=poly(x,_SMALL["Sa"])
    return [
        [qvv,qvp,qvS,qva],
        [qvp,qpp,qpS,qpa],
        [qvS,qpS,qSS,qSa],
        [qva,qpa,qSa,qaa],
    ]


def _large_scaled_matrix(x: Interval):
    """Exact source >=1e-2 branch followed by the physical D_h similarity."""
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


def certified_cell_rho(x: Interval) -> float:
    if not (0.0 < x.lo <= x.hi <= VT.MAX_ABS_ARGUMENT):
        raise ValueError("x=h/tau outside validated range")
    if x.hi < BRANCH_X:
        rb=certified_rho(small_normalized_matrix(x))
        return down(x.lo*rb) if rb>0.0 else 0.0
    if x.lo >= BRANCH_X:
        return certified_rho(_large_scaled_matrix(x))
    left=Interval(x.lo,math.nextafter(BRANCH_X,-math.inf))
    right=Interval(BRANCH_X,x.hi)
    return min(certified_cell_rho(left),certified_cell_rho(right))


def split_x_cell(x: Interval, depth: int=0):
    rho=certified_cell_rho(x)
    if rho>0.0:
        return [(x,rho)]
    if depth>=18:
        raise RuntimeError(f"cannot certify dependency-preserving scaled OU process cell {x.as_list()}")
    # Never bisect across the implementation branch when it can be split exactly.
    if x.lo < BRANCH_X < x.hi:
        return split_x_cell(Interval(x.lo,math.nextafter(BRANCH_X,-math.inf)),depth+1)+split_x_cell(Interval(BRANCH_X,x.hi),depth+1)
    mid=math.sqrt(x.lo*x.hi)
    return split_x_cell(Interval.outward_bounds(x.lo,mid),depth+1)+split_x_cell(Interval.outward_bounds(mid,x.hi),depth+1)


def validate_range(xlo: float, xhi: float, cells: int=24):
    ratio=(xhi/xlo)**(1.0/cells)
    edges=[xlo]
    for _ in range(cells-1): edges.append(edges[-1]*ratio)
    edges.append(xhi)
    if xlo<BRANCH_X<xhi: edges=sorted(set(edges+[BRANCH_X]))
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
        "source_series_algebraically_identical":True,
        "common_positive_x_factor_preserved_before_LDLT":True,
        "ordinary_floating_eigensolver_used":False,
        "x_range":[a.x_lo,a.x_hi],
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
