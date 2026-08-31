#!/usr/bin/env python3
"""Conditioned process-matrix overlay for source-reachable P3.

The deployed integrated-OU covariance is extremely ill-conditioned in the raw
monomial comparison basis even after the known powers of x=h/tau are cancelled.
This overlay keeps the source formulas exact but certifies the process shape in
an exact rational congruence basis

    C = R L^{-1},

where L is the exact LDL factor of the x->0 process-shape Gramian and
R=diag(1,10,100,2).  The limiting transformed Gramian is exactly

diag(2/3, 5/8, 200/567, 1/2).

For every source x-cell we prove

    C (Q_scaled/x) C' >= rho_tilde I,

with outward-rounded interval LDLT.  Since ||C||_2^2 <= ||C||_1||C||_inf =
37310 exactly, this implies in the original scaled coordinates

    Q_scaled >= x_lo * rho_tilde / 37310 I.

No filter, source domain, comparison metric, or P3 usefulness gate is changed.
The small-x source polynomial is evaluated after symbolic power cancellation;
the >=0.01 source branch uses exact rational Taylor coefficients plus an
explicit validated exponential remainder after the same cancellation.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import functools
import json
import math
from pathlib import Path

from ou3_interval import Interval, hull
import ou3_source_reachable_matrix_p3 as BASE
import ou3_validated_transcendentals as VT

BRANCH_X = BASE.BRANCH_X
DEFAULT_DOMAIN = BASE.DEFAULT_DOMAIN
SCHEMA = BASE.SCHEMA
MIN_USEFUL_DELTA = BASE.MIN_USEFUL_DELTA
EXACT_SERIES_ORDER = 22
F = Fraction


def _FI(q: Fraction | int | float) -> Interval:
    q = q if isinstance(q, Fraction) else Fraction(q)
    f = float(q)
    return Interval.outward_bounds(f, f)


def _poly_value(x: Interval, coeffs: tuple[Fraction | int | float, ...]) -> Interval:
    """Outward Horner evaluation of c0+c1*x+... ."""
    y = _FI(coeffs[-1])
    for c in reversed(coeffs[:-1]):
        y = _FI(c) + x * y
    return y


def _factored_source_series(x: Interval, coeffs: tuple[Fraction | int | float, ...]) -> Interval:
    """Evaluate x*(c0+c1*x+...) after cancelling the known raw powers."""
    return x * _poly_value(x, coeffs)


_SMALL_COEFFS = {
    "vv": (F(2,3),-F(1,2),F(7,30),-F(1,12),F(31,1260),-F(1,160),F(127,90720)),
    "vp": (F(1,4),-F(1,6),F(5,72),-F(1,45),F(17,2880),-F(41,30240)),
    "vS": (F(1,15),-F(1,24),F(41,2520),-F(7,1440),F(109,90720)),
    "va": (F(1),-F(1),F(7,12),-F(1,4),F(31,360),-F(1,40),F(127,20160),-F(17,12096)),
    "pp": (F(1,10),-F(1,18),F(5,252),-F(1,180),F(17,12960)),
    "pS": (F(1,36),-F(1,72),F(13,2880),-F(1,864)),
    "pa": (F(1,3),-F(1,3),F(11,60),-F(13,180),F(19,840),-F(1,168),F(247,181440)),
    "SS": (F(1,126),-F(1,288),F(13,12960)),
    "Sa": (F(1,12),-F(1,12),F(2,45),-F(1,60),F(11,2240),-F(73,60480)),
    "aa": (F(2),-F(2),F(4,3),-F(2,3),F(4,15),-F(4,45),F(8,315),-F(2,315),F(4,2835)),
}


def _small_scaled_q(x: Interval):
    if not (x.lo > 0.0 and x.hi < BRANCH_X):
        raise ValueError("small source-polynomial helper requires x<0.01")
    p = {k:_factored_source_series(x,v) for k,v in _SMALL_COEFFS.items()}
    return [
        [p["vv"],p["vp"],p["vS"],p["va"]],
        [p["vp"],p["pp"],p["pS"],p["pa"]],
        [p["vS"],p["pS"],p["SS"],p["Sa"]],
        [p["va"],p["pa"],p["Sa"],p["aa"]],
    ]


def _small_process_shape(x: Interval):
    """Return Q_scaled/x directly, with x cancelled before interval arithmetic."""
    if not (x.lo > 0.0 and x.hi < BRANCH_X):
        raise ValueError("small process-shape helper requires x<0.01")
    p = {k:_poly_value(x,v) for k,v in _SMALL_COEFFS.items()}
    return [
        [p["vv"],p["vp"],p["vS"],p["va"]],
        [p["vp"],p["pp"],p["pS"],p["pa"]],
        [p["vS"],p["pS"],p["SS"],p["Sa"]],
        [p["va"],p["pa"],p["Sa"],p["aa"]],
    ]


@functools.lru_cache(maxsize=512)
def _exact_series_coefficients(denominator_power: int, exp_terms, polynomial):
    """Exact Taylor coefficients of ``numerator/x^denominator_power``.

    They depend only on the source expression, not on the evaluation cell, so
    the exact ``Fraction`` arithmetic is done once.  The word-horizon read of
    this family evaluates it on two orders of magnitude more cells than the
    per-step read did, and rebuilding these rationals on every cell was the
    whole cost of that.
    """
    N = EXACT_SERIES_ORDER
    coeff = [F(0) for _ in range(N+1)]
    for degree,value in polynomial:
        if degree <= N:
            coeff[degree] += value
    for lam,p in exp_terms:
        for j,pj in p:
            for n in range(j,N+1):
                k=n-j
                coeff[n] += pj * F((-lam)**k, math.factorial(k))
    if any(coeff[n] != 0 for n in range(min(denominator_power,N+1))):
        raise RuntimeError("exponential source cancellation order mismatch")
    # Carry the float enclosures alongside the exact values: the Horner pass
    # below converts every coefficient on every cell otherwise.
    return tuple((c, _FI(c)) for c in coeff[denominator_power:])


def _exact_scaled_entry(
    x: Interval,
    denominator_power: int,
    exp_terms: tuple[tuple[int, dict[int, Fraction]], ...],
    polynomial: dict[int, Fraction],
) -> Interval:
    """Validated exact-series enclosure of numerator/x^denominator_power.

    ``exp_terms`` encodes p(x) exp(-lambda*x).  Exact Fraction arithmetic forms
    all Taylor coefficients.  The only truncation is the exponential tail,
    bounded by exp(lambda*x)<2 on the reachable source domain.
    """
    if x.lo <= 0.0:
        raise ValueError("positive x required")
    N = EXACT_SERIES_ORDER
    scaled = _exact_series_coefficients(
        denominator_power,
        tuple((lam, tuple(sorted(pp.items()))) for lam, pp in exp_terms),
        tuple(sorted(polynomial.items())),
    )
    y=scaled[-1][1]
    for _, ci in reversed(scaled[:-1]):
        y=ci+x*y

    rem=0.0
    for lam,p in exp_terms:
        # The tail of exp(-lam*x) after order k is bounded by its next term times
        # exp(lam*x.hi).  The literal 2.0 this used to carry is that factor at
        # lam*x < ln 2, which held for the per-step read but silently expires the
        # moment the same family is read at a word horizon.
        growth=BASE.up(math.exp(lam*x.hi)*(1.0+1.0e-12))
        for j,pj in p.items():
            k=N-j+1
            if k<=0:
                continue
            term=(
                growth*abs(float(pj))*(float(lam)**k)*(x.hi**(N+1))
                /(math.factorial(k)*(x.lo**denominator_power))
            )
            rem=BASE.up(rem+term)
    return Interval(BASE.down(y.lo-rem),BASE.up(y.hi+rem))


# Each tuple is (raw comparison denominator power, exponential terms, polynomial).
_LARGE_EXPR = {
    "vv": (2,((2,{0:-F(1)}),(1,{0:F(4)})),{1:F(2),0:-F(3)}),
    "vp": (3,((2,{0:F(1)}),(1,{1:F(2),0:-F(2)})),{2:F(1),1:-F(2),0:F(1)}),
    "vS": (4,((2,{0:-F(1)}),(1,{2:F(1),0:F(4)})),{3:F(1,3),2:-F(1),1:F(2),0:-F(3)}),
    "va": (1,((2,{0:F(1)}),(1,{0:-F(2)})),{0:F(1)}),
    "pp": (4,((2,{0:-F(1)}),(1,{1:-F(4)})),{3:F(2,3),2:-F(2),1:F(2),0:F(1)}),
    "pS": (5,((2,{0:F(1)}),(1,{2:-F(1),1:F(2),0:-F(2)})),{4:F(1,4),3:-F(1),2:F(2),1:-F(2),0:F(1)}),
    "pa": (2,((2,{0:-F(1)}),(1,{1:-F(2)})),{0:F(1)}),
    "SS": (6,((2,{0:-F(1)}),(1,{2:F(2),0:F(4)})),{5:F(1,10),4:-F(1,2),3:F(4,3),2:-F(2),1:F(2),0:-F(3)}),
    "Sa": (3,((2,{0:F(1)}),(1,{2:-F(1),0:-F(2)})),{0:F(1)}),
    "aa": (0,((2,{0:-F(1)}),),{0:F(1)}),
}


def _large_family(x: Interval, extra_x_power: int):
    if x.lo < BRANCH_X:
        raise ValueError("exact exponential helper requires x>=0.01")
    p={}
    for k,(shift,exp_terms,polynomial) in _LARGE_EXPR.items():
        p[k]=_exact_scaled_entry(x,shift+extra_x_power,exp_terms,polynomial)
    return [
        [p["vv"],p["vp"],p["vS"],p["va"]],
        [p["vp"],p["pp"],p["pS"],p["pa"]],
        [p["vS"],p["pS"],p["SS"],p["Sa"]],
        [p["va"],p["pa"],p["Sa"],p["aa"]],
    ]


def _large_scaled_q(x: Interval):
    return _large_family(x,0)


def _large_process_shape(x: Interval):
    # One additional exact cancelled power returns Q_scaled/x directly.
    return _large_family(x,1)




def _branch_hull(x: Interval, small_fn, large_fn):
    if x.hi < BRANCH_X:
        return small_fn(x)
    if x.lo >= BRANCH_X:
        return large_fn(x)
    left_hi=math.nextafter(BRANCH_X,-math.inf)
    families=[]
    if x.lo<=left_hi:
        families.append(small_fn(Interval(x.lo,left_hi)))
    if x.hi>=BRANCH_X:
        families.append(large_fn(Interval(BRANCH_X,x.hi)))
    return [[hull(*(A[i][j] for A in families)) for j in range(4)] for i in range(4)]


def step_scaled_q(x: Interval):
    return _branch_hull(x,_small_scaled_q,_large_scaled_q)


def process_shape_q(x: Interval):
    return _branch_hull(x,_small_process_shape,_large_process_shape)


# Exact x->0 shape Gramian in [v,p,S,a_w] comparison coordinates.
_LIMIT_SHAPE = (
    (F(2,3),F(1,4),F(1,15),F(1)),
    (F(1,4),F(1,10),F(1,36),F(1,3)),
    (F(1,15),F(1,36),F(1,126),F(1,12)),
    (F(1),F(1,3),F(1,12),F(2)),
)
_L_INV = (
    (F(1),F(0),F(0),F(0)),
    (-F(3,8),F(1),F(0),F(0)),
    (F(1,15),-F(4,9),F(1),F(0)),
    (-F(15,2),F(30),-F(105,2),F(1)),
)
_R_DIAG=(F(1),F(10),F(100),F(2))
_C_RATIONAL=tuple(tuple(_R_DIAG[i]*_L_INV[i][j] for j in range(4)) for i in range(4))
_LIMIT_CONGRUENCE_DIAG=(F(2,3),F(5,8),F(200,567),F(1,2))


def _qmatmul(A,B):
    return tuple(tuple(sum((A[i][k]*B[k][j] for k in range(len(B))),F(0)) for j in range(len(B[0]))) for i in range(len(A)))


def _qtranspose(A):
    return tuple(tuple(A[j][i] for j in range(len(A))) for i in range(len(A[0])))


def _verify_exact_congruence():
    transformed=_qmatmul(_qmatmul(_C_RATIONAL,_LIMIT_SHAPE),_qtranspose(_C_RATIONAL))
    target=tuple(tuple(_LIMIT_CONGRUENCE_DIAG[i] if i==j else F(0) for j in range(4)) for i in range(4))
    if transformed != target:
        raise RuntimeError("exact RL^-1 limiting congruence identity failed")


_verify_exact_congruence()
_C_NORM_INF=max(sum((abs(v) for v in row),F(0)) for row in _C_RATIONAL)
_C_NORM_ONE=max(sum((abs(_C_RATIONAL[i][j]) for i in range(4)),F(0)) for j in range(4))
_C_NORM2_SQ_UPPER=_C_NORM_INF*_C_NORM_ONE
if (_C_NORM_INF,_C_NORM_ONE,_C_NORM2_SQ_UPPER)!=(F(182),F(205),F(37310)):
    raise RuntimeError("unexpected exact congruence norm bound")
_C_INTERVAL=[[_FI(v) for v in row] for row in _C_RATIONAL]


def _congruence_process_rho(x: Interval) -> tuple[float,float]:
    """Return (original-scaled rho, transformed-shape rho) for one x-cell."""
    shape=process_shape_q(x)
    transformed=BASE.matrix_symmetric_hull(
        BASE.matrix_mul(BASE.matrix_mul(_C_INTERVAL,shape),BASE.matrix_transpose(_C_INTERVAL))
    )
    rho_tilde=BASE.certified_rho(transformed)
    if rho_tilde<=0.0:
        return 0.0,0.0
    # C A C' >= rho_tilde I => A >= rho_tilde/||C||_2^2 I.
    # Q_scaled=x*A and x>=x.lo throughout the source cell.
    numerator=BASE.down(x.lo*rho_tilde)
    rho=BASE.down(numerator/float(_C_NORM2_SQ_UPPER))
    return rho,rho_tilde


def split_x_cell(x: Interval, depth: int=0):
    rho,_=_congruence_process_rho(x)
    if rho>0.0:
        return [(x,rho)]
    if depth>=14:
        raise RuntimeError(f"cannot certify RL^-1-preconditioned OU process cell {x.as_list()}")
    mid=math.sqrt(x.lo*x.hi)
    return split_x_cell(Interval.outward_bounds(x.lo,mid),depth+1)+split_x_cell(Interval.outward_bounds(mid,x.hi),depth+1)


# Patch only the process representation/certifier before BASE builds source cells.
BASE.step_scaled_q=step_scaled_q
BASE.split_x_cell=split_x_cell
BASE._build_cached.cache_clear()


def _frac_text(q: Fraction) -> str:
    return str(q.numerator) if q.denominator==1 else f"{q.numerator}/{q.denominator}"


def build(domain_path: Path=DEFAULT_DOMAIN) -> dict:
    d=BASE.build(domain_path)
    out=dict(d)
    out["small_x_process_representation"]="FACTORED_SOURCE_POLYNOMIAL_AFTER_EXACT_X_CANCELLATION"
    out["large_x_process_representation"]="EXACT_RATIONAL_TAYLOR_AFTER_CANCELLATION_PLUS_VALIDATED_EXPONENTIAL_REMAINDER"
    out["process_representation_algebraically_source_equivalent"]=True
    out["process_congruence_preconditioner"]={
        "form":"C=R*L_inverse",
        "exact_rational":True,
        "R_diagonal":[_frac_text(v) for v in _R_DIAG],
        "L_inverse":[[_frac_text(v) for v in row] for row in _L_INV],
        "C_RL_inverse":[[_frac_text(v) for v in row] for row in _C_RATIONAL],
        "limiting_transformed_diagonal":[_frac_text(v) for v in _LIMIT_CONGRUENCE_DIAG],
        "norm_inf_exact":int(_C_NORM_INF),
        "norm_1_exact":int(_C_NORM_ONE),
        "norm_2_squared_upper_exact":int(_C_NORM2_SQ_UPPER),
        "translation_formula":"rho_Qscaled = x_lo * rho_tilde / 37310",
        "certified_object":"C*(Q_scaled/x)*C^T",
    }
    return out


def validate(d: dict) -> list[str]:
    failures=BASE.validate(d)
    if d.get("small_x_process_representation")!="FACTORED_SOURCE_POLYNOMIAL_AFTER_EXACT_X_CANCELLATION":
        failures.append("small-x process representation is not factored")
    if d.get("large_x_process_representation")!="EXACT_RATIONAL_TAYLOR_AFTER_CANCELLATION_PLUS_VALIDATED_EXPONENTIAL_REMAINDER":
        failures.append("large-x process representation is not conditioned")
    if d.get("process_representation_algebraically_source_equivalent") is not True:
        failures.append("conditioned process representation not source equivalent")
    p=d.get("process_congruence_preconditioner",{})
    if p.get("form")!="C=R*L_inverse" or p.get("exact_rational") is not True:
        failures.append("exact RL^-1 congruence preconditioner missing")
    if p.get("norm_inf_exact")!=182 or p.get("norm_1_exact")!=205 or p.get("norm_2_squared_upper_exact")!=37310:
        failures.append("RL^-1 congruence norm translation changed")
    if p.get("limiting_transformed_diagonal") != ["2/3","5/8","200/567","1/2"]:
        failures.append("RL^-1 limiting congruence diagonal changed")
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
        "process_congruence_preconditioner":d["process_congruence_preconditioner"],
        "cells":d["cell_partition"],
        "H":d["modes"]["H"],
        "A":d["modes"]["A"],
        "failures":failures,
    },indent=2,sort_keys=True))
    return 0 if not failures else 2


if __name__=="__main__":
    raise SystemExit(main())
