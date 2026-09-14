"""Exact binary32 branch decision for shipping IntegratedOUChain Qaxis formulas.

The shipping float specialization computes

    tau_eff = max(tau, float(1e-7))
    inv     = float(1) / tau_eff
    x       = h * inv
    small   = abs(x) < float(1e-2)

before selecting the polynomial or general process-covariance formula.  The
real-arithmetic ALT graph must not select this branch from exact ``h/tau`` near
the threshold because that can disagree with the deployed float comparison.

This module closes only that control decision for finite positive operands.  It
does not claim that the carried tuner ``tau`` has already been proved equal to
the deployed stored float, nor does it enclose the arithmetic inside either
coefficient formula or Eigen PSD hygiene.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as B32

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_common/KalmanOUCoreMath.h'
TAU_FLOOR=B32.binary32_positive(F(1,10_000_000))
SMALL_THRESHOLD=B32.binary32_positive(F(1,100))
ONE=B32.binary32_positive(F(1))
QUALIFICATION='OU3_ALT_QAXIS_BINARY32_BRANCH_V1'


def rn32_pos(x):
    return B32.binary32_positive(F(x))


def div32_pos(a,b):
    a=rn32_pos(a); b=rn32_pos(b)
    if b<=0: raise ValueError('positive binary32 divisor required')
    return rn32_pos(a/b)


def mul32_pos(a,b):
    a=rn32_pos(a); b=rn32_pos(b)
    return rn32_pos(a*b)


@dataclass(frozen=True)
class Branch:
    tau_input:F
    h_input:F
    tau_float:F
    h_float:F
    tau_eff:F
    inv:F
    x:F
    small:bool
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if self.qualification!=QUALIFICATION: raise ValueError('wrong Qaxis branch qualification')
        if self.tau_input<=0 or self.h_input<=0: raise ValueError('positive tau/h required')
        if self.tau_float!=rn32_pos(self.tau_input) or self.h_float!=rn32_pos(self.h_input):
            raise ValueError('Qaxis branch operands detached from rounded source operands')
        if self.tau_eff!=max(self.tau_float,TAU_FLOOR):
            raise ValueError('Qaxis tau_eff detached from shipping max clamp')
        if self.inv!=div32_pos(ONE,self.tau_eff):
            raise ValueError('Qaxis inv detached from shipping binary32 division')
        if self.x!=mul32_pos(self.h_float,self.inv):
            raise ValueError('Qaxis x detached from shipping binary32 multiply')
        if self.small != (abs(self.x)<SMALL_THRESHOLD):
            raise ValueError('Qaxis formula branch detached from shipping binary32 comparison')


def branch(tau,h):
    ti=F(tau); hi=F(h)
    if ti<=0 or hi<=0: raise ValueError('positive tau/h required')
    tf=rn32_pos(ti); hf=rn32_pos(hi); te=max(tf,TAU_FLOOR)
    inv=div32_pos(ONE,te); x=mul32_pos(hf,inv)
    return Branch(ti,hi,tf,hf,te,inv,x,abs(x)<SMALL_THRESHOLD)


def _source_shape_matches():
    s=SOURCE.read_text()
    needles=(
        'const T tau_eff = std::max(tau, T(1e-7));',
        'const T inv = T(1) / tau_eff;',
        'const T x = h * inv;',
        'if (std::abs(x) < T(1e-2))',
    )
    return all(s.count(n)>=2 for n in needles)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_nested_and_final_Qaxis_share_literal_tau_h_branch_shape':_source_shape_matches(),
      'tau_and_h_rounded_to_binary32_before_branch_decision':True,
      'tau_floor_binary32_attached':True,
      'reciprocal_binary32_rounding_attached':True,
      'x_multiply_binary32_rounding_attached':True,
      'small_general_comparison_binary32_attached':True,
      'source_owned_tau_deployment_commit_correspondence_closed':False,
      'coefficient_formula_binary32_roundoff_closed':False,
      'Eigen_PSD_hygiene_deployment_closed':False,
      'ALT_LIVE_PASS':False,
    }
