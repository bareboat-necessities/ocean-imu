"""Source-locked binary32 tau-target correspondence for the default shipping tuner.

The exact shadow used for deployment roundoff must start from the *compiled
binary32 operands*, not from ideal decimal spellings. Shipping evaluates

    f = clamp(f_source, 0.03f, 1.2f)
    tau_raw = 1.0f * 0.5f / f
    tau_target = clamp(tau_raw, 0.02f, 12.0f)

This module interprets those binary32 operands as exact rationals, evaluates the
same expression once in exact arithmetic and once with RNE after each shipping
operation. On the only unclamped target range, 0.5/f is in (0.4,12), hence below
16 and one binary32 division has absolute RNE error <=2^-21. Above 12 both
branches share the 12f clamp except within that same final half-ulp cell. The
multiplication 1.0f*0.5f is exact.

This is independent of WPE/libm correctness: ``f_source`` is merely required to
be an already-stored finite positive binary32 value. Upstream production of
that value remains a separate obligation.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
FLOOR=B.rn32(F(3,100)); CEIL=B.rn32(F(6,5))
TAU_COEFF=B.rn32(1); HALF=B.rn32(F(1,2))
TAU_MIN=B.rn32(F(1,50)); TAU_MAX=B.rn32(12)
TARGET_ERROR_MAX=F(1,1<<21)
QUALIFICATION='OU3_ALT_SHIPPING_TAU_TARGET_BINARY32_V2'


def clamp(x,lo,hi): return min(max(x,lo),hi)


@dataclass(frozen=True)
class TargetPair:
    source_frequency:F
    clamped_frequency:F
    exact_target:F
    binary32_target:F
    error:F
    def __post_init__(self):
        names=('source_frequency','clamped_frequency','exact_target','binary32_target','error')
        vals=[F(getattr(self,n)) for n in names]
        for n,v in zip(names,vals): object.__setattr__(self,n,v)
        if not B.is_binary32(self.source_frequency) or not B.is_binary32(self.clamped_frequency):
            raise ValueError('source/clamped frequency must be deployed binary32')
        if not B.is_binary32(self.binary32_target):
            raise ValueError('shipping tau target must be deployed binary32')
        if self.error != self.binary32_target-self.exact_target:
            raise ValueError('tau target error detached from same source expression')
        if abs(self.error)>TARGET_ERROR_MAX:
            raise ValueError('shipping tau target exceeds certified binary32 rounding cell')


def evaluate(source_frequency):
    f0=F(source_frequency)
    if not B.is_binary32(f0) or f0<=0:
        raise ValueError('positive stored binary32 source frequency required')
    f=clamp(f0,FLOOR,CEIL)
    exact_raw=(TAU_COEFF*HALF)/f
    exact_target=clamp(exact_raw,TAU_MIN,TAU_MAX)
    binary_raw=B.div(B.mul(TAU_COEFF,HALF),f)
    binary_target=clamp(binary_raw,TAU_MIN,TAU_MAX)
    return TargetPair(f0,f,exact_target,binary_target,binary_target-exact_target)


def _source_shape_matches():
    s=SOURCE.read_text()
    needles=(
      'constexpr float MIN_TUNE_FREQ_HZ = 0.03f;',
      'constexpr float MAX_TUNE_FREQ_HZ = 1.2f;',
      'constexpr float MIN_TAU_S   = 0.02f;',
      'constexpr float MAX_TAU_S   = 12.0f;',
      'float tau_coeff_    = 1.0f;',
      'float tau_raw = tau_coeff_ * 0.5f / f_tune;',
      'tau_target_   = std::min(std::max(tau_raw,  min_tau_s_), max_tau_s_);')
    return all(n in s for n in needles)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_tau_frequency_and_target_constant_source_shape_matches':_source_shape_matches(),
      'compiled_binary32_constants_used_as_exact_shadow_operands':True,
      'shipping_tau_target_exact_vs_binary32_cell_bound_source_locked':_source_shape_matches(),
      'target_error_abs_bound_seconds':TARGET_ERROR_MAX,
      'upstream_WPE_binary32_frequency_production_closed':False,
      'shipping_tau_predecessor_domain_inductively_closed':False,
      'source_uniform_tau_roundoff_supply_bound_closed':False,
      'complete_word_finite_identity':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
