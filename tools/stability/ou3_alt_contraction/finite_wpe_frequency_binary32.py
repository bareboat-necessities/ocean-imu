"""Binary32 WPE getter topology and tuner-frequency handoff for ALT.

The exact-real WPE shadow legitimately uses one reciprocal period/frequency
pair.  Shipping does something more specific at deployment: it stores one
binary32 ``log_period_sec_`` and evaluates two distinct libm calls

    std::exp( log_period_sec_)
    std::exp(-log_period_sec_)

for period and frequency.  Their returned floats are therefore retained as two
independent binary32 witnesses.  No bit-level reciprocal identity is assumed.

This module closes only the ancestry/topology edge

    stored binary32 log_period -> exp(-log_period) result
      -> SeaStateAutoTuner binary32 clamp/store.

It does NOT prove how the binary32 log-period state was produced and it does NOT
prove target-libm correctness for either exp call.  Those remain deployment
supplies/obligations before the complete 600-step word can be promoted.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_frequency_binary32 as STORE
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as WPE

SOURCE=Path(__file__).resolve().parents[3]/'src/tuner/WavePeriodEstimator.h'
QUALIFICATION='OU3_ALT_WPE_BINARY32_GETTER_TO_TUNER_STORE_V1'


@dataclass(frozen=True)
class StoredLogPeriod:
    """Actual binary32 WPE log state paired with its exact-real shadow value."""
    shadow_log_period:F
    stored_log_period:F
    residual:F
    def __post_init__(self):
        s=F(self.shadow_log_period); q=F(self.stored_log_period); r=F(self.residual)
        if not B.is_binary32(q):
            raise ValueError('WPE deployment log_period must be an actual binary32 value')
        if r != q-s:
            raise ValueError('WPE log-period deployment residual detached from same shadow state')
        object.__setattr__(self,'shadow_log_period',s)
        object.__setattr__(self,'stored_log_period',q)
        object.__setattr__(self,'residual',r)


@dataclass(frozen=True)
class GetterResult:
    log:StoredLogPeriod
    period_argument:F
    frequency_argument:F
    period_result:F
    frequency_result:F
    def __post_init__(self):
        if not isinstance(self.log,StoredLogPeriod):
            raise TypeError('StoredLogPeriod required')
        pa,fa,pr,fr=map(F,(self.period_argument,self.frequency_argument,
                           self.period_result,self.frequency_result))
        if pa != self.log.stored_log_period or fa != -self.log.stored_log_period:
            raise ValueError('WPE exp getter arguments detached from stored binary32 log-period state')
        if not B.is_binary32(pr) or not B.is_binary32(fr) or pr<=0 or fr<=0:
            raise ValueError('WPE exp getter results must be positive finite binary32 witnesses')
        object.__setattr__(self,'period_argument',pa)
        object.__setattr__(self,'frequency_argument',fa)
        object.__setattr__(self,'period_result',pr)
        object.__setattr__(self,'frequency_result',fr)


def bind_log_state(shadow:WPE.WPEState, stored_log_period):
    if not isinstance(shadow,WPE.WPEState) or shadow.log_period is None:
        raise TypeError('WPE state with canonical log period required')
    q=F(stored_log_period)
    return StoredLogPeriod(F(shadow.log_period),q,q-F(shadow.log_period))


def getters(log:StoredLogPeriod, *, period_exp, frequency_exp):
    """Materialize the two distinct shipping std::exp calls on one log state."""
    if not isinstance(log,StoredLogPeriod): raise TypeError('StoredLogPeriod required')
    return GetterResult(log,log.stored_log_period,-log.stored_log_period,
                        F(period_exp),F(frequency_exp))


def store_frequency(result:GetterResult,min_hz,max_hz):
    """Feed exactly the shipping frequency getter result into tuner storage."""
    if not isinstance(result,GetterResult): raise TypeError('GetterResult required')
    return STORE.store(result.frequency_result,min_hz,max_hz)


def _source_shape_matches():
    s=SOURCE.read_text()
    return all(n in s for n in (
      'return std::isfinite(log_period_sec_) ? std::exp(log_period_sec_) : NAN;',
      'return std::isfinite(log_period_sec_) ? std::exp(-log_period_sec_) : NAN;',
      'float log_period_sec_ = NAN;'))


def readiness():
    st=STORE.readiness()
    return {
      'qualification':QUALIFICATION,
      'shipping_WPE_dual_exp_getter_source_shape_matches':_source_shape_matches(),
      'period_and_frequency_getters_share_same_stored_binary32_log_state':True,
      'period_and_frequency_libm_results_retained_separately':True,
      'binary32_period_frequency_bit_reciprocity_assumed':False,
      'WPE_frequency_getter_to_tuner_binary32_store_topology_closed': bool(
          st['frequency_clamp_and_store_exact_binary32'] and
          st['getFrequencyHz_is_identity_on_stored_binary32']),
      'WPE_binary32_log_period_production_closed':False,
      'WPE_period_exp_target_libm_correspondence_closed':False,
      'WPE_frequency_exp_target_libm_correspondence_closed':False,
      'upstream_WPE_binary32_frequency_production_closed':False,
      'complete_word_finite_identity':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
