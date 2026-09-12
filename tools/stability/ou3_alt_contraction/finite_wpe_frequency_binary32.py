"""Binary32 WPE getter topology and tuner-frequency handoff for ALT.

The exact-real WPE shadow legitimately uses one reciprocal period/frequency
pair. Shipping stores one binary32 ``log_period_sec_`` and evaluates two
distinct libm calls

    std::exp( log_period_sec_)
    std::exp(-log_period_sec_)

for period and frequency. Their returned floats are retained as independent
binary32 witnesses; no bit-level reciprocal identity is assumed.

The SeaState wrapper consumes the WPE frequency from the SAMPLE-ENTRY state.
Until ``hasUsablePeriod()`` latches it uses the literal 0.2f prior; afterwards
it calls ``getFrequencyHz()`` and passes that float to ``SeaStateAutoTuner``.
The machine frequency need not equal the exact-real shadow frequency.  Their
post-clamp difference is retained explicitly as deployment supply instead of
being silently identified.

This module closes topology/ancestry only.  It does NOT prove how the binary32
log-period state was produced and does NOT prove target-libm correctness.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_frequency_binary32 as STORE
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as WPE

SOURCE=Path(__file__).resolve().parents[3]/'src/tuner/WavePeriodEstimator.h'
WRAPPER=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
PRIOR_EXACT=F(1,5); PRIOR=B.rn32(PRIOR_EXACT)
QUALIFICATION='OU3_ALT_WPE_BINARY32_GETTER_TO_TUNER_STORE_V3'


def clamp(x,lo,hi): return min(max(x,lo),hi)


@dataclass(frozen=True)
class StoredLogPeriod:
    shadow_log_period:F
    stored_log_period:F
    residual:F
    def __post_init__(self):
        s=F(self.shadow_log_period); q=F(self.stored_log_period); r=F(self.residual)
        if not B.is_binary32(q): raise ValueError('WPE deployment log_period must be an actual binary32 value')
        if r != q-s: raise ValueError('WPE log-period deployment residual detached from same shadow state')
        object.__setattr__(self,'shadow_log_period',s); object.__setattr__(self,'stored_log_period',q); object.__setattr__(self,'residual',r)


@dataclass(frozen=True)
class GetterResult:
    log:StoredLogPeriod
    period_argument:F
    frequency_argument:F
    period_result:F
    frequency_result:F
    def __post_init__(self):
        if not isinstance(self.log,StoredLogPeriod): raise TypeError('StoredLogPeriod required')
        pa,fa,pr,fr=map(F,(self.period_argument,self.frequency_argument,self.period_result,self.frequency_result))
        if pa != self.log.stored_log_period or fa != -self.log.stored_log_period:
            raise ValueError('WPE exp getter arguments detached from stored binary32 log-period state')
        if not B.is_binary32(pr) or not B.is_binary32(fr) or pr<=0 or fr<=0:
            raise ValueError('WPE exp getter results must be positive finite binary32 witnesses')
        object.__setattr__(self,'period_argument',pa); object.__setattr__(self,'frequency_argument',fa)
        object.__setattr__(self,'period_result',pr); object.__setattr__(self,'frequency_result',fr)


@dataclass(frozen=True)
class TunerFrequencyResult:
    """Exact shadow and deployed frequency consumed by one tuner sample."""
    shadow:WPE.WPEState
    branch:str
    getter:GetterResult|None
    exact_shadow_frequency:F
    exact_clamped_frequency:F
    stored:STORE.StoredFrequency
    machine_minus_shadow:F
    def __post_init__(self):
        if not isinstance(self.shadow,WPE.WPEState) or not isinstance(self.stored,STORE.StoredFrequency):
            raise TypeError('preupdate WPE state and StoredFrequency required')
        ef,ec,r=map(F,(self.exact_shadow_frequency,self.exact_clamped_frequency,self.machine_minus_shadow))
        if ef<=0: raise ValueError('positive exact shadow tuner frequency required')
        expected=clamp(ef,F(self.stored.min_hz),F(self.stored.max_hz))
        if ec!=expected: raise ValueError('exact tuner frequency clamp detached from same shadow input')
        if r!=F(self.stored.stored_hz)-ec:
            raise ValueError('machine-minus-shadow tuner frequency supply detached')
        if self.branch not in ('prior','wpe'): raise ValueError('unknown tuner-frequency source branch')
        if self.branch=='prior':
            if self.shadow.usable_period: raise ValueError('usable WPE state cannot take fixed-prior tuner branch')
            if self.getter is not None: raise ValueError('fixed-prior tuner branch consumes no WPE exp getter')
            if ef!=PRIOR_EXACT or self.stored.input_hz!=PRIOR:
                raise ValueError('fixed-prior tuner branch detached from shipping 0.2f prior')
        else:
            if not self.shadow.usable_period or self.shadow.log_period is None:
                raise ValueError('WPE tuner branch requires usable canonical log-period state')
            if not isinstance(self.getter,GetterResult): raise TypeError('WPE tuner branch requires binary32 getter result')
            if self.getter.log.shadow_log_period != F(self.shadow.log_period):
                raise ValueError('WPE getter detached from sample-entry canonical log-period shadow')
            if self.stored.input_hz != self.getter.frequency_result:
                raise ValueError('tuner store detached from same WPE frequency getter result')
        object.__setattr__(self,'exact_shadow_frequency',ef); object.__setattr__(self,'exact_clamped_frequency',ec)
        object.__setattr__(self,'machine_minus_shadow',r)


def bind_log_state(shadow:WPE.WPEState, stored_log_period):
    if not isinstance(shadow,WPE.WPEState) or shadow.log_period is None:
        raise TypeError('WPE state with canonical log period required')
    q=F(stored_log_period)
    return StoredLogPeriod(F(shadow.log_period),q,q-F(shadow.log_period))


def getters(log:StoredLogPeriod, *, period_exp, frequency_exp):
    if not isinstance(log,StoredLogPeriod): raise TypeError('StoredLogPeriod required')
    return GetterResult(log,log.stored_log_period,-log.stored_log_period,F(period_exp),F(frequency_exp))


def store_frequency(result:GetterResult,min_hz,max_hz):
    if not isinstance(result,GetterResult): raise TypeError('GetterResult required')
    return STORE.store(result.frequency_result,min_hz,max_hz)


def tuner_frequency(shadow:WPE.WPEState, *, min_hz, max_hz, getter:GetterResult|None=None,
                    shadow_frequency=None):
    """Literal sample-entry branch with explicit machine-minus-exact frequency supply."""
    if not isinstance(shadow,WPE.WPEState): raise TypeError('preupdate WPE state required')
    lo,hi=F(min_hz),F(max_hz)
    if shadow.usable_period:
        if getter is None: raise TypeError('usable WPE tuner branch requires getter witness')
        if shadow_frequency is None: raise TypeError('usable WPE tuner branch requires exact shadow frequency')
        if shadow.log_period is None or getter.log.shadow_log_period!=F(shadow.log_period):
            raise ValueError('WPE getter detached from preupdate canonical log-period state')
        exact=F(shadow_frequency)
        stored=STORE.store(getter.frequency_result,min_hz,max_hz)
        ec=clamp(exact,lo,hi)
        return TunerFrequencyResult(shadow,'wpe',getter,exact,ec,stored,F(stored.stored_hz)-ec)
    if getter is not None or shadow_frequency is not None:
        raise ValueError('pre-usable WPE branch consumes no getter/shadow-frequency witness')
    stored=STORE.store(PRIOR,min_hz,max_hz); ec=clamp(PRIOR_EXACT,lo,hi)
    return TunerFrequencyResult(shadow,'prior',None,PRIOR_EXACT,ec,stored,F(stored.stored_hz)-ec)


def _source_shape_matches():
    s=SOURCE.read_text(); w=WRAPPER.read_text()
    return all(n in s for n in (
      'return std::isfinite(log_period_sec_) ? std::exp(log_period_sec_) : NAN;',
      'return std::isfinite(log_period_sec_) ? std::exp(-log_period_sec_) : NAN;',
      'float log_period_sec_ = NAN;')) and all(n in w for n in (
      'const float wave_hz = wave_period_.getFrequencyHz();','if (wave_period_.hasUsablePeriod() &&',
      'return wave_hz;','return tune_freq_prior_hz_;','float tune_freq_prior_hz_     = TUNE_FREQ_PRIOR_HZ;',
      'constexpr float TUNE_FREQ_PRIOR_HZ = 0.2f;'))


def readiness():
    st=STORE.readiness()
    return {
      'qualification':QUALIFICATION,
      'shipping_WPE_dual_exp_and_tuner_source_shape_matches':_source_shape_matches(),
      'period_and_frequency_getters_share_same_stored_binary32_log_state':True,
      'period_and_frequency_libm_results_retained_separately':True,
      'binary32_period_frequency_bit_reciprocity_assumed':False,
      'preusable_WPE_uses_literal_binary32_0p2_prior_without_exp':True,
      'usable_WPE_frequency_getter_bound_to_sample_entry_log_state':True,
      'machine_minus_exact_tuner_frequency_supply_exposed':True,
      'WPE_frequency_getter_to_tuner_binary32_store_topology_closed': bool(
          st['frequency_clamp_and_store_exact_binary32'] and st['getFrequencyHz_is_identity_on_stored_binary32']),
      'WPE_binary32_log_period_production_closed':False,
      'WPE_period_exp_target_libm_correspondence_closed':False,
      'WPE_frequency_exp_target_libm_correspondence_closed':False,
      'source_uniform_WPE_frequency_supply_bound_closed':False,
      'upstream_WPE_binary32_frequency_production_closed':False,
      'complete_word_finite_identity':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
