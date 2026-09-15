"""Binary32 WPE getter topology and tuner-frequency handoff for ALT.

The exact-real WPE shadow legitimately uses one reciprocal period/frequency
pair. Shipping stores one binary32 ``log_period_sec_`` and evaluates two
distinct libm calls

    std::exp( log_period_sec_)
    std::exp(-log_period_sec_)

for period and frequency. Their returned floats are retained as independent
binary32 witnesses; no bit-level reciprocal identity is assumed.

The SeaState wrapper consumes the WPE frequency from the SAMPLE-ENTRY state.
It calls ``getFrequencyHz()`` before testing the usable latch and getter validity.
The selected value is that getter only when both tests pass; otherwise it uses
the literal 0.2f prior. The strong machine-entry relation retains eager getter
execution and invalid post-latch fallback separately from the exact shadow's selection. The legacy
component entry retains its explicit matched-branch preconditions.
The machine frequency need not equal the exact-real shadow frequency.  Their
post-clamp difference is retained explicitly as deployment supply instead of
being silently identified.

This module bounds the final clamped frequency discrepancy and retains the
represented topology/ancestry. It does NOT prove how the binary32
log-period state was produced and does NOT prove target-libm correctness.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_frequency_binary32 as STORE
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as WPE

SOURCE=Path(__file__).resolve().parents[3]/'src/tuner/WavePeriodEstimator.h'
WRAPPER=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
AUDITED_WRAPPER_SHA256='fabd03e9c3eb6069df107c7413ffb4b33fbdcd1ce06d3923b0c1ceeb3bcd7359'
PRIOR_EXACT=F(1,5); PRIOR=B.rn32(PRIOR_EXACT)
QUALIFICATION='OU3_ALT_WPE_BINARY32_GETTER_TO_TUNER_STORE_V3'


def clamp(x,lo,hi): return min(max(x,lo),hi)


@dataclass(frozen=True)
class FrequencySupplyBound:
    exact_interval:tuple
    machine_interval:tuple

    def __post_init__(self):
        for name in ('exact_interval','machine_interval'):
            interval=tuple(map(F,getattr(self,name)))
            if len(interval)!=2 or interval[0]>interval[1]:
                raise ValueError('ordered frequency interval required')
            object.__setattr__(self,name,interval)

    @property
    def residual_interval(self):
        return (self.machine_interval[0]-self.exact_interval[1],
                self.machine_interval[1]-self.exact_interval[0])

    def check(self,exact,machine):
        e,m=F(exact),F(machine)
        if not self.exact_interval[0]<=e<=self.exact_interval[1]:
            raise ValueError('exact frequency outside source-owned clamp range')
        if not self.machine_interval[0]<=m<=self.machine_interval[1]:
            raise ValueError('machine frequency outside source-owned clamp range')
        lo,hi=self.residual_interval
        if not lo<=m-e<=hi: raise AssertionError('ordered-interval subtraction failed')
        return m-e


def outer_supply_bound(*,exact_min_hz,exact_max_hz):
    """All executions reaching the final tuning-frequency assignment.

    Shipping first replaces a nonfinite/below-floor value by the floor, then
    caps above-ceiling values. Thus even a reset NaN or a retained statistics
    state lies in the outer interval. This statement does not require a WPE
    getter to succeed or a statistics update to be accepted.
    """
    lo,hi=F(exact_min_hz),F(exact_max_hz)
    if lo<=0 or hi<lo: raise ValueError('positive ordered tuning bounds required')
    return FrequencySupplyBound((lo,hi),(B.rn32(lo),B.rn32(hi)))


def final_tuning_clamp(value,*,min_hz,max_hz):
    """Literal outer clamp; None denotes any nonfinite input class."""
    lo,hi=F(min_hz),F(max_hz)
    if lo<=0 or hi<lo: raise ValueError('positive ordered tuning bounds required')
    f=lo if value is None or F(value)<lo else F(value)
    return min(f,hi)


def statistics_supply_bound(stats_cfg,*,exact_min_hz,exact_max_hz):
    """Uniform discrepancy bound after the two actual clamps.

    Clamp is monotone. The statistics clamp puts every accepted frequency in
    [s_min,s_max]; applying the outer clamp maps that whole interval to its
    clamped endpoints. This holds separately for real and rounded endpoints.
    Subtracting the two ranges bounds machine-minus-shadow for ALL inputs and
    branch choices, without asserting exp accuracy or branch agreement. The
    exact same-history residual is retained; this range never replaces it.
    """
    from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as R
    if not isinstance(stats_cfg,R.StatsConfig): raise TypeError('carried StatsConfig required')
    lo,hi=F(exact_min_hz),F(exact_max_hz)
    if lo<=0 or hi<lo: raise ValueError('positive ordered tuning bounds required')
    ml,mh=B.rn32(lo),B.rn32(hi)
    exact=tuple(clamp(F(s),lo,hi) for s in (stats_cfg.f_min,stats_cfg.f_max))
    machine=tuple(clamp(B.rn32(s),ml,mh) for s in (stats_cfg.f_min,stats_cfg.f_max))
    return FrequencySupplyBound(exact,machine)


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
class FrequencyExp:
    """The eager frequency getter; None denotes any nonfinite result."""
    argument:F
    result:F|None
    def __post_init__(self):
        a=F(self.argument)
        if not STORE.is_finite_input(a): raise ValueError('binary32 frequency-exp argument required')
        object.__setattr__(self,'argument',a)
        if self.result is not None:
            r=F(self.result)
            if not STORE.is_finite_input(r): raise ValueError('binary32 frequency-exp result required')
            object.__setattr__(self,'result',r)


@dataclass(frozen=True)
class MachineRead:
    log_period:F|None
    usable:bool
    exp:FrequencyExp|None
    def __post_init__(self):
        if type(self.usable) is not bool: raise TypeError('literal machine usable latch required')
        if self.log_period is None:
            if self.exp is not None: raise ValueError('nonfinite log executes no frequency exp')
        else:
            lp=F(self.log_period)
            if not STORE.is_finite_input(lp): raise ValueError('binary32 machine log required')
            object.__setattr__(self,'log_period',lp)
            if not isinstance(self.exp,FrequencyExp) or self.exp.argument!=-lp:
                raise ValueError('eager frequency exp detached from sample-entry machine log')
    @property
    def branch(self):
        good=self.exp is not None and self.exp.result is not None and self.exp.result>0
        return 'wpe' if self.usable and good else 'prior'
    @property
    def frequency(self): return self.exp.result if self.branch=='wpe' else PRIOR


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
    machine_read:MachineRead|None=None
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
        if self.machine_read is not None:
            if not isinstance(self.machine_read,MachineRead): raise TypeError('machine frequency read required')
            if self.getter is not None: raise ValueError('machine read owns the eager getter')
            if self.branch!=self.machine_read.branch or self.stored.input_hz!=self.machine_read.frequency:
                raise ValueError('tuner frequency detached from actual machine latch/getter decision')
            if not self.shadow.usable_period and ef!=PRIOR_EXACT:
                raise ValueError('exact preusable source must retain the exact prior')
            if self.shadow.usable_period and self.shadow.log_period is None:
                raise ValueError('exact usable source requires its own initialized log')
        elif self.branch=='prior':
            if self.shadow.usable_period: raise ValueError('usable WPE state cannot take fixed-prior tuner branch')
            if self.getter is not None: raise ValueError('fixed-prior tuner branch consumes no WPE exp getter result')
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


def machine_frequency(shadow,*,log_period,usable,getter,min_hz,max_hz,shadow_frequency=None):
    """Independent exact and machine choices, with no log-initialization equality."""
    if not isinstance(shadow,WPE.WPEState): raise TypeError('same exact WPE entry required')
    if shadow.usable_period:
        if shadow_frequency is None: raise TypeError('exact usable branch requires exact frequency')
        exact=F(shadow_frequency)
    else:
        if shadow_frequency is not None: raise ValueError('exact prior consumes no frequency witness')
        exact=PRIOR_EXACT
    if isinstance(getter,GetterResult):
        if getter.log.stored_log_period!=log_period or getter.log.shadow_log_period!=shadow.log_period:
            raise ValueError('legacy getter detached from machine/exact entry logs')
        getter=FrequencyExp(getter.frequency_argument,getter.frequency_result)
    read=MachineRead(log_period,usable,getter)
    stored=STORE.store(read.frequency,min_hz,max_hz)
    ec=clamp(exact,F(min_hz),F(max_hz))
    return TunerFrequencyResult(shadow,read.branch,None,exact,ec,stored,stored.stored_hz-ec,read)


def machine_frequencies(shadow,entry,*,logs,separate_getter,fma_getter,shadow_frequency,
                        stats_cfg,exact_min_hz,exact_max_hz):
    """Projection of the persistent full WPE predecessor into the lower tuner.

    The strong startup/Live wrapper supplies entry directly from its carried
    state. The lower log ledger must be identical. Current-sample production is
    checked against the full moment update by that same wrapper on return.
    """
    from tools.stability.ou3_alt_contraction import finite_wpe_machine_binary32 as M
    if not isinstance(entry,M.State) or entry.logs!=logs:
        raise ValueError('full machine WPE entry detached from lower carried log ledger')
    lo,hi=F(exact_min_hz),F(exact_max_hz)
    def one(track,usable,getter):
        q=machine_frequency(shadow,log_period=track.log_period,usable=usable,getter=getter,
            min_hz=B.rn32(lo),max_hz=B.rn32(hi),shadow_frequency=shadow_frequency)
        return through_statistics(q,stats_cfg,exact_min_hz=lo,exact_max_hz=hi)
    return (one(logs.separate,entry.separate_usable,separate_getter),
            one(logs.fma,entry.fma_usable,fma_getter))


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


@dataclass(frozen=True)
class StatisticsFrequencyResult:
    """WPE/prior -> statistics store -> outer tuning clamp (two operations)."""
    external:TunerFrequencyResult
    stats_cfg:object
    exact_tune_bounds:tuple
    stats_stored:STORE.StoredFrequency
    stored:STORE.StoredFrequency
    exact_clamped_frequency:F
    machine_minus_shadow:F
    def __post_init__(self):
        from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as R
        if not isinstance(self.external,TunerFrequencyResult) or not isinstance(self.stats_cfg,R.StatsConfig):
            raise TypeError('WPE frequency relation and carried StatsConfig required')
        q=self.external; c=self.stats_cfg
        lo,hi=map(F,self.exact_tune_bounds)
        object.__setattr__(self,'exact_tune_bounds',(lo,hi))
        if lo<=0 or hi<lo or (B.rn32(lo),B.rn32(hi))!=(q.stored.min_hz,q.stored.max_hz):
            raise ValueError('exact and machine outer tuning bounds detached')
        st=STORE.store(q.stored.input_hz,B.rn32(c.f_min),B.rn32(c.f_max))
        final=STORE.store(st.stored_hz,q.stored.min_hz,q.stored.max_hz)
        exact=clamp(clamp(q.exact_shadow_frequency,c.f_min,c.f_max),lo,hi)
        if self.stats_stored!=st or self.stored!=final:
            raise ValueError('frequency path detached from ordered statistics and tuning clamps')
        if F(self.exact_clamped_frequency)!=exact or F(self.machine_minus_shadow)!=final.stored_hz-exact:
            raise ValueError('two-clamp frequency supply detached from same shadow/machine inputs')
        self.supply_bound.check(exact,final.stored_hz)

    @property
    def supply_bound(self):
        return statistics_supply_bound(self.stats_cfg,exact_min_hz=self.exact_tune_bounds[0],
            exact_max_hz=self.exact_tune_bounds[1])


def through_statistics(external:TunerFrequencyResult,stats_cfg,*,exact_min_hz,exact_max_hz):
    """Do not collapse the tuner statistics bounds into the outer tune bounds."""
    from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as R
    if not isinstance(external,TunerFrequencyResult) or not isinstance(stats_cfg,R.StatsConfig):
        raise TypeError('WPE frequency relation and carried StatsConfig required')
    st=STORE.store(external.stored.input_hz,B.rn32(stats_cfg.f_min),B.rn32(stats_cfg.f_max))
    out=STORE.store(st.stored_hz,external.stored.min_hz,external.stored.max_hz)
    lo,hi=F(exact_min_hz),F(exact_max_hz)
    exact=clamp(clamp(external.exact_shadow_frequency,stats_cfg.f_min,stats_cfg.f_max),lo,hi)
    return StatisticsFrequencyResult(external,stats_cfg,(lo,hi),st,out,exact,out.stored_hz-exact)


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
      'preusable_WPE_uses_literal_binary32_0p2_prior_without_consuming_getter_result':True,
      'eager_frequency_getter_and_invalid_post_latch_fallback_composed':True,
      'machine_frequency_branch_independent_of_exact_usable_latch':True,
      'usable_WPE_frequency_getter_bound_to_sample_entry_log_state':True,
      'machine_minus_exact_tuner_frequency_supply_exposed':True,
      'statistics_store_and_outer_tuning_clamps_retained_in_order':True,
      'exact_outer_bounds_not_identified_with_compiled_binary32_bounds':True,
      'WPE_frequency_getter_to_tuner_binary32_store_topology_closed': bool(
          st['frequency_clamp_and_store_exact_binary32'] and st['getFrequencyHz_is_identity_on_stored_binary32']),
      'WPE_binary32_log_period_production_closed':False,
      'WPE_period_exp_target_libm_correspondence_closed':False,
      'WPE_frequency_exp_target_libm_correspondence_closed':False,
      'source_uniform_WPE_frequency_supply_bound_closed': bool(
          _source_shape_matches() and st['shipping_frequency_store_source_shape_matches'] and
          sha256(WRAPPER.read_bytes()).hexdigest()==AUDITED_WRAPPER_SHA256),
      'frequency_supply_bound_requires_no_libm_accuracy_or_branch_agreement':True,
      'outer_frequency_bound_includes_nonfinite_fallback_and_retained_statistics':True,
      'frequency_supply_bound_proves_execution_totality':False,
      'upstream_WPE_binary32_frequency_production_closed':False,
      'complete_word_finite_identity':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
