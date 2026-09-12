"""Binary32 tau-target and tau-EMA graph for shipping SeaStateAutoTuner.

This isolates the deployment-critical scalar that selects the Integrated-OU
Qaxis coefficient branch.  Shipping computes, in float,

    f_tune      = clamp(f_source, min_freq, max_freq)
    tau_raw     = tau_coeff * 0.5f / f_tune
    tau_target  = clamp(tau_raw, min_tau, max_tau)
    sea_time    = 0.5f / f_tune
    adapt_sec   = clamp(adapt_periods * clamp(sea_time,.5,6), max(dt,.05),35)
    alpha       = 1.0f - exp(-dt/adapt_sec)
    tau_applied += alpha * (tau_target - tau_applied)

All ordinary arithmetic and clamp decisions are exact RNE-binary32 here.  The
``std::exp`` result remains an explicit binary32 witness, constrained to the
same rounded argument by a rigorous real enclosure.  The final EMA supports
both separate mul/add and contracted FMA evaluation, because compiler FP
contraction has not yet been qualified for the shipping MCU build.

The strongest entry consumes ``StoredFrequency`` from
``finite_tuner_frequency_binary32``.  Consequently no theorem caller can
silently quantize an exact-real frequency at the tau edge: the value must
already be the actual float stored by ``SeaStateAutoTuner``.  Production of the
upstream WPE float remains open.

This module therefore closes the arithmetic *shape* and downstream stored
frequency/commit identity, not WPE/libm correctness or compiler contraction
selection.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as CAND
from tools.stability.ou3_alt_contraction import finite_tuner_commit as COMMIT
from tools.stability.ou3_alt_contraction import finite_tuner_frequency_binary32 as FREQ

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
HALF=B.rn32(F(1,2)); ONE=B.rn32(1)
SEA_MIN=B.rn32(F(1,2)); SEA_MAX=B.rn32(6)
HORIZON_MIN=B.rn32(F(1,20)); HORIZON_MAX=B.rn32(35)
MEKF_TAU_FLOOR=B.rn32(F(1,1000))
QUALIFICATION='OU3_ALT_TUNER_TAU_BINARY32_V2'


def clamp(x,lo,hi): return min(max(x,lo),hi)


@dataclass(frozen=True)
class TauStep:
    previous:F
    frequency:F
    tau_target:F
    sea_time:F
    adapt_sec:F
    exp_decay:F
    alpha:F
    next_separate:F
    next_fma:F
    def __post_init__(self):
        vals=[F(getattr(self,n)) for n in ('previous','frequency','tau_target','sea_time','adapt_sec','exp_decay','alpha','next_separate','next_fma')]
        for n,v in zip(('previous','frequency','tau_target','sea_time','adapt_sec','exp_decay','alpha','next_separate','next_fma'),vals): object.__setattr__(self,n,v)
        if not all(B.is_binary32(v) for v in vals): raise ValueError('TauStep stores deployed binary32 values only')
        if self.previous<=0 or self.frequency<=0 or self.tau_target<=0 or self.adapt_sec<=0: raise ValueError('positive tuner tau operands required')


def _floats_from_frequency(frequency,cfg:CAND.CandidateConfig,dt):
    if not isinstance(cfg,CAND.CandidateConfig): raise TypeError('finite tuner config required')
    if not cfg.clamp_enabled: raise ValueError('shipping default clamped tuner branch required')
    f0=F(frequency)
    if not B.is_binary32(f0): raise ValueError('tau edge requires already-stored binary32 tuner frequency')
    f=clamp(f0,B.rn32(cfg.min_freq),B.rn32(cfg.max_freq))
    coeff=B.rn32(cfg.tau_coeff)
    tau_raw=B.div(B.mul(coeff,HALF),f)
    tau_target=clamp(tau_raw,B.rn32(cfg.min_tau),B.rn32(cfg.max_tau))
    sea=B.div(HALF,f)
    periods=B.rn32(cfg.adapt_tau_sea_periods)
    if periods>0 and sea>0:
        safe=clamp(sea,SEA_MIN,SEA_MAX)
        requested=B.mul(periods,safe)
        dtf=B.rn32(dt)
        lo=max(dtf,HORIZON_MIN)
        adapt=clamp(requested,lo,HORIZON_MAX)
    else:
        adapt=B.rn32(cfg.adapt_tau_sec)
    return f,tau_target,sea,adapt


def _step_from_frequency(previous,frequency,cfg:CAND.CandidateConfig,*,dt,exp_decay):
    prev=F(previous)
    if not B.is_binary32(prev): raise ValueError('previous tau_applied is not an actual binary32 stored value')
    f,target,sea,adapt=_floats_from_frequency(frequency,cfg,dt)
    dtf=B.rn32(dt); x=B.div(dtf,adapt)
    e=F(exp_decay)
    if not B.is_binary32(e) or not 0<e<=1: raise ValueError('binary32 std::exp result witness required')
    lo,hi,_,_=EXP.enclosure(x)
    if not lo<=e<=hi: raise ValueError('tuner exp witness detached from SAME rounded -dt/adapt_sec argument')
    alpha=B.sub(ONE,e)
    sep=B.ema(prev,target,alpha,contracted=False)
    fused=B.ema(prev,target,alpha,contracted=True)
    return TauStep(prev,f,target,sea,adapt,e,alpha,sep,fused)


def step(previous,sample:CAND.WaveBandSample,cfg:CAND.CandidateConfig,*,dt,exp_decay):
    """Legacy/local entry; explicitly quantizes the exact-real sample frequency."""
    if not isinstance(sample,CAND.WaveBandSample): raise TypeError('finite tuner sample required')
    return _step_from_frequency(previous,B.rn32(sample.frequency_hz),cfg,dt=dt,exp_decay=exp_decay)


def step_from_stored_frequency(previous,stored:FREQ.StoredFrequency,cfg:CAND.CandidateConfig,*,dt,exp_decay):
    """Strong theorem entry consuming the actual SeaStateAutoTuner stored float."""
    if not isinstance(stored,FREQ.StoredFrequency): raise TypeError('StoredFrequency required')
    return _step_from_frequency(previous,FREQ.get_frequency_hz(stored),cfg,dt=dt,exp_decay=exp_decay)


def committed_tau(result:TauStep,*,contracted:bool):
    """Value passed to set_aw_time_constant at the next pending commit."""
    if not isinstance(result,TauStep): raise TypeError('TauStep required')
    if not isinstance(contracted,bool): raise TypeError('literal compiler contraction branch required')
    tau=result.next_fma if contracted else result.next_separate
    return max(MEKF_TAU_FLOOR,tau)


def _source_shape_matches():
    s=SOURCE.read_text()
    needles=(
      'float tau_raw = tau_coeff_ * 0.5f / f_tune;',
      'const float sea_time_sec = 0.5f / f_tune;',
      'const float alpha = 1.0f - std::exp(-dt / adapt_sec);',
      'tune_.tau_applied   += alpha    * (tau_t   - tune_.tau_applied);',
      'mekf_->set_aw_time_constant(tune_.tau_applied);',
    )
    return all(n in s for n in needles)


def readiness():
    b=B.readiness(); fs=FREQ.readiness()
    return {
      'qualification':QUALIFICATION,
      'shipping_tau_target_EMA_commit_source_shape_matches':_source_shape_matches(),
      'tau_target_float_clamp_graph_materialized':True,
      'dynamic_tau_EMA_horizon_float_graph_materialized':True,
      'tau_EMA_exp_result_bound_to_same_rounded_argument':True,
      'tau_EMA_separate_mul_add_result_materialized':b['separate_mul_add_EMA_shape_materialized'],
      'tau_EMA_contracted_fma_result_materialized':b['contracted_fma_EMA_shape_materialized'],
      'pending_commit_passes_stored_tau_directly_to_MEKF_setter':True,
      'source_tuner_frequency_binary32_store_to_tau_edge_closed': bool(
          fs['shipping_frequency_store_source_shape_matches'] and
          fs['frequency_clamp_and_store_exact_binary32'] and
          fs['getFrequencyHz_is_identity_on_stored_binary32']),
      'tuner_exp_libm_binary32_correspondence_closed':False,
      'shipping_compiler_FP_contraction_mode_qualified':False,
      'upstream_WPE_binary32_frequency_production_closed':False,
      'source_frontend_frequency_binary32_storage_correspondence_closed':True,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
