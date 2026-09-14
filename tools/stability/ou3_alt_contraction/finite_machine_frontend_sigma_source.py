"""Same-sample machine frontend source for the binary32 sigma target.

Shipping order is:

  adaptive band step -> tuner statistics step -> variance/noise readout ->
  stillness projection -> sigma target.

This module closes the first three arrows at machine level.  It carries one
persistent actual-machine AdaptiveWaveBandPass state and one persistent actual-
machine SeaStateAutoTuner moment state.  The tuner statistics input is required
to be the SAME newly stored machine band output.  The resulting acceleration
variance and band-noise sigma are therefore source-owned binary32 descendants,
not free inputs to the sigma target.

The external WPE frequency used by the statistics coefficient producer and the
band reference frequency used by the adaptive band are retained separately,
because shipping may source them from different carried states.  This module
does not identify them.  Stillness machine state, WPE/libm correspondence and
target compiler contraction correspondence remain open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_band_machine_ledger as BAND
from tools.stability.ou3_alt_contraction import finite_band_coefficients_binary32 as BC
from tools.stability.ou3_alt_contraction import finite_band_binary32_contraction as BR
from tools.stability.ou3_alt_contraction import finite_band_noise_floor_binary32 as NOISE
from tools.stability.ou3_alt_contraction import finite_stats_binary32_runtime as STATS

QUALIFICATION='OU3_ALT_MACHINE_FRONTEND_SIGMA_SOURCE_V1'
MAX_SAMPLES=min(BAND.MAX_SAMPLES,STATS.MAX_SAMPLES)


@dataclass(frozen=True)
class State:
    band:BAND.State
    stats:STATS.State
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.band,BAND.State) or not isinstance(self.stats,STATS.State):
            raise TypeError('persistent machine band and stats states required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine frontend qualification')
        if self.band.samples!=self.stats.samples:
            raise ValueError('machine band and statistics sample counts detached')
    @property
    def samples(self): return self.band.samples


@dataclass(frozen=True)
class Result:
    before:State
    state:State
    band:BAND.StepResult
    stats_envelope:STATS.Envelope
    noise:NOISE.Result
    accel_variance:F
    band_noise_sigma:F
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.before,State) or not isinstance(self.state,State): raise TypeError('machine frontend states required')
        if not isinstance(self.band,BAND.StepResult) or not isinstance(self.stats_envelope,STATS.Envelope) or not isinstance(self.noise,NOISE.Result):
            raise TypeError('machine frontend component results required')
        object.__setattr__(self,'accel_variance',F(self.accel_variance)); object.__setattr__(self,'band_noise_sigma',F(self.band_noise_sigma))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine frontend result qualification')
        if self.band.before!=self.before.band or self.band.state!=self.state.band:
            raise ValueError('band step detached from machine frontend predecessor/successor')
        if self.stats_envelope.before!=self.before.stats or not self.stats_envelope.accepts(self.state.stats):
            raise ValueError('statistics successor detached from same machine frontend predecessor')
        if self.stats_envelope.accel!=self.state.band.machine.band:
            raise ValueError('statistics input detached from SAME new machine band output')
        if self.noise.band!=self.state.band:
            raise ValueError('band-noise floor detached from SAME new machine band state')
        if self.accel_variance not in STATS.variance_outcomes(self.state.stats):
            raise ValueError('machine acceleration variance detached from carried statistics successor')
        if self.band_noise_sigma!=self.noise.noise_sigma:
            raise ValueError('machine band-noise sigma detached from carried band-noise result')
        if self.state.samples!=self.before.samples+1:
            raise ValueError('machine frontend did not advance exactly one sample')


def initial(): return State(BAND.initial(),STATS.State())


def step(state:State,*,band_coefficients:BC.Coefficients,band_input,
         band_successor:BR.State|None,stats_coefficients:STATS.Coefficients,
         stats_successor:STATS.State,bench_noise_sigma,noise_sqrt_gain=None,accel_variance=None):
    if not isinstance(state,State): raise TypeError('machine frontend State required')
    if state.samples>=MAX_SAMPLES: raise ValueError('machine frontend exceeded bounded startup+word horizon')
    b=BAND.step(state.band,band_coefficients,x=band_input,successor=band_successor)
    # SeaStateAutoTuner::update consumes the current adaptive-band output, not
    # the pre-step band state and not the raw vertical acceleration.
    accel=b.state.machine.band
    snext,senv=STATS.step(state.stats,stats_coefficients,accel=accel,successor=stats_successor)
    nxt=State(b.state,snext)
    noise=NOISE.evaluate(nxt.band,bench_sigma=bench_noise_sigma,sqrt_gain=noise_sqrt_gain)
    av=STATS.variance(snext) if accel_variance is None else F(accel_variance)
    return Result(state,nxt,b,senv,noise,av,noise.noise_sigma)


@dataclass(frozen=True)
class Pair:
    separate:State
    fma:State
    def __post_init__(self):
        if not isinstance(self.separate,State) or not isinstance(self.fma,State):
            raise TypeError('two persistent machine frontend histories required')
        if self.separate.samples!=self.fma.samples:
            raise ValueError('compiler frontend sample counts differ')
    @property
    def samples(self): return self.separate.samples


def initial_pair(): return Pair(initial(),initial())


def bind_step(previous:State,result:Result,*,frequency,band_cfg,stats_cfg,dt,bench_noise_sigma):
    """Re-execute the local graph with the carried WPE entry/config/history.

    Band corners read the previous statistics frequency; the current statistics
    read the WPE/prior input. Never use the current frequency for a ready band's
    lagged corner. The vertical machine input remains explicit and unqualified.
    """
    from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as W
    if not isinstance(previous,State) or not isinstance(result,Result) or not isinstance(frequency,W.StatisticsFrequencyResult):
        raise TypeError('persistent frontend, executed source result and two-clamp WPE frequency required')
    if result.before!=previous: raise ValueError('machine frontend history detached from predecessor')
    if frequency.stats_cfg!=stats_cfg: raise ValueError('machine statistics frequency config detached')
    external=frequency.external.stored.input_hz
    fref=previous.stats.frequency
    if fref is None or fref<=0: fref=external
    fref=min(max(fref,B.rn32(band_cfg.tune_freq_floor)),B.rn32(band_cfg.tune_freq_ceil))
    old=result.band.coefficients
    bc=BC.produce(band_cfg,f_ref=fref,dt=B.rn32(dt),exp_low=old.exp_low,exp_high=old.exp_high)
    if bc!=old: raise ValueError('band coefficients detached from lagged statistics/WPE/config')
    sc=STATS.coefficients(stats_cfg,frequency=external,dt=B.rn32(dt),exp_decay=result.stats_envelope.coefficients.exp_decay)
    if sc!=result.stats_envelope.coefficients or sc.frequency!=frequency.stats_stored.stored_hz:
        raise ValueError('statistics coefficients detached from current WPE entry/config')
    checked=step(previous,band_coefficients=bc,
        band_input=result.band.envelope.x if bc.active else B.rn32(0),
        band_successor=result.state.band.machine if bc.active else None,
        stats_coefficients=sc,stats_successor=result.state.stats,
        bench_noise_sigma=B.rn32(bench_noise_sigma),noise_sqrt_gain=result.noise.sqrt_gain,
        accel_variance=result.accel_variance)
    if checked!=result: raise ValueError('machine frontend result detached from re-executed source graph')
    return result


def require_sigma(result:Result,target):
    """Sigma's variance/readiness/noise must be outputs of this same sample."""
    from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as T
    if not isinstance(result,Result) or not isinstance(target,T.Target):
        raise TypeError('machine frontend result and sigma target required')
    if (target.var_ready,target.accel_variance,target.band_noise_sigma)!=(
            result.state.stats.var_ready,result.accel_variance,result.band_noise_sigma):
        raise ValueError('sigma variance/readiness/noise detached from persistent machine frontend')


def boundary_floors(pair:Pair,*,bench_noise_sigma,required,separate_sqrt_gain=None,fma_sqrt_gain=None):
    """Read the pre-sample band covariance without advancing either history."""
    if not isinstance(pair,Pair) or not isinstance(required,bool): raise TypeError('frontend pair and literal boundary selector required')
    if not required:
        if separate_sqrt_gain is not None or fma_sqrt_gain is not None:
            raise ValueError('unconsumed boundary requires no machine noise sqrt witnesses')
        return None,None
    bench=B.rn32(bench_noise_sigma)
    return (NOISE.evaluate(pair.separate.band,bench_sigma=bench,sqrt_gain=separate_sqrt_gain),
            NOISE.evaluate(pair.fma.band,bench_sigma=bench,sqrt_gain=fma_sqrt_gain))



def require_boundary(pair:Pair,floors,transaction):
    """Reject detached readouts even when a result is reconstructed directly."""
    from tools.stability.ou3_alt_contraction import finite_tuner_machine_boundary as T
    if not isinstance(pair,Pair) or not isinstance(transaction,T.Boundary):
        raise TypeError('persistent frontend pair and binary32 boundary required')
    if not isinstance(floors,tuple) or len(floors)!=2:
        raise TypeError('both retained boundary noise readouts required')
    if not transaction.consumed:
        if floors!=(None,None): raise ValueError('unconsumed boundary must retain no noise readouts')
        return
    for source,floor,commit in zip((pair.separate,pair.fma),floors,
            (transaction.arithmetic.separate,transaction.arithmetic.fma)):
        if not isinstance(floor,NOISE.Result) or floor.band!=source.band:
            raise ValueError('boundary noise readout detached from pre-sample machine band')
        if floor.noise_sigma!=commit.band_noise_floor_sigma:
            raise ValueError('machine commit noise operand detached from same band readout')
    if floors[0].bench_sigma!=floors[1].bench_sigma:
        raise ValueError('compiler noise readouts must share the configured bench sigma')

def readiness():
    return {
      'qualification':QUALIFICATION,
      'persistent_machine_band_and_stats_states_composed':True,
      'statistics_consumes_same_new_machine_band_output':True,
      'machine_accel_variance_derived_from_same_stats_successor':True,
      'machine_band_noise_floor_derived_from_same_band_successor':True,
      'band_and_stats_sample_counts_advance_together':True,
      'external_stats_frequency_not_falsely_identified_with_band_reference_frequency':True,
      'binary32_sigma_target_can_consume_machine_variance_and_noise':True,
      'stillness_machine_state_attached':False,
      'WPE_and_band_stats_libm_correspondence_closed':False,
      'target_compiler_contraction_membership_closed':False,
      'startup_frontend_machine_history_attached':False,
      'Live_600_step_machine_history_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
