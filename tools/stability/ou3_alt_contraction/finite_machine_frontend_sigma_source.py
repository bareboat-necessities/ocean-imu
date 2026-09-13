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
        if self.accel_variance!=STATS.variance(self.state.stats):
            raise ValueError('machine acceleration variance detached from carried statistics successor')
        if self.band_noise_sigma!=self.noise.noise_sigma:
            raise ValueError('machine band-noise sigma detached from carried band-noise result')
        if self.state.samples!=self.before.samples+1:
            raise ValueError('machine frontend did not advance exactly one sample')


def initial(): return State(BAND.initial(),STATS.State())


def step(state:State,*,band_coefficients:BC.Coefficients,band_input,
         band_successor:BR.State|None,stats_coefficients:STATS.Coefficients,
         stats_successor:STATS.State,bench_noise_sigma,noise_sqrt_gain=None):
    if not isinstance(state,State): raise TypeError('machine frontend State required')
    if state.samples>=MAX_SAMPLES: raise ValueError('machine frontend exceeded bounded startup+word horizon')
    b=BAND.step(state.band,band_coefficients,x=band_input,successor=band_successor)
    # SeaStateAutoTuner::update consumes the current adaptive-band output, not
    # the pre-step band state and not the raw vertical acceleration.
    accel=b.state.machine.band
    snext,senv=STATS.step(state.stats,stats_coefficients,accel=accel,successor=stats_successor)
    nxt=State(b.state,snext)
    noise=NOISE.evaluate(nxt.band,bench_sigma=bench_noise_sigma,sqrt_gain=noise_sqrt_gain)
    return Result(state,nxt,b,senv,noise,STATS.variance(snext),noise.noise_sigma)


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
