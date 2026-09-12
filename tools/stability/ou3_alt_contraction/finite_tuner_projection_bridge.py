"""Bridge tracker-free tuner-stillness projection into finite tuner candidate.

This module removes the dominant-frequency tracker from the OU tuning source
graph exactly.  The adaptive band/statistics successor supplies frequency,
variance and propagated noise.  The projected StillnessAdapter state supplies
only the shipping fields ``isStill`` and ``still_time`` plus the same variance
attenuation.  No tracker frequency/state is an operand.
"""
from __future__ import annotations

from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as CAND
from tools.stability.ou3_alt_contraction import finite_tuner_stillness_projection as STILL
from tools.stability.ou3_alt_contraction.finite_tuner_commit import TuneState


def sample_from_projection(front:BAND.FrontendResult,still:STILL.Result,*,sigma_wave_sqrt):
    if not isinstance(front,BAND.FrontendResult): raise TypeError('finite band/statistics successor required')
    if not isinstance(still,STILL.Result): raise TypeError('projected tuner stillness successor required')
    return CAND.WaveBandSample(
        front.current_tuner_frequency,front.variance_ready,front.accel_variance,
        front.band_noise_sigma,still.state.last_is_still,still.state.still_time,
        still.variance_attenuation,sigma_wave_sqrt)


def step(previous:TuneState,front:BAND.FrontendResult,still:STILL.Result,cfg:CAND.CandidateConfig,*,
         sigma_wave_sqrt,dt,time,last_adapt_time,
         spectral:CAND.SpectralWitness,ema:CAND.EmaWitness):
    sample=sample_from_projection(front,still,sigma_wave_sqrt=sigma_wave_sqrt)
    return CAND.step(previous,sample,cfg,dt=dt,time=time,last_adapt_time=last_adapt_time,
                     spectral=spectral,ema=ema)


def readiness():
    return {
      'band_statistics_successor_attached':True,
      'tracker_free_stillness_projection_attached':True,
      'no_tracker_frequency_operand_in_tuner_candidate':True,
      'dominant_frequency_tracker_removed_from_OU_stability_interconnection':True,
      'sigma_wave_sqrt_transcendental_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
