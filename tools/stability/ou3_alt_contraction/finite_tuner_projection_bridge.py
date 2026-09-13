"""Bridge tracker-free tuner-stillness projection into finite tuner candidates.

The adaptive band/statistics successor supplies frequency, variance and
propagated noise. The projected StillnessAdapter state supplies only shipping
``isStill``/``still_time`` plus the same variance attenuation. No dominant-
frequency tracker state is an operand.

Two theorem surfaces are exposed:

* ``step`` keeps the older exact-rational SpectralWitness path for algebra cells
  whose nonlinear roots happen to be rational;
* ``step_interval`` uses the exact rational root enclosure and interval R_S EMA
  for general shipping cells, including the real 0.2-Hz startup prior.

The second path does not approximate or choose an R_S representative; it is an
outer relation containing the exact real-arithmetic shipping candidate.
"""
from __future__ import annotations

from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as CAND
from tools.stability.ou3_alt_contraction import finite_tuner_candidate_interval as ICAND
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


def step_interval(previous:ICAND.IntervalTuneState,front:BAND.FrontendResult,still:STILL.Result,
                  cfg:CAND.CandidateConfig,*,sigma_wave_sqrt,dt,time,last_adapt_time,
                  ema:CAND.EmaWitness,bits:int=96):
    """Same shipping projection with general real SpectralMSE enclosure."""
    if not isinstance(previous,ICAND.IntervalTuneState):
        raise TypeError('interval tuner state required')
    sample=sample_from_projection(front,still,sigma_wave_sqrt=sigma_wave_sqrt)
    return ICAND.step(previous,sample,cfg,dt=dt,time=time,last_adapt_time=last_adapt_time,
                      ema=ema,bits=bits)


def readiness():
    interval=ICAND.readiness()
    return {
      'band_statistics_successor_attached':True,
      'tracker_free_stillness_projection_attached':True,
      'no_tracker_frequency_operand_in_tuner_candidate':True,
      'dominant_frequency_tracker_removed_from_OU_stability_interconnection':True,
      'general_SpectralMSE_interval_candidate_projection_available': bool(
          interval['SpectralMSE_irrational_target_enclosure_consumed'] and
          interval['RS_interval_propagated_through_literal_nonexpansive_EMA']),
      'actual_0p2Hz_prior_requires_no_fake_rational_spectral_witness':True,
      'sigma_wave_sqrt_transcendental_attached':False,
      'binary32_sqrt_pow_exp_correspondence_closed':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
