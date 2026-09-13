"""Bridge tracker-free tuner-stillness projection into finite tuner candidates.

The adaptive band/statistics successor supplies frequency, variance and
propagated noise. The projected StillnessAdapter state supplies only shipping
``isStill``/``still_time`` plus the same variance attenuation. No dominant-
frequency tracker state is an operand.

Three surfaces are retained: legacy exact-rational component algebra, legacy
interval component algebra, and the deployment-faithful interval path whose
q_eff exact root is distinct from the cached binary32 libm value.
"""
from __future__ import annotations

from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as CAND
from tools.stability.ou3_alt_contraction import finite_tuner_candidate_interval as ICAND
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
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
                  cfg,*,sigma_wave_sqrt,dt,time,last_adapt_time,
                  ema:CAND.EmaWitness,bits:int=96):
    """General exact-real interval path, dispatching by explicit config type."""
    if not isinstance(previous,ICAND.IntervalTuneState):
        raise TypeError('interval tuner state required')
    sample=sample_from_projection(front,still,sigma_wave_sqrt=sigma_wave_sqrt)
    if isinstance(cfg,D.DeploymentConfig):
        return ICAND.step_deployment(previous,sample,cfg,dt=dt,time=time,
            last_adapt_time=last_adapt_time,ema=ema,bits=bits)
    if isinstance(cfg,CAND.CandidateConfig):
        return ICAND.step(previous,sample,cfg,dt=dt,time=time,
            last_adapt_time=last_adapt_time,ema=ema,bits=bits)
    raise TypeError('CandidateConfig or DeploymentConfig required')


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
      'deployment_qeff_exact_root_interval_projection_available':interval['deployment_exact_qeff_root_interval_consumed'],
      'deployment_cached_qeff_not_substituted_into_exact_real_projection':interval['deployment_cached_qeff_machine_value_not_used_as_exact_real_coefficient'],
      'actual_0p2Hz_prior_requires_no_fake_rational_spectral_witness':True,
      'sigma_wave_sqrt_transcendental_attached':False,
      # Backward-compatible fail-closed aggregate retained for existing theorem
      # ledgers/tests.  The more precise key below records that cached q_eff is
      # now also a distinct machine-libm obligation.
      'binary32_sqrt_pow_exp_correspondence_closed':False,
      'binary32_qeff_sqrt_pow_exp_correspondence_closed':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
