"""Strong binary32 SpectralMSE target path with source-produced q_eff cache.

``finite_tuner_spectral_binary32`` intentionally permits a stored qeff cache as
an input so its local arithmetic can be tested in isolation.  The theorem path
must be stronger: the cache consumed by the per-sample SpectralMSE target must
be exactly the cache produced from the carried acceleration-noise density by
``finite_tuner_qeff_cache_binary32``.

This wrapper performs that join.  It does not close any of the three target
libm calls; it only prevents a free qeff cache from being spliced into an
otherwise source-bound per-sample target.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_qeff_cache_binary32 as Q
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_binary32 as S

QUALIFICATION='OU3_ALT_SPECTRAL_CACHED_BINARY32_V1'


@dataclass(frozen=True)
class Result:
    cache:Q.CacheWitness
    spectral:S.Step
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.cache,Q.CacheWitness) or not isinstance(self.spectral,S.Step):
            raise TypeError('produced qeff cache and SpectralMSE step required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong cached SpectralMSE qualification')
        if self.spectral.rs_qeff_pow!=self.cache.result:
            raise ValueError('SpectralMSE target consumed a different qeff cache')


def step(cfg:C.CandidateConfig,cache:Q.CacheWitness,*,tau,sigma,pow_result,sqrt_result,bits=96):
    if not isinstance(cfg,C.CandidateConfig): raise TypeError('CandidateConfig required')
    Q.qualify_candidate_cache(cfg.qeff_pow,cache)
    out=S.step_from_config(cfg,tau=tau,sigma=sigma,pow_result=pow_result,
                           sqrt_result=sqrt_result,bits=bits)
    return Result(cache,out)


def readiness():
    q=Q.readiness(); s=S.readiness()
    return {
      'qualification':QUALIFICATION,
      'source_produced_qeff_cache_consumed_by_same_SpectralMSE_target':True,
      'free_qeff_cache_splice_forbidden':True,
      'ordinary_binary32_SpectralMSE_graph_materialized': bool(
          s['pseudo_cadence_binary32_mul_and_clamp_materialized'] and
          s['sigma_div_floor_tau2_u_binary32_graph_materialized'] and
          s['final_coefficient_multiply_multiply_divide_graph_materialized']),
      'all_three_SpectralMSE_libm_sites_have_same_argument_supply_coordinates': bool(
          q['cached_pow_error_interval_against_exact_fourteenth_root_exposed'] and
          s['pow_witness_error_interval_against_exact_root_of_same_machine_u_exposed'] and
          s['sqrt_witness_error_interval_against_exact_root_of_same_machine_TS_exposed']),
      'qeff_cache_target_libm_correspondence_closed':False,
      'per_sample_pow_target_libm_correspondence_closed':False,
      'per_sample_sqrt_target_libm_correspondence_closed':False,
      'binary32_SpectralMSE_target_correspondence_closed':False,
      'source_uniform_complete_startup_reachability_closed':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
