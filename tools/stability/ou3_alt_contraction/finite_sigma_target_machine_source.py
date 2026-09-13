"""Source-owned binary32 sigma target from machine frontend + stillness state.

This is the strongest local entry to ``finite_tuner_sigma_binary32.target``.
It accepts one machine frontend successor and one sigma-relevant StillnessAdapter
successor from the SAME sample ordinal, then derives every non-libm sigma input:

* readiness from the carried machine DebiasedEMA weights;
* acceleration variance from the same machine statistics successor;
* band-noise sigma from the same machine adaptive-band successor;
* still predicate/time/attenuation from the same machine stillness successor.

Only the final ``sqrt(var_wave)`` machine result remains an explicit libm
witness.  Separate/FMA tuner histories are intentionally allowed to consume
different instances of this source relation: that is a conservative
cross-product over unknown compiler contraction choices and cannot exclude an
actual shipping execution.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_machine_frontend_sigma_source as F
from tools.stability.ou3_alt_contraction import finite_stillness_sigma_binary32 as S
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as T

QUALIFICATION='OU3_ALT_SIGMA_TARGET_MACHINE_SOURCE_V1'


@dataclass(frozen=True)
class Result:
    frontend:F.Result
    stillness:S.Step
    target:T.Target
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.frontend,F.Result) or not isinstance(self.stillness,S.Step) or not isinstance(self.target,T.Target):
            raise TypeError('machine frontend, stillness and sigma target relations required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong sigma machine-source qualification')
        if self.frontend.state.samples!=self.stillness.state.samples:
            raise ValueError('machine frontend and stillness successors come from different sample ordinals')
        if self.target.var_ready!=self.frontend.state.stats.var_ready:
            raise ValueError('sigma readiness detached from machine statistics weights')
        if self.target.accel_variance!=self.frontend.accel_variance:
            raise ValueError('sigma acceleration variance detached from machine frontend')
        if self.target.band_noise_sigma!=self.frontend.band_noise_sigma:
            raise ValueError('sigma band-noise input detached from machine frontend')
        if self.target.still!=self.stillness.state.is_still or self.target.still_time!=self.stillness.state.still_time:
            raise ValueError('sigma stillness branch detached from machine stillness successor')
        if self.target.attenuation!=self.stillness.attenuation:
            raise ValueError('sigma attenuation detached from same machine stillness exp')


def derive(frontend:F.Result,stillness:S.Step,cfg:D.DeploymentConfig,*,sqrt_result):
    if not isinstance(frontend,F.Result) or not isinstance(stillness,S.Step) or not isinstance(cfg,D.DeploymentConfig):
        raise TypeError('machine frontend, stillness step and DeploymentConfig required')
    if frontend.state.samples!=stillness.state.samples:
        raise ValueError('machine frontend and stillness must be same sample ordinal')
    still=stillness.state.is_still
    target=T.target(cfg,var_ready=frontend.state.stats.var_ready,
                    accel_variance=frontend.accel_variance,
                    band_noise_sigma=frontend.band_noise_sigma,
                    still=still,
                    still_time=stillness.state.still_time,
                    still_exp_result=stillness.exp_result if still else None,
                    sqrt_result=sqrt_result)
    return Result(frontend,stillness,target)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'same_sample_machine_frontend_and_stillness_required':True,
      'sigma_var_ready_derived_from_machine_stats_weights':True,
      'sigma_accel_variance_derived_from_machine_stats_successor':True,
      'sigma_band_noise_derived_from_machine_band_successor':True,
      'sigma_still_flag_time_and_attenuation_derived_from_machine_stillness_successor':True,
      'no_exact_real_sigma_source_inputs_substituted':True,
      'separate_and_FMA_tuner_tracks_may_consume_distinct_machine_frontend_sources':True,
      'final_sigma_sqrt_target_libm_correspondence_closed':False,
      'upstream_vertical_LP_and_WPE_machine_ancestry_closed':False,
      'target_compiler_contraction_membership_closed':False,
      'startup_frontend_machine_history_attached':False,
      'Live_600_step_machine_history_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
