"""Same-source exact-real / binary32 join for the shipping sigma target.

The exact frontend theorem produces one ``WaveBandSample``.  Deployment executes
the sigma target in binary32, so its variance/noise/stillness operands must not be
silently identified with that exact-real sample.  This module attaches a
``finite_tuner_sigma_binary32.Target`` to the SAME exact sample and exposes the
machine-minus-exact operand and target discrepancies as explicit supplies.

This is ancestry, not a smallness theorem: no source-uniform bounds are imposed
on the supplies here.  The later finite master must bound them from actual
frontend binary32/libm correspondence.  Consequently this join cannot be used
to prune admitted physical histories or authorize storage.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_common_alpha_qualification as A
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as M

QUALIFICATION='OU3_ALT_SIGMA_MACHINE_REAL_JOIN_V1'


def configs_match(candidate:C.CandidateConfig, deployment:D.DeploymentConfig):
    """Common EMA ancestry plus the sigma TARGET's own compiled constants."""
    return A._configs_match(candidate,deployment) and all(
        B.rn32(getattr(candidate,n))==F(getattr(deployment,n))
        for n in ('sigma_coeff','max_sigma'))


@dataclass(frozen=True)
class Supply:
    accel_variance:F
    band_noise_sigma:F
    still_time:F
    attenuation:F
    variance_wave:F
    sigma_target:F
    def __post_init__(self):
        for n in ('accel_variance','band_noise_sigma','still_time','attenuation','variance_wave','sigma_target'):
            object.__setattr__(self,n,F(getattr(self,n)))


@dataclass(frozen=True)
class Join:
    exact_sample:C.WaveBandSample
    exact_target:C.TargetState
    candidate_cfg:C.CandidateConfig
    deployment_cfg:D.DeploymentConfig
    machine:M.Target
    supply:Supply
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.exact_sample,C.WaveBandSample) or not isinstance(self.exact_target,C.TargetState):
            raise TypeError('exact WaveBandSample and target required')
        if not isinstance(self.candidate_cfg,C.CandidateConfig) or not isinstance(self.deployment_cfg,D.DeploymentConfig):
            raise TypeError('candidate and deployment configs required')
        if not isinstance(self.machine,M.Target) or not isinstance(self.supply,Supply):
            raise TypeError('machine sigma target and explicit supply required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong sigma machine-real qualification')
        if not configs_match(self.candidate_cfg,self.deployment_cfg):
            raise ValueError('sigma exact/machine configs disagree after binary32 compilation')
        if self.machine.cfg!=self.deployment_cfg:
            raise ValueError('machine sigma target detached from joined deployment config')
        if self.machine.var_ready!=self.exact_sample.variance_ready or self.machine.still!=self.exact_sample.still:
            raise ValueError('machine sigma branch detached from exact frontend branch')
        expected=C.targets(self.exact_sample,self.candidate_cfg)
        if self.exact_target!=expected:
            raise ValueError('exact sigma target detached from same frontend sample')
        s=self.supply
        if s.accel_variance!=self.machine.accel_variance-self.exact_sample.accel_variance:
            raise ValueError('accel-variance supply detached from same source operands')
        if s.band_noise_sigma!=self.machine.band_noise_sigma-self.exact_sample.band_noise_sigma:
            raise ValueError('band-noise supply detached from same source operands')
        if s.still_time!=self.machine.still_time-self.exact_sample.still_time:
            raise ValueError('still-time supply detached from same source operands')
        if s.attenuation!=self.machine.attenuation-self.exact_sample.still_attenuation:
            raise ValueError('attenuation supply detached from same source operands')
        if s.variance_wave!=self.machine.var_wave-self.exact_target.variance_wave:
            raise ValueError('wave-variance supply detached from same target pair')
        if s.sigma_target!=self.machine.sigma_target-self.exact_target.sigma_target:
            raise ValueError('sigma-target supply detached from same target pair')


def join(exact_sample:C.WaveBandSample,candidate_cfg:C.CandidateConfig,
         deployment_cfg:D.DeploymentConfig,machine:M.Target):
    if not isinstance(exact_sample,C.WaveBandSample) or not isinstance(machine,M.Target):
        raise TypeError('exact WaveBandSample and binary32 Target required')
    exact=C.targets(exact_sample,candidate_cfg)
    supply=Supply(
      machine.accel_variance-exact_sample.accel_variance,
      machine.band_noise_sigma-exact_sample.band_noise_sigma,
      machine.still_time-exact_sample.still_time,
      machine.attenuation-exact_sample.still_attenuation,
      machine.var_wave-exact.variance_wave,
      machine.sigma_target-exact.sigma_target)
    return Join(exact_sample,exact,candidate_cfg,deployment_cfg,machine,supply)


def readiness():
    return {
      'exact_frontend_sample_and_binary32_sigma_branch_joined':True,
      'candidate_and_deployment_common_scalars_compared_after_binary32_compilation':True,
      'sigma_target_scale_and_maximum_compared_after_binary32_compilation':True,
      'variance_ready_and_stillness_branch_identity_enforced':True,
      'machine_minus_exact_variance_noise_time_attenuation_supplies_exposed':True,
      'machine_minus_exact_wave_variance_and_sigma_target_supplies_exposed':True,
      'sigma_sqrt_and_still_exp_remain_bound_to_machine_operands':True,
      'upstream_frontend_binary32_correspondence_closed':False,
      'sigma_sqrt_and_still_exp_target_libm_correspondence_closed':False,
      'source_uniform_sigma_input_supply_bounds_closed':False,
      'startup_frontend_machine_TuneState_product_attached':False,
      'Live_600_step_machine_TuneState_product_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
