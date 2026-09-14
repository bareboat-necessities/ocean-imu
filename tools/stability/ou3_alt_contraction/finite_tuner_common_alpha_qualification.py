"""Qualify the common shipping tau/sigma EMA alpha for deployment use.

``adapt_mekf`` computes ONE ``alpha = 1-exp(-dt/adapt_sec)`` and uses it for
both ``tau_applied`` and ``sigma_applied``.  A raw ``TauStep`` is convenient but
constructible, so the sigma theorem must not trust its alpha field by type alone.

This module re-derives the frequency-dependent horizon from the declared tuner
configuration, verifies the TauStep's target/sea-time/horizon/exp/alpha and BOTH
compiler successors, and checks that the deployment SpectralMSE config carries
the same scalars that actually determine this common alpha.  The legacy
exact-rational ``CandidateConfig`` is compared to the deployment configuration
only after each shared scalar is compiled to binary32; exact-rational equality
to a compiled float is not a shipping invariant.

Important separation: ``sigma_coeff`` and ``max_sigma`` determine the sigma
TARGET, not the common tau/sigma EMA alpha.  They are therefore intentionally
NOT prerequisites here.  Their same-source consistency is enforced by the
sigma machine/real join and the coherent whole-TuneState product.  Mixing those
obligations here would reject valid alpha ancestry for an unrelated target
configuration difference and obscure which theorem actually failed.

Target ``std::exp`` platform correspondence remains open; this proves ancestry
and arithmetic shape, not libm implementation correctness.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_tau_binary32 as T

QUALIFICATION='OU3_ALT_COMMON_ALPHA_QUALIFICATION_V3'


@dataclass(frozen=True)
class Qualified:
    step:T.TauStep
    candidate_cfg:C.CandidateConfig
    deployment_cfg:D.DeploymentConfig
    dt:F
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.step,T.TauStep) or not isinstance(self.candidate_cfg,C.CandidateConfig) or not isinstance(self.deployment_cfg,D.DeploymentConfig):
            raise TypeError('TauStep, CandidateConfig and DeploymentConfig required')
        object.__setattr__(self,'dt',F(self.dt))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong common-alpha qualification')
        if not B.is_binary32(self.dt) or self.dt<=0: raise ValueError('common-alpha dt must be positive binary32')

    @property
    def alpha(self): return self.step.alpha


def _configs_match(c:C.CandidateConfig,d:D.DeploymentConfig):
    # Exactly the shipping fields that determine the frequency clamp, tau target,
    # sea-period horizon and hence the ONE alpha shared by tau and sigma EMAs.
    # Sigma-target scale/clamp fields are checked in SIGJOIN, not here.
    names=('min_freq','max_freq','tau_coeff','min_tau','max_tau',
           'adapt_tau_sec','adapt_tau_sea_periods')
    return all(B.rn32(getattr(c,n))==F(getattr(d,n)) for n in names) and c.clamp_enabled==d.clamp_enabled


def qualify(step:T.TauStep,candidate_cfg:C.CandidateConfig,deployment_cfg:D.DeploymentConfig,*,dt):
    if not isinstance(step,T.TauStep) or not isinstance(candidate_cfg,C.CandidateConfig) or not isinstance(deployment_cfg,D.DeploymentConfig):
        raise TypeError('TauStep, CandidateConfig and DeploymentConfig required')
    h=F(dt)
    if not B.is_binary32(h) or h<=0: raise ValueError('common-alpha dt must be positive binary32')
    if not _configs_match(candidate_cfg,deployment_cfg):
        raise ValueError('tau and deployment configs disagree after binary32 compilation of common-alpha scalars')
    f,target,sea,adapt=T._floats_from_frequency(step.frequency,candidate_cfg,h)
    if (step.frequency,step.tau_target,step.sea_time,step.adapt_sec)!=(f,target,sea,adapt):
        raise ValueError('TauStep detached from declared common-alpha frequency/config horizon')
    x=B.div(h,adapt); elo,ehi,_,_=EXP.enclosure(x)
    if not elo<=step.exp_decay<=ehi:
        raise ValueError('TauStep exp witness detached from same rounded common-alpha argument')
    alpha=B.sub(B.rn32(1),step.exp_decay)
    if step.alpha!=alpha: raise ValueError('TauStep alpha detached from 1-exp common coefficient')
    if step.next_separate!=B.ema(step.previous,step.tau_target,alpha,contracted=False):
        raise ValueError('TauStep separate successor detached from common alpha')
    if step.next_fma!=B.ema(step.previous,step.tau_target,alpha,contracted=True):
        raise ValueError('TauStep FMA successor detached from common alpha')
    return Qualified(step,candidate_cfg,deployment_cfg,h)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'candidate_common_alpha_scalars_compiled_to_binary32_before_deployment_comparison':True,
      'sigma_target_scale_and_clamp_deliberately_deferred_to_sigma_join':True,
      'tau_step_rederived_from_declared_frequency_config_and_dt':True,
      'candidate_and_deployment_configs_share_common_tau_sigma_adaptation_scalars':True,
      'same_exp_decay_and_alpha_drive_verified_tau_successors':True,
      'qualified_alpha_available_for_sigma_machine_recurrence':True,
      'target_libm_exp_correspondence_closed':False,
      'sigma_machine_ledger_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
