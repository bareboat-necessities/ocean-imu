"""Local deployment join between exact-real and binary32 SpectralMSE targets.

For one SAME deployment configuration and one SAME candidate target whose tau
and sigma are actual binary32 values, compose:

* the exact-real SpectralMSE enclosure, including the exact q_eff^(1/14) root;
* the literal binary32 target graph, including the produced cached q_eff value
  and explicit per-call pow/sqrt witnesses;
* the shipping R_S clamp.

The result is an exact rational interval for

    R_S_machine - R_S_exact_real.

This interval deliberately includes every ordinary binary32 rounding difference
inside the target graph plus the three libm supplies (cached q_eff pow,
per-candidate pow, sqrt).  It is witness-local.  No source-uniform target-libm
bound or platform correspondence is claimed here.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_binary32 as M
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as R

QUALIFICATION='OU3_ALT_SPECTRAL_MACHINE_REAL_JOIN_V1'


@dataclass(frozen=True)
class Join:
    cfg:D.DeploymentConfig
    target:C.TargetState
    exact:R.Enclosure
    machine:M.Step
    machine_target_RS:F
    residual_lo:F
    residual_hi:F
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.cfg,D.DeploymentConfig) or not isinstance(self.target,C.TargetState):
            raise TypeError('deployment config and target required')
        if not isinstance(self.exact,R.Enclosure) or not isinstance(self.machine,M.Step):
            raise TypeError('exact and machine target relations required')
        for n in ('machine_target_RS','residual_lo','residual_hi'):
            object.__setattr__(self,n,F(getattr(self,n)))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine-real join qualification')
        if self.exact.target!=self.target: raise ValueError('exact target detached from joined candidate')
        if self.machine.tau!=self.target.tau_target or self.machine.sigma!=self.target.sigma_target:
            raise ValueError('machine target detached from same tau/sigma candidate')
        if self.machine.rs_qeff_pow!=self.cfg.qeff_cache.result:
            raise ValueError('machine graph detached from joined qeff cache')
        if self.residual_lo!=self.machine_target_RS-self.exact.target_RS_hi or self.residual_hi!=self.machine_target_RS-self.exact.target_RS_lo:
            raise ValueError('machine-real residual detached from exact target interval')
        if self.residual_lo>self.residual_hi: raise ValueError('invalid machine-real residual interval')


def join(cfg:D.DeploymentConfig,target:C.TargetState,*,pow_result,sqrt_result,bits:int=R.DEFAULT_BITS):
    if not isinstance(cfg,D.DeploymentConfig) or not isinstance(target,C.TargetState):
        raise TypeError('DeploymentConfig and TargetState required')
    tau=F(target.tau_target); sigma=F(target.sigma_target)
    if not B.is_binary32(tau): raise ValueError('joined tau target must be actual binary32')
    if not B.is_binary32(sigma): raise ValueError('joined sigma target must be actual binary32')
    exact=R.enclose_deployment(cfg,target,bits=bits)
    machine=M.step_from_deployment_config(cfg,tau=tau,sigma=sigma,
             pow_result=pow_result,sqrt_result=sqrt_result,bits=bits)
    machine_target=C.clamp(machine.raw_RS,cfg.min_RS,cfg.max_RS) if cfg.clamp_enabled else machine.raw_RS
    lo=machine_target-exact.target_RS_hi; hi=machine_target-exact.target_RS_lo
    return Join(cfg,target,exact,machine,machine_target,lo,hi)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'same_deployment_config_feeds_exact_real_and_binary32_targets':True,
      'same_binary32_tau_sigma_candidate_feeds_both_relations':True,
      'shipping_RS_clamp_composed_on_machine_target':True,
      'machine_minus_exact_real_RS_residual_interval_exposed':True,
      'residual_includes_qeff_cache_pow_and_per_candidate_pow_sqrt_supplies':True,
      'target_libm_correspondence_closed':False,
      'source_uniform_machine_real_RS_residual_bound_closed':False,
      'RS_commit_correspondence_closed':False,
      'source_uniform_complete_startup_reachability_closed':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
