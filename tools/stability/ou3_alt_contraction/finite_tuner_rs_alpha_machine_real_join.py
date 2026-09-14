"""Same-source exact-real/binary32 join for the SpectralMSE R_S EMA alpha.

The older exact tuner accepted a free ``decay_RS`` witness.  That is useful for
component algebra but is not sufficient for a deployment theorem.  This module
roots BOTH sides at the same candidate tau, adaptation multiplier and dt:

Exact-real shadow:
    safe_tau = clamp(tau, .5, 6)
    horizon  = clamp(mult*safe_tau, max(dt,.05), 35)
    x        = dt/horizon
    exp(-x) lies in the rigorous exact-rational enclosure
    alpha    = 1-exp(-x), hence alpha lies in [1-exp_hi,1-exp_lo].

Machine path:
    ``finite_tuner_rs_alpha_binary32.Step`` carries the literal binary32 clamps,
    multiply/divide, explicit std::exp result and ``1-exp`` subtraction.

The result exposes machine-minus-exact alpha as an interval.  It does NOT claim
target-libm correctness or a source-uniform tight libm error bound.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_binary32 as M

QUALIFICATION='OU3_ALT_RS_ALPHA_MACHINE_REAL_JOIN_V1'
REAL_TIME_MIN=F(1,2); REAL_TIME_MAX=F(6)
REAL_HORIZON_MIN=F(1,20); REAL_HORIZON_MAX=F(35)


def clamp(x,lo,hi): return min(max(F(x),F(lo)),F(hi))


@dataclass(frozen=True)
class Join:
    cfg:D.DeploymentConfig
    tau_target:F
    dt:F
    machine:M.Step
    exact_safe_tau:F
    exact_requested_horizon:F
    exact_RS_sec:F
    exact_exp_argument:F
    exact_decay_lo:F
    exact_decay_hi:F
    exact_alpha_lo:F
    exact_alpha_hi:F
    machine_minus_exact_alpha_lo:F
    machine_minus_exact_alpha_hi:F
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.cfg,D.DeploymentConfig) or not isinstance(self.machine,M.Step):
            raise TypeError('deployment config and machine RS-alpha step required')
        for n in ('tau_target','dt','exact_safe_tau','exact_requested_horizon','exact_RS_sec','exact_exp_argument',
                  'exact_decay_lo','exact_decay_hi','exact_alpha_lo','exact_alpha_hi',
                  'machine_minus_exact_alpha_lo','machine_minus_exact_alpha_hi'):
            object.__setattr__(self,n,F(getattr(self,n)))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong RS alpha machine-real qualification')
        if self.machine.mult!=self.cfg.adapt_RS_mult or self.machine.tau_target!=self.tau_target or self.machine.dt!=self.dt:
            raise ValueError('machine alpha path detached from same config/tau/dt')
        if not 0<=self.exact_alpha_lo<=self.exact_alpha_hi<=1:
            raise ValueError('invalid exact alpha interval')
        if self.machine_minus_exact_alpha_lo!=self.machine.alpha-self.exact_alpha_hi or self.machine_minus_exact_alpha_hi!=self.machine.alpha-self.exact_alpha_lo:
            raise ValueError('alpha residual interval detached')


def join(cfg:D.DeploymentConfig,*,tau_target,dt,exp_decay,bits=None):
    if not isinstance(cfg,D.DeploymentConfig): raise TypeError('DeploymentConfig required')
    tau=F(tau_target); h=F(dt)
    if not B.is_binary32(tau) or not B.is_binary32(h):
        raise ValueError('joined deployment tau and dt must be actual binary32 values')
    machine=M.step(mult=cfg.adapt_RS_mult,tau_target=tau,dt=h,exp_decay=exp_decay,
                   slew_log=cfg.adapt_RS_slew_log)
    # Exact-real source shadow uses the mathematical values of the same compiled
    # binary32 constants/operands, but no intermediate float rounding.
    safe=clamp(tau,REAL_TIME_MIN,REAL_TIME_MAX)
    requested=F(cfg.adapt_RS_mult)*safe
    lo=min(max(h,REAL_HORIZON_MIN),REAL_HORIZON_MAX)
    rssec=clamp(requested,lo,REAL_HORIZON_MAX)
    x=h/rssec
    dlo,dhi,_,_=EXP.enclosure(x)
    alo=F(1)-dhi; ahi=F(1)-dlo
    return Join(cfg,tau,h,machine,safe,requested,rssec,x,dlo,dhi,alo,ahi,
                machine.alpha-ahi,machine.alpha-alo)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'exact_real_RS_horizon_rooted_at_same_config_tau_dt':True,
      'exact_real_exp_decay_rigorously_enclosed_at_same_horizon_argument':True,
      'exact_real_alpha_interval_not_free_witness':True,
      'machine_alpha_path_rooted_at_same_config_tau_dt':True,
      'machine_minus_exact_alpha_interval_exposed':True,
      'target_libm_exp_correspondence_closed':False,
      'source_uniform_tight_alpha_supply_bound_closed':False,
      'RS_EMA_interval_consumes_alpha_interval':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
