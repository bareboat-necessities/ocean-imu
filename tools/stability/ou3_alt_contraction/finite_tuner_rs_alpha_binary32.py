"""Binary32 production graph for shipping SpectralMSE ``alpha_RS``.

For the deployed default ``adapt_RS_slew_log_ == 0`` the optional logarithmic
slew-shortening branch in ``adaptiveSmoothingHorizonSec`` is dead.  Shipping
therefore evaluates, in float,

    safe_tau = clamp(tau_t, 0.5f, 6.0f)
    requested = adapt_RS_mult_ * safe_tau
    RS_sec = clamp(requested, max(dt, 0.05f), 35.0f)
    decay = std::exp(-dt / RS_sec)
    alpha_RS = 1.0f - decay

The ordinary binary32 operations and both clamps are materialized exactly here.
``std::exp`` remains an explicit positive binary32 witness tied to the SAME
rounded ``dt / RS_sec`` argument through the rigorous exact-real exp enclosure.
No platform-libm correctness or source-uniform error bound is assumed.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_tau_binary32 as TAU

COMMON=Path(__file__).resolve().parents[3]/'src/kalman_common/SeaStateFusionFilterCommon.h'
FILTER=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
ONE=B.rn32(1)
TIME_MIN=B.rn32(F(1,2)); TIME_MAX=B.rn32(6)
HORIZON_MIN=B.rn32(F(1,20)); HORIZON_MAX=B.rn32(35)
QUALIFICATION='OU3_ALT_RS_ALPHA_BINARY32_V1'


def clamp(x,lo,hi): return min(max(x,lo),hi)

def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(f'{name} must be actual binary32')
    return q


@dataclass(frozen=True)
class Step:
    mult:F
    tau_target:F
    dt:F
    safe_tau:F
    requested_horizon:F
    lower_horizon:F
    RS_sec:F
    exp_argument:F
    exp_decay:F
    alpha:F
    qualification:str=QUALIFICATION
    exp_profile:str='legacy-enclosure'
    def __post_init__(self):
        for n in ('mult','tau_target','dt','safe_tau','requested_horizon','lower_horizon',
                  'RS_sec','exp_argument','exp_decay','alpha'):
            object.__setattr__(self,n,F(getattr(self,n)))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong RS alpha qualification')
        if not all(B.is_binary32(getattr(self,n)) for n in
                   ('mult','tau_target','dt','safe_tau','requested_horizon','lower_horizon','RS_sec','exp_argument','exp_decay','alpha')):
            raise ValueError('RS alpha graph stores binary32 values only')
        if self.mult<=0 or self.tau_target<=0 or self.dt<=0 or self.RS_sec<=0:
            raise ValueError('positive RS smoothing operands required')
        if self.safe_tau!=clamp(self.tau_target,TIME_MIN,TIME_MAX):
            raise ValueError('safe tau detached from shipping dynamic-time clamp')
        if self.requested_horizon!=B.mul(self.mult,self.safe_tau):
            raise ValueError('requested horizon detached from source-order multiply')
        if self.lower_horizon!=min(max(self.dt,HORIZON_MIN),HORIZON_MAX):
            raise ValueError('lower horizon detached from dt guard')
        if self.RS_sec!=clamp(self.requested_horizon,self.lower_horizon,HORIZON_MAX):
            raise ValueError('RS_sec detached from shipping final horizon clamp')
        if self.exp_argument!=B.div(self.dt,self.RS_sec):
            raise ValueError('exp argument detached from rounded dt/RS_sec')
        if self.alpha!=B.sub(ONE,self.exp_decay):
            raise ValueError('alpha_RS detached from 1-exp source operation')


def step(*,mult,tau_target,dt,exp_decay,slew_log=0,exp_profile='legacy-enclosure'):
    """Materialize deployed slew-disabled alpha_RS graph.

    ``target`` and ``applied`` are deliberately absent because source does not
    read them when ``slew_log <= 0``.  This prevents a ghost log/ratio witness
    from entering the default theorem path.
    """
    m=_q(mult,'adapt_RS_mult'); tau=_q(tau_target,'tau target'); h=_q(dt,'dt')
    sl=_q(slew_log,'adapt_RS_slew_log')
    if sl!=0: raise ValueError('current deployment theorem requires literal slew_log=0 branch')
    if m<=0 or tau<=0 or h<=0: raise ValueError('positive RS smoothing operands required')
    safe=clamp(tau,TIME_MIN,TIME_MAX)
    requested=B.mul(m,safe)
    lo=min(max(h,HORIZON_MIN),HORIZON_MAX)
    rssec=clamp(requested,lo,HORIZON_MAX)
    x=B.div(h,rssec)
    e=_q(exp_decay,'RS exp result')
    if not 0<e<=1: raise ValueError('RS exp result must lie in (0,1]')
    try:
        TAU.check_exp_result(x,e,exp_profile)
    except ValueError as exc:
        raise ValueError('RS exp witness detached from SAME rounded -dt/RS_sec argument/profile') from exc
    a=B.sub(ONE,e)
    return Step(m,tau,h,safe,requested,lo,rssec,x,e,a,exp_profile=exp_profile)


def _source_shape_matches():
    c=COMMON.read_text(); f=FILTER.read_text()
    common_needles=(
      'const float safe_time_scale =',
      'clampDynamicEmaTimeScaleSec(tau_sec);',
      'float horizon = mult * safe_time_scale;',
      'if (slew_log > 0.0f && target > 0.0f && applied > 0.0f &&',
      'return seastate::tuner::limits::clampDynamicEmaHorizonSec(horizon, dt);')
    filter_needles=(
      'const float RS_sec = seastate::common::adaptiveSmoothingHorizonSec(',
      'adapt_RS_mult_, tau_t, RS_t, tune_.RS_applied, adapt_RS_slew_log_, dt);',
      'const float alpha_RS = 1.0f - std::exp(-dt / RS_sec);')
    return all(x in c for x in common_needles) and all(x in f for x in filter_needles)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_RS_horizon_alpha_source_shape_matches':_source_shape_matches(),
      'deployed_slew_log_zero_removes_log_ratio_branch':True,
      'dynamic_tau_time_scale_binary32_clamp_materialized':True,
      'RS_horizon_binary32_multiply_and_final_clamp_materialized':True,
      'rounded_dt_over_RS_horizon_argument_materialized':True,
      'RS_exp_witness_bound_to_same_rounded_argument':True,
      'one_minus_exp_binary32_alpha_materialized':True,
      'target_libm_exp_correspondence_closed':False,
      'source_uniform_alpha_RS_supply_bound_closed':False,
      'persistent_RS_machine_real_recurrence_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
