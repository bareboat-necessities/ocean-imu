"""Finite shipping R_acc selection immediately before the accelerometer event.

Shipping calls ``apply_racc_vibration_inflation_`` after the current
AccelVibrationGuard step and before ``measurement_update_acc_only``.  This
module materializes that exact control relation in real arithmetic.  The
shipping default vibration gain is 0.75; quiet operation is nevertheless an
exact no-write path because guard ``excessRms`` is zero below engagement.

Optional vessel-RAO weighting is Live-only and uses the previous measurement-
only schedule.  The RAO hypot and sqrt evaluations and the final effective-std
sqrt are exact witnesses; deployed binary32 ancestry remains open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as GUARD
from tools.stability.ou3_alt_contraction.finite_tuner_commit import TuneState


def R(x): return M.rational(x)
def V(x): return tuple(M.vec(x,3))
def clamp(x,lo,hi): return min(max(x,lo),hi)
ONE=(F(1),F(1),F(1))


@dataclass(frozen=True)
class RaoConfig:
    max_std_scale:F=F(1)
    transition_snr:F=F(2)
    horizontal_x_tau:F=F(7,10)
    horizontal_y_tau:F=F(1)
    heave_period:F=F(12,5)
    heave_damping:F=F(9,20)
    two_pi:F=F(6283185307179586,10**15)
    def __post_init__(self):
        for n in ('max_std_scale','transition_snr','horizontal_x_tau','horizontal_y_tau','heave_period','heave_damping','two_pi'):
            object.__setattr__(self,n,R(getattr(self,n)))
        if self.max_std_scale<0 or self.transition_snr<=0 or self.horizontal_x_tau<0 or self.horizontal_y_tau<0 or self.heave_period<=0 or self.heave_damping<=0 or self.two_pi<=0:
            raise ValueError('invalid RAO noise-weight configuration')


@dataclass(frozen=True)
class Config:
    vibration_gain:F=F(3,4)
    sigma_coeff:F=F(9,10)
    rao:RaoConfig=RaoConfig()
    def __post_init__(self):
        vg,sc=R(self.vibration_gain),R(self.sigma_coeff)
        if vg<0 or sc<=0: raise ValueError('nonnegative vibration gain and positive sigma coefficient required')
        object.__setattr__(self,'vibration_gain',vg); object.__setattr__(self,'sigma_coeff',sc)
        if not isinstance(self.rao,RaoConfig): raise TypeError('RaoConfig required')


@dataclass(frozen=True)
class State:
    inflated:bool=False
    effective_std:tuple=(F(0),F(0),F(0))
    def __post_init__(self):
        if not isinstance(self.inflated,bool): raise TypeError('literal Racc inflated flag required')
        e=V(self.effective_std)
        if any(x<0 for x in e): raise ValueError('Racc std must be nonnegative')
        object.__setattr__(self,'effective_std',e)


@dataclass(frozen=True)
class RaoWitness:
    inverse_heave:F
    scale_x:F
    scale_y:F
    def __post_init__(self):
        for n in ('inverse_heave','scale_x','scale_y'):
            x=R(getattr(self,n)); object.__setattr__(self,n,x)
            if x<0: raise ValueError('nonnegative RAO witness required')


@dataclass(frozen=True)
class EffectiveSqrtWitness:
    values:tuple
    def __post_init__(self):
        v=V(self.values)
        if any(x<0 for x in v): raise ValueError('nonnegative effective-std sqrt witness required')
        object.__setattr__(self,'values',v)


@dataclass(frozen=True)
class Result:
    state:State
    effective_std:tuple
    covariance:tuple
    scales:tuple
    excess_rms:F
    wrote_Racc:bool
    restored_nominal:bool


def _diag_sq(std):
    s=V(std)
    return tuple(tuple(s[i]*s[i] if i==j else F(0) for j in range(3)) for i in range(3))


def _rao_scales(cfg:RaoConfig,*,live,vertical_std,frequency,nominal,witness:RaoWitness|None):
    if not isinstance(live,bool): raise TypeError('literal Live branch required')
    vertical_std,frequency=R(vertical_std),R(frequency); nominal=V(nominal)
    armed=(cfg.max_std_scale>1 and cfg.transition_snr>0 and vertical_std>=0 and frequency>0)
    if not live or not armed:
        if witness is not None: raise ValueError('inactive RAO branch consumes no hypot/sqrt witness')
        return ONE
    if witness is None: raise ValueError('active RAO branch requires same-expression witnesses')
    fT=frequency*cfg.heave_period
    inv2=(1-fT*fT)**2+(2*cfg.heave_damping*fT)**2
    if witness.inverse_heave*witness.inverse_heave != inv2:
        raise ValueError('inverse-heave hypot witness detached from same frequency/response')
    maximum=min(cfg.max_std_scale,F(4))
    out=[F(1),F(1),F(1)]
    for i,(tau,root) in enumerate(((cfg.horizontal_x_tau,witness.scale_x),(cfg.horizontal_y_tau,witness.scale_y))):
        if nominal[i]<=0: continue
        wt=cfg.two_pi*frequency*tau
        snr=vertical_std*witness.inverse_heave/((1+wt*wt)*nominal[i])
        u=clamp(snr/cfg.transition_snr-1,F(0),F(1))
        weight=1-u*u*(3-2*u)
        rad=1+(maximum*maximum-1)*weight
        if root*root != rad: raise ValueError('RAO scale sqrt witness detached from same SNR smoothstep')
        out[i]=root
    return tuple(out)


def step(state:State,cfg:Config,guard:GUARD.Result,*,nominal_std,tune:TuneState,
         preupdate_frequency,live,rao_witness:RaoWitness|None=None,
         effective_sqrt:EffectiveSqrtWitness|None=None):
    if not isinstance(state,State) or not isinstance(cfg,Config) or not isinstance(guard,GUARD.Result) or not isinstance(tune,TuneState):
        raise TypeError('Racc state/config, same guard result and previous TuneState required')
    base=V(nominal_std); f=R(preupdate_frequency)
    if cfg.vibration_gain<=0 and cfg.rao.max_std_scale<=1 and not state.inflated:
        if rao_witness is not None or effective_sqrt is not None: raise ValueError('fully dormant Racc branch consumes no witnesses')
        return Result(state,state.effective_std,_diag_sq(state.effective_std),ONE,guard.excess_rms,False,False)
    if min(base)<=0:
        if rao_witness is not None or effective_sqrt is not None: raise ValueError('unknown nominal Racc branch consumes no scale witnesses')
        return Result(state,state.effective_std,_diag_sq(state.effective_std),ONE,guard.excess_rms,False,False)

    vertical_std=tune.sigma_applied/cfg.sigma_coeff
    scales=_rao_scales(cfg.rao,live=live,vertical_std=vertical_std,frequency=f,nominal=base,witness=rao_witness)
    excess=R(guard.excess_rms)
    if excess<=0 and max(scales)<=1:
        if effective_sqrt is not None: raise ValueError('dormant restore branch consumes no effective sqrt witness')
        if state.inflated:
            nxt=State(False,base)
            return Result(nxt,base,_diag_sq(base),scales,excess,True,True)
        return Result(state,state.effective_std,_diag_sq(state.effective_std),scales,excess,False,False)

    if effective_sqrt is None: raise ValueError('inflating Racc branch requires effective sqrt witness')
    added=cfg.vibration_gain*excess
    vals=[]
    for i in range(3):
        rad=(base[i]*scales[i])**2+added**2
        root=effective_sqrt.values[i]
        if root*root != rad: raise ValueError('effective Racc sqrt witness detached from same base/RAO/vibration expression')
        vals.append(root)
    eff=tuple(vals); nxt=State(True,eff)
    return Result(nxt,eff,_diag_sq(eff),scales,excess,True,False)


def readiness():
    return {
      'shipping_default_vibration_gain_three_quarters':True,
      'same_guard_excess_RMS_drives_Racc':True,
      'previous_TuneState_drives_low_wave_scale':True,
      'preupdate_frequency_drives_low_wave_scale':True,
      'RAO_disabled_default_branch_materialized':True,
      'RAO_live_horizontal_smoothstep_branch_materialized':True,
      'one_time_nominal_restore_materialized':True,
      'effective_std_to_diagonal_covariance_matches_set_Racc_std':True,
      'hypot_and_sqrt_binary32_ancestry_attached':False,
      'nominal_Racc_stage_ancestry_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
