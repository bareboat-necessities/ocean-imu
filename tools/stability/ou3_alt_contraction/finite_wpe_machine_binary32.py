"""Persistent dual-compiler binary32 WPE moment + log machine product.

Each global compiler history carries its own full moment state and canonical
log-period state from reset. Both histories consume the SAME machine vertical
sample and static WPE configuration, but their stored arithmetic may differ.
The moment-horizon ``exp(log_period)`` witness is bound to the same stored log
state and, on a valid smoothed update, must be exactly the same exp result used
again by ``update_log_period_`` in that compiler history.

A valid raw period is therefore no longer detached from the log/tuner machine
history. Target libm/compiler qualification and source-uniform bounds remain
explicitly open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as SHADOW
from tools.stability.ou3_alt_contraction import finite_wpe_moment_binary32 as MOM
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as LOG

QUALIFICATION='OU3_ALT_WPE_MACHINE_BINARY32_PRODUCT_V1'


def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(name+' must be binary32')
    return q


@dataclass(frozen=True)
class HorizonPeriodWitness:
    log_period:F
    period_exp:F
    def __post_init__(self):
        object.__setattr__(self,'log_period',_q(self.log_period,'WPE horizon log period'))
        object.__setattr__(self,'period_exp',_q(self.period_exp,'WPE horizon period exp'))
        if self.period_exp<=0: raise ValueError('positive WPE horizon period required')


@dataclass(frozen=True)
class RawLogBinding:
    raw_period:F
    log_raw:F
    def __post_init__(self):
        object.__setattr__(self,'raw_period',_q(self.raw_period,'WPE raw period binding'))
        object.__setattr__(self,'log_raw',_q(self.log_raw,'WPE log(raw) binding'))
        if self.raw_period<=0: raise ValueError('positive WPE raw period binding required')


@dataclass(frozen=True)
class ModeWitnesses:
    moment:dict
    horizon:HorizonPeriodWitness|None=None
    raw_log:RawLogBinding|None=None
    log:LOG.InitWitness|LOG.SmoothWitness|None=None
    def __post_init__(self):
        if not isinstance(self.moment,dict): raise TypeError('WPE moment witness dictionary required')


@dataclass(frozen=True)
class State:
    cfg:MOM.Config
    separate:MOM.State=MOM.State()
    fma:MOM.State=MOM.State()
    logs:LOG.State=LOG.State()
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.cfg,MOM.Config) or not isinstance(self.separate,MOM.State) or not isinstance(self.fma,MOM.State) or not isinstance(self.logs,LOG.State):
            raise TypeError('machine WPE config/moment/log states required')
        if not (self.separate.samples==self.fma.samples==self.logs.samples):
            raise ValueError('dual WPE machine sample counts detached')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong WPE machine product qualification')


def initial(shadow_cfg:SHADOW.WPEConfig):
    return State(MOM.Config.from_shadow(shadow_cfg))


def _period(track:LOG.Track,w:HorizonPeriodWitness|None):
    if track.log_period is None:
        if w is not None: raise ValueError('uninitialized WPE log consumes no horizon exp witness')
        return None
    if not isinstance(w,HorizonPeriodWitness): raise TypeError('initialized WPE log requires horizon exp witness')
    if w.log_period!=track.log_period: raise ValueError('WPE horizon exp detached from same stored log period')
    return w.period_exp


def _bind_log(moment:MOM.StepResult,track:LOG.Track,w:ModeWitnesses):
    produced=moment.branch=='valid-period'
    if not produced:
        if w.raw_log is not None or w.log is not None:
            raise ValueError('nonproducing machine WPE branch consumes no raw-log/log-update witness')
        return False,None
    if not isinstance(w.raw_log,RawLogBinding): raise TypeError('valid machine WPE period requires raw/log binding')
    if w.raw_log.raw_period!=moment.raw_period: raise ValueError('WPE std::log input detached from same machine raw period')
    if track.log_period is None:
        if not isinstance(w.log,LOG.InitWitness): raise TypeError('first machine WPE period requires InitWitness')
    else:
        if not isinstance(w.log,LOG.SmoothWitness): raise TypeError('initialized machine WPE period requires SmoothWitness')
        if not isinstance(w.horizon,HorizonPeriodWitness): raise TypeError('smoothed machine WPE period requires horizon witness')
        if w.log.sea_period_exp!=w.horizon.period_exp:
            raise ValueError('same exp(log_period) call value not reused across moment horizon/log smoothing ancestry')
    if w.log.log_raw!=w.raw_log.log_raw:
        raise ValueError('machine WPE log witness detached from std::log result for same raw period')
    return True,w.log


def step(state:State,*,dt,vertical_accel,separate:ModeWitnesses,fma:ModeWitnesses):
    if not isinstance(state,State): raise TypeError('WPE machine State required')
    if not isinstance(separate,ModeWitnesses) or not isinstance(fma,ModeWitnesses):
        raise TypeError('both compiler WPE witnesses required')
    x=_q(vertical_accel,'common machine WPE vertical input'); h=_q(dt,'common machine WPE dt')
    sp=_period(state.logs.separate,separate.horizon); fp=_period(state.logs.fma,fma.horizon)
    sr=MOM.step(state.separate,state.cfg,dt=h,vertical_accel=x,canonical_period=sp,**separate.moment)
    fr=MOM.step(state.fma,state.cfg,dt=h,vertical_accel=x,canonical_period=fp,**fma.moment)
    sb,sw=_bind_log(sr,state.logs.separate,separate); fb,fw=_bind_log(fr,state.logs.fma,fma)
    logs=LOG.advance_modes(state.logs,separate_produced=sb,fma_produced=fb,
                           separate_witness=sw,fma_witness=fw)
    nxt=State(state.cfg,sr.state,fr.state,logs.state)
    return nxt,sr,fr,logs


def readiness():
    m=MOM.readiness(); l=LOG.readiness()
    return {
      'qualification':QUALIFICATION,
      'two_persistent_machine_moment_histories_rooted_at_reset':True,
      'same_machine_vertical_sample_drives_both_compiler_WPE_histories':True,
      'mode_specific_raw_period_derived_from_same_mode_moments':m['raw_period_derived_from_same_machine_moments_and_sqrt_result'],
      'raw_period_bound_to_same_mode_std_log_input':True,
      'moment_horizon_period_bound_to_same_stored_log_state':True,
      'same_exp_log_period_value_bound_across_horizon_and_log_smoothing':True,
      'per_compiler_period_branch_divergence_retained':l['per_compiler_period_branch_divergence_representable'],
      'dual_compiler_persistent_moment_history_attached':True,
      'WPE_raw_period_binary32_production_attached_to_log_history':True,
      'target_exp_log_sqrt_libm_correspondence_closed':False,
      'compiler_profile_selection_closed':False,
      'source_uniform_WPE_machine_supply_bounds_closed':False,
      'startup_frontend_vertical_ancestry_attached':False,
      'Live_600_step_WPE_machine_history_attached':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,
    }
