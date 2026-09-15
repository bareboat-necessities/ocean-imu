"""Persistent dual-compiler binary32 WPE moment, log and usability machine product.

Each global compiler history carries its own full moment state and canonical
log-period state and one-way usable latch from reset. Both histories consume the SAME machine vertical
sample and static WPE configuration, but their stored arithmetic may differ.
The moment-horizon ``exp(log_period)`` witness is bound to the same stored log
state and, on a valid smoothed update, must be exactly the same exp result used
again by ``update_log_period_`` in that compiler history.

The usable gate consumes each mode's post-update log getter and rounded
elapsed/history comparisons. It is skipped on every nonproducing update.

A valid raw period is therefore no longer detached from the log/tuner machine
history. Target libm/compiler qualification and source-uniform bounds remain
explicitly open.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as SHADOW
from tools.stability.ou3_alt_contraction import finite_wpe_moment_binary32 as MOM
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as LOG
from tools.stability.ou3_alt_contraction import finite_wpe_usable_binary32 as USABLE
from tools.stability.ou3_alt_contraction import finite_wpe_uniform_bounds as BOUNDS

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
    usable:USABLE.PeriodWitness|None=None
    def __post_init__(self):
        if not isinstance(self.moment,dict): raise TypeError('WPE moment witness dictionary required')


@dataclass(frozen=True)
class State:
    cfg:MOM.Config
    separate:MOM.State=MOM.State()
    fma:MOM.State=MOM.State()
    logs:LOG.State=LOG.State()
    qualification:str=QUALIFICATION
    separate_usable:bool=False
    fma_usable:bool=False
    bounded_profile:bool=False
    supply_scale_exponent:int=0
    libm_profile:str=BOUNDS.RNE_PROFILE
    def __post_init__(self):
        if type(self.separate_usable) is not bool or type(self.fma_usable) is not bool:
            raise TypeError('literal per-compiler WPE usability latches required')
        if not isinstance(self.cfg,MOM.Config) or not isinstance(self.separate,MOM.State) or not isinstance(self.fma,MOM.State) or not isinstance(self.logs,LOG.State):
            raise TypeError('machine WPE config/moment/log states required')
        if not (self.separate.samples==self.fma.samples==self.logs.samples):
            raise ValueError('dual WPE machine sample counts detached')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong WPE machine product qualification')
        if type(self.bounded_profile) is not bool: raise TypeError('literal WPE bounds qualification flag required')
        if type(self.supply_scale_exponent) is not int or self.supply_scale_exponent<0:
            raise ValueError('nonnegative integer WPE amplitude envelope required')
        if not self.bounded_profile and self.supply_scale_exponent:
            raise ValueError('unqualified WPE state cannot carry an amplitude certificate')
        BOUNDS._profile_errors(self.libm_profile)
        if self.bounded_profile: BOUNDS.check_state(self)


def initial(shadow_cfg:SHADOW.WPEConfig,*,bounded_profile=False,libm_profile=BOUNDS.RNE_PROFILE):
    return State(MOM.Config.from_shadow(shadow_cfg),bounded_profile=bounded_profile,libm_profile=libm_profile)


def widen_supply(state:State,*,input_abs_upper):
    """Weaken the quantitative envelope without changing any machine state.

    This consumes a source-domain bound, not a trace-derived admission claim.
    It is useful at goLive, whose sensor residual domain can exceed startup's.
    """
    if not isinstance(state,State) or not state.bounded_profile:
        raise TypeError('previously qualified WPE history required for envelope widening')
    exponent=max(state.supply_scale_exponent,BOUNDS.scale_for_input(input_abs_upper,state.libm_profile))
    if exponent==state.supply_scale_exponent: return state
    return replace(state,supply_scale_exponent=exponent)


def same_machine_history(before:State,after:State):
    """Only monotone proof-envelope metadata may change at the Live join."""
    return (isinstance(before,State) and isinstance(after,State)
            and all(getattr(before,k) is getattr(after,k) for k in ('cfg','separate','fma','logs'))
            and all(getattr(before,k)==getattr(after,k) for k in
                    ('separate_usable','fma_usable','bounded_profile','qualification','libm_profile'))
            and after.supply_scale_exponent>=before.supply_scale_exponent)


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
    if state.bounded_profile:
        BOUNDS.check_state(state)
        if h!=BOUNDS.DT or abs(x)>BOUNDS.build(state.supply_scale_exponent,state.libm_profile)['vertical_input_abs_upper']:
            raise ValueError('WPE source input exceeds bounded profile')
        BOUNDS.check_mode_before(state.logs.separate,separate,state.supply_scale_exponent,state.libm_profile)
        BOUNDS.check_mode_before(state.logs.fma,fma,state.supply_scale_exponent,state.libm_profile)
    sp=_period(state.logs.separate,separate.horizon); fp=_period(state.logs.fma,fma.horizon)
    sr=MOM.step(state.separate,state.cfg,dt=h,vertical_accel=x,canonical_period=sp,
                libm_profile=state.libm_profile,**separate.moment)
    fr=MOM.step(state.fma,state.cfg,dt=h,vertical_accel=x,canonical_period=fp,
                libm_profile=state.libm_profile,**fma.moment)
    sb,sw=_bind_log(sr,state.logs.separate,separate); fb,fw=_bind_log(fr,state.logs.fma,fma)
    logs=LOG.advance_modes(state.logs,separate_produced=sb,fma_produced=fb,
                           separate_witness=sw,fma_witness=fw)
    sg=USABLE.update(state.separate_usable,produced_period=sb,
        log_period=logs.state.separate.log_period,elapsed=sr.state.elapsed,
        lambda_=state.cfg.lambda_,witness=separate.usable)
    fg=USABLE.update(state.fma_usable,produced_period=fb,
        log_period=logs.state.fma.log_period,elapsed=fr.state.elapsed,
        lambda_=state.cfg.lambda_,witness=fma.usable)
    nxt=State(state.cfg,sr.state,fr.state,logs.state,
        separate_usable=sg.after,fma_usable=fg.after,bounded_profile=state.bounded_profile,
        supply_scale_exponent=state.supply_scale_exponent,libm_profile=state.libm_profile)
    return nxt,sr,fr,logs


def readiness():
    m=MOM.readiness(); l=LOG.readiness()
    u=USABLE.readiness()
    return {
      'bounded_input_uniform_moment_raw_log_supplies_closed':BOUNDS.build()['reset_to_every_finite_prefix_bounded_input_induction_closed'],
      'bounded_profile_checks_same_log_exp_arguments_without_target_promotion':True,
      'qualification':QUALIFICATION,
      'two_persistent_machine_moment_histories_rooted_at_reset':True,
      'same_machine_vertical_sample_drives_both_compiler_WPE_histories':True,
      'mode_specific_raw_period_derived_from_same_mode_moments':m['raw_period_derived_from_same_machine_moments_and_sqrt_result'],
      'raw_period_bound_to_same_mode_std_log_input':True,
      'moment_horizon_period_bound_to_same_stored_log_state':True,
      'same_exp_log_period_value_bound_across_horizon_and_log_smoothing':True,
      'per_compiler_period_branch_divergence_retained':l['per_compiler_period_branch_divergence_representable'],
      'source_produced_per_compiler_usable_latches_retained':u['post_update_binary32_usability_predicates_materialized'],
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
