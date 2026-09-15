"""Attach canonical binary64 a_w-sync clock qualification to the strongest word.

The lower joined product already owns the single admitted physical/filter event,
including the machine guard -> private Mahony -> frontend/tuner -> Racc ->
accelerometer histories. This wrapper executes no second event. It inspects the
already-produced per-compiler a_w-sync predecessor/successor and proves the
shipping double ``time_ - last_aw_cov_sync_sec_ > adapt_every_secs_`` decision
matches the exact 5 ms graph on every represented IMU edge.

The certificate is deliberately limited to the current default 0.1 s cadence and
the bounded startup+600-edge horizon. The cadence is read from the persistent
runtime/deployment tuner configuration; event-local overrides are forbidden.
Indefinite clock lifetime remains fail-closed. No storage search is authorized.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_admitted_machine_joined_frontend_racc_interleaved_prefix as LOWER
from tools.stability.ou3_alt_contraction import finite_admitted_machine_measurement_supply_interleaved_prefix as MEAS
from tools.stability.ou3_alt_contraction import finite_aw_sync_clock_binary64 as CLOCK
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B

QUALIFICATION='OU3_ALT_ADMITTED_MACHINE_CLOCK_QUALIFIED_INTERLEAVER_V1'


@dataclass(frozen=True)
class State:
    base:LOWER.State
    entry_source_steps:int
    clock_steps:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.base,LOWER.State): raise TypeError('joined frontend/Racc state required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong clock-qualified word')
        if not isinstance(self.entry_source_steps,int) or not isinstance(self.clock_steps,int):
            raise TypeError('integer clock qualification counters required')
        if self.entry_source_steps<0 or self.clock_steps<0: raise ValueError('nonnegative counters required')
        if self.base.source_steps-self.entry_source_steps!=self.clock_steps:
            raise ValueError('clock qualification count detached from joined IMU count')


@dataclass(frozen=True)
class ModeClock:
    mode:str
    now_step:int
    last_step:int
    deployed_elapsed:F
    deployed_cadence:F
    due:bool
    def __post_init__(self):
        if self.mode not in ('separate','fma'): raise ValueError('literal compiler mode required')
        if not isinstance(self.due,bool): raise TypeError('literal due branch required')


@dataclass(frozen=True)
class ImuResult:
    state:State
    lower:LOWER.ImuResult
    separate:ModeClock
    fma:ModeClock
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.ImuResult):
            raise TypeError('clock-qualified IMU result malformed')
        if self.lower.state!=self.state.base: raise ValueError('clock-qualified successor detached from same lower event')


@dataclass(frozen=True)
class CompleteWord:
    state:State
    lower:LOWER.CompleteWord
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.CompleteWord):
            raise TypeError('clock-qualified complete word requires lower complete word')
        if self.lower.state!=self.state.base: raise ValueError('clock-qualified complete word detached')
        if self.state.clock_steps!=SOURCE.TRANSITIONS:
            raise ValueError('binary64 a_w-sync clock not qualified on all 600 IMU edges')


def _carried_cadence(base:LOWER.State):
    if not isinstance(base,LOWER.State): raise TypeError('joined frontend/Racc state required')
    mt=LOWER._mtune_state(base.base)
    runtime=LOWER._runtime(base.base)
    exact=F(runtime.candidate_cfg.adapt_every_sec)
    deployed=F(mt.deployment_cfg.adapt_every_sec)
    if B.rn32(exact)!=deployed:
        raise ValueError('a_w-sync cadence detached from carried runtime/deployment tuner config')
    if exact!=CLOCK.ADAPT_REAL or deployed!=CLOCK.ADAPT_FLOAT:
        raise ValueError('current bounded clock certificate requires carried shipping-default cadence')
    return exact,deployed


def begin(base:LOWER.State):
    if not isinstance(base,LOWER.State): raise TypeError('joined frontend/Racc state required')
    _carried_cadence(base)
    return State(base,base.source_steps,0)


def _mode_clock(mode,event,time):
    q=CLOCK.qualify_exact_tick(time,event.aw_after_prediction.last_sync_time,CLOCK.ADAPT_REAL)
    due=q['exact_due']
    if due:
        if event.aw_after_tick.last_sync_time!=F(time):
            raise ValueError(mode+' a_w-sync due edge failed to snapshot current exact time')
    elif event.aw_after_tick.last_sync_time!=event.aw_after_prediction.last_sync_time:
        raise ValueError(mode+' a_w-sync not-due edge changed last-sync clock')
    return ModeClock(mode,q['now_step'],q['last_step'],q['deployed_elapsed'],q['deployed_cadence'],due)


def imu_step(state:State,**kwargs):
    if not isinstance(state,State): raise TypeError('clock-qualified State required')
    if 'aw_sync_adapt_every' in kwargs:
        raise TypeError('a_w-sync cadence is carried by runtime configuration, not an event-local operand')
    exact_cadence,_=_carried_cadence(state.base)
    lower=LOWER.imu_step(state.base,aw_sync_adapt_every=exact_cadence,**kwargs)
    measurement=lower.lower.lower
    if not isinstance(measurement,MEAS.ImuResult):
        raise TypeError('joined word lost same-event machine measurement supply')
    live=MEAS._live_result(measurement.lower)
    t=live.tuner_suffix.state.time
    sep=_mode_clock('separate',measurement.separate,t)
    fma=_mode_clock('fma',measurement.fma,t)
    nxt=State(lower.state,state.entry_source_steps,state.clock_steps+1)
    return ImuResult(nxt,lower,sep,fma)


def mag_step(state:State,**kwargs):
    if not isinstance(state,State): raise TypeError('clock-qualified State required')
    base,event=LOWER.mag_step(state.base,**kwargs)
    return State(base,state.entry_source_steps,state.clock_steps),event


def set_hold(state:State,*,hold):
    if not isinstance(state,State): raise TypeError('clock-qualified State required')
    base,event=LOWER.set_hold(state.base,hold=hold)
    return State(base,state.entry_source_steps,state.clock_steps),event


def complete(state:State): return CompleteWord(state,LOWER.complete(state.base))


def readiness():
    c=CLOCK.readiness(); low=LOWER.readiness()
    return {
      'strong_joined_guard_frontend_Racc_word_consumed':low['machine_guard_excess_drives_machine_Racc_recurrence'],
      'shipping_inner_binary64_aw_sync_source_shape_bound':c['shipping_inner_time_and_aw_sync_source_shape_matches'],
      'default_aw_sync_due_partition_binary64_closed':c['canonical_aw_sync_binary64_predicate_closed'],
      'same_executed_event_aw_sync_clock_qualified':True,
      'complete_word_requires_clock_qualification_on_all_600_IMU_edges':True,
      'goLive_aw_sync_clock_reset_modeled_below':True,
      'adapt_every_runtime_value_ancestry_closed_for_current_word':True,
      'per_event_adapt_every_override_forbidden':True,
      'canonical_5ms_source_dt_ancestry_closed_for_current_word':True,
      'mutable_adapt_every_setter_history_before_runtime_root_needed':False,
      'arbitrary_dt_inner_clock_closed':False,
      'indefinite_inner_clock_lifetime_closed':False,
      'all_event_arithmetic_witnesses_source_uniformly_qualified':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
