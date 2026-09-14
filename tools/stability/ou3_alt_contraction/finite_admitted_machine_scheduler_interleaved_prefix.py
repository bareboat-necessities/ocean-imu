"""Admitted machine-prediction word with persistent per-compiler S schedulers.

A machine pseudo-period displacement is hybrid, not merely additive: at a
sample-entry pending boundary shipping retargets the persistent S scheduler
before prediction, and after prediction the scheduler may choose a different
due/not-due branch. Therefore a coefficient-only ``delta period`` supply is
insufficient.

This layer strengthens ``finite_admitted_machine_prediction_interleaved_prefix``
with TWO persistent scheduler histories, one for each global compiler mode.
Each scheduler:

* is rooted at goLive by retargeting the SAME startup predecessor scheduler as
  the exact shadow;
* is retargeted iff the same IMU event consumed the machine pending boundary;
* then advances by the SAME admitted physical step duration;
* retains its literal due/not-due decision as theorem data;
* is unchanged by asynchronous MAG/HOLD events.

No second filter transition is executed. The machine S measurement itself is
still open: this module stops at the scheduler branch and leaves R_S/LDLT
composition to the next layer. Binary32 ``nextafter`` correspondence remains
an explicit witness obligation on overdue retargets.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_machine_prediction_interleaved_prefix as LOWER
from tools.stability.ou3_alt_contraction import finite_admitted_machine_tunestate_interleaved_prefix as MTUNE
from tools.stability.ou3_alt_contraction import finite_startup_live_machine_tunestate_bridge as GO
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import finite_post_prediction as POST
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE

QUALIFICATION='OU3_ALT_ADMITTED_MACHINE_SCHEDULER_INTERLEAVER_V1'


@dataclass(frozen=True)
class State:
    base:LOWER.State
    separate_scheduler:POST.Scheduler
    fma_scheduler:POST.Scheduler
    entry_imu_steps:int
    scheduler_steps:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.base,LOWER.State): raise TypeError('admitted machine-prediction state required')
        if not isinstance(self.separate_scheduler,POST.Scheduler) or not isinstance(self.fma_scheduler,POST.Scheduler):
            raise TypeError('both persistent machine S schedulers required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine-scheduler qualification')
        if not isinstance(self.entry_imu_steps,int) or not isinstance(self.scheduler_steps,int) or self.entry_imu_steps<0 or self.scheduler_steps<0:
            raise ValueError('nonnegative integer scheduler counters required')
        self.base.base.separate_active.require_scheduler(self.separate_scheduler)
        self.base.base.fma_active.require_scheduler(self.fma_scheduler)
        if LOWER._imu_steps(self.base.base)-self.entry_imu_steps!=self.scheduler_steps:
            raise ValueError('machine scheduler-step count detached from admitted IMU count')


@dataclass(frozen=True)
class SchedulerEvent:
    before:POST.Scheduler
    after_retarget:POST.Scheduler
    after_step:POST.Scheduler
    due:bool
    retargeted:bool
    def __post_init__(self):
        if not all(isinstance(x,POST.Scheduler) for x in (self.before,self.after_retarget,self.after_step)):
            raise TypeError('scheduler event requires literal scheduler states')
        if not isinstance(self.due,bool) or not isinstance(self.retargeted,bool):
            raise TypeError('literal scheduler branch flags required')


@dataclass(frozen=True)
class ImuResult:
    state:State
    lower:LOWER.ImuResult
    separate:SchedulerEvent
    fma:SchedulerEvent
    exact_due:bool
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.ImuResult):
            raise TypeError('machine-scheduler IMU result malformed')
        if not isinstance(self.separate,SchedulerEvent) or not isinstance(self.fma,SchedulerEvent) or not isinstance(self.exact_due,bool):
            raise TypeError('machine/exact scheduler decisions required')
        if self.state.separate_scheduler!=self.separate.after_step or self.state.fma_scheduler!=self.fma.after_step:
            raise ValueError('carried machine schedulers detached from same-event successors')


@dataclass(frozen=True)
class CompleteWord:
    state:State
    lower:LOWER.CompleteWord
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.CompleteWord):
            raise TypeError('machine-scheduler complete word requires lower complete word')
        if self.lower.state!=self.state.base:
            raise ValueError('machine-scheduler complete word detached from machine-prediction word')
        if self.state.scheduler_steps!=SOURCE.TRANSITIONS:
            raise ValueError('machine scheduler recurrence not attached to all 600 admitted IMU edges')


def begin(base:LOWER.State,separate_scheduler:POST.Scheduler,fma_scheduler:POST.Scheduler):
    if not isinstance(base,LOWER.State): raise TypeError('admitted machine-prediction state required')
    return State(base,separate_scheduler,fma_scheduler,LOWER._imu_steps(base.base),0)


def begin_from_goLive(base:LOWER.State,go:GO.Result,*,scheduler_before:POST.Scheduler,
                      exact_scheduler_park:ACTIVE.NextafterParkWitness|None=None,
                      separate_scheduler_park:ACTIVE.NextafterParkWitness|None=None,
                      fma_scheduler_park:ACTIVE.NextafterParkWitness|None=None):
    """Root exact and both machine schedulers at one startup predecessor state."""
    if not isinstance(base,LOWER.State) or not isinstance(go,GO.Result) or not isinstance(scheduler_before,POST.Scheduler):
        raise TypeError('machine-prediction Live state, whole-machine goLive result and predecessor scheduler required')
    exact=ACTIVE.retarget_scheduler(go.live.active,scheduler_before,park=exact_scheduler_park)
    if exact!=go.live.state.scheduler:
        raise ValueError('declared startup scheduler predecessor detached from exact goLive scheduler')
    sep=ACTIVE.retarget_scheduler(go.separate_active,scheduler_before,park=separate_scheduler_park)
    fma=ACTIVE.retarget_scheduler(go.fma_active,scheduler_before,park=fma_scheduler_park)
    return begin(base,sep,fma)


def _advance_one(before:POST.Scheduler,active:ACTIVE.ActiveParameters,*,boundary_consumed:bool,h,
                 park:ACTIVE.NextafterParkWitness|None):
    if not isinstance(boundary_consumed,bool): raise TypeError('literal boundary-consumed flag required')
    if boundary_consumed:
        rooted=ACTIVE.retarget_scheduler(active,before,park=park)
    else:
        if park is not None: raise ValueError('no machine boundary retarget consumes no nextafter witness')
        active.require_scheduler(before)
        rooted=before
    due,after=rooted.step(h)
    return SchedulerEvent(before,rooted,after,due,boundary_consumed)


def imu_step(state:State,*,separate_scheduler_park=None,fma_scheduler_park=None,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted machine-scheduler State required')
    restricted=kwargs.get('restricted')
    if restricted is None: raise TypeError('same admitted physical restriction required')
    h=restricted.segment.h
    lower=LOWER.imu_step(state.base,**kwargs)
    mb=lower.lower.machine_boundary
    sep=_advance_one(state.separate_scheduler,lower.state.base.separate_active,
                     boundary_consumed=mb.consumed,h=h,park=separate_scheduler_park)
    fma=_advance_one(state.fma_scheduler,lower.state.base.fma_active,
                     boundary_consumed=mb.consumed,h=h,park=fma_scheduler_park)
    live=MTUNE._live_result(lower.lower.lower)
    exact_due=live.post_prediction.S_service_due
    nxt=State(lower.state,sep.after_step,fma.after_step,state.entry_imu_steps,state.scheduler_steps+1)
    return ImuResult(nxt,lower,sep,fma,exact_due)


def mag_step(state:State,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted machine-scheduler State required')
    base,event=LOWER.mag_step(state.base,**kwargs)
    return State(base,state.separate_scheduler,state.fma_scheduler,state.entry_imu_steps,state.scheduler_steps),event


def set_hold(state:State,*,hold):
    if not isinstance(state,State): raise TypeError('admitted machine-scheduler State required')
    base,event=LOWER.set_hold(state.base,hold=hold)
    return State(base,state.separate_scheduler,state.fma_scheduler,state.entry_imu_steps,state.scheduler_steps),event


def complete(state:State): return CompleteWord(state,LOWER.complete(state.base))


def readiness():
    low=LOWER.readiness()
    return {
      'admitted_machine_prediction_word_consumed':low['machine_prediction_root_relation_attached_to_admitted_event'],
      'separate_and_FMA_pseudoS_scheduler_states_carried_persistently':True,
      'goLive_machine_schedulers_rooted_from_same_startup_predecessor_as_exact_shadow':True,
      'pending_boundary_retargets_each_machine_scheduler_before_same_sample_due_decision':True,
      'same_admitted_physical_dt_advances_exact_and_both_machine_scheduler_histories':True,
      'machine_due_not_due_branch_retained_per_compiler_history':True,
      'MAG_and_HOLD_preserve_both_machine_schedulers_by_identity':True,
      'complete_word_requires_both_scheduler_recurrences_on_all_600_IMU_edges':True,
      'machine_pseudo_period_scheduler_effect_attached':True,
      'machine_RS_measurement_effect_attached':False,
      'scheduler_nextafter_binary32_correspondence_closed':False,
      'machine_due_branch_S_measurement_composed':False,
      'source_uniform_machine_coefficient_supply_bound_closed':False,
      'all_event_arithmetic_witnesses_source_uniformly_qualified':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
