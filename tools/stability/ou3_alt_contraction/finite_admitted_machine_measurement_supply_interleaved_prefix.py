"""Admitted machine finite-event supply through S service and accelerometer.

This strengthens ``finite_admitted_machine_prediction_supply_interleaved_prefix``.
The lower word already retains, on each admitted IMU edge, both global compiler
machine ActiveParameters, their persistent pseudo-S scheduler histories and the
full joint24/21x21 prediction displacement relative to the one executed exact
shipping prediction.

This layer propagates each local machine-coefficient perturbation through the
remaining measurement part of that SAME IMU event:

  machine prediction
    -> that compiler history's persistent queued a_w covariance floor
    -> covariance symmetry hygiene
    -> that compiler history's literal due/not-due S scheduler branch
    -> S=0 service with that compiler history's actual applied R_S
    -> held guarded accelerometer correction
    -> end-of-sample a_w floor request/snapshot using that machine Sigma_aw.

The two machine a_w-sync histories share the shipping request clock/predicate
with the exact event but retain their own historical target matrices.  A target
queued before the current tuner boundary therefore cannot be rewritten by that
boundary and cannot be borrowed from the exact shadow after the compiler tracks
have diverged.  Separate floor/eigensolver and measurement LDLT witnesses are
kept because a covariance/R_S perturbation may change those branches.

The exact side is not recomputed or refitted. It is the already executed lower
shipping event. This is a finite same-event arithmetic-supply relation, not a
second physical source history. The machine Racc path, guard/private-Mahony
floating arithmetic, native clock/Eigen/libm correspondence and source-uniform
supply bounds remain open. No storage search is authorized.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_admitted_machine_prediction_supply_interleaved_prefix as LOWER
from tools.stability.ou3_alt_contraction import finite_admitted_machine_tunestate_interleaved_prefix as MTUNE
from tools.stability.ou3_alt_contraction import finite_machine_prediction_displacement as DISP
from tools.stability.ou3_alt_contraction import finite_active_runtime_word as WORD
from tools.stability.ou3_alt_contraction import finite_post_prediction as POST
from tools.stability.ou3_alt_contraction import finite_measurement_runtime as MEAS
from tools.stability.ou3_alt_contraction import finite_periodic_aw_sync as AWSYNC
from tools.stability.ou3_alt_contraction import finite_source_bound_prediction_word as EXACTROOT
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE

QUALIFICATION='OU3_ALT_ADMITTED_MACHINE_MEASUREMENT_SUPPLY_INTERLEAVER_V3'


def _imu_steps(state):
    return LOWER._imu_steps(state.base.base)


def _preword(state):
    return LOWER._preword(state.base.base)


def _exact_aw(state):
    return _preword(state).live.live.live.aw_sync


def _same_aw_control(machine:AWSYNC.State,exact:AWSYNC.State):
    return machine.pending==exact.pending and machine.last_sync_time==exact.last_sync_time


def _subvec(a,b): return tuple(F(x)-F(y) for x,y in zip(a,b))
def _submat(a,b): return tuple(tuple(F(x)-F(y) for x,y in zip(ra,rb)) for ra,rb in zip(a,b))
def _supply(machine,exact): return DISP.Supply(_subvec(machine.z,exact.z),_submat(machine.covariance,exact.covariance))


@dataclass(frozen=True)
class ModeEvent:
    active:ACTIVE.ActiveParameters
    prediction:DISP.Relation
    scheduler_event:object
    aw_before:AWSYNC.State
    post_prediction:POST.PostPrediction
    aw_after_prediction:AWSYNC.State
    S_service:POST.ServicedPostPrediction
    accelerometer:MEAS.MeasurementRuntimeResult
    aw_after_tick:AWSYNC.State
    post_supply:DISP.Supply
    S_supply:DISP.Supply
    accel_supply:DISP.Supply
    mode:str
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.active,ACTIVE.ActiveParameters) or not isinstance(self.prediction,DISP.Relation):
            raise TypeError('machine active parameters and prediction relation required')
        if self.mode not in ('separate','fma'): raise ValueError('invalid compiler mode')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine measurement-supply qualification')
        if self.prediction.machine_roots.active!=self.active or self.prediction.machine_roots.active_join.mode!=self.mode:
            raise ValueError('measurement supply detached from same compiler active/prediction roots')
        if not all(isinstance(x,AWSYNC.State) for x in (self.aw_before,self.aw_after_prediction,self.aw_after_tick)):
            raise TypeError('persistent machine aw-sync ancestry required')
        if not isinstance(self.post_prediction,POST.PostPrediction) or not isinstance(self.S_service,POST.ServicedPostPrediction):
            raise TypeError('machine post-prediction and S-service results required')
        if not isinstance(self.accelerometer,MEAS.MeasurementRuntimeResult):
            raise TypeError('machine accelerometer result required')
        if self.post_prediction.scheduler!=self.scheduler_event.after_step:
            raise ValueError('machine post-prediction scheduler successor detached from persistent scheduler history')
        if self.post_prediction.S_service_due!=self.scheduler_event.due:
            raise ValueError('machine S due/not-due branch detached from persistent scheduler history')
        if self.S_service.prefix!=self.post_prediction:
            raise ValueError('machine S service detached from same post-prediction prefix')
        if self.accelerometer.state.reference!=self.prediction.machine.reference:
            raise ValueError('machine accelerometer correction changed the physical endpoint')
        if self.aw_after_prediction!=AWSYNC.consume_at_prediction(self.aw_before):
            raise ValueError('machine aw-floor consume detached from same pre-event snapshot')


@dataclass(frozen=True)
class State:
    base:LOWER.State
    separate_aw_sync:AWSYNC.State
    fma_aw_sync:AWSYNC.State
    entry_imu_steps:int
    measurement_steps:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.base,LOWER.State): raise TypeError('admitted machine prediction-supply state required')
        if not isinstance(self.separate_aw_sync,AWSYNC.State) or not isinstance(self.fma_aw_sync,AWSYNC.State):
            raise TypeError('both persistent compiler aw-sync states required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong admitted machine measurement-supply qualification')
        exact=_exact_aw(self)
        if not _same_aw_control(self.separate_aw_sync,exact) or not _same_aw_control(self.fma_aw_sync,exact):
            raise ValueError('machine aw-sync request control detached from same shipping event clock')
        if not isinstance(self.entry_imu_steps,int) or not isinstance(self.measurement_steps,int):
            raise TypeError('integer measurement-supply counters required')
        if self.entry_imu_steps<0 or self.measurement_steps<0:
            raise ValueError('nonnegative measurement-supply counters required')
        if _imu_steps(self)-self.entry_imu_steps!=self.measurement_steps:
            raise ValueError('machine measurement-supply count detached from admitted IMU count')


@dataclass(frozen=True)
class ImuResult:
    state:State
    lower:LOWER.ImuResult
    separate:ModeEvent
    fma:ModeEvent
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.ImuResult):
            raise TypeError('admitted machine measurement-supply result malformed')
        if not isinstance(self.separate,ModeEvent) or not isinstance(self.fma,ModeEvent):
            raise TypeError('both compiler measurement-supply events required')
        if self.separate.prediction!=self.lower.separate or self.fma.prediction!=self.lower.fma:
            raise ValueError('measurement supplies detached from lower same-event prediction supplies')
        if self.separate.mode!='separate' or self.fma.mode!='fma':
            raise ValueError('compiler measurement histories crossed')
        if self.state.separate_aw_sync!=self.separate.aw_after_tick or self.state.fma_aw_sync!=self.fma.aw_after_tick:
            raise ValueError('persistent machine aw-sync successor detached from same-event results')


@dataclass(frozen=True)
class CompleteWord:
    state:State
    lower:LOWER.CompleteWord
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.CompleteWord):
            raise TypeError('measurement-supply complete word requires lower complete word')
        if self.lower.state!=self.state.base:
            raise ValueError('measurement-supply complete word detached from lower prediction-supply word')
        if self.state.measurement_steps!=SOURCE.TRANSITIONS:
            raise ValueError('machine measurement supply not attached to all 600 admitted IMU edges')


def begin(base:LOWER.State):
    if not isinstance(base,LOWER.State): raise TypeError('admitted machine prediction-supply state required')
    _probe=type('_Probe',(),{'base':base,'separate_aw_sync':AWSYNC.State(),'fma_aw_sync':AWSYNC.State()})()
    exact=LOWER._preword(base.base).live.live.live.aw_sync
    return State(base,exact,exact,LOWER._imu_steps(base.base),0)


def _live_result(lower:LOWER.ImuResult):
    # prediction-supply -> scheduler -> machine-prediction -> machine-TuneState
    mtune=lower.lower.lower.lower
    return MTUNE._live_result(mtune.lower)


def _mode_event(*,mode,relation,scheduler_event,aw_before,live,segment,
                floor_solver_success=None,floor_eigenvectors=None,floor_eigenvalues=None,
                S_ldlt=None,accel_ldlt=None,S_alpha=1,S_radius=None,
                accel_alpha=1,accel_radius=F(2,5),temperature_c,
                aw_sync_adapt_every=F(1,10),aw_sync_enabled=True):
    if mode not in ('separate','fma'): raise ValueError('invalid compiler mode')
    if not isinstance(aw_sync_enabled,bool): raise TypeError('literal machine aw-sync enable branch required')
    active=relation.machine_roots.active
    target=AWSYNC.floor_target(aw_before)
    if target is None: target=active.Sigma_aw
    predicted=WORD.Predicted(active,relation.machine)
    post=WORD.post_prediction_from_active(predicted,h=segment.h,pending_aw_floor=aw_before.pending,
        aw_floor_target=target,scheduler=scheduler_event.after_retarget,
        floor_solver_success=floor_solver_success,floor_eigenvectors=floor_eigenvectors,
        floor_eigenvalues=floor_eigenvalues)
    aw_after_prediction=AWSYNC.consume_at_prediction(aw_before)
    if post.scheduler!=scheduler_event.after_step or post.S_service_due!=scheduler_event.due:
        raise ValueError('re-executed machine scheduler branch detached from carried scheduler event')
    serviced=WORD.service_S_from_active(active,post,ldlt=S_ldlt,alpha=S_alpha,radius=S_radius)
    conditioning=EXACTROOT._accel_conditioning(temperature_c)
    accel=MEAS.accelerometer_from_held_guarded_racc(serviced.state,segment,live.guarded,
        conditioning,live.racc,ldlt=accel_ldlt,alpha=accel_alpha,radius=accel_radius)
    t=live.tuner_suffix.state.time; cadence=F(aw_sync_adapt_every)
    due=aw_sync_enabled and (t-aw_after_prediction.last_sync_time)>cadence
    sync=AWSYNC.tick(aw_after_prediction,time=t,adapt_every=cadence,live=True,
        active_sigma=active.Sigma_aw if due else None,enabled=aw_sync_enabled,
        congruent=False,legacy=False)
    if sync.requested_now!=live.aw_sync.requested_now:
        raise ValueError('machine aw-sync request predicate detached from exact shipping event control')
    return ModeEvent(active,relation,scheduler_event,aw_before,post,aw_after_prediction,
        serviced,accel,sync.state,
        _supply(post.state,live.post_prediction.state),
        _supply(serviced.state,live.S_service.state),
        _supply(accel.state,live.accelerometer.state),mode)


def imu_step(state:State,*,
             separate_floor_solver_success=None,separate_floor_eigenvectors=None,separate_floor_eigenvalues=None,
             fma_floor_solver_success=None,fma_floor_eigenvectors=None,fma_floor_eigenvalues=None,
             separate_S_ldlt=None,fma_S_ldlt=None,
             separate_accel_ldlt:MEAS.SafeLDLT,fma_accel_ldlt:MEAS.SafeLDLT,
             **kwargs):
    if not isinstance(state,State): raise TypeError('admitted machine measurement-supply State required')
    if not isinstance(separate_accel_ldlt,MEAS.SafeLDLT) or not isinstance(fma_accel_ldlt,MEAS.SafeLDLT):
        raise TypeError('both compiler machine accelerometer LDLT branches required')
    restricted=kwargs.get('restricted'); temperature_c=kwargs.get('temperature_c')
    if restricted is None or temperature_c is None:
        raise TypeError('same admitted physical restriction and source-owned temperature required')
    segment=restricted.segment
    lower=LOWER.imu_step(state.base,**kwargs)
    live=_live_result(lower)
    sched=lower.lower
    common=dict(live=live,segment=segment,
                S_alpha=kwargs.get('S_alpha',1),S_radius=kwargs.get('S_radius'),
                accel_alpha=kwargs.get('accel_alpha',1),accel_radius=kwargs.get('accel_radius',F(2,5)),
                temperature_c=temperature_c,
                aw_sync_adapt_every=kwargs.get('aw_sync_adapt_every',F(1,10)),
                aw_sync_enabled=kwargs.get('aw_sync_enabled',True))
    sep=_mode_event(mode='separate',relation=lower.separate,scheduler_event=sched.separate,
        aw_before=state.separate_aw_sync,
        floor_solver_success=separate_floor_solver_success,
        floor_eigenvectors=separate_floor_eigenvectors,floor_eigenvalues=separate_floor_eigenvalues,
        S_ldlt=separate_S_ldlt,accel_ldlt=separate_accel_ldlt,**common)
    fma=_mode_event(mode='fma',relation=lower.fma,scheduler_event=sched.fma,
        aw_before=state.fma_aw_sync,
        floor_solver_success=fma_floor_solver_success,
        floor_eigenvectors=fma_floor_eigenvectors,floor_eigenvalues=fma_floor_eigenvalues,
        S_ldlt=fma_S_ldlt,accel_ldlt=fma_accel_ldlt,**common)
    nxt=State(lower.state,sep.aw_after_tick,fma.aw_after_tick,
              state.entry_imu_steps,state.measurement_steps+1)
    return ImuResult(nxt,lower,sep,fma)


def mag_step(state:State,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted machine measurement-supply State required')
    base,event=LOWER.mag_step(state.base,**kwargs)
    return State(base,state.separate_aw_sync,state.fma_aw_sync,state.entry_imu_steps,state.measurement_steps),event


def set_hold(state:State,*,hold):
    if not isinstance(state,State): raise TypeError('admitted machine measurement-supply State required')
    base,event=LOWER.set_hold(state.base,hold=hold)
    return State(base,state.separate_aw_sync,state.fma_aw_sync,state.entry_imu_steps,state.measurement_steps),event


def complete(state:State): return CompleteWord(state,LOWER.complete(state.base))


def readiness():
    low=LOWER.readiness()
    return {
      'admitted_machine_prediction_supply_word_consumed':low['machine_root_effect_injected_into_joint24_prediction_relation'],
      'persistent_separate_and_FMA_aw_sync_snapshot_histories_attached':True,
      'queued_machine_aw_floor_target_is_historical_same_mode_state':True,
      'pending_machine_aw_floor_branch_reexecuted_from_same_mode_snapshot':True,
      'machine_post_prediction_floor_and_hygiene_reexecuted_from_perturbed_covariance':True,
      'persistent_per_compiler_scheduler_due_not_due_branch_reused_and_checked':True,
      'machine_due_S_service_uses_same_compiler_applied_RS':True,
      'machine_RS_measurement_effect_attached':True,
      'held_guarded_accelerometer_reexecuted_from_machine_post_S_state':True,
      'accelerometer_measurement_propagates_machine_prediction_supply':True,
      'full_joint24_and_21x21_supply_retained_after_post_S_and_accelerometer':True,
      'machine_Racc_coefficient_displacement_attached':False,
      'machine_guard_private_Mahony_and_conditioning_float_correspondence_closed':False,
      'machine_aw_sync_clock_binary64_correspondence_closed':False,
      'scheduler_nextafter_binary32_correspondence_closed':False,
      'machine_floor_eigensolver_and_measurement_LDLT_finite_precision_closed':False,
      'source_uniform_machine_event_supply_bound_closed':False,
      'all_target_libm_and_Eigen_correspondence_closed':False,
      'all_event_arithmetic_witnesses_source_uniformly_qualified':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
