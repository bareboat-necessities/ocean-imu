"""Admitted whole-machine Live word with mandatory machine prediction roots.

This strengthens ``finite_admitted_machine_tunestate_interleaved_prefix``.  The
lower product already carries one admitted COMPLETE-BRMM history, one admitted
BIAS history, bounded IMU ISS forcing, coherent WPE/tau/sigma/R_S machine
histories, common pending-boundary semantics, and the applied machine
ActiveParameters for both global compiler histories.

The remaining prediction gap was that the machine tau/Sigma values were merely
joined to the exact shadow state; their actual OU/Q-axis coefficient roots were
not mandatory on the admitted IMU edge.  Here every IMU successor must also
materialize BOTH machine prediction-root relations on the SAME pre-event source
segment/raw packet.  Crucially, the exact side of each join is the
``ActiveParameters`` actually consumed by the executed prediction *after* any
sample-entry pending boundary, not the stale predecessor active state.

No second physical transition or filter successor is executed.  These are
coefficient relations attached to the one admitted event.  Pseudo-period
scheduler effects and R_S/S-service effects remain separate open obligations,
as do source-uniform bounds and target libm/Eigen correspondence.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_machine_tunestate_interleaved_prefix as LOWER
from tools.stability.ou3_alt_contraction import finite_machine_active_prediction_roots as MPRED
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE

QUALIFICATION='OU3_ALT_ADMITTED_MACHINE_PREDICTION_INTERLEAVER_V1'


def _word(state:LOWER.State):
    try:
        return state.base.base.prefix.prefix.live.live_word
    except AttributeError as exc:
        raise TypeError('whole-machine admitted state lost source-owning Live word') from exc


def _imu_steps(state:LOWER.State):
    try:
        return state.base.base.prefix.imu_steps
    except AttributeError as exc:
        raise TypeError('whole-machine admitted state lost IMU/source ordinal') from exc


@dataclass(frozen=True)
class State:
    base:LOWER.State
    entry_imu_steps:int
    prediction_steps:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.base,LOWER.State):
            raise TypeError('admitted whole-machine TuneState state required')
        if self.qualification!=QUALIFICATION:
            raise ValueError('wrong admitted machine-prediction qualification')
        if not isinstance(self.entry_imu_steps,int) or not isinstance(self.prediction_steps,int):
            raise TypeError('integer machine-prediction counters required')
        if self.entry_imu_steps<0 or self.prediction_steps<0:
            raise ValueError('nonnegative machine-prediction counters required')
        if _imu_steps(self.base)-self.entry_imu_steps!=self.prediction_steps:
            raise ValueError('machine prediction-root count detached from admitted IMU count')


@dataclass(frozen=True)
class ImuResult:
    state:State
    lower:LOWER.ImuResult
    separate_roots:MPRED.Roots
    fma_roots:MPRED.Roots
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.ImuResult):
            raise TypeError('admitted machine-prediction IMU result malformed')
        if not isinstance(self.separate_roots,MPRED.Roots) or not isinstance(self.fma_roots,MPRED.Roots):
            raise TypeError('both global compiler machine prediction roots required')
        if self.separate_roots.active!=self.state.base.separate_active:
            raise ValueError('separate prediction roots detached from applied machine state')
        if self.fma_roots.active!=self.state.base.fma_active:
            raise ValueError('FMA prediction roots detached from applied machine state')
        if self.separate_roots.active_join.mode!='separate' or self.fma_roots.active_join.mode!='fma':
            raise ValueError('machine prediction compiler histories crossed')


@dataclass(frozen=True)
class CompleteWord:
    state:State
    lower:LOWER.CompleteWord
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.CompleteWord):
            raise TypeError('machine-prediction complete word requires lower complete word')
        if self.lower.state!=self.state.base:
            raise ValueError('machine-prediction complete word detached from lower whole-machine word')
        if self.state.prediction_steps!=SOURCE.TRANSITIONS:
            raise ValueError('machine prediction roots were not attached to all 600 admitted IMU edges')


def begin(base:LOWER.State):
    if not isinstance(base,LOWER.State): raise TypeError('admitted whole-machine state required')
    return State(base,_imu_steps(base),0)


def imu_step(state:State,*,separate_machine_root_kwargs,fma_machine_root_kwargs,**kwargs):
    """Execute one lower admitted IMU edge and attach both machine root families.

    ``*_machine_root_kwargs`` contain only machine-coefficient arithmetic
    witnesses (OU exp/expm1, Q-axis branch exp/PSD and attitude witnesses).  The
    physical segment, raw packet and exact executed ActiveParameters are derived
    from the one lower admitted event and cannot be replaced through those dicts.
    """
    if not isinstance(state,State): raise TypeError('admitted machine-prediction State required')
    if not isinstance(separate_machine_root_kwargs,dict) or not isinstance(fma_machine_root_kwargs,dict):
        raise TypeError('both compiler machine-root witness dictionaries required')
    forbidden={'state','physical','raw','machine_active','mode','exact_active'}
    if forbidden & set(separate_machine_root_kwargs) or forbidden & set(fma_machine_root_kwargs):
        raise TypeError('machine-root dictionaries cannot override source/event ancestry')

    preword=_word(state.base)
    restricted=kwargs.get('restricted'); witness=kwargs.get('witness'); raw=kwargs.get('raw')
    if restricted is None or not isinstance(witness,SOURCE.StepWitness) or raw is None:
        raise TypeError('strong admitted IMU restriction, source witness and raw packet required')
    segment=restricted.segment
    physical=SOURCE.QualifiedPhysicalSegment(preword.source.root,witness,segment)

    lower=LOWER.imu_step(state.base,**kwargs)
    live=LOWER._live_result(lower.lower)
    exact_active=live.prediction.active
    sep=MPRED.build(preword,physical,raw,lower.state.separate_active,
                    mode='separate',exact_active=exact_active,**separate_machine_root_kwargs)
    fma=MPRED.build(preword,physical,raw,lower.state.fma_active,
                    mode='fma',exact_active=exact_active,**fma_machine_root_kwargs)
    nxt=State(lower.state,state.entry_imu_steps,state.prediction_steps+1)
    return ImuResult(nxt,lower,sep,fma)


def mag_step(state:State,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted machine-prediction State required')
    base,event=LOWER.mag_step(state.base,**kwargs)
    return State(base,state.entry_imu_steps,state.prediction_steps),event


def set_hold(state:State,*,hold):
    if not isinstance(state,State): raise TypeError('admitted machine-prediction State required')
    base,event=LOWER.set_hold(state.base,hold=hold)
    return State(base,state.entry_imu_steps,state.prediction_steps),event


def complete(state:State):
    return CompleteWord(state,LOWER.complete(state.base))


def readiness():
    low=LOWER.readiness(); roots=MPRED.readiness()
    return {
      'admitted_whole_machine_TuneState_word_consumed':low['Live_600_step_machine_TuneState_product_attached'],
      'same_admitted_source_segment_and_raw_packet_drive_both_machine_root_families':True,
      'post_boundary_exact_ActiveParameters_taken_from_executed_prediction':True,
      'separate_and_FMA_machine_OU_Qaxis_roots_required_on_every_strong_IMU_edge':True,
      'machine_prediction_roots_counted_against_same_admitted_IMU_ordinal':True,
      'complete_word_requires_machine_prediction_roots_on_all_600_IMU_edges':True,
      'machine_root_relation_uses_recomputed_machine_tau_Qaxis_branch':roots['small_general_Qaxis_branch_is_recomputed_from_machine_tau_not_shadow_tau'],
      'machine_prediction_root_relation_attached_to_admitted_event':True,
      'machine_pseudo_period_scheduler_effect_attached':False,
      'machine_RS_measurement_effect_attached':False,
      'OU_Qaxis_target_libm_and_Eigen_correspondence_closed':False,
      'source_uniform_machine_coefficient_supply_bound_closed':False,
      'all_event_arithmetic_witnesses_source_uniformly_qualified':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
