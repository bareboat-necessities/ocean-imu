"""Admitted machine-scheduler word with same-event full prediction supply.

This layer attaches the exact ``finite_machine_prediction_displacement`` relation
to every represented admitted IMU edge.  It does not accept a detached exact
prediction snapshot.  Instead it rebuilds the source-owned exact roots from the
same pre-event source word, physical segment, raw packet and exact runtime
witnesses, evaluates the paired prediction, and requires that result to equal
the prediction state actually produced inside the lower shipping IMU event.

Both global compiler machine-root families are then compared against that one
executed exact prediction.  The result retains full joint24 and 21x21 covariance
prediction displacement for separate and FMA histories.  Scheduler branching is
already carried by the lower layer.  Post-prediction floor/S-service and the
accelerometer correction are not re-evaluated from the perturbed prediction yet,
so this is still an intermediate finite-event supply, not the complete machine
shipping word and not a storage certificate.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_machine_scheduler_interleaved_prefix as LOWER
from tools.stability.ou3_alt_contraction import finite_admitted_machine_prediction_interleaved_prefix as MPREFIX
from tools.stability.ou3_alt_contraction import finite_admitted_machine_tunestate_interleaved_prefix as MTUNE
from tools.stability.ou3_alt_contraction import finite_machine_prediction_displacement as DISP
from tools.stability.ou3_alt_contraction import finite_source_bound_prediction_word as EXACTROOT
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE

QUALIFICATION='OU3_ALT_ADMITTED_MACHINE_PREDICTION_SUPPLY_INTERLEAVER_V1'
ROOT_KEYS=(
 'ou_alpha','ou_em1','boundary_noise_sqrt','rs_sqrt_scale',
 'bias_phi','bias_em1_2','angular_full','angular_half',
 'qaxis_marginal_exp','qaxis_final_exp','qaxis_marginal_psd','qaxis_final_psd')


def _preword(state:LOWER.State):
    # scheduler -> machine-prediction -> whole-machine TuneState -> WPE/tau...
    return MPREFIX._word(state.base.base)


def _imu_steps(state:LOWER.State):
    return MPREFIX._imu_steps(state.base.base)


@dataclass(frozen=True)
class State:
    base:LOWER.State
    entry_imu_steps:int
    supply_steps:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.base,LOWER.State): raise TypeError('admitted machine-scheduler state required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong admitted machine-prediction-supply qualification')
        if not isinstance(self.entry_imu_steps,int) or not isinstance(self.supply_steps,int) or self.entry_imu_steps<0 or self.supply_steps<0:
            raise ValueError('nonnegative integer prediction-supply counters required')
        if _imu_steps(self.base)-self.entry_imu_steps!=self.supply_steps:
            raise ValueError('machine prediction-supply count detached from admitted IMU count')


@dataclass(frozen=True)
class ImuResult:
    state:State
    lower:LOWER.ImuResult
    separate:DISP.Relation
    fma:DISP.Relation
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.ImuResult):
            raise TypeError('admitted machine prediction-supply result malformed')
        if not isinstance(self.separate,DISP.Relation) or not isinstance(self.fma,DISP.Relation):
            raise TypeError('both compiler prediction-displacement relations required')
        if self.separate.machine_roots!=self.lower.lower.separate_roots:
            raise ValueError('separate prediction supply detached from same-event machine roots')
        if self.fma.machine_roots!=self.lower.lower.fma_roots:
            raise ValueError('FMA prediction supply detached from same-event machine roots')
        if self.separate.exact!=self.fma.exact:
            raise ValueError('compiler histories detached from one executed exact prediction')


@dataclass(frozen=True)
class CompleteWord:
    state:State
    lower:LOWER.CompleteWord
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.CompleteWord):
            raise TypeError('prediction-supply complete word requires lower complete word')
        if self.lower.state!=self.state.base:
            raise ValueError('prediction-supply complete word detached from machine-scheduler word')
        if self.state.supply_steps!=SOURCE.TRANSITIONS:
            raise ValueError('full machine prediction supply not attached to all 600 admitted IMU edges')


def begin(base:LOWER.State):
    if not isinstance(base,LOWER.State): raise TypeError('admitted machine-scheduler state required')
    return State(base,_imu_steps(base),0)


def _exact_root_kwargs(kwargs):
    missing=[k for k in ('ou_alpha','ou_em1','qaxis_marginal_psd','qaxis_final_psd') if k not in kwargs]
    if missing: raise TypeError('exact source-owned prediction witnesses missing '+repr(missing))
    return {k:kwargs.get(k) for k in ROOT_KEYS}


def imu_step(state:State,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted machine prediction-supply State required')
    restricted=kwargs.get('restricted'); witness=kwargs.get('witness'); raw=kwargs.get('raw')
    if restricted is None or not isinstance(witness,SOURCE.StepWitness) or raw is None:
        raise TypeError('same admitted restriction, source witness and raw packet required')
    preword=_preword(state.base); segment=restricted.segment
    physical=SOURCE.QualifiedPhysicalSegment(preword.source.root,witness,segment)
    exact_roots=EXACTROOT.build(preword,physical,raw,**_exact_root_kwargs(kwargs))

    lower=LOWER.imu_step(state.base,**kwargs)
    live=MTUNE._live_result(lower.lower.lower.lower)
    pre_core=preword.live.live.live.mekf
    common=dict(Qbase=preword.runtime.Qbase,
        use_exact_attitude_Q=kwargs.get('use_exact_attitude_Q',True),
        attitude_first_ldlt_success=kwargs.get('attitude_first_ldlt_success',True),
        attitude_second_ldlt_success=kwargs.get('attitude_second_ldlt_success'))
    sep=DISP.compare(pre_core,segment,raw,exact_roots,lower.lower.separate_roots,
        machine_predecessor=kwargs.get('separate_machine_predecessor'),**common)
    fma=DISP.compare(pre_core,segment,raw,exact_roots,lower.lower.fma_roots,
        machine_predecessor=kwargs.get('fma_machine_predecessor'),**common)
    if sep.exact!=live.prediction.state or fma.exact!=live.prediction.state:
        raise ValueError('rebuilt exact prediction detached from executed admitted IMU event')
    nxt=State(lower.state,state.entry_imu_steps,state.supply_steps+1)
    return ImuResult(nxt,lower,sep,fma)


def mag_step(state:State,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted machine prediction-supply State required')
    base,event=LOWER.mag_step(state.base,**kwargs)
    return State(base,state.entry_imu_steps,state.supply_steps),event


def set_hold(state:State,*,hold):
    if not isinstance(state,State): raise TypeError('admitted machine prediction-supply State required')
    base,event=LOWER.set_hold(state.base,hold=hold)
    return State(base,state.entry_imu_steps,state.supply_steps),event


def complete(state:State): return CompleteWord(state,LOWER.complete(state.base))


def readiness():
    low=LOWER.readiness(); disp=DISP.readiness()
    return {
      'admitted_machine_scheduler_word_consumed':low['machine_pseudo_period_scheduler_effect_attached'],
      'source_owned_exact_prediction_rebuilt_from_same_pre_event_word':True,
      'rebuilt_exact_prediction_required_equal_executed_shipping_prediction':True,
      'separate_and_FMA_full_joint24_covariance_prediction_supplies_attached_same_event':True,
      'machine_prediction_supply_uses_full_24_and_21x21_relation':disp['full_prediction_coefficient_displacement_relation_attached'],
      'complete_word_requires_full_prediction_supply_on_all_600_IMU_edges':True,
      'machine_root_effect_injected_into_joint24_prediction_relation':True,
      'post_prediction_floor_and_S_service_propagate_machine_prediction_supply':False,
      'machine_RS_measurement_effect_attached':False,
      'accelerometer_measurement_propagates_machine_prediction_supply':False,
      'source_uniform_machine_prediction_supply_bound_closed':False,
      'all_target_libm_and_Eigen_correspondence_closed':False,
      'all_event_arithmetic_witnesses_source_uniformly_qualified':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
