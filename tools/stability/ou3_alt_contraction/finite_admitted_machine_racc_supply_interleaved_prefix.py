"""Admitted machine Racc coefficient supply on the same finite IMU event.

The lower machine-measurement product already reexecutes each compiler history's
prediction, historical a_w floor, S scheduler/service and accelerometer using the
exact executed Racc object. This wrapper removes that last tuner-coefficient
identification without inventing a second physical history.

Shipping Racc reads the SAMPLE-ENTRY TuneState ``sigma_applied`` on EVERY
IMU call.  Despite its name, this field advances in the tuner candidate step
on every post-Cold sample; it is not held until the next pending MEKF commit.
The two sigma readouts below are the last values consumed by Racc, not a second
held parameter state.  Each next event reads the carried machine TuneState
again, before that event's candidate update.

Racc likewise consumes ``tuner_frequency_hz_()`` directly: the same preupdate
WPE getter/prior BEFORE statistics or outer tuning clamps.  The later band/tuner
path retains its own two clamps.  Each compiler's Racc state persists, and its
coefficient displacement is propagated through the same post-S accelerometer
relation. Guard/libm displacement remains a separate deployment obligation.

The exact shipping event remains the comparison shadow. Separate/FMA Racc
hypot/sqrt witnesses and accelerometer LDLT branches remain distinct. Native
binary32 ancestry for the Racc arithmetic, guard arithmetic and RAO libm paths
is still open, as are source-uniform supplies. No storage search is authorized.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_admitted_machine_measurement_supply_interleaved_prefix as LOWER
from tools.stability.ou3_alt_contraction import finite_admitted_machine_tunestate_interleaved_prefix as MTUNE
from tools.stability.ou3_alt_contraction import finite_racc_runtime as RACC
from tools.stability.ou3_alt_contraction import finite_measurement_runtime as MEAS
from tools.stability.ou3_alt_contraction import finite_source_bound_prediction_word as EXACTROOT
from tools.stability.ou3_alt_contraction import finite_machine_prediction_displacement as DISP
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE

QUALIFICATION='OU3_ALT_ADMITTED_MACHINE_RACC_SUPPLY_INTERLEAVER_V2'


def _exact_live_state(base:LOWER.State):
    if not isinstance(base,LOWER.State): raise TypeError('machine measurement-supply state required')
    return LOWER._preword(base).live.live.live


def _runtime(base:LOWER.State):
    if not isinstance(base,LOWER.State): raise TypeError('machine measurement-supply state required')
    runtime=LOWER._preword(base).runtime
    if not isinstance(runtime.racc_cfg,RACC.Config): raise TypeError('persistent shipping Racc config required')
    return runtime


def _mtune_result(lower:LOWER.ImuResult):
    if not isinstance(lower,LOWER.ImuResult): raise TypeError('machine measurement-supply result required')
    # measurement -> prediction-supply -> scheduler -> machine-prediction -> TuneState
    out=lower.lower.lower.lower.lower
    if not isinstance(out,MTUNE.ImuResult): raise TypeError('lower word lost same-event machine TuneState result')
    return out


def _subvec(a,b): return tuple(F(x)-F(y) for x,y in zip(a,b))
def _submat(a,b): return tuple(tuple(F(x)-F(y) for x,y in zip(ra,rb)) for ra,rb in zip(a,b))
def _supply(machine,exact): return DISP.Supply(_subvec(machine.z,exact.z),_submat(machine.covariance,exact.covariance))


@dataclass(frozen=True)
class State:
    base:LOWER.State
    separate_applied_sigma:F
    fma_applied_sigma:F
    separate_racc:RACC.State
    fma_racc:RACC.State
    entry_imu_steps:int
    racc_steps:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.base,LOWER.State): raise TypeError('machine measurement-supply state required')
        for n in ('separate_applied_sigma','fma_applied_sigma'):
            q=F(getattr(self,n)); object.__setattr__(self,n,q)
            if q<0: raise ValueError('nonnegative machine applied sigma required')
        if not isinstance(self.separate_racc,RACC.State) or not isinstance(self.fma_racc,RACC.State):
            raise TypeError('both persistent machine Racc states required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine Racc-supply qualification')
        if not isinstance(self.entry_imu_steps,int) or not isinstance(self.racc_steps,int) or self.entry_imu_steps<0 or self.racc_steps<0:
            raise ValueError('nonnegative integer Racc counters required')
        if self.base.measurement_steps-self.entry_imu_steps!=self.racc_steps:
            raise ValueError('machine Racc recurrence detached from lower measurement-event count')


@dataclass(frozen=True)
class ModeEvent:
    mode:str
    applied_sigma_before:F
    applied_sigma:F
    frequency:F
    racc:RACC.Result
    accelerometer:MEAS.MeasurementRuntimeResult
    accel_supply:DISP.Supply
    def __post_init__(self):
        if self.mode not in ('separate','fma'): raise ValueError('literal compiler mode required')
        object.__setattr__(self,'applied_sigma_before',F(self.applied_sigma_before))
        object.__setattr__(self,'applied_sigma',F(self.applied_sigma))
        object.__setattr__(self,'frequency',F(self.frequency))
        if self.applied_sigma_before<0 or self.applied_sigma<0 or self.frequency<=0:
            raise ValueError('valid machine Racc scalar ancestry required')
        if not isinstance(self.racc,RACC.Result) or not isinstance(self.accelerometer,MEAS.MeasurementRuntimeResult):
            raise TypeError('machine Racc and accelerometer results required')
        if not isinstance(self.accel_supply,DISP.Supply): raise TypeError('machine Racc accelerometer supply required')


@dataclass(frozen=True)
class ImuResult:
    state:State
    lower:LOWER.ImuResult
    separate:ModeEvent
    fma:ModeEvent
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.ImuResult):
            raise TypeError('machine Racc admitted result malformed')
        if self.separate.mode!='separate' or self.fma.mode!='fma': raise ValueError('machine Racc compiler histories crossed')
        if self.state.separate_applied_sigma!=self.separate.applied_sigma or self.state.fma_applied_sigma!=self.fma.applied_sigma:
            raise ValueError('persistent applied sigma detached from same-event Racc results')
        if self.state.separate_racc!=self.separate.racc.state or self.state.fma_racc!=self.fma.racc.state:
            raise ValueError('persistent Racc state detached from same-event result')


@dataclass(frozen=True)
class CompleteWord:
    state:State
    lower:LOWER.CompleteWord
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.CompleteWord):
            raise TypeError('machine Racc complete word requires lower complete word')
        if self.lower.state!=self.state.base: raise ValueError('machine Racc complete word detached from lower word')
        if self.state.racc_steps!=SOURCE.TRANSITIONS:
            raise ValueError('machine Racc recurrence not attached to all 600 admitted IMU edges')


def _mtune_state(base:LOWER.State):
    if not isinstance(base,LOWER.State): raise TypeError('machine measurement-supply state required')
    out=base.base.base.base.base
    if not isinstance(out,MTUNE.State): raise TypeError('measurement word lost machine TuneState predecessor')
    return out


def begin(base:LOWER.State):
    exact=_exact_live_state(base)
    machine=_mtune_state(base).machine
    return State(base,machine.sigma.separate,machine.sigma.fma,
                 exact.racc,exact.racc,base.measurement_steps,0)


def _applied_sigmas(state:State,mtune:MTUNE.ImuResult):
    """Read the pre-candidate memory, regardless of the pending commit bit."""
    previous=_mtune_state(state.base).machine
    if mtune.machine_boundary.arithmetic.before!=previous:
        raise ValueError('Racc sigma source detached from same sample-entry TuneState')
    return F(previous.sigma.separate),F(previous.sigma.fma)


def _racc_frequency(frequency):
    from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as WPEF
    if not isinstance(frequency,WPEF.StatisticsFrequencyResult):
        raise TypeError('same preupdate WPE/statistics frequency relation required')
    return F(frequency.external.stored.input_hz)


def _mode(*,mode,state_racc,old_sigma,new_sigma,frequency,lower_mode,live,segment,
          cfg,nominal_std,rao_witness,racc_sqrt,accel_ldlt,temperature_c,
          accel_alpha=1,accel_radius=F(2,5)):
    rr=RACC.step_from_applied_sigma(state_racc,cfg,live.guarded.guard,
        nominal_std=nominal_std,sigma_applied=new_sigma,
        preupdate_frequency=frequency,live=True,rao_witness=rao_witness,
        effective_sqrt=racc_sqrt)
    conditioning=EXACTROOT._accel_conditioning(temperature_c)
    accel=MEAS.accelerometer_from_held_guarded_racc(lower_mode.S_service.state,segment,
        live.guarded,conditioning,rr,ldlt=accel_ldlt,alpha=accel_alpha,radius=accel_radius)
    return ModeEvent(mode,F(old_sigma),F(new_sigma),F(frequency),rr,accel,
                     _supply(accel.state,live.accelerometer.state))


def imu_step(state:State,*,
             separate_racc_accel_ldlt:MEAS.SafeLDLT,fma_racc_accel_ldlt:MEAS.SafeLDLT,
             separate_rao_witness=None,separate_racc_sqrt=None,
             fma_rao_witness=None,fma_racc_sqrt=None,
             **kwargs):
    if not isinstance(state,State): raise TypeError('machine Racc State required')
    if not isinstance(separate_racc_accel_ldlt,MEAS.SafeLDLT) or not isinstance(fma_racc_accel_ldlt,MEAS.SafeLDLT):
        raise TypeError('both machine-Racc accelerometer LDLT branches required')
    restricted=kwargs.get('restricted'); temperature_c=kwargs.get('temperature_c')
    if restricted is None or temperature_c is None:
        raise TypeError('same admitted source and source-owned temperature required')
    runtime=_runtime(state.base)
    cfg=runtime.racc_cfg; nominal=runtime.nominal_racc_std
    # Persistent runtime configuration is already injected by the admitted lower
    # word. Accepting event-local Racc config here would create a detached second
    # execution surface, so the wrapper deliberately has no such parameters.
    lower=LOWER.imu_step(state.base,**kwargs)
    mtune=_mtune_result(lower); live=LOWER._live_result(lower.lower)
    sep_sigma,fma_sigma=_applied_sigmas(state,mtune)
    freq=mtune.lower
    common=dict(live=live,segment=restricted.segment,cfg=cfg,nominal_std=nominal,
                temperature_c=temperature_c,accel_alpha=kwargs.get('accel_alpha',1),
                accel_radius=kwargs.get('accel_radius',F(2,5)))
    sep=_mode(mode='separate',state_racc=state.separate_racc,
        old_sigma=state.separate_applied_sigma,new_sigma=sep_sigma,
        frequency=_racc_frequency(freq.separate_frequency),lower_mode=lower.separate,
        rao_witness=separate_rao_witness,racc_sqrt=separate_racc_sqrt,
        accel_ldlt=separate_racc_accel_ldlt,**common)
    fma=_mode(mode='fma',state_racc=state.fma_racc,
        old_sigma=state.fma_applied_sigma,new_sigma=fma_sigma,
        frequency=_racc_frequency(freq.fma_frequency),lower_mode=lower.fma,
        rao_witness=fma_rao_witness,racc_sqrt=fma_racc_sqrt,
        accel_ldlt=fma_racc_accel_ldlt,**common)
    nxt=State(lower.state,sep_sigma,fma_sigma,sep.racc.state,fma.racc.state,
              state.entry_imu_steps,state.racc_steps+1)
    return ImuResult(nxt,lower,sep,fma)


def mag_step(state:State,**kwargs):
    if not isinstance(state,State): raise TypeError('machine Racc State required')
    base,event=LOWER.mag_step(state.base,**kwargs)
    return State(base,state.separate_applied_sigma,state.fma_applied_sigma,
                 state.separate_racc,state.fma_racc,state.entry_imu_steps,state.racc_steps),event


def set_hold(state:State,*,hold):
    if not isinstance(state,State): raise TypeError('machine Racc State required')
    base,event=LOWER.set_hold(state.base,hold=hold)
    return State(base,state.separate_applied_sigma,state.fma_applied_sigma,
                 state.separate_racc,state.fma_racc,state.entry_imu_steps,state.racc_steps),event


def complete(state:State): return CompleteWord(state,LOWER.complete(state.base))


def readiness():
    low=LOWER.readiness(); rr=RACC.readiness()
    return {
      'machine_measurement_supply_word_consumed':low['accelerometer_measurement_propagates_machine_prediction_supply'],
      'persistent_Racc_config_and_nominal_std_consumed_from_admitted_runtime':True,
      'raw_applied_sigma_carried_separately_from_stationary_Sigma_aw':True,
      'pending_machine_boundary_updates_same_mode_applied_sigma_before_Racc':True,
      'nonpending_boundary_preserves_machine_applied_sigma':False,
      'every_Racc_event_reads_sample_entry_machine_TuneState_sigma':True,
      'Racc_frequency_bypasses_statistics_and_outer_tuning_clamps':True,
      'persistent_separate_and_FMA_Racc_states_attached':True,
      'same_guard_drives_exact_and_both_machine_Racc_histories':rr['same_guard_excess_RMS_drives_Racc'],
      'same_mode_preupdate_WPE_frequency_drives_machine_Racc':True,
      'machine_Racc_coefficient_displacement_attached':True,
      'machine_Racc_effect_propagated_through_accelerometer':True,
      'full_joint24_and_21x21_supply_retained_after_machine_Racc_accelerometer':True,
      'machine_Racc_binary32_hypot_sqrt_correspondence_closed':False,
      'machine_guard_private_Mahony_and_conditioning_float_correspondence_closed':False,
      'source_uniform_machine_event_supply_bound_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
