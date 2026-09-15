"""Persistent compiler CORE/control continuation for the conditional finite word.

The strongest wrapper injects each preceding CORE into the prediction roots and
propagates it through covariance floor, S service, Racc and accelerometer. This
module then executes that mode's finite watchdog/gravity suffix and its MAG/HOLD
runtime algorithms. Full state, covariance, chart, magnetic memory and control
feed the next event. No generic callback or exact-shadow restart is used.

These are conditional finite-real event graphs with machine tuner/frontend
operands. Source-uniform totality, target rounding/solver correspondence and
machine startup CORE qualification remain separate. ``observe_*`` preserve the
older diagnostic-only relations so their failure cannot masquerade as this
concrete continuation; they never qualify completion.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_live_interleave as LIVE
from tools.stability.ou3_alt_contraction import finite_tilt_watchdog as WATCH
from tools.stability.ou3_alt_contraction import finite_tilt_reset_runtime as RESET
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_machine_accel_guard_binary32 as GUARD


@dataclass(frozen=True)
class Observation:
    kind: str
    imu_ordinal: int
    separate_predecessor_connected: bool | None
    fma_predecessor_connected: bool | None
    separate_event: object | None = None
    fma_event: object | None = None


@dataclass(frozen=True)
class State:
    separate: CORE.State
    fma: CORE.State
    observations: tuple[Observation, ...] = ()
    missing: tuple[str, ...] = ()
    separate_runtime: LIVE.State | None = None
    fma_runtime: LIVE.State | None = None

    def __post_init__(self):
        if not isinstance(self.separate, CORE.State) or not isinstance(self.fma, CORE.State):
            raise TypeError('both full compiler CORE states required')
        if self.separate.reference != self.fma.reference:
            raise ValueError('compiler CORE states must share the physical/bias endpoint')
        if (self.separate_runtime is None)!=(self.fma_runtime is None):
            raise ValueError('both compiler runtime control histories required')
        for mode in ('separate','fma'):
            runtime=getattr(self,mode+'_runtime')
            if runtime is not None:
                if not isinstance(runtime,LIVE.State) or runtime.live.live.mekf!=getattr(self,mode):
                    raise ValueError('compiler CORE detached from persistent runtime control')
        object.__setattr__(self, 'observations', tuple(self.observations))
        object.__setattr__(self, 'missing', tuple(sorted(set(self.missing))))
        n = 0
        for event in self.observations:
            if not isinstance(event, Observation) or event.kind not in ('imu', 'mag', 'hold'):
                raise TypeError('literal IMU/MAG/HOLD continuation observation required')
            n += event.kind == 'imu'
            if event.imu_ordinal != n:
                raise ValueError('machine CORE event ledger detached from source ordinal')
        if n > SOURCE.TRANSITIONS:
            raise ValueError('machine CORE continuation exceeds 600 IMU events')

    @property
    def imu_steps(self):
        return sum(event.kind == 'imu' for event in self.observations)


def begin(separate: CORE.State, fma: CORE.State):
    """Conditional root only; startup machine-state qualification is external."""
    return State(separate, fma)


def begin_from_interleaved(runtime:LIVE.State):
    """Conditional startup attachment; do not infer target seed correspondence."""
    if not isinstance(runtime,LIVE.State):
        raise TypeError('startup-produced interleaved runtime required')
    core=runtime.live.live.mekf
    return State(core,core,separate_runtime=runtime,fma_runtime=runtime)


def _machine_guarded(exact:SENSOR.GuardedImuSample,machine:GUARD.Result):
    """Bind the already-executed machine conditioning to the reset operand.

    Only the guarded output is read by the reused tilt/reset algorithm. The
    copied exact guard metadata is not a second guard recurrence or certificate.
    """
    if not isinstance(exact,SENSOR.GuardedImuSample) or not isinstance(machine,GUARD.Result):
        raise TypeError('same-event exact packet and executed machine guard required')
    GUARD.api_vec(exact.raw_accel_body,machine.raw_acc,'same reset packet')
    return SENSOR.GuardedImuSample(exact.raw,replace(exact.guard,output=machine.conditioned_acc))


def _imu_suffix(runtime,template,core,proxy,sample,dt,*,reset,gravity):
    """Execute the literal per-compiler watchdog and ungauged gravity suffix."""
    if not isinstance(runtime,LIVE.State) or not isinstance(template,LIVE.State):
        raise TypeError('persistent compiler control and same-event runtime template required')
    if core.reference!=template.live.live.mekf.reference:
        raise ValueError('machine suffix detached from the source-owned endpoint')
    if not isinstance(reset,dict): raise TypeError('same-mode reset witness dictionary required')
    wd=WATCH.step_over_limit(runtime.live.watchdog,dt=dt,
                            over_limit=RESET.watchdog_over_limit(core))
    if wd.fired:
        operands=dict(reset)
        sigmas={k:operands.pop(k) for k in ('tilt_sigma','yaw_sigma') if k in operands}
        witness=RESET.preserve_yaw_witness(core,sample,**operands)
        core=WATCH.preserve_yaw_reset(core,sample,witness,**sigmas)
    elif reset:
        raise ValueError('nonfiring machine watchdog consumes no reset witnesses')
    magnetic=runtime.magnetic
    if isinstance(magnetic,LIVE.MAG.UngaugedLiveState):
        if gravity is None:
            raise ValueError('ungauged machine IMU requires its own gravity-gate witnesses')
        gate=LIVE.MAG.GRAVITY.imu_step(magnetic.startup.word.gate,magnetic.memory.cfg.gravity,
            q_proxy_bw=LIVE.P.quat_conj(core.q_hat),acc_body=sample.raw_accel_body,
            gyro_body=sample.raw_gyro_body,dt=dt,**gravity)
        prefix=replace(magnetic.startup.word.mag,proxy=proxy)
        word=replace(magnetic.startup.word,gate=gate.state,mag=prefix)
        magnetic=replace(magnetic,startup=replace(magnetic.startup,word=word))
    elif gravity is not None:
        raise ValueError('gauged machine word consumes no initial gravity witnesses')
    # MAG/HOLD read the CORE, proxy and these persistent control fields only.
    # Other template fields continue to be owned by the outer tuner graph.
    imu=replace(template.live.live,mekf=core,
                tuner=replace(template.live.live.tuner,vertical=proxy))
    return LIVE.State(LIVE.LIVE.State(imu,wd.state),magnetic,runtime.clock,runtime.schedule)


def advance_imu(state:State,*,separate_predecessor,fma_predecessor,
                separate_after_accel,fma_after_accel,template,exact_guarded,machine_guard,
                separate_proxy,fma_proxy,dt,separate_reset=None,fma_reset=None,
                separate_gravity=None,fma_gravity=None,
                separate_measurement,fma_measurement,separate_racc,fma_racc):
    """Substitute actual compiler predecessors through the complete core suffix.

    The lower producers execute prediction/floor/S/accelerometer from these
    predecessors. This function executes the remaining control-dependent core
    maps and persists their results for the next IMU/MAG/HOLD operation.
    Arithmetic/solver witnesses remain conditional; none are inferred here.
    """
    if state.separate_runtime is None or state.missing:
        raise ValueError('actual machine continuation cannot extend a diagnostic-only prefix')
    if (separate_predecessor,fma_predecessor)!=(state.separate,state.fma):
        raise ValueError('machine prediction discarded the preceding full CORE successor')
    for mode,pre,acc,meas,racc in (
        ('separate',separate_predecessor,separate_after_accel,separate_measurement,separate_racc),
        ('fma',fma_predecessor,fma_after_accel,fma_measurement,fma_racc)):
        if meas.mode!=mode or racc.mode!=mode or meas.prediction.machine_predecessor!=pre or racc.accelerometer.state!=acc:
            raise ValueError('same-mode prediction/measurement/Racc event chain required')
    sample=_machine_guarded(exact_guarded,machine_guard)
    common=dict(template=template,sample=sample,dt=dt)
    sr=_imu_suffix(state.separate_runtime,core=separate_after_accel,proxy=separate_proxy,
                   reset={} if separate_reset is None else separate_reset,
                   gravity=separate_gravity,**common)
    fr=_imu_suffix(state.fma_runtime,core=fma_after_accel,proxy=fma_proxy,
                   reset={} if fma_reset is None else fma_reset,
                   gravity=fma_gravity,**common)
    record=Observation('imu',state.imu_steps+1,True,True,
                       (separate_measurement,separate_racc,sr),
                       (fma_measurement,fma_racc,fr))
    return State(sr.live.live.mekf,fr.live.live.mekf,state.observations+(record,),(),sr,fr)


def advance_hold(state:State,*,hold):
    if state.separate_runtime is None or state.missing:
        raise ValueError('machine HOLD requires an actual compiler control prefix')
    sr=LIVE.set_hold(state.separate_runtime,hold=hold)
    fr=LIVE.set_hold(state.fma_runtime,hold=hold)
    record=Observation('hold',state.imu_steps,True,True,sr,fr)
    return State(sr.state.live.live.mekf,fr.state.live.live.mekf,
                 state.observations+(record,),(),sr.state,fr.state),(sr,fr)


def advance_mag(state:State,*,separate=None,fma=None,**kwargs):
    """Execute both full magnetic/control words at the same physical endpoint."""
    if state.separate_runtime is None or state.missing:
        raise ValueError('machine MAG requires an actual compiler control prefix')
    sw={} if separate is None else dict(separate)
    fw={} if fma is None else dict(fma)
    forbidden={'residual_body','packet_id','begun','have_last_imu','origin'}
    if forbidden & (set(sw)|set(fw)):
        raise TypeError('compiler magnetic witnesses cannot override common source operands')
    sr=LIVE.mag_step(state.separate_runtime,**dict(kwargs,**sw))
    fr=LIVE.mag_step(state.fma_runtime,**dict(kwargs,**fw))
    record=Observation('mag',state.imu_steps,True,True,sr,fr)
    return State(sr.state.live.live.mekf,fr.state.live.live.mekf,
                 state.observations+(record,),(),sr.state,fr.state),(sr,fr)


def observe_imu(state: State, *, separate_predecessor: CORE.State,
                fma_predecessor: CORE.State, separate_after_accel: CORE.State,
                fma_after_accel: CORE.State):
    """Retain outputs of the actual local producer and audit the next input.

    This records failures instead of discarding a legal source history.  In
    particular, a mismatch is not a new source-admission restriction.  Final
    completion is refused until actual machine event algorithms close it.
    """
    if not isinstance(state, State):
        raise TypeError('persistent machine CORE continuation required')
    inputs = (separate_predecessor, fma_predecessor)
    outputs = (separate_after_accel, fma_after_accel)
    if not all(isinstance(x, CORE.State) for x in inputs + outputs):
        raise TypeError('full predecessor and post-accelerometer CORE states required')
    before = state.separate.reference
    after = separate_after_accel.reference
    if any(x.reference != before for x in inputs):
        raise ValueError('machine prediction precursor detached from current physical endpoint')
    if fma_after_accel.reference != after:
        raise ValueError('machine successors use different physical/bias histories')
    if after.time - before.time != SOURCE.DT:
        raise ValueError('machine IMU continuation must consume the canonical source step')
    for name in ('history_id', 'bias_root', 'bias_family', 'live_origin'):
        if getattr(before, name) != getattr(after, name):
            raise ValueError('machine CORE continuation restarted ' + name)
    connected = (state.separate == separate_predecessor, state.fma == fma_predecessor)
    missing = set(state.missing)
    for mode, ok in zip(('separate', 'fma'), connected):
        if not ok:
            missing.add(mode + '_machine_CORE_successor_not_consumed')
    missing.add('machine_post_accelerometer_event_and_control_suffix_not_attached')
    observation = Observation('imu', state.imu_steps + 1, *connected)
    return State(separate_after_accel, fma_after_accel,
                 state.observations + (observation,), tuple(missing))


def observe_async(state: State, *, kind: str):
    if not isinstance(state, State) or kind not in ('mag', 'hold'):
        raise TypeError('persistent machine CORE and literal MAG/HOLD kind required')
    # Retain the last known compiler CORE; do not replace it with the shadow's
    # result or falsely turn the unrepresented machine event into an identity.
    observation = Observation(kind, state.imu_steps, None, None)
    return State(state.separate, state.fma, state.observations + (observation,),
                 state.missing + (kind + '_machine_CORE_and_control_successor_not_attached',))


def require_complete(state: State):
    if not isinstance(state, State):
        raise TypeError('persistent full machine CORE continuation required')
    if state.imu_steps != SOURCE.TRANSITIONS:
        raise ValueError('machine CORE continuation lacks all 600 IMU events')
    if state.missing:
        raise ValueError('machine CORE word is incomplete: ' + ', '.join(state.missing))
    if state.separate_runtime is None:
        raise ValueError('machine CORE event/control successor algorithms not attached')
    if any(event.separate_event is None or event.fma_event is None for event in state.observations):
        raise ValueError('machine CORE event/control successor algorithms not attached')
    return True


def readiness():
    return {
        'compiler_specific_full_CORE_local_successors_persist': True,
        'next_local_prediction_predecessor_checked_against_carried_CORE': True,
        'joint24_full_21_covariance_nominal_attitude_atlas_and_reference_checked': True,
        'missing_machine_MAG_HOLD_and_IMU_suffix_explicit': True,
        'conditional_IMU_watchdog_MAG_and_HOLD_successors_execute_from_carried_CORE':True,
        'source_histories_restricted_to_zero_machine_discrepancy': False,
        'full_machine_event_and_control_successor_algorithms_attached': True,
        'source_uniform_complete_600_step_word_qualified': False,
        'storage_search_allowed': False,
    }
