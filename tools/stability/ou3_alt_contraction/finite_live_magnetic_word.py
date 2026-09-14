"""Same-history startup/Live magnetic events, including default calibration.

The one qualified raw packet enters continuous accumulation, reference learning,
coupled hard-iron/reference application, and finally the MEKF in shipping order.
The actual default continuous estimator is NOT disabled to simplify the word.
This is a conditional finite-real composition, not source-uniform stability or
binary32/libm/solver qualification. See docs/ou3-alt-live-magnetic-word.md.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_mag_tilt_frame as TILT
from tools.stability.ou3_alt_contraction import finite_mag_tuner_default as TUNER
from tools.stability.ou3_alt_contraction import finite_mag_startup_word as START
from tools.stability.ou3_alt_contraction import finite_mag_startup_ready as READY
from tools.stability.ou3_alt_contraction import finite_mag_startup_source as SOURCE
from tools.stability.ou3_alt_contraction import finite_mag_startup_physical_bridge as BRIDGE
from tools.stability.ou3_alt_contraction import finite_mag_source_qualification as QUAL
from tools.stability.ou3_alt_contraction import finite_mag_gravity_gate as GRAVITY
from tools.stability.ou3_alt_contraction import finite_mag_runtime as MAG
from tools.stability.ou3_alt_contraction import finite_mag_reference_runtime as REF
from tools.stability.ou3_alt_contraction import finite_mag_reference_events as EVENTS
from tools.stability.ou3_alt_contraction import finite_mag_bias_gate as GATE
from tools.stability.ou3_alt_contraction import finite_live_async_mag as ASYNC
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as VERT
from tools.stability.ou3_alt_contraction import finite_continuous_mag_runtime as HI

R = M.rational


@dataclass(frozen=True)
class Config:
    gate: GATE.Config = GATE.Config()
    gravity: GRAVITY.Config = GRAVITY.Config()
    tuner: TUNER.Config = TUNER.Config(min_samples=128, min_window=F(15))
    sigma_internal: tuple = (F(3, 10),)*3
    refinement_enabled: bool = True
    refinement_start: F = F(90)
    refinement_window: F = F(30)
    sample_dt: F = F(1, 200)
    min_reference_norm: F = F(1, 1000)
    continuous_enabled: bool = True
    continuous: HI.Config = HI.Config()
    apply_fraction: F = F(1)
    slew_tau: F = F(45)

    def __post_init__(self):
        if not isinstance(self.gate, GATE.Config) or not isinstance(self.gravity, GRAVITY.Config):
            raise TypeError('persistent shipping magnetic/gravity configuration required')
        if not isinstance(self.tuner, TUNER.Config) or not isinstance(self.continuous, HI.Config):
            raise TypeError('persistent acquisition/continuous calibration configuration required')
        if (self.gate.with_mag, self.gate.mag_delay) != (self.gravity.with_mag, self.gravity.mag_delay):
            raise ValueError('inner and outer magnetic configuration must share delay/enable')
        if self.tuner.quality_weighting or self.tuner.estimate_hard_iron:
            raise NotImplementedError('optional startup weighting/hard-iron solve remains open')
        for name in ('refinement_enabled', 'continuous_enabled'):
            if not isinstance(getattr(self, name), bool):
                raise TypeError('literal wrapper magnetic branch required')
        for name in ('refinement_start', 'refinement_window', 'sample_dt',
                     'min_reference_norm', 'apply_fraction', 'slew_tau'):
            object.__setattr__(self, name, R(getattr(self, name)))
        if min(self.refinement_start, self.refinement_window, self.min_reference_norm) < 0 or self.sample_dt <= 0:
            raise ValueError('valid magnetic clock and reference configuration required')
        sig = tuple(M.vec(self.sigma_internal, 3))
        if min(sig) < 0:
            raise ValueError('nonnegative constructor magnetic sigma required')
        object.__setattr__(self, 'sigma_internal', sig)


@dataclass(frozen=True)
class Memory:
    cfg: Config
    model: SOURCE.Model
    history_id: str
    producer_root: str
    continuous: HI.State = HI.State()
    last_hi_time: F | None = None
    applied: HI.Applied = HI.Applied()

    def __post_init__(self):
        if not isinstance(self.cfg, Config) or not isinstance(self.model, SOURCE.Model):
            raise TypeError('bound configuration and fixed physical magnetic model required')
        if any(not isinstance(x, str) or not x for x in (self.history_id, self.producer_root)):
            raise ValueError('persistent physical and reference-producer roots required')
        if not isinstance(self.continuous, HI.State) or not isinstance(self.applied, HI.Applied):
            raise TypeError('persistent continuous magnetic state required')
        if self.applied.startup_bias != HI.ZERO:
            raise NotImplementedError('nondefault startup hard-iron fit remains open')
        if self.last_hi_time is not None:
            t = R(self.last_hi_time)
            if t < 0:
                raise ValueError('nonnegative last hard-iron sample time required')
            object.__setattr__(self, 'last_hi_time', t)


def _source(memory, physical, residual_body, packet_id):
    # CORE.Reference has a physical history token; startup endpoints may be
    # PhysicalKinematics without one, in which case the separately held token
    # is required. Neither token alone qualifies a BRMM history.
    if hasattr(physical, 'history_id') and physical.history_id != memory.history_id:
        raise ValueError('magnetic source detached from current physical history')
    sample = BRIDGE.sample(physical, memory.history_id, memory.model, residual_body, packet_id)
    q = QUAL.qualify_squared(sample)  # fixed named envelope, not a per-call bound
    QUAL.assert_source_qualified(q)
    return q


def _tilt(proxy, *, q_norm, yaw_half):
    if not isinstance(proxy, VERT.State):
        raise TypeError('same persistent private observer required')
    if not isinstance(q_norm, TILT.SqrtWitness) or q_norm.value <= F(1, 10**6):
        raise NotImplementedError('tiny/nonfinite private-quaternion branch remains open')
    return TILT.yaw_removed(proxy.q, q_norm=q_norm, yaw_half=yaw_half)


def _accumulate(memory, source, proxy, tilt, *, decay, eigen_success, spectrum):
    cfg = memory.cfg
    if not cfg.continuous_enabled or not proxy.initialized:
        if any(x is not None for x in (decay, eigen_success, spectrum)):
            raise ValueError('disabled/uninitialized continuous path consumes no witnesses')
        return memory, None
    if tilt is None:
        raise TypeError('continuous calibration needs the same private proxy tilt')
    t = source.physical.time
    dt = t-memory.last_hi_time if memory.last_hi_time is not None and t > memory.last_hi_time else cfg.sample_dt
    out = HI.update(memory.continuous, cfg.continuous, dt=dt,
                    rotation=CORE.rotation(tilt.q_tilt), raw_body=source.raw_body,
                    decay=decay, eigen_success=eigen_success, spectrum=spectrum)
    return replace(memory, continuous=out.state, last_hi_time=t), out


@dataclass(frozen=True)
class StartupState:
    word: START.State
    memory: Memory
    ready: READY.Result | None = None

    def __post_init__(self):
        if not isinstance(self.word, START.State) or not isinstance(self.memory, Memory):
            raise TypeError('persistent startup word and continuous magnetic memory required')
        if self.word.source_history_id not in (None, self.memory.history_id):
            raise ValueError('startup magnetic physical ancestry detached')
        if self.word.source_model_root not in (None, self.memory.model.model_root):
            raise ValueError('startup magnetic source-model ancestry detached')
        if self.ready is not None:
            if not isinstance(self.ready, READY.Result) or self.ready.tuner_state != self.word.mag.tuner:
                raise ValueError('startup ready state detached from carried acquisition')
            if self.ready.active_reference.root != self.memory.producer_root:
                raise ValueError('startup active magnetic producer detached')


def begin_startup(word, cfg, model, history_id, producer_root):
    """Initialize the magnetic supplement before ANY magnetic call, not at Live."""
    if not isinstance(word, START.State):
        raise TypeError('initial startup magnetic/gravity word required')
    if (word.mag.tuner != TUNER.State() or word.mag.last_mag_time is not None
            or word.source_model_root is not None or word.source_history_id is not None):
        raise ValueError('continuous magnetic statistics cannot be restarted mid-startup')
    return StartupState(word, Memory(cfg, model, history_id, producer_root))


@dataclass(frozen=True)
class StartupResult:
    state: StartupState
    startup: START.MagResult | None
    continuous: HI.Result | None
    qualification: QUAL.Qualification | None


def startup_call(state, physical, *, residual_body=None, packet_id=None,
                 begun=True, have_last_imu=True, proxy_q_norm=None, proxy_yaw_half=None,
                 hi_decay=None, hi_eigen_success=None, hi_spectrum=None,
                 mag_norm=None, mean_norm=None, horizontal_sqrt=None, ready_yaw_half=None):
    """Outer startup updateMag, including accumulation AHEAD of startup gates."""
    if not isinstance(state, StartupState):
        raise TypeError('persistent startup magnetic product state required')
    cfg = state.memory.cfg
    if not isinstance(begun, bool) or not isinstance(have_last_imu, bool):
        raise TypeError('literal startup wrapper branches required')
    operands = (residual_body, packet_id, proxy_q_norm, proxy_yaw_half, hi_decay,
                hi_eigen_success, hi_spectrum, mag_norm, mean_norm, horizontal_sqrt, ready_yaw_half)
    if not begun or not cfg.gate.with_mag or physical.time < cfg.gate.mag_delay:
        if any(x is not None for x in operands):
            raise ValueError('outer-gated startup call consumes no source/arithmetic operands')
        return StartupResult(state, None, None, None)
    qualified = _source(state.memory, physical, residual_body, packet_id)
    source = qualified.sample
    proxy = state.word.mag.proxy
    need_tilt = cfg.continuous_enabled and proxy.initialized
    # Startup admission is computed independently of continuous accumulation;
    # computing its predicate early has no side effect on either operand graph.
    admission = GRAVITY.startup_mag_admission(state.word.gate, cfg.gravity,
        wrapper_time=physical.time, begun=begun, have_last_imu=have_last_imu,
        mag_ref_set=state.ready is not None)
    tilt = _tilt(proxy, q_norm=proxy_q_norm, yaw_half=proxy_yaw_half) if need_tilt or admission.admitted else None
    if tilt is None and (proxy_q_norm is not None or proxy_yaw_half is not None):
        raise ValueError('unread private tilt consumes no quaternion witnesses')
    memory, continuous = _accumulate(state.memory, source, proxy, tilt,
        decay=hi_decay, eigen_success=hi_eigen_success, spectrum=hi_spectrum)
    out = START.update_mag_source_call(state.word, cfg.gravity, cfg.tuner, source,
        begun=begun, have_last_imu=have_last_imu, mag_ref_set=state.ready is not None,
        sample_dt=cfg.sample_dt, boat_q_norm=proxy_q_norm if admission.admitted else None,
        yaw_half=proxy_yaw_half if admission.admitted else None,
        mag_norm=mag_norm, mean_norm=mean_norm, horizontal_sqrt=horizontal_sqrt)
    out = replace(out, qualification=qualified)
    ready = state.ready
    due = (ready is None and out.magnetic is not None and out.magnetic.tuner_step.returned_ready
           and M.dot(out.state.mag.tuner.world_reference, out.state.mag.tuner.world_reference) > cfg.min_reference_norm**2)
    if due:
        ready = READY.transition(out.state.mag.tuner, READY.Config(cfg.sigma_internal, memory.producer_root),
                                 yaw_half=ready_yaw_half)
    elif ready_yaw_half is not None:
        raise ValueError('no provisional reference write consumes no gauge witness')
    return StartupResult(StartupState(out.state, memory, ready), out, continuous, qualified)


@dataclass(frozen=True)
class LiveState:
    memory: Memory
    tuner: TUNER.State
    last_mag_time: F | None
    active: REF.State
    control: GATE.State
    refinement_started: bool = False
    refinement_done: bool = False
    refinement_time: F | None = None

    def __post_init__(self):
        if not isinstance(self.memory, Memory) or not isinstance(self.tuner, TUNER.State):
            raise TypeError('carried magnetic memory and tuner required')
        if not isinstance(self.active, REF.State) or not isinstance(self.control, GATE.State):
            raise TypeError('persistent active reference and bias gate required')
        if self.active.root != self.memory.producer_root or self.active.model.sigma_internal != self.memory.cfg.sigma_internal:
            raise ValueError('active reference detached from same producer/constructor Rmag')
        for name in ('last_mag_time', 'refinement_time'):
            if getattr(self, name) is not None:
                object.__setattr__(self, name, R(getattr(self, name)))
        if any(not isinstance(x, bool) for x in (self.refinement_started, self.refinement_done)):
            raise TypeError('literal refinement state required')
        if self.last_mag_time is not None and self.last_mag_time < 0:
            raise ValueError('nonnegative magnetic acquisition clock required')
        if self.refinement_time is not None and self.refinement_time < 0:
            raise ValueError('nonnegative refinement completion clock required')
        if self.refinement_done != (self.refinement_time is not None):
            raise ValueError('refinement completion and its timestamp must agree')
        if self.refinement_done and not self.refinement_started:
            raise ValueError('refinement cannot finish before starting')


def enter_live(state: StartupState):
    if not isinstance(state, StartupState) or state.ready is None:
        raise NotImplementedError('nongauged Live entry remains a separate shipping branch')
    # begin() calls setAccBiasHold(mag_refine_enabled). There have been no inner
    # MEKF updateMag calls before the outer handoff.
    control = GATE.State(0, None, True,
                         state.memory.cfg.gate.with_mag and state.memory.cfg.refinement_enabled)
    return LiveState(state.memory, state.word.mag.tuner, state.word.mag.last_mag_time,
                     state.ready.active_reference, control)


def _overwrite_absolute_yaw(core, mean, gauge_half, q_norm, yaw_half):
    if not isinstance(gauge_half, TILT.YawHalfWitness) or (gauge_half.c, gauge_half.s) != tuple(mean[:2]):
        raise ValueError('refinement yaw gauge detached from same accepted mean')
    tilt = TILT.yaw_removed(P.quat_conj(core.q_hat), q_norm=q_norm, yaw_half=yaw_half)
    qa = (gauge_half.cos_half, F(0), F(0), -gauge_half.sin_half)
    q_hat = tuple(P.quat_conj(P.quat_mul(qa, tilt.q_tilt)))
    z = list(core.z)
    z[:3] = CORE.cayley(P.quat_mul(core.reference.q_world_to_body, P.quat_conj(q_hat)))
    # set_quaternion_boat writes qref and zero local attitude slots ONLY. It
    # neither resets P nor reanchors any physical/motion/bias state.
    return CORE.State(core.mode, tuple(z), core.covariance, q_hat, core.reference)


@dataclass(frozen=True)
class LiveResult:
    state: LiveState
    filter: CORE.State
    qualification: QUAL.Qualification | None
    continuous: HI.Result | None
    refinement: TUNER.StepResult | None
    refinement_filter: CORE.State
    hold_release: GATE.GateResult | None
    continuous_apply: HI.ApplyResult | None
    measurement: ASYNC.Result | None
    effective_residual: tuple | None


def live_call(state, core, proxy, *, residual_body=None, packet_id=None,
              proxy_q_norm=None, proxy_yaw_half=None,
              hi_decay=None, hi_eigen_success=None, hi_spectrum=None,
              mag_norm=None, mean_norm=None, horizontal_sqrt=None,
              gauge_half=None, mekf_q_norm=None, mekf_yaw_half=None,
              apply_decay=None, apply_new_norm=None, apply_anchor_norm=None,
              ldlt=None, alpha=1, radius=None):
    """One Live call at the CURRENT physical endpoint; never advances IMU time."""
    if not isinstance(state, LiveState) or not isinstance(core, CORE.State) or not isinstance(proxy, VERT.State):
        raise TypeError('persistent Live magnetic/core/private observer states required')
    cfg = state.memory.cfg
    t = core.reference.time
    operands = (residual_body, packet_id, proxy_q_norm, proxy_yaw_half, hi_decay,
        hi_eigen_success, hi_spectrum, mag_norm, mean_norm, horizontal_sqrt,
        gauge_half, mekf_q_norm, mekf_yaw_half, apply_decay, apply_new_norm, apply_anchor_norm, ldlt, radius)
    if not cfg.gate.with_mag or t < cfg.gate.mag_delay:
        if any(x is not None for x in operands) or alpha != 1:
            raise ValueError('outer-gated Live call consumes no source/arithmetic operands')
        return LiveResult(state, core, None, None, None, core, None, None, None, None)
    qualified = _source(state.memory, core.reference, residual_body, packet_id)
    source = qualified.sample
    refining = cfg.refinement_enabled and not state.refinement_done and t >= cfg.refinement_start
    need_tilt = refining or (cfg.continuous_enabled and proxy.initialized)
    tilt = _tilt(proxy, q_norm=proxy_q_norm, yaw_half=proxy_yaw_half) if need_tilt else None
    if tilt is None and (proxy_q_norm is not None or proxy_yaw_half is not None):
        raise ValueError('unread private tilt consumes no quaternion witnesses')
    memory, continuous = _accumulate(state.memory, source, proxy, tilt,
        decay=hi_decay, eigen_success=hi_eigen_success, spectrum=hi_spectrum)
    active, control = state.active, state.control
    tuner, last = state.tuner, state.last_mag_time
    started, done, done_time = state.refinement_started, state.refinement_done, state.refinement_time
    refinement = release = None
    if refining:
        if not started:
            tuner, last, started = TUNER.State(), None, True
        dt = t-last if last is not None and t > last else cfg.sample_dt
        last = t  # persists on acquisition rejection
        tc = replace(cfg.tuner, min_window=cfg.refinement_window)
        corrected = tuple(source.raw_body[i]-memory.applied.total_bias[i] for i in range(3))
        refinement = TUNER.step(tuner, tc, q_tilt_bw=tilt.q_tilt, mag_body=corrected, dt=dt,
            q_norm=None if tuner.ready else TUNER.SqrtWitness(1, 1), mag_norm=mag_norm,
            mean_norm=mean_norm, horizontal_sqrt=horizontal_sqrt)
        tuner = refinement.state
        if refinement.returned_ready and M.dot(tuner.world_reference, tuner.world_reference) > cfg.min_reference_norm**2:
            active = EVENTS.rewrite(active, EVENTS.WriteWitness('live_refinement',
                MAG.Model(tuner.world_reference, cfg.sigma_internal), memory.producer_root))
            core = _overwrite_absolute_yaw(core, tuner.mean, gauge_half, mekf_q_norm, mekf_yaw_half)
            done, done_time = True, t
            release = GATE.set_hold(control, core, cfg.gate, hold=False, live=True)
            control, core = release.state, release.filter_state
        elif any(x is not None for x in (gauge_half, mekf_q_norm, mekf_yaw_half)):
            raise ValueError('unfinished refinement consumes no MEKF yaw-write witnesses')
    elif any(x is not None for x in (mag_norm, mean_norm, horizontal_sqrt, gauge_half, mekf_q_norm, mekf_yaw_half)):
        raise ValueError('inactive refinement consumes no acquisition/yaw operands')
    after_refinement = core
    applied = HI.apply(memory.applied, memory.continuous, active.model.world_reference,
        time=t, sample_dt=cfg.sample_dt, fraction=cfg.apply_fraction, slew_tau=cfg.slew_tau,
        enabled=cfg.continuous_enabled, refinement_enabled=cfg.refinement_enabled,
        refinement_done=done, proxy_min_norm=cfg.min_reference_norm,
        decay=apply_decay, new_norm=apply_new_norm, anchor_norm=apply_anchor_norm)
    memory = replace(memory, applied=applied.state)
    if applied.wrote_reference:
        active = EVENTS.rewrite(active, EVENTS.WriteWitness('continuous_hard_iron',
            MAG.Model(applied.reference, cfg.sigma_internal), memory.producer_root))
    corrected = tuple(source.raw_body[i]-memory.applied.total_bias[i] for i in range(3))
    field_difference = tuple(source.model.world_field[i]-active.model.world_reference[i] for i in range(3))
    rotated_difference = MAG.SENSOR.q_rotate(core.reference.q_world_to_body, field_difference)
    nu = tuple(rotated_difference[i]+source.model.hard_iron_body[i]
               -memory.applied.total_bias[i]+source.residual_body[i] for i in range(3))
    packet = MAG.Sample(core.reference, corrected, nu, active.model)
    measured = ASYNC.update_mag_call(ASYNC.State(core, control, active), cfg.gate,
        time=t, live=True, sample=packet, ldlt=ldlt, alpha=alpha, radius=radius)
    nxt = LiveState(memory, tuner, last, active, measured.state.control, started, done, done_time)
    return LiveResult(nxt, measured.state.filter, qualified, continuous, refinement,
                      after_refinement, release, applied, measured, nu)


def readiness():
    return {
        'same_physical_packet_to_startup_continuous_and_Live_magnetic_paths': True,
        'continuous_statistics_carried_from_startup_not_restarted_at_Live': True,
        'refinement_reset_clock_private_tilt_reference_yaw_hold_order_composed': True,
        'continuous_application_precedes_same_packet_measurement': True,
        'physical_field_reference_and_applied_offset_discrepancies_retained': True,
        'reference_generation_and_constructor_Rmag_persist': True,
        'rejected_measurements_still_advance_attempt_gate': True,
        'startup_full_history_and_entry_qualified': False,
        'source_uniform_complete_magnetic_word_qualified': False,
        'deployment_finite_precision_closed': False,
        'complete_word_finite_identity': False,
        'ALT_STARTUP_PASS': False,
        'ALT_LIVE_PASS': False,
        'ALT_END_TO_END_PASS': False,
    }
