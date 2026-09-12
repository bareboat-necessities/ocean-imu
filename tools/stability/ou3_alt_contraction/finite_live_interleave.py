"""Startup-rooted successive Live IMU/magnetic/hold event composition for ALT.

A magnetic call reads the private observer at the current IMU endpoint and
replaces ONLY the same product state's MEKF/control/calibration components.
The next IMU consumes that exact successor, including all 21 covariance rows.
No covariance, reference, physical origin, model or tuner is restarted.

The theorem-facing IMU edge consumes one source-qualified object containing BOTH
the exact physical transition and raw IMU packet. The theorem-facing magnetic
edge consumes a physical endpoint obtained as the before/after endpoint of an
admitted O^601_BRMM/BIAS transition, rather than relying only on a matching
history string. Lower-level ``imu_step``/``mag_step`` remain conditional finite
algebra and are not source admission. The Live magnetic edge uses the stronger
dual-clock composer: physical/inner-MEKF time and the outer binary32 wrapper
clock remain distinct all the way through delay/refinement/continuous-calibration
and same-packet measurement. This module cannot enable storage by itself.
"""
from __future__ import annotations
from dataclasses import dataclass, replace

from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_live_tilt_prefix as LIVE
from tools.stability.ou3_alt_contraction import finite_tilt_watchdog as WATCH
from tools.stability.ou3_alt_contraction import finite_startup_live_runtime_bridge as START
from tools.stability.ou3_alt_contraction import finite_startup_handoff_seed as SEED
from tools.stability.ou3_alt_contraction import finite_live_magnetic_word as MAG
from tools.stability.ou3_alt_contraction import finite_live_magnetic_dual_clock as MAGCLOCK
from tools.stability.ou3_alt_contraction import finite_mag_bias_gate as GATE
from tools.stability.ou3_alt_contraction import finite_mag_call_schedule as SCHEDULE
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE


@dataclass(frozen=True)
class State:
    live: LIVE.State
    magnetic: MAG.LiveState
    clock: SCHEDULE.Clock
    schedule: SCHEDULE.Schedule

    def __post_init__(self):
        if not isinstance(self.live, LIVE.State) or not isinstance(self.magnetic, MAG.LiveState):
            raise TypeError('same Live/watchdog and magnetic product state required')
        if not isinstance(self.clock, SCHEDULE.Clock) or not isinstance(self.schedule, SCHEDULE.Schedule):
            raise TypeError('persistent declared magnetic schedule required')
        core = self.live.live.mekf
        if core.reference.history_id != self.magnetic.memory.history_id:
            raise ValueError('magnetic and IMU physical history detached')
        if core.reference.live_origin != self.clock.live_time:
            raise ValueError('magnetic clock cannot restart the one-time Live origin')
        if self.magnetic.control.updates != self.clock.calls:
            raise ValueError('counted magnetic calls detached from literal control state')
        c = self.magnetic.control
        if (c.first_time is None) != (self.clock.calls == 0):
            raise ValueError('first magnetic attempt timestamp detached from count')
        if c.first_time is not None and not self.clock.live_time <= c.first_time <= self.clock.last_time:
            raise ValueError('first magnetic attempt clock detached from Live call history')
        # These three timestamps are outer-wrapper clock coordinates.  Comparing
        # them with the current physical endpoint is only a causal upper check;
        # they are deliberately not identified with the physical clock.
        for time in (self.magnetic.last_mag_time, self.magnetic.memory.last_hi_time,
                     self.magnetic.memory.applied.last_time):
            if time is not None and time > core.reference.time:
                raise ValueError('magnetic memory comes from a future physical endpoint')
        if (c.locked or c.hold) != (core.mode == 'H'):
            raise ValueError('H18/A21 state detached from bias lock/hold control')
        # The declared call schedule is indexed by physical source endpoints;
        # shipping magnetic branch timing itself is handled by MAGCLOCK.
        SCHEDULE.check_prefix(self.clock, self.schedule, time=core.reference.time)


@dataclass(frozen=True)
class Result:
    state: State
    event: object


def from_startup(bridge: START.Result, magnetic: MAG.StartupState, *,
                 proxy_q_norm, proxy_yaw_half, schedule=None):
    if not isinstance(bridge, START.Result) or not isinstance(magnetic, MAG.StartupState):
        raise TypeError('startup Live bridge and startup magnetic history required')
    if magnetic.ready is None:
        raise NotImplementedError('nongauged startup/Live branch remains open')
    live = bridge.state
    if magnetic.word.mag.proxy != bridge.frontend_before.tuner.vertical:
        raise ValueError('magnetic handoff detached from the same startup private observer')
    if live.mekf.reference.history_id != magnetic.memory.history_id:
        raise ValueError('magnetic handoff detached from startup physical history')
    for time in (magnetic.word.mag.last_mag_time, magnetic.memory.last_hi_time):
        if time is not None and time > live.mekf.reference.time:
            raise ValueError('startup magnetic state comes from a future sample')
    seed = SEED.seed(magnetic.word.mag.proxy, magnetic.ready.pending_yaw,
                     proxy_q_norm=proxy_q_norm, proxy_yaw_half=proxy_yaw_half)
    if tuple(P.quat_conj(seed.q_seed)) != live.mekf.q_hat:
        raise ValueError('Live attitude detached from the actual pending magnetic yaw gauge')
    if live.mekf.mode != 'H' or live.mekf.reference.live_origin != live.mekf.reference.time:
        raise ValueError('fresh gauged Live begins in H18 at its one-time physical origin')
    return State(LIVE.State(live, WATCH.State()), MAG.enter_live(magnetic),
                 SCHEDULE.Clock(live.mekf.reference.time),
                 SCHEDULE.default_schedule() if schedule is None else schedule)


def imu_step(state: State, raw, segment, **kwargs):
    """Conditional finite-algebra IMU step; not source admission by itself."""
    if not isinstance(state, State):
        raise TypeError('startup-rooted interleaved state required')
    if not isinstance(raw,SENSOR.RawImuSample):
        raise TypeError('same physical raw IMU packet required')
    if raw.deheel_body_to_internal != SENSOR.IDENTITY3:
        raise ValueError('nonzero deheel map violates existing ALT zero-wind-heel scope')
    if 'tilt_deg' in kwargs or 'reset_witness' in kwargs:
        raise TypeError('interleaved theorem word accepts no free tilt/reset output')
    out = LIVE.step_from_shipping_operands(state.live, raw, segment, **kwargs)
    return Result(State(out.state, state.magnetic, state.clock, state.schedule), out)


def imu_step_source_qualified(state: State, packet: SOURCE.QualifiedRawImuSample, **kwargs):
    if not isinstance(state,State) or not isinstance(packet,SOURCE.QualifiedRawImuSample):
        raise TypeError('interleaved state and source-qualified raw IMU packet required')
    qualified=packet.physical; core=state.live.live.mekf
    if qualified.root.history_id != core.reference.history_id:
        raise ValueError('qualified COMPLETE-BRMM history detached from current Live state')
    if qualified.root.live_origin != core.reference.live_origin:
        raise ValueError('qualified source restarted the one-time Live S origin')
    if qualified.segment.before != core.reference:
        raise ValueError('qualified source segment does not start at current physical endpoint')
    return imu_step(state,packet.raw,qualified.segment,**kwargs)


def mag_step(state: State, **kwargs):
    """Conditional dual-clock magnetic algebra; not source admission by itself."""
    if not isinstance(state, State):
        raise TypeError('startup-rooted interleaved state required')
    out = MAGCLOCK.live_call(state.magnetic, state.live.live.mekf,
                             state.live.live.tuner.vertical, **kwargs)
    clock = state.clock
    if out.measurement is not None and out.measurement.wrapper_attempted:
        # This ledger constrains the physical endpoint cadence.  The event's
        # internal outer-wrapper timestamps were already derived by MAGCLOCK.
        clock = SCHEDULE.record_call(clock, state.schedule, time=out.filter.reference.time)
    live = replace(state.live, live=replace(state.live.live, mekf=out.filter))
    return Result(State(live, out.state, clock, state.schedule), out)


def mag_step_source_qualified(state: State, endpoint, **kwargs):
    """Theorem-facing async magnetic edge at a checked outer physical endpoint.

    The endpoint may be a represented transition endpoint or the fresh Live
    origin.  Neither form upgrades necessary outer checks to full source
    admission.  Once qualified, the literal magnetic successor is evaluated by
    the dual-clock composer rather than identifying wrapper and physical time.
    """
    endpoint_types=(SOURCE.QualifiedPhysicalEndpoint,SOURCE.QualifiedPhysicalOrigin)
    if not isinstance(state,State) or not isinstance(endpoint,endpoint_types):
        raise TypeError('interleaved state and source-checked physical endpoint required')
    core=state.live.live.mekf; root=endpoint.root
    if root.history_id != core.reference.history_id or root.history_id != state.magnetic.memory.history_id:
        raise ValueError('qualified magnetic endpoint detached from persistent physical history')
    if root.live_origin != core.reference.live_origin:
        raise ValueError('qualified magnetic endpoint restarted Live origin')
    if endpoint.endpoint != core.reference:
        raise ValueError('asynchronous magnetic call not attached to current admitted endpoint')
    return mag_step(state,**kwargs)


def set_hold(state: State, *, hold):
    if not isinstance(state, State):
        raise TypeError('startup-rooted interleaved state required')
    event = GATE.set_hold(state.magnetic.control, state.live.live.mekf,
                          state.magnetic.memory.cfg.gate, hold=hold, live=True)
    magnetic = replace(state.magnetic, control=event.state)
    live = replace(state.live, live=replace(state.live.live, mekf=event.filter_state))
    return Result(State(live, magnetic, state.clock, state.schedule), event)


def readiness():
    src=SOURCE.readiness(); magclock=MAGCLOCK.readiness()
    return {
        'startup_gauge_and_private_observer_attached_at_Live_entry': True,
        'successive_IMU_mag_IMU_events_share_full_state_covariance': True,
        'continuous_magnetic_memory_not_restarted_at_Live': True,
        'magnetic_refinement_and_continuous_application_composed': True,
        'live_magnetic_outer_inner_dual_clock_composed':magclock['dual_clock_magnetic_word_composed'],
        'live_magnetic_wrapper_clock_prefix_arithmetic_closed':magclock['canonical_prefix_wrapper_clock_arithmetic_closed'],
        'interleaved_IMU_uses_same_operand_tilt_reset_entry': True,
        'free_watchdog_angle_and_reset_quaternion_forbidden': True,
        'source_qualified_finite_IMU_entry_available': True,
        'source_qualified_raw_IMU_packet_owned_by_same_physical_step': src['raw_IMU_packet_bound_to_same_qualified_physical_predecessor'],
        'persistent_sensor_residual_histories_required': src['persistent_gyro_and_accel_residual_history_tokens_required'],
        'source_qualified_async_magnetic_endpoint_entry_available': src['qualified_async_endpoint_comes_from_admitted_transition'],
        'source_checked_async_magnetic_endpoint_entry_available': bool(
            src['qualified_async_endpoint_comes_from_checked_outer_transition']
            and src['sample_zero_checked_outer_endpoint_available_without_fake_transition']),
        'sample_zero_full_source_membership_proved': src['sample_zero_full_source_membership_proved'],
        'correlated_COMPLETE_BRMM_left_inclusion_consumed': src['correlated_COMPLETE_BRMM_left_inclusion_consumed'],
        'persistent_BIAS_parameter_token_available': src['one_bias_family_parameter_token_over_word_required'],
        'external_hold_and_count_release_feed_next_IMU_mode': True,
        'same_history_every_represented_event_successor_exposed': True,
        'finite_prefix_mag_call_deadlines_checked': True,
        'quantitative_sensor_residual_ISS_envelope_attached': False,
        'finite_estimator_coefficients_bound_to_same_source_continuation': False,
        # The qualified entry exists, but the complete master has not yet been
        # rewritten to forbid every lower-level conditional magnetic call.
        'finite_magnetic_source_bound_to_same_COMPLETE_BRMM_history': False,
        'startup_magnetic_dual_clock_history_required_at_handoff': False,
        'infinite_schedule_qualified_by_finite_prefix': False,
        'source_uniform_complete_600_step_word_qualified': False,
        'storage_search_allowed': False,
        'deployment_clock_and_integer_lifetime_qualified': False,
        'ALT_STARTUP_PASS': False,
        'ALT_LIVE_PASS': False,
        'ALT_END_TO_END_PASS': False,
    }
