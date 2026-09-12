"""Startup-rooted successive Live IMU/magnetic/hold event composition for ALT.

A magnetic call reads the private observer at the current IMU endpoint and
replaces ONLY the same product state's MEKF/control/calibration components.
The next IMU consumes that exact successor, including all 21 covariance rows.
No covariance, reference, physical origin, model or tuner is restarted.

Finite-prefix timing checks are not qualification of an infinite schedule.
The supplying component arithmetic/source premises remain open, and this
module cannot enable the source-uniform master/storage gate.
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
from tools.stability.ou3_alt_contraction import finite_mag_bias_gate as GATE
from tools.stability.ou3_alt_contraction import finite_mag_call_schedule as SCHEDULE


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
        for time in (self.magnetic.last_mag_time, self.magnetic.memory.last_hi_time,
                     self.magnetic.memory.applied.last_time):
            if time is not None and time > core.reference.time:
                raise ValueError('magnetic memory comes from a future physical endpoint')
        if (c.locked or c.hold) != (core.mode == 'H'):
            raise ValueError('H18/A21 state detached from bias lock/hold control')
        SCHEDULE.check_prefix(self.clock, self.schedule, time=core.reference.time)


@dataclass(frozen=True)
class Result:
    state: State
    event: object


def from_startup(bridge: START.Result, magnetic: MAG.StartupState, *,
                 proxy_q_norm, proxy_yaw_half, schedule=None):
    """Use carried startup magnetic state and the SAME gauged handoff, not a fake Live root."""
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
    if not isinstance(state, State):
        raise TypeError('startup-rooted interleaved state required')
    if not isinstance(raw,SENSOR.RawImuSample):
        raise TypeError('same physical raw IMU packet required')
    if raw.deheel_body_to_internal != SENSOR.IDENTITY3:
        raise ValueError('nonzero deheel map violates existing ALT zero-wind-heel scope')
    out = LIVE.step(state.live, raw, segment, **kwargs)
    return Result(State(out.state, state.magnetic, state.clock, state.schedule), out)


def mag_step(state: State, **kwargs):
    if not isinstance(state, State):
        raise TypeError('startup-rooted interleaved state required')
    # Neither a new physical endpoint nor a second proxy can be supplied by a
    # caller. The exact IMU predecessor supplies BOTH.
    out = MAG.live_call(state.magnetic, state.live.live.mekf,
                        state.live.live.tuner.vertical, **kwargs)
    clock = state.clock
    if out.measurement is not None and out.measurement.wrapper_attempted:
        clock = SCHEDULE.record_call(clock, state.schedule, time=out.filter.reference.time)
    live = replace(state.live, live=replace(state.live.live, mekf=out.filter))
    return Result(State(live, out.state, clock, state.schedule), out)


def set_hold(state: State, *, hold):
    """Compose a literal external setAccBiasHold event without inventing eventual release."""
    if not isinstance(state, State):
        raise TypeError('startup-rooted interleaved state required')
    event = GATE.set_hold(state.magnetic.control, state.live.live.mekf,
                          state.magnetic.memory.cfg.gate, hold=hold, live=True)
    magnetic = replace(state.magnetic, control=event.state)
    live = replace(state.live, live=replace(state.live.live, mekf=event.filter_state))
    return Result(State(live, magnetic, state.clock, state.schedule), event)


def readiness():
    return {
        'startup_gauge_and_private_observer_attached_at_Live_entry': True,
        'successive_IMU_mag_IMU_events_share_full_state_covariance': True,
        'continuous_magnetic_memory_not_restarted_at_Live': True,
        'magnetic_refinement_and_continuous_application_composed': True,
        'external_hold_and_count_release_feed_next_IMU_mode': True,
        'same_history_every_represented_event_successor_exposed': True,
        'finite_prefix_mag_call_deadlines_checked': True,
        'infinite_schedule_qualified_by_finite_prefix': False,
        'source_uniform_complete_600_step_word_qualified': False,
        'storage_search_allowed': False,
        'deployment_clock_and_integer_lifetime_qualified': False,
        'ALT_STARTUP_PASS': False,
        'ALT_LIVE_PASS': False,
        'ALT_END_TO_END_PASS': False,
    }
