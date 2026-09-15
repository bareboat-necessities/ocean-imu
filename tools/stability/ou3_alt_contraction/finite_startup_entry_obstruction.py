"""Exact physical obstruction to universal fresh entry in one Cayley chart.

The family is a level boat at rest at any constant heading, with zero wave
translation, zero sensor biases and zero disturbances. COMPLETE-BRMM explicitly
admits quiet zero waves. Before any magnetic call the same IMU measurements
drive the same startup state at every heading. Shipping allows ungauged timeout
handoff. Its scalar identity attitude has Cayley denominator zero for a boat
facing south, and arbitrarily small denominator at nearby headings.

This is a representation failure, not a filter-instability counterexample.
Native public-API correspondence is tested separately; no replay supplies a
source-uniform bound. See docs/ou3-alt-startup-pre-rho.md for the induction.
"""
from __future__ import annotations

from fractions import Fraction as F
from dataclasses import replace

import ou3_brmm_physical_wave_condition as PHYSICS
from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as CLOCK
from tools.stability.ou3_alt_contraction import finite_startup_handoff_control as CONTROL
from tools.stability.ou3_alt_contraction import finite_binary32_mahony_startup as SEED
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as VERTICAL

GRAVITY = F(196133, 20000)
ZERO = (F(0),) * 3
IDENTITY = (F(1), F(0), F(0), F(0))
SOUTH = (F(0), F(0), F(0), F(1))


def quiet_reference(q_world_to_body=SOUTH, *, bias_family='BIAS0'):
    """Endpoint of the explicitly defined all-time constant physical history."""
    t = CLOCK.STARTUP_TIMEOUT_STEPS * CLOCK.DT_REAL
    return CORE.Reference(t, q_world_to_body, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO,
                          t, 'quiet-constant-heading', 'zero-bias', bias_family)


def quiet_packet(reference):
    return SENSOR.RawImuSample(reference, ZERO, ZERO, ZERO, ZERO,
                               (F(0), F(0), -GRAVITY))


def near_south(n):
    """Unit rational quaternion with Cayley z = n - 1/n, for every n >= 1."""
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError('positive integer family index required')
    return (F(2*n, n*n+1), F(0), F(0), F(n*n-1, n*n+1))


def quiet_observer_invariant():
    """Exact reachable fixed point, not sampled extrema or a fitted bound.

    The recurrence for q/integral/up does not read the elapsed clock. After the
    second literal update those coordinates are a fixed point. Checking its
    exact successor therefore closes induction for every later constant input
    on the named scalar arithmetic profile; elapsed time remains separately
    bounded by CLOCK. Normalizing a nonzero scalar quaternion is identity.
    """
    cfg = VERTICAL.Config(SEED.rn(F(1,5)), SEED.rn(F(1,50)), SEED.rn(GRAVITY), 20)
    kwargs = dict(dt=CLOCK.DT_FLOAT, gyro=ZERO, acc=(F(0), F(0), -cfg.gravity))
    first = SEED.step(VERTICAL.State(), cfg, **kwargs).vertical.state
    second = SEED.step(first, cfg, **kwargs).vertical.state
    third = SEED.step(second, cfg, **kwargs).vertical.state
    if replace(third, elapsed=second.elapsed) != second:
        raise AssertionError('quiet private observer is not an exact fixed point')
    for state in (first, second):
        if state.q[0] <= 0 or state.q[1:] != ZERO or state.integral != ZERO:
            raise AssertionError('quiet proxy lost scalar identity attitude')
    return second


def build():
    physical = PHYSICS.physical_condition()
    if not physical['quiet_zero_wave_is_admissible']:
        raise RuntimeError('quiet source admission changed; re-audit entry obstruction')
    for family in ('BIAS0', 'BIAS1', 'BIAS2'):
        quiet_packet(quiet_reference(bias_family=family))
    fixed = quiet_observer_invariant()
    control = CONTROL.Decision(CONTROL.Config(), CLOCK.TIMEOUT_CROSSING,
                               True, True, 'TunerReady', True, F(10), False)
    if not control.by_timeout or control.by_quality:
        raise AssertionError('ungauged timeout no longer follows literal control')
    try:
        CORE.cayley(SOUTH)
    except ValueError as exc:
        error = str(exc)
    else:
        raise AssertionError('Cayley pole unexpectedly admitted')
    return {
        'qualification': 'OU3_ALT_UNGAUGED_ENTRY_CHART_OBSTRUCTION_V1',
        'physical_history': 'p=v=a=S=0; constant yaw; beta=b_g=noise=0',
        'same_IMU_history_for_every_constant_heading': True,
        'quiet_proxy_fixed_point_scalar_binary32': fixed.q,
        'quiet_proxy_fixed_point_after_second_update_proved': True,
        'no_preLive_mag_calls_required_by_MAG_CALL_SCHEDULE_v1': True,
        'physical_zero_wave_and_zero_bias_history_used': True,
        'ungauged_timeout_control_permitted': control.by_timeout,
        'default_timeout_sample': CLOCK.STARTUP_TIMEOUT_STEPS,
        'relative_quaternion_at_south': SOUTH,
        'Cayley_denominator_at_south': SOUTH[0],
        'Cayley_constructor_error': error,
        'near_south_Cayley_z_identity': 'n - 1/n for every positive integer n',
        'universal_single_Cayley_fresh_entry_possible': False,
        'source_uniform_bounded_fresh_Cayley_radius_possible': False,
        'shipping_filter_instability_claimed': False,
        'multi_chart_or_quotient_entry_closed': False,
        'storage_search_allowed': False,
        'ALT_STARTUP_PASS': False,
        'ALT_LIVE_PASS': False,
        'ALT_END_TO_END_PASS': False,
    }


def validate(report):
    return [key + ' differs from the entry obstruction'
            for key, value in build().items() if report.get(key) != value]
