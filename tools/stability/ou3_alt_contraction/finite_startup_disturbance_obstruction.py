"""A finite residual bound alone cannot certify shipping startup reachability.

The exact quiet physical family below leaves the measured gravity direction
unchanged but its magnitude below the literal first-sample seed threshold.
This refutes extending universal startup to arbitrary bounded IMU residuals.
It does not infer hardware admission from a finite bound and does not refute
a separately specified small-disturbance startup theorem or post-Live ISS.
"""
from __future__ import annotations

from dataclasses import replace
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

import ou3_brmm_physical_wave_condition as PHYSICS
from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_binary32_mahony_startup as SEED
from tools.stability.ou3_alt_contraction import finite_machine_accel_guard_binary32 as GUARD
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as VERTICAL
from tools.stability.ou3_alt_contraction import finite_startup_handoff_control as CONTROL
from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as CLOCK

QUALIFICATION='OU3_ALT_STARTUP_BOUNDED_RESIDUAL_OBSTRUCTION_V1'
GRAVITY=F(196133,20000)
EPSILON=F(1,2048)
ZERO=(F(0),)*3
IDENTITY=(F(1),F(0),F(0),F(0))
RAW_ACCEL=(F(0),F(0),-EPSILON)
RESIDUAL=(F(0),F(0),GRAVITY-EPSILON)
ROOT=Path(__file__).resolve().parents[3]
SOURCES={
    'src/tuner/VerticalAccelComplementary.h':'7b2e38abc8d3260e3ce0f0d1111a59b2855932c3196298072783532292ae2c0a',
    'src/tuner/AccelVibrationGuard.h':'09ec79b607a4c7e2cc7961930813a3c0cf66d6f82c7bddde4b7dc4084a01398e',
    'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h':'fabd03e9c3eb6069df107c7413ffb4b33fbdcd1ce06d3923b0c1ceeb3bcd7359',
}


def quiet_packet(time=F(0),*,bias_family='BIAS0'):
    """Restriction of p=v=a=beta=b_g=0, q=identity and one constant residual.

    Zero S is independent of its integration origin. The endpoint's origin=0
    is only the existing physical container's anchor; no Live event is assumed.
    Zero bias/driver solves every BIAS0/1/2 recurrence for every allowed phi.
    """
    reference=CORE.Reference(F(time),IDENTITY,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
                             F(0),'quiet-small-measured-gravity','zero-bias',bias_family)
    return SENSOR.RawImuSample(reference,ZERO,ZERO,RESIDUAL,ZERO,RAW_ACCEL,
                               gravity_world=(F(0),F(0),GRAVITY))


def _unique_exp(magnitude):
    lo,hi=GUARD.EXP.exp_minus_enclosure(magnitude)
    a,b=B.rn32(lo),B.rn32(hi)
    if a!=b: raise ArithmeticError('default guard exp enclosure crosses an RNE cell')
    return a


def guard_successor(before):
    """Same constant-input guard step under the named finite RNE relation."""
    c=GUARD.Config(); h=CLOCK.DT_FLOAT
    if not before.initialized:
        return GUARD.step(before,c,raw_gyro=ZERO,raw_acc=RAW_ACCEL,dt=h)
    def decay(f):
        return _unique_exp(B.mul(B.mul(B.mul(GUARD.TWO,GUARD.PI_F),f),h))
    alpha,gamma=decay(c.cutoff_hz),decay(c.detect_hz)
    for coefficient in (alpha,gamma):
        if GUARD._lin_values(coefficient,-EPSILON,-EPSILON)!=(-EPSILON,):
            raise AssertionError('constant guard value not preserved by every named contraction')
    return GUARD.step(before,c,raw_gyro=ZERO,raw_acc=RAW_ACCEL,dt=h,
        lp_alpha_exp=alpha,lp_successors=(RAW_ACCEL,)*c.poles,
        detect_gamma_exp=gamma,detect_successors=(RAW_ACCEL,ZERO),
        removed_beta_exp=decay(c.removed_rms_hz),removed_ms_successor=ZERO,
        removed_rms_sqrt=F(0),slew_exp=_unique_exp(B.div(h,c.slew_tau)),
        weight_successor=F(0))


def observer_fixed_point():
    """One exact successor equals reset; induction handles every later input."""
    cfg=VERTICAL.Config(SEED.rn(F(1,5)),SEED.rn(F(1,50)),SEED.rn(GRAVITY),20)
    before=VERTICAL.State()
    result=SEED.step(before,cfg,dt=CLOCK.DT_FLOAT,gyro=ZERO,acc=RAW_ACCEL)
    result.startup.validate()
    if result.startup.branch!='acc-norm-too-small' or result.vertical.state!=before:
        raise AssertionError('literal small-norm startup no longer preserves reset')
    return result


def seed_magnitude_margin(*,physical_acceleration_cap,bias_cap,residual_cap,
                          conversion_guard_and_norm_error_cap):
    """Sufficient norm margin by reverse triangle inequality, not admission.

    The last cap must cover API conversion, conditioning and norm evaluation.
    All four caps must be proved for the same startup history before this
    implication can qualify shipping. No covariance is used as a noise bound.
    """
    caps=tuple(F(x) for x in (physical_acceleration_cap,bias_cap,residual_cap,
                              conversion_guard_and_norm_error_cap))
    if any(x<0 for x in caps): raise ValueError('nonnegative source/error caps required')
    return GRAVITY-sum(caps)-SEED.rn(F(1,1000))


def build():
    for path,expected in SOURCES.items():
        if sha256((ROOT/path).read_bytes()).hexdigest()!=expected:
            raise RuntimeError('startup obstruction needs shipping source re-audit: '+path)
    if not PHYSICS.physical_condition()['quiet_zero_wave_is_admissible']:
        raise RuntimeError('quiet physical source admission changed')
    for family in ('BIAS0','BIAS1','BIAS2'):
        packet=quiet_packet(bias_family=family)
        if packet.internal_accel!=RAW_ACCEL or packet.physical.beta!=ZERO:
            raise AssertionError('constant residual detached from zero-bias physical source')
    seed=observer_fixed_point()
    first=guard_successor(GUARD.State()); second=guard_successor(first.state)
    if replace(second.state,samples=first.state.samples)!=first.state:
        raise AssertionError('constant-input guard is not a fixed point modulo bookkeeping')
    if first.conditioned_acc!=RAW_ACCEL or second.conditioned_acc!=RAW_ACCEL:
        raise AssertionError('dormant guard altered the seed operand')
    # Grant every other favorable handoff predicate. Neither path can overcome
    # a false proxy_initialized, at any representable clock or tuner stage.
    decision=CONTROL.Decision(CONTROL.Config(),CLOCK.TIMEOUT_CROSSING,
                              True,False,'TunerReady',True,F(10),True)
    if decision.handoff_due: raise AssertionError('uninitialized proxy can hand off')
    return {
        'qualification':QUALIFICATION,
        'physical_history':'p=v=a=S=0; q=identity; beta=b_g=0',
        'raw_acceleration':RAW_ACCEL,
        'accel_residual':RESIDUAL,
        'residual_norm':GRAVITY-EPSILON,
        'measured_gravity_direction_error':F(0),
        'raw_norm_binary32':seed.startup.norm,
        'seed_threshold_binary32':SEED.rn(F(1,1000)),
        'seed_norm_margin':seed.startup.norm-SEED.rn(F(1,1000)),
        'same_physical_and_residual_history_retained':True,
        'zero_wave_and_zero_BIAS0_1_2_history':True,
        'guard_constant_input_fixed_point':True,
        'guard_remains_dormant':second.state.weight==0,
        'proxy_reset_is_exact_fixed_point':True,
        'quality_and_timeout_both_require_initialized_proxy':not decision.handoff_due,
        'finite_prefix_no_initialization_induction':True,
        'arbitrary_bounded_residuals_imply_universal_startup':False,
        'post_Live_ISS_bound_supplies_startup_raw_norm_premise':False,
        'specified_small_disturbance_startup_theorem_falsified':False,
        'hardware_residual_admission_claimed':False,
        'target_libm_compiler_qualification_claimed':False,
        'startup_sensor_contract_selected':False,
        'storage_search_allowed':False,
        'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,'ALT_END_TO_END_PASS':False,
    }


def validate(report):
    return [key+' differs from startup residual obstruction'
            for key,value in build().items() if report.get(key)!=value]
