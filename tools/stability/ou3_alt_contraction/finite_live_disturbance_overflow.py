"""Bounded raw residuals do not establish deployed Live arithmetic totality.

This is an input-domain obstruction, not an admitted finite runtime word. The
same quiet physical history has ordinary startup followed by one representable
gyro residual. Its amplitude is bounded but the shipping squared norms overflow.
The native companion executes construction, actual goLive and 600 subsequent
calls; it never installs a filter state or edits the implementation.

The commissioned Live MEMS contract now excludes this raw pulse before execution;
this remains a regression against dropping that input premise. The existing
BoundedHistory bounds an executed finite forcing. No such object
is claimed for the nonfinite execution. Using failure to construct that object
to reject raw inputs would make the desired totality implication circular.
"""
from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path

import ou3_brmm_physical_wave_condition as PHYSICS
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_binary32_mahony as BITS
from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR

QUALIFICATION='OU3_ALT_LIVE_BOUNDED_RAW_RESIDUAL_OVERFLOW_V1'
ROOT=Path(__file__).resolve().parents[3]
ZERO=(F(0),)*3
IDENTITY=(F(1),F(0),F(0),F(0))
GRAVITY=F(196133,20000)
DT=B.rn32(F(1,200))
PULSE=F(1<<80)
RAW_ACCEL=(F(0),F(0),-B.rn32(GRAVITY))
RAW_ACCEL_RESIDUAL=(F(0),F(0),GRAVITY-B.rn32(GRAVITY))
MAX_FINITE=BITS.value(BITS.MAX_FINITE)
OVERFLOW_THRESHOLD=(MAX_FINITE+F(1<<128))/2


def quiet_packet(time=F(0),*,pulse=False,bias_family='BIAS0'):
    """One unchanged physical history, with a declared source-owned residual.

    ``pulse`` identifies the first post-Live sample, not a physical angular
    velocity. The physical body stays level and stationary on every sample.
    The ordinary startup accelerometer discrepancy is only float conversion.
    """
    if not isinstance(pulse,bool): raise TypeError('literal pulse event required')
    reference=CORE.Reference(F(time),IDENTITY,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
        F(0),'quiet-live-overflow-source','zero-bias',bias_family)
    gyro=(PULSE,F(0),F(0)) if pulse else ZERO
    return SENSOR.RawImuSample(reference,ZERO,gyro,RAW_ACCEL_RESIDUAL,gyro,
                               RAW_ACCEL,gravity_world=(F(0),F(0),GRAVITY))


def build():
    from tools.stability.ou3_alt_contraction import finite_live_input_contract as INPUT
    try: INPUT.check_packet(quiet_packet(pulse=True))
    except ValueError as error:
        if 'gyro exceeds' not in str(error): raise
    else: raise AssertionError('unphysical pulse entered the commissioned MEMS domain')
    if not PHYSICS.physical_condition()['quiet_zero_wave_is_admissible']:
        raise RuntimeError('quiet physical history admission changed')
    vertical=(ROOT/'src/tuner/VerticalAccelComplementary.h').read_text()
    mahony=(ROOT/'src/ahrs/Mahony_AHRS.h').read_text()
    if not all(x in vertical for x in ('!gyro.allFinite() || !acc.allFinite()',
        'ahrs_.update(gyro.x(), gyro.y(), gyro.z(),',
        'return Eigen::Quaternionf::Identity();',
        'up_ms2_ = -(down_row.dot(acc) + gravity_ms2);')):
        raise RuntimeError('vertical overflow/source/getter path needs re-audit')
    if 'invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3)' not in mahony:
        raise RuntimeError('Mahony normalization path needs re-audit')
    for family in ('BIAS0','BIAS1','BIAS2'):
        raw=quiet_packet(pulse=True,bias_family=family)
        if raw.internal_gyro!=(PULSE,0,0) or raw.physical.gyro_bias!=ZERO:
            raise AssertionError('gyro residual detached from quiet physical source')
    # Quiet aligned Mahony has q=(q0,0,0,0), zero integral, with q0 near one.
    # Even granting the much weaker q0>=1/2, literal float Euler integration
    # creates a component whose square exceeds the RNE overflow threshold.
    half_step=B.mul(F(1,2),DT)
    gyro_step=B.mul(PULSE,half_step)
    component_lower=B.mul(F(1,2),gyro_step)
    if component_lower**2<=OVERFLOW_THRESHOLD:
        raise AssertionError('selected pulse no longer forces normalization overflow')
    try: B.mul(component_lower,component_lower)
    except ValueError as error:
        if 'overflow' not in str(error): raise
    else: raise AssertionError('finite kernel accepted a nonfinite norm product')
    # The source-only amplitude is finite. This is not a produced post-event
    # RestrictedForcing and must never be attached to a CompleteWord object.
    raw_bound=PULSE+1
    if PULSE**2+RAW_ACCEL_RESIDUAL[2]**2>raw_bound**2:
        raise AssertionError('declared raw residual ceiling violated')
    return {
        'qualification':QUALIFICATION,
        'physical_history':'p=v=a=S=0; q=identity; beta=b_g=0',
        'pulse_gyro_residual':(PULSE,F(0),F(0)),
        'raw_residual_norm_ceiling':raw_bound,
        'all_raw_packet_components_finite_binary32':all(B.is_binary32(x) for x in (PULSE,*RAW_ACCEL)),
        'physical_quiet_and_zero_BIAS0_1_2_ancestry':True,
        'startup_gyro_residual_zero':True,
        'startup_accel_residual_is_only_API_rounding':RAW_ACCEL_RESIDUAL,
        'gyro_half_step_binary32':gyro_step,
        'quiet_q0_at_least_one_half_component_lower':component_lower,
        'squared_component_lower':component_lower**2,
        'RNE_finite_overflow_threshold':OVERFLOW_THRESHOLD,
        'squared_component_exceeds_overflow_threshold':True,
        'bounded_raw_residual_implies_finite_execution':False,
        'selected_MEMS_input_profile_excludes_pulse':True,
        'selected_commissioned_Live_input_contract_falsified':False,
        'nonfinite_execution_has_produced_RestrictedForcing':False,
        'source_uniform_deployment_supplies_closed':False,
        'target_hardware_counterexample_claimed':False,
        'shipping_filter_changed':False,
        'commissioned_startup_contract_falsified':False,
        'finite_storage_master_ready':False,
        'storage_search_allowed':False,
        'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,'ALT_END_TO_END_PASS':False,
    }


def validate(report):
    return [key+' differs from Live raw-residual overflow analysis'
            for key,value in build().items() if report.get(key)!=value]
