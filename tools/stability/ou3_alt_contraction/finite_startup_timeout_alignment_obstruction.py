"""Admitted physical/source witness for the failed 150-second startup deadline.

Source membership is analytic. Native execution is a finite counterexample to
that deadline on the tested compiler, not an all-target arithmetic certificate
or a counterexample to eventual startup. The same witness later reaches Live.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from math import factorial

from tools.stability.ou3_alt_contraction import finite_startup_sensor_contract as SENSOR
from tools.stability.ou3_alt_contraction import finite_startup_direction_sampling as TRIG
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT

QUALIFICATION='OU3_ALT_TIMEOUT_ALIGNMENT_PHYSICAL_WITNESS_V1'
A=F('3.313664')
OMEGA=F('0.64')
GYRO_RESIDUAL_NOMINAL=F('0.019')
ACCEL_RESIDUAL_AUDIT=F('0.0001')
PI_LO=F('3.141592653589')
PI_HI=F('3.141592653590')
# At a breakpoint the right derivative is the physical sampled body rate.
YAW_PIECES=((F(0),F(5),F('.2')), (F(5),F(125),F('.61')),
            (F(125),F(130),F('.41')), (F(130),F(135),F('-.2')))
FINAL_RATE=F('.2')


def yaw_and_rate(time):
    t=F(time)
    if t<0: raise ValueError('nonnegative startup time required')
    yaw=F(0)
    for start,end,rate in YAW_PIECES:
        if t<end: return yaw+rate*(t-start),rate
        yaw+=rate*(end-start)
    return yaw+FINAL_RATE*(t-YAW_PIECES[-1][1]),FINAL_RATE


@lru_cache(maxsize=1)
def build():
    domain=SENSOR.domain()
    d_lo,d_hi=ROOT.sqrt_enclosure(SENSOR.G**2+A**2)
    mean_ideal=1-SENSOR.G/d_hi
    # The audited finite API accelerometer values are restrictions of one
    # continuous bounded residual interpolation. Put its direction difference
    # into m(t); no derivative or frequency assumption on that residual is used.
    direction_residual=2*ACCEL_RESIDUAL_AUDIT/(d_lo-ACCEL_RESIDUAL_AUDIT)
    mean=mean_ideal+direction_residual
    primitive=A/(OMEGA*d_lo)
    physical={'position_norm_upper':A/OMEGA**2,
              'velocity_norm_upper':A/OMEGA,
              'acceleration_norm_upper':A,
              'centered_position_primitive_upper':2*A/OMEGA**3,
              'body_rate_norm_upper':max(abs(x[2]) for x in YAW_PIECES),
              'direction_mean_norm_upper':mean,
              'direction_primitive_norm_upper':primitive,
              'accel_residual_norm_upper':ACCEL_RESIDUAL_AUDIT}
    lo,hi=TRIG.pi_interval()
    if not PI_LO<lo<hi<PI_HI: raise ArithmeticError('pi range-reduction endpoints failed')
    assert physical['position_norm_upper']<F('8.1')
    assert physical['velocity_norm_upper']<F('5.5')
    assert A<SENSOR.A and physical['centered_position_primitive_upper']<1100
    assert physical['body_rate_norm_upper']<35*PI_LO/180
    assert F('.018')<OMEGA/(2*PI_HI)<OMEGA/(2*PI_LO)<F('.88')
    assert mean<F(domain['total_direction_mean_chord_norm_upper'])
    assert primitive<F(domain['total_direction_primitive_norm_upper_s'])
    for p in domain['profiles'].values():
        assert ACCEL_RESIDUAL_AUDIT<F(p['accel_residual_budget']['fast_error_vector_upper_mps2'])
        assert GYRO_RESIDUAL_NOMINAL<F(p['gyro_residual_vector_norm_upper_rad_s'])
    return {'qualification':QUALIFICATION,'exact':{k:str(v) for k,v in physical.items()},
            'physical_kinematic_caps_and_centered_primitive_closed':True,
            'same_history_continuous_direction_budget_closed':True,
            'zero_bias_and_zero_driver_solve_all_BIAS_families':True,
            'both_commissioned_sensor_profiles_covered_by_packet_audit':True,
            'first_sample_requires_no_installed_filter_state':True,
            'pre_Live_magnetic_calls_required_by_current_schedule':False,
            'native_timeout_branch_failure_must_be_checked_by_execution':True,
            'all_target_compiler_correspondence_closed':False,
            'eventual_finite_startup_refuted':False,
            'later_uniform_startup_deadline_proved':False,
            'storage_search_allowed':False}


def _sin_cos_interval(angle):
    """Exact rational Taylor/range-reduction enclosure, no platform libm."""
    x=F(angle); middle=(PI_LO+PI_HI)/2
    turns=(x+middle)//(2*middle)
    x-=2*turns*middle
    reduction=abs(turns)*(PI_HI-PI_LO)
    sine=sum(((-1)**k*x**(2*k+1)/factorial(2*k+1) for k in range(9)),F(0))
    cosine=sum(((-1)**k*x**(2*k)/factorial(2*k) for k in range(9)),F(0))
    sine_error=abs(x)**19/factorial(19)+reduction
    cosine_error=abs(x)**18/factorial(18)+reduction
    return (sine-sine_error,sine+sine_error),(cosine-cosine_error,cosine+cosine_error)


def audit_native_packets(lines):
    """Check every emitted API packet against ONE analytic physical history.

    The hexadecimal values are the actual binary32 API operands. The physical
    trajectory remains the same circular wave/yaw history for the entire file.
    Peak residuals supply a continuous bounded extension between samples; this
    is not an inference of BRMM membership from sampled extrema.
    """
    build(); count=0; accel_max2=F(0); gyro_max2=F(0)
    smallest_gyro_cap=min(F(p['gyro_residual_vector_norm_upper_rad_s'])
                          for p in SENSOR.domain()['profiles'].values())
    for line in lines:
        fields=line.split()
        if len(fields)!=7: raise ValueError('native source packet must have ordinal and six API components')
        ordinal=int(fields[0]); count+=1
        if ordinal!=count: raise ValueError('native source packet sequence restarted or skipped')
        values=tuple(F(float.fromhex(x)) for x in fields[1:])
        gyro,acc=values[:3],values[3:]
        t=F(ordinal-1)*SENSOR.DT; yaw,rate=yaw_and_rate(t)
        sine,cosine=_sin_cos_interval(OMEGA*t-yaw)
        expected=((A*cosine[0],A*cosine[1]),(A*sine[0],A*sine[1]),(-SENSOR.G,-SENSOR.G))
        acc_error2=sum((max(abs(a-lo),abs(a-hi))**2 for a,(lo,hi) in zip(acc,expected)),F(0))
        gyro_error2=gyro[0]**2+gyro[1]**2+(gyro[2]-rate)**2
        if acc_error2>ACCEL_RESIDUAL_AUDIT**2:
            raise ValueError('native accelerometer detached from the admitted circular-wave packet')
        if gyro_error2>smallest_gyro_cap**2:
            raise ValueError('native gyro exceeds commissioned source residual')
        accel_max2=max(accel_max2,acc_error2); gyro_max2=max(gyro_max2,gyro_error2)
    if count<30602: raise ValueError('native source audit did not cover the disputed startup+word horizon')
    return {'samples':count,'accel_residual_norm2_upper':accel_max2,
            'gyro_residual_norm2_upper':gyro_max2,
            'same_analytic_physical_history_for_every_packet':True,
            'both_commissioned_profiles_pass':True}


def validate(report):
    expected=build()
    return [k+' differs from timeout witness scope' for k,v in expected.items() if report.get(k)!=v]
