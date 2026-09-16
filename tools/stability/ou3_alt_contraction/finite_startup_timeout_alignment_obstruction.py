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
DOCKING_A=F('2.9')
DOCKING_OMEGA=F('.6')
DOCKING_ACCEL_RESIDUAL=F('.01')
DOCKING_TIME=F('150.01')
# A fixed constructive source, not an optimizer dependency or a fitted rho.
DOCKING_RATES=tuple(map(F,('''
.181930600531 .515399372842 .540507226271 .552654215095 .560708466752
.566385833303 .589351557302 .550697071356 .584179369245 .576549738968
.579617345047 .580925573491 .581591761612 .582454083615 .582985293928
.583455328412 .583430350242 .583885522476 .583319266784 .584607984154
.575530312911 .586086065677 .588020650327 .586860570603 .591085103126
.584217236131 -.140746326636 .210382990552 .040083210406 .111985855580
''').split()))


def docking_yaw_and_rate(time):
    t=F(time)
    if t<0: raise ValueError('nonnegative startup time required')
    yaw=F(0)
    for j,rate in enumerate(DOCKING_RATES):
        start=F(5*j); end=DOCKING_TIME if j==29 else F(5*(j+1))
        if t<end: return yaw+rate*(t-start),rate
        yaw+=rate*(end-start)
    return yaw+DOCKING_OMEGA*(t-DOCKING_TIME),DOCKING_OMEGA


def docking_source_bounds():
    """All-time physical circle; sensor-tail bounds still require induction.

    Unlike the timeout witness, the indefinite residual controller is only
    finite-execution evidence. The mean bound below is conditional on its
    stated all-time residual bound, not a substitute for proving that bound.
    """
    lo,hi=ROOT.sqrt_enclosure(SENSOR.G**2+DOCKING_A**2)
    bounds={'position_norm':DOCKING_A/DOCKING_OMEGA**2,
            'velocity_norm':DOCKING_A/DOCKING_OMEGA,
            'acceleration_norm':DOCKING_A,
            'centered_position_primitive':2*DOCKING_A/DOCKING_OMEGA**3,
            'yaw_rate_norm':max(DOCKING_OMEGA,*map(abs,DOCKING_RATES)),
            'direction_mean_if_residual_bounded':1-SENSOR.G/hi+
                2*DOCKING_ACCEL_RESIDUAL/(lo-DOCKING_ACCEL_RESIDUAL),
            'direction_primitive':DOCKING_A/(DOCKING_OMEGA*lo),
            'accel_residual_required':DOCKING_ACCEL_RESIDUAL}
    assert bounds['position_norm']<F('8.1')
    assert bounds['velocity_norm']<F('5.5')
    assert bounds['acceleration_norm']<SENSOR.A
    assert bounds['yaw_rate_norm']<35*PI_LO/180
    assert bounds['centered_position_primitive']<1100
    assert F('.018')<DOCKING_OMEGA/(2*PI_HI)<DOCKING_OMEGA/(2*PI_LO)<F('.88')
    assert bounds['direction_mean_if_residual_bounded']<F('.1')
    assert bounds['direction_primitive']<F('1.5')
    for profile in SENSOR.domain()['profiles'].values():
        assert DOCKING_ACCEL_RESIDUAL<F(profile['accel_residual_budget']['fast_error_vector_upper_mps2'])
    return {'bounds':bounds,'same_all_time_physical_circle':True,
            'zero_bias_and_driver_solve_all_BIAS_families':True,
            'finite_prefix_requires_actual_packet_audit':True,
            'indefinite_tail_residual_bound_proved':False,
            'indefinite_tail_rounding_tube_proved':False,
            'eventual_finite_startup_refuted':False,
            'later_uniform_startup_deadline_proved':False,
            'storage_search_allowed':False}


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


def audit_native_packets(lines, *, docking=False):
    """Check every emitted API packet against ONE analytic physical history.

    The hexadecimal values are the actual binary32 API operands. The physical
    trajectory remains the same circular wave/yaw history for the entire file.
    Peak residuals supply a continuous bounded extension between samples; this
    is not an inference of BRMM membership from sampled extrema.
    """
    build(); count=0; accel_max2=F(0); gyro_max2=F(0)
    amplitude=DOCKING_A if docking else A
    omega=DOCKING_OMEGA if docking else OMEGA
    accel_cap=DOCKING_ACCEL_RESIDUAL if docking else ACCEL_RESIDUAL_AUDIT
    smallest_gyro_cap=min(F(p['gyro_residual_vector_norm_upper_rad_s'])
                          for p in SENSOR.domain()['profiles'].values())
    for line in lines:
        fields=line.split()
        if len(fields)!=7: raise ValueError('native source packet must have ordinal and six API components')
        ordinal=int(fields[0]); count+=1
        if ordinal!=count: raise ValueError('native source packet sequence restarted or skipped')
        values=tuple(F(float.fromhex(x)) for x in fields[1:])
        gyro,acc=values[:3],values[3:]
        t=F(ordinal-1)*SENSOR.DT
        yaw,rate=(docking_yaw_and_rate if docking else yaw_and_rate)(t)
        sine,cosine=_sin_cos_interval(omega*t-yaw)
        expected=((amplitude*cosine[0],amplitude*cosine[1]),(amplitude*sine[0],amplitude*sine[1]),(-SENSOR.G,-SENSOR.G))
        acc_error2=sum((max(abs(a-lo),abs(a-hi))**2 for a,(lo,hi) in zip(acc,expected)),F(0))
        gyro_error2=gyro[0]**2+gyro[1]**2+(gyro[2]-rate)**2
        if acc_error2>accel_cap**2:
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


def conditional_inverted_axis_fixed_point():
    """Exact bad observer invariant; its reachability is deliberately UNPROVED.

    A true level vessel may yaw at RN32(.2). With I_z=-RN32(.2), zero
    horizontal I, and the inverted x-axis quaternion, every feedback cross
    product and corrected gyro component is exactly zero. This is useful for
    an eventual-capture architecture review, never an installed startup state.
    """
    from dataclasses import replace
    from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
    from tools.stability.ou3_alt_contraction import finite_binary32_mahony as MAH
    from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as VERTICAL
    c=F(16748919,16777216)
    spin=B.rn32(F(1,5)); zero=F(0)
    cfg=VERTICAL.Config(B.rn32(F(1,5)),B.rn32(F(1,50)),B.rn32(SENSOR.G),F(20))
    before=VERTICAL.State(q=(zero,c,zero,zero),integral=(zero,zero,-spin),initialized=True)
    one=MAH.step_initialized(before,cfg,dt=B.rn32(SENSOR.DT),
                            gyro=(zero,zero,spin),acc=(zero,zero,-B.rn32(SENSOR.G)))
    MAH.Arithmetic(list(one.operations),list(one.normalizations)).verify()
    seated=one.vertical.state
    two=MAH.step_initialized(seated,cfg,dt=B.rn32(SENSOR.DT),
                            gyro=(zero,zero,spin),acc=(zero,zero,-B.rn32(SENSOR.G)))
    MAH.Arithmetic(list(two.operations),list(two.normalizations)).verify()
    if replace(two.vertical.state,elapsed=seated.elapsed)!=seated:
        raise AssertionError('inverted axis observer is not a fixed point modulo elapsed clock')
    if seated.q!=before.q or seated.integral!=before.integral:
        raise AssertionError('inverted fixed point changed quaternion/integral')
    return {'quaternion':seated.q,'integral':seated.integral,
            'physical_yaw_rate':spin,'up_output':seated.up,
            'exact_named_binary32_observer_fixed_point':True,
            'normalized_projected_accel_z_is_positive_gravity':True,
            'source_first_sample_seed_can_be_replaced_by_this_state':False,
            'admitted_startup_reaches_inverted_fixed_point':False,
            'eventual_finite_startup_refuted':False,
            'target_compiler_equivalence_qualified':False}


def conditional_circular_antialignment():
    """Exact-real invariant enlarging the possible docking target.

    This is a continuous observer/source construction. Neither its discrete
    preservation nor entry from the real startup seed is asserted.
    """
    h=(F(20,101),F(0),F(-99,101))
    down=tuple(-x for x in h)
    yaw_rate=F(-3,5); axial_rate=F(73,100)
    integral=(axial_rate*h[0],F(0),axial_rate*h[2]-yaw_rate)
    gyro=(F(0),F(0),yaw_rate)
    acc=(SENSOR.G*F(20,99),F(0),-SENSOR.G)
    magnitude=SENSOR.G*F(101,99)
    def cross(a,b):
        return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
    corrected=tuple(a+b for a,b in zip(gyro,integral))
    assert sum(x*x for x in h)==1
    assert tuple(-a/magnitude for a in acc)==down
    assert cross(down,h)==(0,0,0)
    assert cross(corrected,h)==(0,0,0)
    projection=sum(x*y for x,y in zip(h,acc))
    assert projection==magnitude>0
    amplitude=acc[0]; frequency=abs(yaw_rate)
    caps={'position_norm':amplitude/frequency**2,
          'velocity_norm':amplitude/frequency,
          'acceleration_norm':amplitude,
          'centered_position_primitive':2*amplitude/frequency**3,
          'direction_mean_norm':F(2,101),
          'direction_primitive_norm':F(20,101)/frequency}
    assert caps['position_norm']<F('8.1') and caps['velocity_norm']<F('5.5')
    assert amplitude<SENSOR.A and caps['centered_position_primitive']<1100
    assert caps['direction_mean_norm']<F('.10') and caps['direction_primitive_norm']<F('1.5')
    assert abs(yaw_rate)<35*PI_LO/180
    return {'observer_down':h,'measured_down':down,'integral':integral,
            'raw_gyro':gyro,'raw_accel':acc,'projected_force':projection,
            'physical_and_temporal_bounds':caps,
            'continuous_observer_and_source_invariant_closed':True,
            'finite_binary32_preservation_qualified':False,
            'admitted_startup_reaches_this_invariant':False,
            'physical_transition_into_new_circle_qualified':False,
            'eventual_finite_startup_refuted':False}
