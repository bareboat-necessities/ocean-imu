"""Declared startup IMU profiles and exact consequences on the same raw packet.

The profile is a theorem hypothesis selected from engineering/device evidence.
It is not hardware qualification, nor an oracle for complete BRMM membership.
Peak checks do not establish the separately declared temporal direction budget.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
import json
from functools import lru_cache
from hashlib import sha256

import ou3_fast_inv_sqrt_interval as INV

from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT

DOMAIN=Path(__file__).resolve().parents[1]/'ou3_alt_startup_sensor_domain.json'
QUALIFICATION='OU3_ALT_COMMISSIONED_IMU_STARTUP_V1'
G=F(196133,20000)
A=F(44,5)
BIAS_COMPONENT=F(13,100)
BIAS_NORM=F(113,500)  # strictly outward of sqrt(3)*0.13
DT=F(1,200)
U=F(1,2**24)
ETA=F(1,2**150)
SOURCES={
    'src/ahrs/Mahony_AHRS.h':'ce0e93bf9c2b6213594985e3c9d3cf06edfd7418ba2bbfcace666f6ac93cddf3',
    'src/tuner/VerticalAccelComplementary.h':'7b2e38abc8d3260e3ce0f0d1111a59b2855932c3196298072783532292ae2c0a',
}


def domain():
    d=json.loads(DOMAIN.read_text())
    if d['qualification']!=QUALIFICATION: raise ValueError('wrong startup sensor qualification')
    expected={'gravity_mps2':G,'physical_acceleration_norm_upper_mps2':A,
              'true_accel_bias_component_upper_mps2':BIAS_COMPONENT,'imu_dt_s':DT}
    if any(F(d[k])!=v for k,v in expected.items()):
        raise ValueError('startup sensor domain detached from physical/BIAS scope')
    for name,p in d['profiles'].items():
        b=p['accel_residual_budget']
        total=F(b['commissioned_scale_cross_axis_nonlinearity_operator_upper'])*(G+A)+F(b['fast_error_vector_upper_mps2'])+F(b['acquisition_conversion_delay_model_vector_upper_mps2'])
        if total>F(p['accel_residual_vector_norm_upper_mps2']):
            raise ValueError(name+' component budget exceeds selected residual cap')
    return d


@dataclass(frozen=True)
class Profile:
    name:str
    accel_residual:F
    gyro_residual:F
    def __post_init__(self):
        p=domain()['profiles'].get(self.name)
        if p is None or (self.accel_residual,self.gyro_residual)!=(
                F(p['accel_residual_vector_norm_upper_mps2']),F(p['gyro_residual_vector_norm_upper_rad_s'])):
            raise ValueError('profile differs from declared commissioned sensor bounds')


def profile(name):
    p=domain()['profiles'][name]
    return Profile(name,F(p['accel_residual_vector_norm_upper_mps2']),F(p['gyro_residual_vector_norm_upper_rad_s']))


def norm2(v): return sum((F(x)**2 for x in v),F(0))


def magnitude_certificate(p:Profile):
    """Reverse triangle inequality, API RNE, nonnegative square/sum and sqrt.

    RN(x) has |RN(x)-x|<=u|x|+eta, including subnormals. Each
    squared component traverses at most three rounded nodes in the scalar
    norm graph; all summands are nonnegative. Monotonic square root plus its
    final rounding gives a source-uniform bound without sampling directions.
    """
    if not isinstance(p,Profile): raise TypeError('declared startup profile required')
    assert BIAS_NORM**2>=3*BIAS_COMPONENT**2
    lo=G-A-BIAS_NORM-p.accel_residual
    hi=G+A+BIAS_NORM+p.accel_residual
    conversion=U*hi+2*ETA
    stored_lo=lo-conversion; stored_hi=hi+conversion
    sq_lo=(1-U)**3*stored_lo**2-6*ETA
    sq_hi=(1+U)**3*stored_hi**2+6*ETA
    root_lo=ROOT.sqrt_enclosure(sq_lo)[0]
    root_hi=ROOT.sqrt_enclosure(sq_hi)[1]
    machine_lo=(1-U)*root_lo-ETA
    machine_hi=(1+U)*root_hi+ETA
    threshold=B.rn32(F(1,1000))
    return {'raw_norm_lower':lo,'raw_norm_upper':hi,
            'API_conversion_norm_error_upper':conversion,
            'stored_norm_lower':stored_lo,'stored_norm_upper':stored_hi,
            'computed_norm_squared_lower':sq_lo,'computed_norm_squared_upper':sq_hi,
            'computed_norm_lower':machine_lo,'computed_norm_upper':machine_hi,
            'seed_threshold':threshold,'seed_margin_lower':machine_lo-threshold,
            'every_admitted_sample_passes_seed_norm_test':machine_lo>threshold,
            'norm_intermediates_finite_and_positive':F(1,100)<sq_lo<sq_hi<400,
            'startup_SVD_or_handoff_qualified':False}


@lru_cache(maxsize=1)
def _normalized_quaternion_bound():
    """Cover every finite nonnegative norm word, INCLUDING zero/subnormals.

    The reviewed interval engine encloses the literal integer seed/Newton
    program, with 16 complete bit cells per exponent. This is exhaustive
    interval coverage, not evaluation at sampled float arguments.
    """
    shell=F(0); ymax=F(0)
    for exponent in range(255):
        for k in range(16):
            start=(exponent<<23)+k*(1<<19); end=start+(1<<19)-1
            x,y=INV._cell_inverse_sqrt(start,end)
            if y.lo<=0: raise ArithmeticError('inverse-sqrt enclosure lost positive output')
            shell=max(shell,F(INV.F32.mul(INV.F32.mul(x,y),y).hi))
            ymax=max(ymax,F(y.hi))
    assert ymax<2**65
    # n is the rounded 4-term sum of squares (at most four nodes per
    # path, seven absolute underflow charges). t=RN(RN(n*y)*y).
    ny2=(shell+ETA*ymax+ETA)/(1-U)**2
    ideal=(ny2+7*ETA*ymax*ymax)/(1-U)**4
    assert ideal<2
    actual=(1+U)**2*ideal+10*ETA  # four component multiplications
    assert actual<F(139,125)
    return actual,shell,ymax


def vertical_supply_certificate(p:Profile):
    """Bound the actual projection AFTER a defined scalar Mahony update.

    Does not assume accurate tilt. It does not prove that the preceding
    seed, feedback, Euler update or norm sum is defined for every source.
    """
    actual,shell,ymax=_normalized_quaternion_bound(); q2=F(139,125)
    # For the exact quadratic down row, ||d(q)||_2=||q||_2^2.
    # Eight relative roundings and sixteen absolute underflow charges
    # dominate each scalar down-row component's actual polynomial tree.
    derr=4*q2*((1+U)**8-1)+16*ETA
    value=(q2+2*derr)*magnitude_certificate(p)['stored_norm_upper']
    for _ in range(8): value=(1+U)*value+ETA  # dot products and sums
    value=(1+U)*(value+B.rn32(G))+ETA
    assert value<32
    return {'computed_quaternion_norm2_upper':actual,
            'normalization_scalar_shell_upper':shell,'inverse_sqrt_abs_upper':ymax,
            'all_nonnegative_finite_norm_words_including_zero_subnormal_covered':True,
            'vertical_abs_upper_after_defined_Mahony_update':value,
            'vertical_WPE_input_abs_32_implication_closed':True,
            'preceding_Mahony_operations_source_uniformly_total':False}


@dataclass(frozen=True)
class History:
    """One declared sensor history beside one physical/BIAS history.

    The quantifier includes the domain's temporal direction budget. Labels
    preserve ancestry; they do not prove that an arbitrary trace has it.
    """
    profile:Profile
    physical_history_id:str
    bias_root:str
    bias_family:str
    gyro_residual_history_id:str
    accel_residual_history_id:str
    def __post_init__(self):
        if not isinstance(self.profile,Profile): raise TypeError('declared startup profile required')
        if self.bias_family not in ('BIAS0','BIAS1','BIAS2'): raise ValueError('declared BIAS family required')
        for x in (self.physical_history_id,self.bias_root,self.gyro_residual_history_id,self.accel_residual_history_id):
            if not isinstance(x,str) or not x: raise ValueError('persistent history identity required')


def check_packet(history:History,raw:SENSOR.RawImuSample,*,ordinal:int):
    """Necessary restrictions of the declared history, never trace admission."""
    if not isinstance(history,History) or not isinstance(raw,SENSOR.RawImuSample):
        raise TypeError('same declared startup history and raw physical packet required')
    if type(ordinal) is not int or ordinal<1: raise ValueError('positive startup ordinal required')
    r=raw.physical
    if (r.history_id,r.bias_root,r.bias_family)!=(history.physical_history_id,history.bias_root,history.bias_family):
        raise ValueError('startup sensor profile detached from physical/BIAS history')
    if r.time!=(ordinal-1)*DT: raise ValueError('startup sensor sample detached from 5ms source clock')
    if raw.deheel_body_to_internal!=SENSOR.IDENTITY3:
        raise ValueError('startup profile requires zero heel')
    if raw.gravity_world!=(F(0),F(0),G): raise ValueError('startup gravity differs from declared source')
    d=domain(); p=history.profile
    if norm2(r.acceleration)>A*A: raise ValueError('physical startup acceleration exceeds existing BRMM cap')
    if any(abs(x)>BIAS_COMPONENT for x in r.beta): raise ValueError('startup true bias exceeds BIAS envelope')
    if norm2(raw.accel_residual_internal)>p.accel_residual**2: raise ValueError('startup accelerometer residual exceeds profile')
    if norm2(raw.gyro_residual_internal)>p.gyro_residual**2: raise ValueError('startup gyro residual exceeds profile')
    # Rational 22/7 exceeds pi, hence also bounds the declared 35 deg/s rate.
    body_rate=F(35)*F(22,7)/180
    if norm2(raw.omega_sample_internal)>body_rate**2: raise ValueError('startup physical body rate exceeds bound')
    b=F(d['true_gyro_bias_initial_norm_upper_rad_s']) if ordinal==1 else F(d['true_gyro_bias_norm_upper_rad_s'])
    if norm2(r.gyro_bias)>b*b: raise ValueError('startup true gyro bias exceeds profile')
    cert=magnitude_certificate(p)
    if not cert['raw_norm_lower']**2<=norm2(raw.internal_accel)<=cert['raw_norm_upper']**2:
        raise AssertionError('same physical sensor identity violated triangle bound')
    return cert


def build():
    root=Path(__file__).resolve().parents[3]
    if any(sha256((root/path).read_bytes()).hexdigest()!=digest for path,digest in SOURCES.items()):
        raise RuntimeError('startup normalization/projection source changed; re-audit sensor supply proof')
    d=domain(); profiles={name:magnitude_certificate(profile(name)) for name in d['profiles']}
    return {'qualification':QUALIFICATION,'profiles':profiles,
            'declared_domain':d,
            'reviewed_source_hashes':dict(SOURCES),
            'vertical_supplies':{name:vertical_supply_certificate(profile(name)) for name in d['profiles']},
            'sensor_profile_selected':True,'same_packet_peak_restrictions_materialized':True,
            'source_uniform_seed_norm_margin_closed':all(x['every_admitted_sample_passes_seed_norm_test'] for x in profiles.values()),
            'temporal_direction_budget_is_additional_source_premise':True,
            'temporal_budget_follows_from_peak_checks':False,
            'old_Mahony_invariant_covers_new_profiles':False,
            'hardware_admission_qualified':False,'target_compiler_libm_qualified':False,
            'universal_startup_reachability_closed':False,'storage_search_allowed':False}


def validate(report):
    return [k+' differs from declared startup sensor qualification' for k,v in build().items() if report.get(k)!=v]
