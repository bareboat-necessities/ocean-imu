"""User-authorized MEMS input domain and finite private-observer induction.

These are theorem premises on actual filter API arguments, not deployed clamps.
The new Live input condition intersects the physical/bias/ISS histories; the
stricter startup residual and temporal contracts continue independently.
The scalar observer induction preserves the startup integral on every finite
prefix using an invariant of the rounding lattice, without replacing Live entry
by a reset. A sharper 30,602-sample bound is retained as an optional consequence.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_binary32_mahony as MAH
from tools.stability.ou3_alt_contraction import finite_mahony_prefix_totality as PREFIX
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SOURCE
from tools.stability.ou3_alt_contraction import finite_startup_sensor_contract as START
from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as CLOCK

QUALIFICATION='OU3_ALT_LIVE_MEMS_INPUT_V1'
PROFILE='BMI270_MPU6886_RAW_API_V1'
DOMAIN=Path(__file__).resolve().parents[1]/'ou3_alt_live_input_domain.json'
GYRO_ABS_CAP=F(35)
ACCEL_ABS_CAP=F(160)
U,ETA=START.U,START.ETA


def source_cell_upper(cap):
    """Upper endpoint of the cap's RNE cell in the pre-API real shadow."""
    q=F(cap)
    if q<=0 or not B.is_binary32(q): raise ValueError('positive representable cap required')
    return q+B._pow2(max(-149,B._floor_log2_positive(q)-23))/2


def domain():
    data=json.loads(DOMAIN.read_text())
    if data['schema']!=QUALIFICATION or data['profile']!=PROFILE:
        raise ValueError('wrong Live MEMS input profile')
    if (F(data['gyro_component_abs_upper_rad_s']),F(data['accel_component_abs_upper_mps2']))!=(GYRO_ABS_CAP,ACCEL_ABS_CAP):
        raise ValueError('Live MEMS input caps changed; rederive source-uniform inequalities')
    return data


@dataclass(frozen=True)
class PacketAdmission:
    raw: SOURCE.RawImuSample
    machine_gyro: tuple
    machine_accel: tuple
    profile: str=PROFILE

    def __post_init__(self):
        if not isinstance(self.raw,SOURCE.RawImuSample) or self.profile!=PROFILE:
            raise TypeError('same physical-source raw packet and fixed MEMS profile required')
        object.__setattr__(self,'machine_gyro',tuple(self.machine_gyro))
        object.__setattr__(self,'machine_accel',tuple(self.machine_accel))
        for label,exact,machine,cap in (
            ('gyro',self.raw.raw_gyro_body,self.machine_gyro,GYRO_ABS_CAP),
            ('accelerometer',self.raw.raw_accel_body,self.machine_accel,ACCEL_ABS_CAP)):
            if any(abs(F(x))>source_cell_upper(cap) for x in exact):
                raise ValueError(label+' exceeds commissioned Live raw API component cap')
            expected=tuple(B.rn32(x) for x in exact)
            if tuple(machine)!=expected:
                raise ValueError(label+' machine values detached from same source API rounding')
            if any(abs(x)>cap for x in expected):
                raise ValueError(label+' exceeds commissioned Live raw API component cap')


def check_packet(raw):
    """Admit the same packet before executing any coefficient or state update."""
    domain()
    if not isinstance(raw,SOURCE.RawImuSample): raise TypeError('source-owned raw packet required')
    # Check the domain before attempting conversion, including enormous finite
    # rational inputs that would otherwise overflow the API rounding graph.
    for label,values,cap in (('gyro',raw.raw_gyro_body,GYRO_ABS_CAP),
                            ('accelerometer',raw.raw_accel_body,ACCEL_ABS_CAP)):
        if any(abs(x)>source_cell_upper(cap) for x in values):
            raise ValueError(label+' exceeds commissioned Live raw API component cap')
    return PacketAdmission(raw,tuple(B.rn32(x) for x in raw.raw_gyro_body),
                            tuple(B.rn32(x) for x in raw.raw_accel_body))


require_raw_packet=check_packet


def conversion_certificate():
    # Rational pi upper bound is sufficient for the real hardware envelope;
    # explicit float compilation matches the app's literal DEG2RAD expression.
    pi_upper=F(314159265358979324,10**17)
    pi_literal=B.rn32(F(314159265358979323846,10**20))
    deg2rad=B.div(pi_literal,180)
    gyro_fullscale=B.mul(2000,deg2rad)
    acc_fullscale=B.mul(16,B.rn32(START.G))
    assert 2000*pi_upper/180<GYRO_ABS_CAP
    assert gyro_fullscale<GYRO_ABS_CAP and acc_fullscale<ACCEL_ABS_CAP
    return {'compiled_gyro_2000dps_rad_s':gyro_fullscale,
            'compiled_accel_16g_mps2':acc_fullscale,
            'gyro_conversion_margin_rad_s':GYRO_ABS_CAP-gyro_fullscale,
            'accel_conversion_margin_mps2':ACCEL_ABS_CAP-acc_fullscale,
            'arbitrary_calibration_preserves_caps_inferred':False}


def default_application_temperature_certificate():
    """The shipped AtomS3R OU3 call site uses the API's literal 35-C default.

    This is a call-site consequence, not a temperature bound imposed on the
    general four-argument filter API or its explicit thermal-history graph.
    """
    root=Path(__file__).resolve().parents[3]
    app=(root/'sensors/full_marine_ins/atomS3R_ins_kalman_ou3/atomS3R_ins_kalman_ou3.ino').read_text()
    wrapper=(root/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h').read_text()
    mekf=(root/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h').read_text()
    if app.count('fusion_.update(')!=1 or 'fusion_.update(dt_, w_cal_, a_cal_);' not in app:
        raise RuntimeError('OU3 application temperature call needs re-audit')
    if 'float tempC = 35.0f)' not in wrapper or 'impl_.updateTime(dt, gyro_body_ned, acc_body_ned, tempC);' not in wrapper:
        raise RuntimeError('OU3 wrapper default temperature propagation needs re-audit')
    if 'static constexpr T tempC_ref = T(35.0);' not in mekf or 'k_a_ * (tempC - tempC_ref)' not in mekf:
        raise RuntimeError('MEKF thermal model reference needs re-audit')
    delta=B.sub(35,35); modeled=B.mul(B.rn32(F(1,500)),delta)
    assert delta==0 and modeled==0
    return {'application_default_temperature_c':F(35),
            'application_default_thermal_model_component':modeled,
            'shipping_three_argument_call_thermal_term_zero':True,
            'general_four_argument_API_temperature_restricted':False,
            'general_thermal_history_totality_qualified':False}


@lru_cache(maxsize=1)
def _finite_horizon_certificate():
    """Scalar source-uniform totality for the initialized shared prefix.

    Assumptions: the source-bound seed enters Q2, integral starts at zero at
    actual construction, at most MAX_STEPS calls, actual h=RN32(1/200), and
    the dormant guard passes the admitted packet unchanged. These are the
    existing finite-word premises; no Live integral reset is introduced.
    """
    domain()
    # Recheck source/gains and both original startup profile proofs. The new
    # bound cannot enlarge a seed region or replace startup reachability.
    inherited=PREFIX.build()
    q2_actual,_,_=START._normalized_quaternion_bound()
    q2,qp=PREFIX.Q2,PREFIX.QP
    assert q2_actual<q2<qp**2
    upper=PREFIX.upper
    # Every raw accel component is at most 160. The squared norm is finite,
    # nonnegative, and the all-finite-word inverse-sqrt lemma supplies Q2 for
    # its normalized 3-vector, even for zero/subnormal norm sums.
    acc_square=upper(ACCEL_ABS_CAP**2)
    accel_norm_sum=upper(upper(2*acc_square)+acc_square)
    assert accel_norm_sum<100000
    half_gravity=(1+U)**4*q2/2+10*ETA
    feedback=upper(2*upper(qp*half_gravity))
    assert feedback<PREFIX.ERROR
    kp,ki=B.rn32(F(1,5)),B.rn32(F(1,50))
    n=CLOCK.MAX_STEPS; dt=CLOCK.DT_FLOAT
    delta=upper(upper(ki*PREFIX.ERROR)*dt)
    integral=n*upper(delta)/(1-n*U)
    assert n*U<1 and integral<F(41,10)
    corrected=upper(upper(GYRO_ABS_CAP+integral)+upper(kp*PREFIX.ERROR))
    assert corrected<40
    scaled=upper(corrected*B.rn32(dt/2))
    term=upper(qp*scaled)
    un=upper(qp+upper(upper(2*term)+term))
    assert un<F(7,5)
    norm=upper(upper(upper(2*upper(un*un))+upper(un*un))+upper(un*un))
    assert norm<8
    # Finite normalized successor returns to Q2. The induction is independent
    # of which admitted gyro axis/direction occurs at any of the 30,602 calls.
    # Project the same raw accel with ||d(q)||2=||q||2^2. sqrt(3)<7/4 avoids
    # irrational witnesses while retaining the vector correlation.
    acc_norm=F(7,4)*ACCEL_ABS_CAP
    derr=4*q2*((1+U)**8-1)+16*ETA
    vertical=(q2+2*derr)*acc_norm
    for _ in range(8): vertical=upper(vertical)
    vertical=upper(vertical+B.rn32(START.G))
    assert vertical<322
    return {
        'arithmetic_profile':MAH.PROFILE,
        'max_steps':n,
        'startup_integral_not_reset_at_Live':True,
        'raw_gyro_component_upper':GYRO_ABS_CAP,
        'raw_accel_component_upper':ACCEL_ABS_CAP,
        'real_shadow_gyro_component_upper':source_cell_upper(GYRO_ABS_CAP),
        'real_shadow_accel_component_upper':source_cell_upper(ACCEL_ABS_CAP),
        'accel_norm_sum_upper':accel_norm_sum,
        'feedback_component_upper':feedback,
        'integral_component_upper':integral,
        'feedback_corrected_rate_component_upper':corrected,
        'Euler_component_upper':un,
        'Euler_norm_sum_upper':norm,
        'normalized_quaternion_norm2_upper':q2_actual,
        'vertical_abs_upper':vertical,
        'initialized_scalar_Mahony_prefix_totality_closed':True,
        'vertical_input_abs_512_closed_under_prefix_premises':True,
        'startup_seed_and_reachability_inferred':False,
        'target_compiler_and_libm_qualified':False,
        'complete_MEKF_covariance_solver_word_totality_closed':False,
        'inherited_startup_source_hashes':inherited['reviewed_source_hashes'],
    }


@lru_cache(maxsize=1)
def prefix_certificate():
    """Every finite initialized prefix under the raw-input/scalar profile.

    This uses a rounding-lattice barrier, not accumulated linear-in-time
    error. The per-call integral increment D is below the upper half-ULP at
    4096. Since RNE is monotone and RN(4096+D)=4096, [-4096,4096] is forward
    invariant. The symmetric negative boundary follows identically. The
    elapsed float has the same invariant argument at 131072 seconds.
    """
    sharp=_finite_horizon_certificate()
    upper=PREFIX.upper
    q2,qp=PREFIX.Q2,PREFIX.QP
    kp,ki=B.rn32(F(1,5)),B.rn32(F(1,50))
    dt=CLOCK.DT_FLOAT
    increment=upper(upper(ki*PREFIX.ERROR)*dt)
    integral=F(4096)
    half_ulp=B._pow2(12-24)
    assert increment<half_ulp
    assert B.add(integral,increment)==integral
    assert B.add(-integral,-increment)==-integral
    corrected=upper(upper(GYRO_ABS_CAP+integral)+upper(kp*PREFIX.ERROR))
    assert corrected<4132
    scaled=upper(corrected*B.rn32(dt/2))
    term=upper(qp*scaled)
    un=upper(qp+upper(upper(2*term)+term))
    norm=upper(upper(upper(2*upper(un*un))+upper(un*un))+upper(un*un))
    assert norm<5000
    elapsed=F(131072)
    assert B.add(elapsed,dt)==elapsed
    # All three-square raw-norm operands and all four-square quaternion-norm
    # operands remain finite. The same exhaustive invsqrt lemma returns the
    # normalized state to Q2; its vertical projection bound is pointwise and
    # has no accumulated horizon dependence.
    return {
        'arithmetic_profile':MAH.PROFILE,
        'max_steps':None,
        'startup_integral_not_reset_at_Live':True,
        'raw_gyro_component_upper':GYRO_ABS_CAP,
        'raw_accel_component_upper':ACCEL_ABS_CAP,
        'accel_norm_sum_upper':sharp['accel_norm_sum_upper'],
        'feedback_component_upper':sharp['feedback_component_upper'],
        'integral_increment_upper':increment,
        'integral_component_upper':integral,
        'integral_barrier_upper_half_ulp':half_ulp,
        'integral_barrier_rounding_margin':half_ulp-increment,
        'feedback_corrected_rate_component_upper':corrected,
        'Euler_component_upper':un,
        'Euler_norm_sum_upper':norm,
        'normalized_quaternion_norm2_upper':sharp['normalized_quaternion_norm2_upper'],
        'elapsed_float_upper':elapsed,
        'vertical_abs_upper':sharp['vertical_abs_upper'],
        'initialized_scalar_Mahony_prefix_totality_closed':True,
        'all_initialized_finite_prefixes_totality_closed':True,
        'integral_rounding_lattice_invariant_closed':True,
        'elapsed_rounding_lattice_invariant_closed':True,
        'vertical_input_abs_512_closed_under_prefix_premises':True,
        'finite_horizon_refinement':sharp,
        'startup_seed_and_reachability_inferred':False,
        'target_compiler_and_libm_qualified':False,
        'complete_MEKF_covariance_solver_word_totality_closed':False,
        'indefinite_clock_progress_inferred':False,
        'inherited_startup_source_hashes':sharp['inherited_startup_source_hashes'],
    }


def build():
    data=domain(); proof=prefix_certificate()
    return {
        'qualification':QUALIFICATION,'profile':PROFILE,
        'raw_gyro_component_upper':GYRO_ABS_CAP,
        'raw_accel_component_upper':ACCEL_ABS_CAP,
        'devices':data['devices'],'conversion':conversion_certificate(),
        'application_thermal':default_application_temperature_certificate(),
        'prefix':proof,
        'input_domain_is_independent_of_execution_success':True,
        'startup_residual_and_temporal_contracts_unchanged':True,
        'arbitrary_bounded_raw_gyro_overflow_history_admitted':False,
        'finite_raw_input_cap_is_full_deployment_certificate':False,
        'storage_search_allowed':False,
        'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,'ALT_END_TO_END_PASS':False,
    }


def validate(report):
    return [key+' differs from commissioned Live input contract'
            for key,value in build().items() if report.get(key)!=value]
