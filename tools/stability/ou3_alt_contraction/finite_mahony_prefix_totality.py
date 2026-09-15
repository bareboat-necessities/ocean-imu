"""Finite scalar Mahony totality on the commissioned startup prefix.

The controlling prerequisite is existence of the same-history startup word,
before its dissipation ratio can be formed.  A projection bound conditional on
an already-defined update did not establish that existence.  This module closes
the scalar state-update part for ordinary FromTwoVectors seeding, including
integral feedback, through the existing timeout-plus-word horizon.

This is a conditional arithmetic theorem, not a timeout/capture theorem.  It
uses the declared RNE/gradual-underflow/no-FMA graph and correctly rounded sqrt.
The near-antiparallel solver and actual compiler/libm remain separate premises.
The discarded Euler-angle outputs are outside this state-update theorem.
"""
from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_mahony as M
from tools.stability.ou3_alt_contraction import finite_binary32_mahony_startup as SEED
from tools.stability.ou3_alt_contraction import finite_startup_sensor_contract as SENSOR
from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as CLOCK

QUALIFICATION = 'OU3_ALT_SCALAR_MAHONY_PREFIX_TOTALITY_V1'
U, ETA = SENSOR.U, SENSOR.ETA
Q2 = F(139, 125)
QP = F(53, 50)
ERROR = F(13, 10)
SOURCE = 'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
SOURCE_HASH = 'fabd03e9c3eb6069df107c7413ffb4b33fbdcd1ce06d3923b0c1ceeb3bcd7359'


def upper(x):
    """RNE magnitude bound, valid including subnormal results."""
    return (1 + U) * F(x) + ETA


def ordinary_seed_certificate(p):
    """Correlated norm bound; never divide independent norm intervals.

    For .1 <= ||a|| <= 21, the scalar square/sum/root graph has
      ((1-u)^5-10000 eta)||a||^2 <= s^2
          <= ((1+u)^5+10000 eta)||a||^2.
    Six sum/square underflow charges and the sqrt cross term consume less
    than 100 eta before division by ||a||^2 >= .01.  The three component
    divisions consume less than 20 eta in squared norm.  Both normalizations
    therefore have the same bound, without multiplying two unrelated shells.

    Since the second vector is exactly e3, c=v0.z and the cross product is
    exactly (v0.y,-v0.x,0).  In the ordinary branch, m=1+cutoff>0.  Writing
    ||v0||^2 <= 1+delta gives the cancellation-preserving seed inequality
       ||q||^2 <= K [1+delta/(2m)] + 100 eta.
    The 1-c^2 term is retained instead of bounding cross/s separately.
    """
    packet = SENSOR.magnitude_certificate(p)
    assert F(1, 10) < packet['stored_norm_lower']
    assert packet['stored_norm_upper'] < 21
    normalized_lo = (1-U)**2 / ((1+U)**5 + 10000*ETA) - 20*ETA
    normalized_hi = (1+U)**2 / ((1-U)**5 - 10000*ETA) + 20*ETA
    delta = F(1, 10**6)
    assert 1-delta < normalized_lo <= normalized_hi < 1+delta
    # The first normalized vector meets the same .1--21 domain for the
    # second normalization.  This is the required induction of the ratio.
    assert F(1, 100) < normalized_lo and normalized_hi < 21**2
    cutoff = SEED.rn(-1 + SEED.rn(F(1, 100000)))
    gap = 1 + cutoff
    assert gap > 0
    k = (1+U)**7 / (1-U)**3
    seed_q2 = k * (1 + delta/(2*gap)) + 100*ETA
    assert seed_q2 < F(53, 50) < Q2
    return {
        'profile': p.name,
        'normalized_norm2_lower': normalized_lo,
        'normalized_norm2_upper': normalized_hi,
        'ordinary_branch_one_plus_c_lower': gap,
        'ordinary_seed_norm2_upper': seed_q2,
        'ordinary_seed_operations_finite_and_denominators_positive': True,
        'ordinary_seed_enters_normalized_state_shell': True,
        'near_antiparallel_branch_qualified': False,
    }


def prefix_certificate(p):
    """Induction on the actual scalar recurrence for at most MAX_STEPS.

    The reviewed inverse-sqrt interval lemma covers every nonnegative finite
    norm word, including zero/subnormals.  Its four-term normalization bound
    also covers three terms (fewer nonnegative rounding/underflow charges).
    Thus both normalized acceleration and quaternion have squared norm < Q2.

    Each exact half-gravity component has magnitude <= ||q||^2/2.  Four
    rounded operations plus ten eta dominate the actual scalar expression.
    Feedback is bounded per component, retaining the integral recurrence.
    For I[n+1] <= (1+u)I[n]+D, I[0]=0,
       I[n] <= n D/(1-nu),
    by (1+u)^n <= 1/(1-nu).  This bound is conditional on the finite horizon;
    it is never used as an all-time integral invariant.
    """
    seed = ordinary_seed_certificate(p)
    q2_actual, _, _ = SENSOR._normalized_quaternion_bound()
    assert q2_actual < Q2 < QP**2
    half_gravity = (1+U)**4 * Q2/2 + 10*ETA
    feedback = upper(2*upper(QP*half_gravity))
    assert feedback < ERROR
    kp, ki = SEED.rn(F(1,5)), SEED.rn(F(1,50))
    dt = CLOCK.DT_FLOAT
    delta_i = upper(upper(ki*ERROR)*dt)
    n = CLOCK.MAX_STEPS
    assert n*U < 1
    # RN(I+delta_i) rounds the increment as well as the previous accumulator.
    integral = n*upper(delta_i)/(1-n*U)
    assert integral < F(41,10)
    # API conversion follows the existing sensor relation.  The larger
    # ongoing true-gyro-bias bound also covers the first-sample bound.
    raw_gyro = F(35)*F(22,7)/180 + F(1,50) + p.gyro_residual
    gyro = upper(raw_gyro)
    assert gyro < F(2,3)
    corrected_rate = upper(upper(gyro+integral)+upper(kp*ERROR))
    assert corrected_rate < 6
    half_dt = SEED.rn(dt/2)
    scaled_rate = upper(corrected_rate*half_dt)
    q_term = upper(QP*scaled_rate)
    un = upper(QP+upper(upper(2*q_term)+q_term))
    assert un < F(28,25)
    norm = upper(upper(upper(2*upper(un*un))+upper(un*un))+upper(un*un))
    assert norm < 6
    # A rounded finite nonnegative norm enters the exhaustive Newton lemma;
    # the normalization successor returns to Q2, closing the induction.
    elapsed = n*upper(dt)/(1-n*U)
    assert elapsed < 154
    vertical = SENSOR.vertical_supply_certificate(p)
    assert vertical['vertical_abs_upper_after_defined_Mahony_update'] < 32
    return {
        'profile': p.name,
        'max_steps': n,
        'ordinary_seed': seed,
        'normalized_quaternion_norm2_upper': q2_actual,
        'half_gravity_component_abs_upper': half_gravity,
        'feedback_component_abs_upper': feedback,
        'integral_component_abs_upper': integral,
        'raw_gyro_component_abs_upper': gyro,
        'feedback_corrected_rate_component_abs_upper': corrected_rate,
        'Euler_component_abs_upper': un,
        'Euler_norm_sum_upper': norm,
        'elapsed_upper': elapsed,
        'vertical_abs_upper': vertical['vertical_abs_upper_after_defined_Mahony_update'],
        'ordinary_seed_and_scalar_state_prefix_totality_closed': True,
        'initialized_prefix_totality_given_Q2_seed_and_zero_integral_closed': True,
        'startup_timeout_reachability_inferred_from_horizon': False,
        'near_antiparallel_seed_totality_closed': False,
        'unused_Euler_angle_library_calls_qualified': False,
        'target_compiler_and_sqrt_qualified': False,
        'indefinite_integral_invariant_closed': False,
    }


def build():
    root = Path(__file__).resolve().parents[3]
    sources = {**SENSOR.SOURCES, SOURCE: SOURCE_HASH}
    if any(sha256((root/path).read_bytes()).hexdigest()!=digest
           for path,digest in sources.items()):
        raise RuntimeError('Mahony source/gains changed; re-audit prefix totality')
    profiles = {name: prefix_certificate(SENSOR.profile(name))
                for name in SENSOR.domain()['profiles']}
    return {
        'qualification': QUALIFICATION,
        'arithmetic_profile': M.PROFILE,
        'reviewed_source_hashes': sources,
        'profiles': profiles,
        'ordinary_seed_scalar_prefix_totality_closed': all(
            p['ordinary_seed_and_scalar_state_prefix_totality_closed'] for p in profiles.values()),
        'every_startup_branch_totality_closed': False,
        'startup_deployment_supply_bounds_closed': False,
        'storage_search_allowed': False,
    }


def validate(report):
    return [key+' differs from scalar-prefix totality theorem'
            for key,value in build().items() if report.get(key)!=value]
