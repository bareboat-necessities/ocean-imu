"""Continuous commissioned direction budgets do not imply sampled budgets.

An exact enclosure certificate for a physical sinusoid plus a continuous 200 Hz
sensor residual. This is a failure of an inference used in a possible
startup proof, not a counterexample to shipping capture or a storage result.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from math import factorial

from tools.stability.ou3_alt_contraction import finite_startup_sensor_contract as SENSOR
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT

QUALIFICATION = "OU3_ALT_CONTINUOUS_TO_SAMPLED_DIRECTION_GAP_V1"
A = F(349, 50)  # horizontal physical acceleration amplitude, 6.98 m/s^2
N = F(3, 25)  # horizontal 200 Hz residual peak, 0.12 m/s^2
T = F(2)
DT = F(1, 200)
RESIDUAL_CYCLES_PER_PERIOD = 400
COUNT = 400
BITS = 80


def _floor(x):
    return F((x.numerator << BITS) // x.denominator, 1 << BITS)


def _ceil(x):
    return -_floor(-x)


def _atan_interval(x, terms):
    value = sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(terms)), F(0))
    following = (-1)**terms*x**(2*terms+1)/F(2*terms+1)
    return min(value, value+following), max(value, value+following)


@lru_cache(maxsize=1)
def pi_interval():
    # Machin identity and the alternating series establish both endpoints.
    a, b = _atan_interval(F(1, 5), 30)
    c, d = _atan_interval(F(1, 239), 10)
    return 16*a-4*d, 16*b-4*c


def _sine_interval(k):
    # Reduce the full period to [-pi,pi], without approximate range reduction.
    k = k if k <= COUNT//2 else k-COUNT
    p, q = pi_interval()
    a, b = sorted((F(2*k, COUNT)*p, F(2*k, COUNT)*q))
    x = (a+b)/2
    value = sum(((-1)**j*x**(2*j+1)/factorial(2*j+1) for j in range(22)), F(0))
    error = abs(x)**45/factorial(45)+(b-a)/2
    return _floor(value-error), _ceil(value+error)


def _direction(x):
    # Exact real direction of (x,0,g); roots are integer/rational enclosures.
    g = SENSOR.G
    lo, hi = ROOT.sqrt_enclosure(g*g+x*x, bits=100)
    if x >= 0:
        rx = (_floor(x/hi), _ceil(x/lo))
    else:
        rx = (_floor(x/lo), _ceil(x/hi))
    return rx, (_floor(g/hi-1), _ceil(g/lo-1))


def _sample_mean(offset):
    sx_lo = sx_hi = sz_lo = sz_hi = F(0)
    for k in range(COUNT):
        a, b = _sine_interval(k)
        a, b = A*a+offset, A*b+offset
        # x/sqrt(g^2+x^2) is monotone; the axial component is even.
        rx_a, rz_a = _direction(a)
        rx_b, rz_b = _direction(b)
        sx_lo += rx_a[0]
        sx_hi += rx_b[1]
        sz_lo += min(rz_a[0], rz_b[0])
        sz_hi += F(0) if a <= 0 <= b else max(rz_a[1], rz_b[1])
    return (sx_lo/COUNT, sx_hi/COUNT), (sz_lo/COUNT, sz_hi/COUNT)


@lru_cache(maxsize=1)
def build():
    domain = SENSOR.domain()
    if F(domain['total_direction_mean_chord_norm_upper']) != F(1, 10):
        raise ValueError('certificate tied to the selected continuous mean budget')
    if F(domain['total_direction_primitive_norm_upper_s']) != F(3, 2):
        raise ValueError('certificate tied to the selected continuous primitive budget')
    if T/DT != COUNT or RESIDUAL_CYCLES_PER_PERIOD != COUNT:
        raise AssertionError('period/sample/residual frequency ancestry is inconsistent')
    pi_lo, pi_hi = pi_interval()
    baseline_x, baseline_z = _sample_mean(F(0))
    sampled_x, sampled_z = _sample_mean(N)

    # The continuous baseline horizontal mean is exactly zero by half-period
    # symmetry. A periodic midpoint rule bounds the axial integral. For
    # f(theta)=(1+(A/g)^2 sin(theta)^2)^(-1/2),
    # |f''| <= (A/g)^2 + 3 (A/g)^4/4.
    ratio2 = (A/SENSOR.G)**2
    second_derivative = ratio2+F(3, 4)*ratio2**2
    quadrature_error = second_derivative*(2*pi_hi/COUNT)**2/24
    baseline_mean_upper = -baseline_z[0]+quadrature_error

    # Let n(t)=N*cos(400*pi*t). Each sample reads +N. For the normalized
    # curve u(x)=(x,0,g)/sqrt(x*x+g*g), ||u''(x)|| <= 1/g^2. Taylor's
    # integral remainder is at most N^2/(2*g^2). The first-order term
    # Du(A*sin(pi*t))*n(t) has zero boundary contribution under integration
    # by parts. Its mean norm is at most N*A/(400*g^2), since |a'|<=A*pi.
    # No approximate cancellation of a sampled quadrature is credited here.
    residual_mean_change_upper = N*A/(RESIDUAL_CYCLES_PER_PERIOD*SENSOR.G**2)+N*N/(2*SENSOR.G**2)
    continuous_mean_upper = baseline_mean_upper+residual_mean_change_upper

    # The 2s-periodic function r-m has zero integral. Its primitive can use the
    # representative integration interval [-T/2,T/2], hence this all-time bound.
    chord_upper = (A+N)/SENSOR.G
    primitive_upper = T/2*(chord_upper+continuous_mean_upper)
    sampled_mean_norm2_lower = sampled_x[0]**2+sampled_z[1]**2

    # Exact-real detector bound, including the first-update transient. This is
    # evidence that the sampled residual does not create high-frequency samples;
    # native compiler/guard rounding qualification is explicitly separate.
    gamma_upper = F(23, 50)
    # gamma=exp(-pi/4). A positive Taylor partial sum is a strict lower
    # bound for exp(pi/4), so this inequality verifies gamma<0.46.
    if sum(((pi_lo/4)**k/factorial(k) for k in range(5)), F(0)) <= 1/gamma_upper:
        raise ArithmeticError('exact-real detector coefficient bound failed')
    phase_step = 2*pi_hi/COUNT
    first_difference = A*phase_step
    second_difference = A*phase_step**2
    detector_upper = gamma_upper**2*first_difference + gamma_upper**2/(1-gamma_upper)**2*second_difference

    values = {
        'continuous_mean_norm_upper': continuous_mean_upper,
        'continuous_primitive_norm_upper_s': primitive_upper,
        'sampled_mean_x_lower': sampled_x[0],
        'sampled_mean_x_upper': sampled_x[1],
        'sampled_mean_z_lower': sampled_z[0],
        'sampled_mean_z_upper': sampled_z[1],
        'sampled_mean_norm2_lower': sampled_mean_norm2_lower,
        'sampled_mean_norm2_excess_lower': sampled_mean_norm2_lower-F(1, 100),
        'baseline_periodic_quadrature_error_upper': quadrature_error,
        'residual_mean_change_upper': residual_mean_change_upper,
        'physical_position_norm_upper_m': A/pi_lo**2,
        'physical_velocity_norm_upper_mps': A/pi_lo,
        'centered_position_primitive_upper_ms': 2*A/pi_lo**3,
        'exact_real_guard_detector_abs_upper': detector_upper,
    }
    if not (continuous_mean_upper < F(1, 10) and primitive_upper < F(3, 2)):
        raise ArithmeticError('continuous source budget witness failed')
    if not (sampled_x[0] > 0 and sampled_z[1] < 0 and sampled_mean_norm2_lower > F(1, 100)):
        raise ArithmeticError('sampled mean separation failed')
    if not (A < SENSOR.A and A/pi_lo < F(11, 2) and A/pi_lo**2 < F(81, 10)):
        raise ArithmeticError('physical waveform exceeded unchanged caps')
    if not (2*A/pi_lo**3 < 1100 and detector_upper < F(3, 100)):
        raise ArithmeticError('physical primitive or exact-real guard margin failed')
    return {
        'qualification': QUALIFICATION,
        'exact': {k: str(v) for k, v in values.items()},
        'display': {k: float(v) for k, v in values.items()},
        'continuous_direction_budget_witness_closed': True,
        'same_numeric_sampled_budget_implication_refuted': True,
        'residual_is_continuous_200Hz_sinusoid': True,
        'sampled_values_are_smooth_wave_plus_constant_residual': True,
        'physical_wave_frequency_hz': '1/2',
        'zero_bias_solves_BIAS0_BIAS1_BIAS2': True,
        'residual_below_both_commissioned_fast_error_allowances': bool(all(
            N <= F(p['accel_residual_budget']['fast_error_vector_upper_mps2'])
            for p in domain['profiles'].values())),
        'full_shipping_guard_binary32_correspondence_closed': False,
        'shipping_startup_nonreachability_proved': False,
        'universal_startup_capture_closed': False,
        'storage_search_allowed': False,
    }


def validate(result):
    failures = []
    if result.get('qualification') != QUALIFICATION:
        failures.append('qualification mismatch')
    for key in ('continuous_direction_budget_witness_closed',
                'same_numeric_sampled_budget_implication_refuted',
                'residual_is_continuous_200Hz_sinusoid',
                'sampled_values_are_smooth_wave_plus_constant_residual',
                'zero_bias_solves_BIAS0_BIAS1_BIAS2',
                'residual_below_both_commissioned_fast_error_allowances'):
        if result.get(key) is not True:
            failures.append(key+' not proved')
    for key in ('full_shipping_guard_binary32_correspondence_closed',
                'shipping_startup_nonreachability_proved',
                'universal_startup_capture_closed', 'storage_search_allowed'):
        if result.get(key) is not False:
            failures.append(key+' promoted by a premise-transfer obstruction')
    return failures
