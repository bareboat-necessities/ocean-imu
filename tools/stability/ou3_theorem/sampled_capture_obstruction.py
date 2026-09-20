"""All-time sampling ambiguity obstructing the declared six-degree capture.

This is an analytic continuous-history construction, not a finite simulation.
The 2x2 certificate discharges the actual magnetic-service premise on its
stationary source record. It does not assert observer divergence or refute
conditional local dissipativity. See docs/ou3-sampled-capture-obstruction.md.
"""
from fractions import Fraction as F
from .matrix_certificates import add, encoded, ldlt, matmul, transpose


def transition(t):
    return [[F(1), t], [F(0), F(1)]]


def process(t):
    sg, qb = F('.00135'), F('1e-10')
    return [[sg*sg*t+qb*t**3/3, qb*t*t/2], [qb*t*t/2, qb*t]]


def predict(p, t):
    f = transition(t)
    return add(matmul(matmul(f, p), transpose(f)), process(t))


def correct(p, r):
    return [[p[i][j]-p[i][0]*p[0][j]/(p[0][0]+r)
             for j in range(2)] for i in range(2)]


def service_certificate():
    # Literal deployed profile: 200 Hz IMU, 25 Hz magnetometer; sigma_m=.8 uT.
    # At zero rate the source Simpson process integral is exact, and its
    # positive-definite 6x6 LDL hygiene does not add a numerical floor.
    spacing, horizon, count, sb = F('.04'), F('1.04'), 25, F('.02')
    r = F('.8')**2/F(75)**2
    upper = [[F('.008'), F(0)], [F(0), F('.00002')]]
    prefix = [[F('.009'), F(0)], [F(0), F('.00003')]]
    # Handoff need not coincide with the 25 Hz magnetic clock. Charge the
    # entire first partial cell before using the periodic posterior invariant.
    seed = [[F('.087')**2, F(0)], [F(0), F('1e-6')]]
    first = predict(seed, spacing)
    ia, ib, ic = upper[0][0]-first[0][0], upper[1][1]-first[1][1], first[0][1]
    assert ia > 0 and ib > 0 and ia*ib > ic*ic
    pre = predict(upper, spacing)
    invariant = add(upper, correct(pre, r), F(-1))
    il, ip = ldlt(invariant)
    # All t in [0,.04], not just the right endpoint. Diagonal upper bounds
    # and the absolute cross bound are monotone on this scalar positive chain.
    a = prefix[0][0]-pre[0][0]
    b = prefix[1][1]-pre[1][1]
    c = pre[0][1]
    assert a > 0 and b > 0 and a*b > c*c
    # Any one-second window has 25 consecutive applied magnetic observations.
    # det(O'O) is invariant under a common shift of observation times.
    determinant = sb**2*count**2*(count**2-1)*spacing**2/12
    trace_ceiling = count*(1+sb**2*horizon**2)
    process_noise_ceiling = count*process(horizon)[0][0]
    noise_ceiling = r+process_noise_ceiling
    raw_information_floor = determinant/(trace_ceiling*noise_ceiling)
    root_covariance_norm = max(prefix[0][0], prefix[1][1]/sb**2)
    loss_floor = 1/(1/raw_information_floor+root_covariance_norm)
    assert loss_floor > 3
    # If heading is injected about true down instead of nominal down, its z
    # component is 3/5. Other invariant axis groups add PSD information.
    true_axis_floor = F(9,25)*loss_floor
    assert true_axis_floor > 1
    return {
        'post_mag_covariance_upper': encoded(upper),
        'every_prefix_covariance_upper': encoded(prefix),
        'handoff_covariance': encoded(seed),
        'first_partial_cell_determinant_lower': str(ia*ib-ic*ic),
        'invariant_residual': encoded(invariant),
        'invariant_residual_ldlt_L': encoded(il),
        'invariant_residual_ldlt_D': [str(x) for x in ip],
        'all_phase_prefix_determinant_lower': str(a*b-c*c),
        'observation_gram_determinant': str(determinant),
        'observation_gram_trace_upper': str(trace_ceiling),
        'observation_noise_operator_upper': str(noise_ceiling),
        'normalized_root_covariance_upper': str(root_covariance_norm),
        'actual_innovation_service_lower': str(loss_floor),
        'true_heading_axis_service_lower': str(true_axis_floor),
        'required_service_floor': '1',
        'service_window_s': '1',
        'all_time_service_verified_in_real_arithmetic': True,
        'float32_all_time_service_verified': False,
    }


def witness_certificate():
    g, h, cosine, sine = F('9.80665'), F('.005'), F(3,5), F(4,5)
    # R maps body to world. The two signs are two separate persistent histories.
    rotation = [[F(1),F(0),F(0)], [F(0),cosine,-sine], [F(0),sine,cosine]]
    a = [F(0), g*sine, g*(1-cosine)]
    specific = [[a[0]], [a[1]], [a[2]-g]]
    body = matmul(transpose(rotation), specific)
    assert body == [[F(0)], [F(0)], [-g]]
    amplitude_squared = sum(x*x for x in a)
    # omega=2*pi/h > 6/h. These rational ceilings hold on the entire history.
    omega_lower = 6/h
    maxima_squared = {
        'acceleration': amplitude_squared,
        'velocity': amplitude_squared/omega_lower**2,
        'displacement': amplitude_squared/omega_lower**4,
        'primitive_difference': 4*amplitude_squared/omega_lower**6,
    }
    limits = {'acceleration': F('8.8'), 'velocity': F('5.5'),
              'displacement': F('8.1'), 'primitive_difference': F(1100)}
    assert all(maxima_squared[k] < limits[k]**2 for k in limits)
    return {
        'qualification': 'OU3_SAMPLED_PHYSICAL_CAPTURE_OBSTRUCTION_V1',
        'verified': True,
        'scope': 'declared six-degree physical tilt capture under regular real-arithmetic shipping execution',
        'sample_period_s': str(h),
        'physical_frequency_hz': '200',
        'body_to_world_rotation_plus': encoded(rotation),
        'acceleration_cosine_amplitude_plus': [str(x) for x in a],
        'tilt_cosine': str(cosine),
        'tilt_angle_exact': 'acos(3/5)',
        'all_time_squared_envelope_upper': {k:str(v) for k,v in maxima_squared.items()},
        'body_accelerometer_at_every_sample': ['0', '0', str(-g)],
        'body_gyro_at_every_sample': ['0', '0', '0'],
        'body_magnetometer_uT': ['75', '0', '0'],
        'physical_bias_and_sensor_error': 'zero',
        'both_histories_nonzero_motion': True,
        'same_physical_history_across_all_words': True,
        'service': service_certificate(),
        'declared_six_degree_eventual_capture_refuted': True,
        'conditional_local_A21_stability_refuted': False,
        'filter_divergence_claimed': False,
        'new_physical_assumption_adopted': False,
    }
