"""Exact constants for physical sampling and joint vector information.

Role in the finite-error proof: sampling defects are physical prediction
supply; vector information is an input to observability, not a replacement
for the full corrected 21-coordinate loss. See ou3-sampling-fidelity.md.
"""
from fractions import Fraction as F
import json
from pathlib import Path
from .lin_path_certificate import small_x_source_defect


def trapezoidal_weights(steps):
    """Nonuniform normalized endpoint weights for a window of complete cells."""
    steps = [F(h) for h in steps]
    if not steps or any(h <= 0 for h in steps):
        raise ValueError("positive cell lengths required")
    total = sum(steps)
    return ([steps[0]/(2*total)]
            + [(a+b)/(2*total) for a,b in zip(steps,steps[1:])]
            + [steps[-1]/(2*total)])


def vector_information_floor(epsilon, magnetic_residual, horizontal=F(1,5)):
    """Joint 3-D Gram floor; inputs use g and B_max normalization."""
    epsilon, magnetic_residual = F(epsilon), F(magnetic_residual)
    if epsilon < 0 or magnetic_residual < 0:
        raise ValueError("nonnegative residual bounds required")
    area = horizontal-epsilon-magnetic_residual*(1+epsilon)
    trace = (1+epsilon)**2+(1+magnetic_residual)**2
    return max(F(0),area)**2/trace


def certificate():
    c = json.loads(Path(__file__).with_name('constants.json').read_text(),parse_float=F)
    motion, sensor, bias, mag = (c[k] for k in
                               ('marine_motion','sensor_model','imu_bias','magnetic_service'))
    g, jerk, h = F('9.80665'), F(motion['J_max_mps3']), F(sensor['sample_period_max_s'])
    velocity = F(motion['V_max_mps'])
    acc_residual = F(bias['B_a_mps2'])+F(sensor['accel_fast_residual_norm_max_mps2'])
    # Both magnetic residual envelopes are charged; neither is discarded.
    mag_residual = F(mag['hard_iron_residual_norm_max_uT'])+F(mag['measurement_residual_norm_max_uT'])
    eta = jerk*h/4
    alias_window, information_window = F(32), F(64)
    alias_error = 2*velocity/alias_window+eta+acc_residual
    # pi > 3.14159; sin(x) >= x-x^3/6 on [0,pi/2]. Tilt 6 deg has half-angle pi/60.
    angle_lower = F('3.14159')/60
    six_degree_chord_lower = 2*g*(angle_lower-angle_lower**3/6)
    assert alias_error < six_degree_chord_lower
    epsilon = (2*velocity/information_window+eta+acc_residual)/g
    beta = mag_residual/F(mag['field_norm_max_uT'])
    horizontal = F(mag['horizontal_field_min_uT'])/F(mag['field_norm_max_uT'])
    floor = vector_information_floor(epsilon,beta,horizontal)
    assert floor > 0
    # Old witness: omega=400*pi > 1200 and ||A||^2=4g^2/5.
    old_jerk_squared_lower = F(1200)**2*4*g*g/5
    assert old_jerk_squared_lower > jerk*jerk
    # Physical OU forcing u=a'+lambda*a; this is a mismatch identity, not
    # an OU law on truth. Keep the source matrix comparison and polynomial
    # transition defect from the independently checked LIN path certificate.
    amplitude=F(motion['A_max_mps2'])
    sigma_min=F('.05')
    lambdas=[1/F(c['a21_schedule'][k]) for k in ('tau_min_s','tau_max_s')]
    action_rate=max((jerk+lam*amplitude)**2/(2*sigma_min**2*lam) for lam in lambdas)
    eps, transition_defect, _, _=small_x_source_defect()
    assert 0 <= eps < 1
    theta=F(1,100)
    source_action_rate=((1+theta)*action_rate+(1+1/theta)*amplitude**2*transition_defect/
                        F(sensor['sample_period_min_s']))/(1-eps)
    return {
        'qualification':'OU3_JERK_SAMPLING_FIDELITY_V1',
        'verified':True,
        'acceleration_locally_absolutely_continuous_required':True,
        'jerk_limit_mps3':str(jerk),
        'sample_period_max_s':str(h),
        'per_cell_acceleration_change_upper_mps2':str(jerk*h),
        'trapezoidal_mean_sampling_defect_upper_mps2':str(eta),
        'left_sample_mean_sampling_defect_upper_mps2':str(jerk*h/2),
        'physical_left_taylor_defects':{
            'velocity_mps':str(jerk*h*h/2),
            'displacement_m':str(jerk*h**3/6),
            'primitive_m_s':str(jerk*h**4/24),
        },
        'old_witness_jerk_squared_lower':str(old_jerk_squared_lower),
        'old_witness_excluded':True,
        'stationary_sample_alias':{
            'constant_true_attitude_required':True,
            'window_s':str(alias_window),
            'hidden_specific_force_upper_mps2':str(alias_error),
            'six_degree_chord_lower_mps2':str(six_degree_chord_lower),
            'six_degree_or_larger_tilt_excluded':True,
        },
        'joint_vector_information':{
            'coordinate_dimension':3,
            'window_min_s':str(information_window),
            'constant_world_magnetic_field_required':True,
            'accelerometer_world_transport_is_true_rotation':True,
            'magnetometer_rows_are_actually_applied_events':True,
            'accelerometer_residual_total_mps2':str(acc_residual),
            'magnetometer_residual_total_uT':str(mag_residual),
            'normalized_acceleration_mean_error_upper':str(epsilon),
            'normalized_magnetic_mean_error_upper':str(beta),
            'normalized_joint_gram_eigenvalue_lower':str(floor),
            'full_21_corrected_loss_certified':False,
        },
        'physical_LIN_prediction_supply':{
            'scope':'regular real-arithmetic default A21; actual realized tau and covariance',
            'input_identity':'u(t)=a_dot(t)+a(t)/tau_k on each physical cell',
            'ideal_action_rate_upper':str(action_rate),
            'source_action_rate_upper':str(source_action_rate),
            'source_small_x_relative_Q_defect_upper':str(eps),
            'source_transition_defect_action_per_acceleration_squared':str(transition_defect),
            'finite_error_prediction_metric_supply_verified':True,
            'physical_OU_prior_assumed':False,
            'whole_word_retained_radius_certified':False,
        },
        'shipping_capture_certified':False,
        'nonlinear_retention_certified':False,
        'quiet_water_admitted':True,
    }
