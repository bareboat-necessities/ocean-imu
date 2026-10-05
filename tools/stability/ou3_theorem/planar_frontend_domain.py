"""Analytical planar frontend domain, with explicit arithmetic qualifications.

This is a subordinate domain lemma for the complete finite-error inequality.
There is no replay, time loop, secant or sampled phase set. The MEKF nominal
domain and its complete-word contraction are NOT consequences of this lemma.
"""
from __future__ import annotations

from fractions import Fraction as F
from math import factorial

from .mahony_raw_normalization import U, certificate as norm_certificate
from .matrix_certificates import add, encoded, ldlt, matmul, transpose
from .planar_service_frontend_binding import rn
from .planar_service_guard import certificate as guard_certificate


def rounding_charge(magnitude):
    """RN error <= u |exact result| + min-subnormal/2."""
    return U*F(magnitude) + F(1, 2**150)


def one_step_charges():
    """Forward error of the planar scalar DAG; common normalization cancels.

    Bounds are on exact represented operands at the start of the operation.
    In particular we do not charge a generic quaternion-norm error as an
    attitude error. The same positive inverse-square-root multiplies w and y.
    """
    h, kp, ki = rn(F(1, 200)), rn(F(1, 5)), rn(F(1, 50))
    w, y = F(1001, 1000), F(14, 1250)
    ax, az = F(203, 10000), F(1001, 1000)
    # ax,az before the final component multiplies in accel normalization.
    eax, eaz = rounding_charge(ax), rounding_charge(az)
    hx, hz = w*y, (w*w+y*y)/2
    ehx = rounding_charge(hx)
    ehz = (rounding_charge(w*w)+rounding_charge(y*y)
            + rounding_charge(w*w+y*y+rounding_charge(w*w)
                              +rounding_charge(y*y)))/2
    # Products az*halfvx and ax*halfvz, then their subtraction. The linked
    # ideal difference is a*(delta-e), bounded by .51*(.00025+.0021).
    ep = eaz*hx+(az+eaz)*ehx+rounding_charge((az+eaz)*(hx+ehx))
    eq = eax*hz+(ax+eax)*ehz+rounding_charge((ax+eax)*(hz+ehz))
    feedback = F(51, 100)*(F(1, 4000)+F(21, 10000))
    ef = ep+eq+rounding_charge(feedback+ep+eq)
    e_ki = ki*ef+rounding_charge(ki*(feedback+ef))
    e_ih = h*e_ki+rounding_charge(h*(ki*feedback+e_ki))
    integral = F(17, 100000)
    ei = e_ih+rounding_charge(integral+ki*h*feedback+e_ih)
    omega = F(11, 1750)
    eg = rounding_charge(omega)
    egi = eg+ei+rounding_charge(omega+eg+integral+ki*h*feedback+ei)
    epf = kp*ef+rounding_charge(kp*(feedback+ef))
    rate = omega+integral+ki*h*feedback+kp*feedback
    erate = egi+epf+rounding_charge(rate+egi+epf)
    assert rate+erate < F(67, 10000)
    v = h*rate/2
    ev = h*erate/2+rounding_charge(h*(rate+erate)/2)
    ewy = y*ev+rounding_charge(y*(v+ev))
    eyw = w*ev+rounding_charge(w*(v+ev))
    ew = ewy+rounding_charge(w+y*v+ewy)
    ey = eyw+rounding_charge(y+w*v+eyw)
    # 2 atan(y/w), with w_after >= .997-y*v > .996. Straight-segment
    # differentiation gives |d pitch| <= 2(|dy|+.012|dw|)/.996.
    assert F(997, 1000)-y*v-ew > F(249, 250)
    assert (y+w*v+ey)/(F(997, 1000)-y*v-ew) < F(3, 250)
    euler_angle = 2*(ey+F(3, 250)*ew)/F(249, 250)
    # The final common positive scalar cancels from y/w. Two final RN
    # multiplies change the ratio by at most 2u/(1-u), plus subnormal charge.
    normalize_angle = 4*F(3, 250)*U/(1-U)+F(1, 2**130)
    h0 = F(1, 200)
    source_curvature = h0*h0*F(1, 50)*F(11, 35)**2/2
    step_mismatch = abs(h-h0)*omega
    arctan_defect = h**3*rate**3/12  # |2 atan(h rate/2)-h rate|
    pitch_charge = euler_angle+normalize_angle+source_curvature+step_mismatch+arctan_defect
    assert pitch_charge < F(46, 10**9)
    assert 10*ei < F(2, 10**9)
    return {'pitch_charge': pitch_charge, 'scaled_integral_charge': 10*ei,
            'feedback_roundoff': ef, 'corrected_rate_upper': rate+erate,
            'euler_direction_charge': euler_angle,
            'normalization_direction_charge': normalize_angle,
            'source_sampling_charge': source_curvature+step_mismatch,
            'arctan_charge': arctan_defect}


def guard_and_seed():
    """Rational induction margins; source derivation is in the appendix.

    sqrt and exp accuracy are explicit implementation qualifications. This
    does not assert that ISO C++ alone specifies a libm error bound.
    """
    physical = guard_certificate()
    # Source expression 2*pi*f*h, including represented pi/f/h and products.
    x_detect_lo = 2*F(314159, 100000)*25*F(1, 200)*(1-U)**5
    exp_detect_lo = sum((x_detect_lo**j/F(factorial(j)) for j in range(9)), F(0))
    assert (1+2*U)/exp_detect_lo < F(1, 2)
    x_rms_lo = 2*F(314159, 100000)*F(1, 20)*F(1, 200)*(1-U)**5
    x_rms_hi = 2*F(22, 7)*F(1, 20)*F(1, 200)*(1+U)**5
    beta_lo = 1-(1-x_rms_lo+x_rms_lo**2/2)*(1+2*U)-2*U
    beta_hi = 1-(1-x_rms_hi)*(1-2*U)+2*U
    assert F(1, 1000) < beta_lo < beta_hi < F(1, 100)
    delta = F(physical['successive_sample_force_difference_upper'])+20*U
    gamma = F(1, 2)
    lp_round = 100*U
    first_error = F(1, 1000)
    assert gamma*(first_error+delta)+lp_round < first_error
    first_output = F(11, 10000)
    assert first_error+rounding_charge(first_error) < first_output
    second_state = F(11, 5000)
    assert (first_output+second_state)/2+lp_round < second_state
    second_output = F(1, 250)
    assert first_output+second_state+rounding_charge(first_output+second_state) < second_output
    variance = F(1, 10000)
    rms_round = 20*U*variance+F(1, 2**140)
    assert F(1, 1000)*(variance-second_output**2) > rms_round
    assert 3*variance*(1+10*U) < F(3, 100)**2
    # Two vector normalizations use <=12u relative component error here;
    # the positive-dot FromTwoVectors branch uses sqrt(2(1+c)), reciprocal,
    # product. gamma_64 dominates its norm/direction error (appendix).
    gamma64 = 64*U/(1-64*U)
    dv = 12*U/(1-12*U)  # two planar vector normalizations
    ds = 12*U/(1-12*U)  # positive-dot quaternion formula
    direction_error = F(22, 1000)*(2*dv/(1-dv)+ds*(1+dv)/(1-dv))+F(1, 2**130)
    norm_error_before_round = (2*dv+dv**2)/F(39, 10)
    norm_error = norm_error_before_round+ds*(1+norm_error_before_round)+F(1, 2**130)
    assert direction_error < gamma64 and norm_error < gamma64
    assert gamma64 < F(1, 100000)
    seed_angle = F(1, 4000)+gamma64
    assert seed_angle < F(1, 1000)
    assert 10*seed_angle**2 < F(3, 500)**2
    assert 1-gamma64 > F(995, 1000) and 1+gamma64 < F(1001, 1000)
    return {'first_detector_unrounded_error_cap': str(first_error),
            'second_detector_component_cap': str(second_output),
            'per_axis_RMS_square_cap': str(variance),
            'seed_pitch_error_cap': str(seed_angle),
            'seed_squared_norm_error_cap': str(gamma64),
            'required_detector_gamma_interval': ['0', '1/2'],
            'required_RMS_beta_interval': ['1/1000', '1/100'],
            'derived_RMS_beta_bounds': list(map(str, (beta_lo, beta_hi))),
            'derived_seed_direction_roundoff': str(direction_error),
            'derived_seed_norm_roundoff': str(norm_error),
            'library_qualification': 'sqrt RN; exp relative error <=2u at the fixed detector/RMS/slew arguments'}


def certificate():
    h, kp, ki = rn(F(1, 200)), rn(F(1, 5)), rn(F(1, 50))
    g, cs, nu = F(196133, 20000), F(39999, 40001), F(11, 35)
    # Input component rounding changes planar direction by <= 4u here.
    acc_offset = F(1, 50)*nu**2/(g*cs)+4*U
    assert acc_offset < F(1, 4000)
    norm = norm_certificate()
    assert F(norm['raw_output_squared_norm_bounds'][0]) > F(995, 1000)
    assert F(norm['pre_component_scalar_squared_norm_bounds'][0]) > F(995, 1000)
    assert F(norm['pre_component_scalar_squared_norm_bounds'][1]) < F(1001, 1000)
    charges = one_step_charges()
    radius, margin = F(3, 500), F(27, 100000)
    G = [[F(10), F(-5)], [F(-5), F(15)]]
    endpoints = []
    for a in (F(49, 100), F(51, 100)):
        A = [[1-a*(kp*h+ki*h*h), h/10], [-10*ki*h*a, F(1)]]
        B = [[a*(kp*h+ki*h*h)], [10*ki*h*a]]
        gap = add(add(G, matmul(matmul(transpose(A), G), A), -1), G, -margin)
        _, pivots = ldlt(gap)
        assert matmul(matmul(transpose(B), G), B)[0][0] < F(1, 500)**2
        endpoints.append({'gain': str(a), 'positive_pivots': list(map(str, pivots))})
    # G^-1=[[3/25,1/25],[1/25,2/25]]. Coordinate bounds, no eigenvector fit.
    assert radius**2*F(3, 25) < F(21, 10000)**2
    assert radius**2*F(2, 25)/100 < F(17, 100000)**2
    # Raw q^2 and pre-component accelerometer normalization factor squared
    # in [.995,1.001]; sqrt(.995)>.995, sqrt(1.001)<1.001.
    gain_low = F(1, 2)*F(995, 1000)**2*(1-(F(21, 10000)+F(1, 4000))**2/6)
    gain_high = F(1, 2)*F(1001, 1000)**2
    assert F(49, 100) < gain_low <= gain_high < F(51, 100)
    charge = F(1, 500)*F(1, 4000)+F(19, 6)*charges['pitch_charge']+F(31, 8)*charges['scaled_integral_charge']
    assert charge < radius*margin/2
    # Raw pre-normalization norm: exact planar Euler action multiplies norm
    # squared by 1+(h*rate/2)^2. A deliberately coarse 1e-5 rounding charge
    # dominates the component-error bound above and leaves [1/4,4].
    pre_lo = F(995, 1000)-F(1, 100000)
    pre_hi = F(1001, 1000)*(1+(h*F(67, 10000)/2)**2)+F(1, 100000)
    assert F(1, 4) < pre_lo < pre_hi < 4
    # After update the reference at current packet phase differs by one
    # physical sample from the pre-update error coordinate.
    reference_angle = F(21, 10000)+F(1, 200)*F(11, 1750)
    reference_cap = F(11, 5000)
    assert reference_angle < reference_cap
    mu_x = 75*(1-reference_cap**2/2)
    return {
        'qualification': 'OU3_PLANAR_FRONTEND_DOMAIN_V1',
        'result_type': 'PROVED — analytical conditional domain implication',
        'scope': 'private frontend and real weighted reference geometry; not complete inherited MEKF domain',
        'arithmetic_premises': ['IEEE binary32 RN, gradual underflow, no unsafe reassociation',
            'same-record planar sensor samples are componentwise RN of the exact physical source',
            'sqrt is correctly rounded; exp relative error <=2u at fixed detector/RMS/slew arguments',
            'Eigen positive-dot FromTwoVectors scalar operation order as derived in the appendix'],
        'coordinates': 'pre-update e=2 atan2(q2,q0)-psi(t_k), z=(e,10 integralFBy)',
        'G': encoded(G), 'root_G_norm_radius': str(radius),
        'literal_h_kp_ki': list(map(str, (h, kp, ki))),
        'source_acc_angle_offset_upper': str(acc_offset),
        'effective_gain_bounds': list(map(str, (gain_low, gain_high))),
        'endpoint_certificates': endpoints, 'squared_norm_decrement': str(margin),
        'source_derived_one_step_charges': {k: str(v) for k, v in charges.items()},
        'G_norm_charge_upper': str(charge), 'strict_radial_margin': str(radius*margin/2-charge),
        'raw_pre_normalization_squared_norm_bounds': list(map(str, (pre_lo, pre_hi))),
        'pitch_error_cap': '21/10000', 'integral_cap': '17/100000',
        'all_prefix_frontend_retention_under_premises': True,
        'reference_geometry_real': {'angle_cap': str(reference_cap), 'Bx_lower': str(mu_x),
            'abs_Bz_upper': str(75*reference_cap), 'norm_upper': '75',
            'scope': 'any nonempty same-history nonnegative weighted mean; default zero hard iron; exact tilt normalization'},
        'guard_real_certificate': guard_certificate()['qualification'],
        'guard_and_seed_rational_margins': guard_and_seed(),
        'initialization_and_guard_transfer_under_arithmetic_premises': True,
        'target_toolchain_libm_qualification_verified': False,
        'reference_float_accumulation_and_acquisition_verified': False,
        'inherited_nominal_mekf_domain_verified': False,
        'complete_word_storage_contraction_verified': False,
        'structures_preserved': ['raw quaternion norm', 'literal new-integral-before-proportional update',
            'binary32 h/kp/ki', 'sample phase and gyro quadrature', 'inherited pitch/integral state'],
        'relaxations_introduced': ['rational directed error majorants of scalar operation DAG',
            'guard and initialization are explicit unmet arithmetic bindings, not erased'],
        'all_time_magnetic_service_verified': False, 'theorem_closed': False,
    }


if __name__ == '__main__':
    import json
    print(json.dumps(certificate(), indent=2, sort_keys=True))
