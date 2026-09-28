"""World-frame factorization of the historical AG rows.

Exact real-arithmetic algebra for regular A21: fixed deployment frame, zero
lever arm, constant committed reference, literal predictions, applied
corrections and resets. It enters V_next<=rho V+c_d|d|^2 only through the
historical reader: the six-column floor (c, a, b0), then B_*, J_AG, the full
covariance upper bound and rho_0<1. The injection budget also bounds the
injection-squared term of the finite-angle reset supply. Nothing here proves
the source-uniform premises; see docs/ou3-moving-six-pivots.md.
"""
from fractions import Fraction as F
from math import isqrt

from .matrix_certificates import add, identity, is_psd, matmul, transpose


def skew(v):
    x, y, z = map(F, v)
    return [[F(0), -z, y], [z, F(0), -x], [-y, x, F(0)]]


def _dot(u, v):
    return sum((F(a)*F(b) for a, b in zip(u, v)), F(0))


def _cross(u, v):
    return [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]


def _column(v):
    return [[F(x)] for x in v]


def quaternion_rotation(q):
    """Exact rotation of a nonzero rational quaternion (w,x,y,z), Eigen convention."""
    w, x, y, z = map(F, q)
    n = w*w+x*x+y*y+z*z
    if n == 0:
        raise ValueError('nonzero quaternion required')
    return [[(w*w+x*x-y*y-z*z)/n, 2*(x*y-w*z)/n, 2*(x*z+w*y)/n],
            [2*(x*y+w*z)/n, (w*w-x*x+y*y-z*z)/n, 2*(y*z-w*x)/n],
            [2*(x*z-w*y)/n, 2*(y*z+w*x)/n, (w*w-x*x-y*y+z*z)/n]]


def world_factorization(root_rotation, word):
    """Literal body rows versus world rows on one supplied word.

    word entries: ('predict', R_next, Rs, Bs), ('reset', R_next, d) or
    ('row', vector) with vector the WORLD nominal force a_hat-g e_z or the
    committed reference. Body AG transport follows the literal source
    ([[Rs,Bs],[0,I]] then G=I+[d]/2); the world transport uses only
    W=R_next' Rs R at predictions and N=R_next' G R at resets. Returns the
    exact residual of O=blockdiag(-R_k) Otilde diag(R_0',I) and every
    discrepancy factor. Predictions are world-neutral exactly when the
    covariance rotation equals the mean increment R_next R'.
    """
    rot = [[F(x) for x in r] for r in root_rotation]
    body_a, body_b = identity(3), [[F(0)]*3 for _ in range(3)]
    world_a, world_b = identity(3), [[F(0)]*3 for _ in range(3)]
    residual, factors, rows = F(0), [], 0
    for item in word:
        kind = item[0]
        if kind == 'predict':
            nxt, rs, bs = ([[F(x) for x in r] for r in m] for m in item[1:])
            body_a, body_b = matmul(rs, body_a), add(matmul(rs, body_b), bs)
            w = matmul(transpose(nxt), matmul(rs, rot))
            world_a = matmul(w, world_a)
            world_b = add(matmul(w, world_b), matmul(transpose(nxt), bs))
            factors.append(('predict', w))
            rot = nxt
        elif kind == 'reset':
            nxt = [[F(x) for x in r] for r in item[1]]
            g = add(identity(3), skew(item[2]), F(1, 2))
            body_a, body_b = matmul(g, body_a), matmul(g, body_b)
            n = matmul(transpose(nxt), matmul(g, rot))
            world_a, world_b = matmul(n, world_a), matmul(n, world_b)
            factors.append(('reset', n, matmul(transpose(rot), _column(item[2]))))
            rot = nxt
        elif kind == 'row':
            rows += 1
            f = [F(x) for x in item[1]]
            body_f = [r[0] for r in matmul(rot, _column(f))]
            h = [[-x for x in r] for r in skew(body_f)]
            lhs = [a + b for a, b in zip(matmul(h, body_a), matmul(h, body_b))]
            m = [[-x for x in r] for r in matmul(rot, skew(f))]
            rhs = [a + b for a, b in zip(matmul(m, matmul(world_a, transpose(root_rotation))),
                                         matmul(m, world_b))]
            residual = max([residual] + [abs(x-y) for p, q in zip(lhs, rhs) for x, y in zip(p, q)])
        else:
            raise ValueError('unsupported world-frame word entry')
    return {'rows': rows, 'row_factorization_residual': residual,
            'world_attitude_transport': world_a, 'world_gyro_transport': world_b,
            'factors': factors}


def reset_discrepancy_gram(x):
    """N'N for N=Exp(-[x])(I+[x]/2): exactly I+(|x|^2 I-xx')/4 >= I."""
    x = [F(v) for v in x]
    outer = matmul(_column(x), [x])
    return add(identity(3), add([[_dot(x, x)*F(i == j) for j in range(3)] for i in range(3)],
                                outer, F(-1)), F(1, 4))


def reset_angle_upper(theta):
    """Rational upper bound on the angle between N(x)^-1 b and b, |x|=theta<4.

    N^-1 rotates the transverse part by theta-atan(theta/2)<=theta/2+theta^3/24
    and scales it by rho=(1+theta^2/4)^(-1/2); the scaling turns a line by at
    most atan((1-rho)/(2 sqrt rho))<=theta^2/(16-theta^2).
    """
    t = F(theta)
    if not 0 <= t < 4:
        raise ValueError('injection norm in [0,4) required')
    return t/2 + t**3/24 + t*t/(16-t*t)


def reset_distance_upper(theta):
    """|N(x)-I|<=theta/2+theta^2+theta^3/4 for |x|=theta<=2."""
    t = F(theta)
    if not 0 <= t <= 2:
        raise ValueError('injection norm in [0,2] required')
    return t/2 + t*t + t**3/4


def _sqrt_lower(value, digits=12):
    """Decimal rational q<=sqrt(value), value>=0."""
    value = F(value)
    if value < 0:
        raise ValueError('nonnegative value required')
    scale = 10**digits
    scaled = value*scale*scale
    return F(isqrt(scaled.numerator//scaled.denominator), scale)


def decimal_lower(value, digits=12):
    """Decimal rational q<=value, retained for readable exact reports."""
    value = F(value)
    scale = 10**digits
    return F((value.numerator*scale)//value.denominator, scale)


def same_cell_world_floor(force_world, field_world, injection_upper):
    """Rational L<=sigma_min(C)^2 for [skew(R f); skew(R_2 b)G], any attitude R.

    C'C is congruent by rotations to skew(f)'skew(f)+N'skew(b)'skew(b)N,
    so only the WORLD nominal force, committed reference and |d| enter.
    With s the sine of the f/b line angle and psi the reset angle bound,
    L=|f|^2|b|^2 (s-psi)_+^2/(|f|^2+|b|^2) (least eigenvalue >= product/sum).
    """
    f, b = [F(x) for x in force_world], [F(x) for x in field_world]
    ff, bb = _dot(f, f), _dot(b, b)
    if ff <= 0 or bb <= 0:
        raise ValueError('nonzero nominal force and reference required')
    sine_squared = _dot(_cross(f, b), _cross(f, b))/(ff*bb)
    margin = _sqrt_lower(sine_squared) - reset_angle_upper(injection_upper)
    return ff*bb*max(F(0), margin)**2/(ff+bb)


def injection_budget(gain_theta, innovation_covariance, innovation, prior_theta=None):
    """Exact Loewner checks for the literal attitude injection d=K_theta r.

    d d' <= NIS K_theta S K_theta' holds for any gain (Cauchy--Schwarz).
    For the actual Kalman gain K=PH'S^-1 with S-HPH'>=0, Joseph gives
    K S K'<=P, hence d d'<=NIS P_theta,theta (prior). Returns NIS and flags.
    """
    k = [[F(x) for x in r] for r in gain_theta]
    s = [[F(x) for x in r] for r in innovation_covariance]
    r = [F(x) for x in innovation]
    from .lin_path_certificate import inverse
    nis = matmul([r], matmul(inverse(s), _column(r)))[0][0]
    d = matmul(k, _column(r))
    dd = matmul(d, transpose(d))
    ksk = matmul(k, matmul(s, transpose(k)))
    out = {'NIS': nis, 'injection': [x[0] for x in d],
           'gain_action_dominates': is_psd(add([[nis*x for x in row] for row in ksk], dd, F(-1)))}
    if prior_theta is not None:
        p = [[F(x) for x in row] for row in prior_theta]
        out['prior_dominates_gain_action'] = is_psd(add(p, ksk, F(-1)))
        out['prior_dominates_injection'] = is_psd(add([[nis*x for x in row] for row in p], dd, F(-1)))
    return out


def nominal_attitude_column_floor(window_s, aw_error_upper, reference_horizontal_fraction,
                                  gravity=F('9.80665'), velocity=F('5.5'), jerk=F(100),
                                  step_max=F('0.006')):
    """Normalized world Gram floor of the NOMINAL attitude columns on a window.

    Rows: acc [(a_hat_i-g e_z)/g]_x with trapezoidal weights, and one applied
    magnetic row [b_hat_w]_x with weight one. |sum alpha_i a_i|<=2V/T+J h/4
    (sampling fidelity) and |a_hat-a|<=aw_error in the estimator world frame.
    No attitude or attitude-error term enters. Reset discrepancies A_tilde
    with |A_tilde-I|<=beta lower sqrt(floor) by at most beta*sqrt(u_max^2+1),
    u_max the normalized row-force ceiling; that perturbation is separate.
    """
    t, err, frac, g = map(F, (window_s, aw_error_upper, reference_horizontal_fraction, gravity))
    if min(t, g) <= 0 or err < 0 or not 0 < frac <= 1:
        raise ValueError('positive window/gravity, nonnegative AW error, fraction in (0,1]')
    physical = 2*F(velocity)/t + F(jerk)*F(step_max)/4
    eps = (physical+err)/g
    gamma = max(F(0), frac-eps)**2/((1+eps)**2+1)
    return {'normalized_mean_force_error': eps, 'transverse_margin': frac-eps,
            'normalized_gram_floor_without_injection': gamma,
            'aw_error_threshold': g*frac-physical, 'positive': frac > eps}


def collinear_cadence_floor(horizontal_fraction, window_s, gravity=F('9.80665'),
                            velocity=F('5.5'), jerk=F(100)):
    """Lower bound on sum(D_k^2)/sum(D_k) for physical force parallel to b.

    If a(t_k)-g e_z is parallel to b at instants with gaps D_k covering a
    window of length L, then c=g(e_z-(e_z.b)b) has c.a(t_k)=|c|, and the jerk
    bound gives c.a>=|c|-J min(t-t_k,t_(k+1)-t) between them. Integrating and
    using |v(t_n)-v(t_0)|<=2V gives |c| L-J sum D_k^2/4<=2V.
    """
    h, length, g, v, j = map(F, (horizontal_fraction, window_s, gravity, velocity, jerk))
    if not 0 < h <= 1 or min(length, g, v, j) <= 0:
        raise ValueError('horizontal fraction in (0,1] and positive limits required')
    return max(F(0), 4*(g*h-2*v/length)/j)


def collinear_same_cell_witness():
    """MARINE MOTION/IMU BIAS history whose force is parallel to B at 1-Hz instants.

    B=(21,0,72) uT, c=g(e_z-(e_z.b)b), a(t)=c cos(2 pi t), zero biases,
    roll (1/100)sin(t/2). Exact rational checks with 3.14159<pi<3.1416.
    One applied magnetic correction per 1-s window is NOT claimed to satisfy
    MAGNETIC SERVICE: its heading/axial-bias information is nearly rank one.
    """
    import json
    from pathlib import Path
    limits = json.loads(Path(__file__).with_name('constants.json').read_text(), parse_float=F)
    motion, mag = limits['marine_motion'], limits['magnetic_service']
    g, pi_lo, pi_hi = F('9.80665'), F('3.14159'), F('3.1416')
    field = [F(21), F(0), F(72)]
    norm = F(75)
    unit = [x/norm for x in field]
    ez = [F(0), F(0), F(1)]
    c = [g*(e-_dot(ez, unit)*u) for e, u in zip(ez, unit)]
    c_norm = g*F(7, 25)
    assert _dot(c, c) == c_norm**2 and _dot(c, unit) == 0
    integer_force = [a-g*e for a, e in zip(c, ez)]
    half_force = [-a-g*e for a, e in zip(c, ez)]
    group = skew(integer_force) + skew(field)
    kernel = matmul(group, _column(unit))
    roll, rate = F(1, 100), F(1, 2)
    checks = {
        'field_norm_uT': norm == 75 and _dot(field, field) == norm**2,
        'horizontal_fraction_admitted': F(21, 75) >= F(mag['horizontal_field_min_uT'])/F(mag['field_norm_max_uT']),
        'acceleration_upper': c_norm <= F(motion['A_max_mps2']),
        'jerk_upper': 2*pi_hi*c_norm <= F(motion['J_max_mps3']),
        'velocity_upper': c_norm/(2*pi_lo) <= F(motion['V_max_mps']),
        'displacement_upper': c_norm/(4*pi_lo**2) <= F(motion['P_max_m']),
        'primitive_upper': c_norm/(8*pi_lo**3) <= F(motion['P_AC_max_m_s']),
        'attitude_rate_upper': roll*rate <= F(motion['Omega_max_rad_s']),
        'integer_time_force_parallel_field': _cross(integer_force, unit) == [0, 0, 0],
        'ideal_same_cell_kernel_along_field': kernel == [[0]]*6,
        'half_integer_transverse_force_squared': _dot(_cross(half_force, unit), _cross(half_force, unit)) == (2*c_norm)**2,
    }
    assert all(checks.values()), checks
    return {'field_world_uT': [str(x) for x in field],
            'acceleration_amplitude_vector_mps2': [str(x) for x in c],
            'acceleration_amplitude_mps2': str(c_norm),
            'jerk_upper_mps3': str(2*pi_hi*c_norm),
            'roll_amplitude_rad': str(roll), 'roll_frequency_rad_s': str(rate),
            'excited_windows': 'complete T_E>=4 pi windows have gravity span 1/50',
            'half_integer_transverse_force_mps2': str(2*c_norm),
            'marine_motion_and_imu_bias_admitted': True,
            'magnetic_service_admission_claimed': False,
            **{k: v for k, v in checks.items()}}


def certificate():
    """Exact supplied algebra audits and the physical witness; no uniform premise."""
    # Supplied word with rational rotations. The covariance prediction uses the
    # exact mean increment, so predictions are world-neutral; resets are not.
    r0 = quaternion_rotation((3, 1, -1, 2))
    step = quaternion_rotation((12, 1, 0, 1))
    r1 = matmul(step, r0)
    bs = [[F(1, 200), F(-1, 5000), F(1, 7000)], [F(1, 5000), F(1, 200), F(0)],
          [F(-1, 7000), F(0), F(1, 200)]]
    d = [F(1, 50), F(-1, 20), F(1, 40)]
    inj = quaternion_rotation((40, F(1, 2), F(-5, 4), F(5, 8)))
    r2 = matmul(inj, r1)
    step2 = quaternion_rotation((25, 0, 1, -1))
    r3 = matmul(step2, r2)
    word = [('row', [F(1, 3), F(-1, 4), F(-49, 5)]), ('row', [F(21), F(0), F(72)]),
            ('predict', r1, step, bs), ('row', [F(-1, 2), F(0), F(-19, 2)]),
            ('reset', r2, d), ('row', [F(21), F(0), F(72)]),
            ('predict', r3, step2, bs), ('row', [F(1), F(1, 5), F(-10)])]
    result = world_factorization(r0, word)
    assert result['row_factorization_residual'] == 0
    predictions = [f[1] for f in result['factors'] if f[0] == 'predict']
    assert all(w == identity(3) for w in predictions)
    reset = [f for f in result['factors'] if f[0] == 'reset'][0]
    x = [row[0] for row in reset[2]]
    gram = matmul(transpose(reset[1]), reset[1])
    # Rational injection rotation differs from Exp([d]); the Gram identity
    # holds for ANY orthogonal mean injection, hence sigma_min(N)>=1.
    assert gram == reset_discrepancy_gram(x)
    assert is_psd(add(gram, identity(3), F(-1)))
    # Same-cell world floor: congruence identity on the supplied acc/mag group.
    f_world, b_world = [F(1, 3), F(-1, 4), F(-49, 5)], [F(21), F(0), F(72)]
    rot, rot2 = r0, matmul(inj, r0)
    body = ([[-v for v in row] for row in skew([r[0] for r in matmul(rot, _column(f_world))])]
            + matmul([[-v for v in row] for row in skew([r[0] for r in matmul(rot2, _column(b_world))])],
                     add(identity(3), skew(d), F(1, 2))))
    n = matmul(transpose(rot2), matmul(add(identity(3), skew(d), F(1, 2)), rot))
    world = skew(f_world) + matmul(skew(b_world), n)
    assert matmul(transpose(body), body) == matmul(rot, matmul(matmul(transpose(world), world), transpose(rot)))
    injection = F(3, 50)
    assert _dot(d, d) <= injection**2
    floor = decimal_lower(same_cell_world_floor(f_world, b_world, injection))
    assert floor > 0 and is_psd(add(matmul(transpose(world), world), identity(3), -floor))
    # Injection budget on a supplied exact Kalman correction.
    p = [[F(4), F(1), F(0), F(1, 2)], [F(1), F(3), F(1, 3), F(0)],
         [F(0), F(1, 3), F(2), F(1, 5)], [F(1, 2), F(0), F(1, 5), F(5)]]
    h = [[F(1), F(0), F(2), F(1)], [F(0), F(1), F(-1), F(1, 2)], [F(1, 3), F(1), F(0), F(2)]]
    rmeas = [[F(1, 2), F(1, 10), F(0)], [F(1, 10), F(1), F(0)], [F(0), F(0), F(3, 2)]]
    s = add(matmul(h, matmul(p, transpose(h))), rmeas)
    from .lin_path_certificate import inverse
    k = matmul(p, matmul(transpose(h), inverse(s)))
    budget = injection_budget(k[:3], s, [F(3), F(-2), F(5, 2)], [row[:3] for row in p[:3]])
    assert budget['gain_action_dominates'] and budget['prior_dominates_gain_action']
    assert budget['prior_dominates_injection']
    # Generalized exact quiet subcase: any constant attitude and any admitted
    # reference. f=-g e_z, sine=horizontal fraction>=1/5, g>=9, |B|>=20.
    quiet_c_squared = F(81*400, 81+400)*F(1, 25)
    quiet_c = F(8, 5)
    assert quiet_c**2 <= quiet_c_squared
    from .moving_pivots import six_pivot_floor_valid, two_group_six_column_floor
    quiet = two_group_six_column_floor(quiet_c, F(4, 125), 1, 0)
    quiet_s = quiet['six_column_singular_lower']
    assert quiet_s == F(16, 635) and six_pivot_floor_valid(quiet_s, 12, quiet_s/4)
    column = nominal_attitude_column_floor(16, 0, F(1, 5))
    threshold16 = column['aw_error_threshold']
    long_window = nominal_attitude_column_floor(10**6, 0, F(1, 5))['aw_error_threshold']
    from .aw_covariance_ceiling import INHERITED_STD, ceiling, storage_radius
    from .lin_path_certificate import small_x_source_defect
    from .nuisance_upper_certificate import bounds as nuisance_bounds
    # |e_aw|^2<=lambda_max(P_aw) V. The isotropic-sync ceiling
    # lambda_max(P_aw)<=(1+eps)16 replaces the inherited 156^2 (still bound here).
    assert nuisance_bounds()[3][3] == INHERITED_STD
    eps = small_x_source_defect()[0]
    witness = collinear_same_cell_witness()
    cadence = collinear_cadence_floor(F(1, 5), 16)
    assert cadence > F(1, 25)
    return {
        'qualification': 'OU3_WORLD_FRAME_HISTORICAL_ROWS_V1',
        'scope': 'real-arithmetic regular A21: fixed deployment frame, zero lever arm, constant committed reference; supplied exact words audit algebra only',
        'world_row_factorization_exact_residual': str(result['row_factorization_residual']),
        'predictions_world_neutral_when_covariance_rotation_matches_mean': True,
        'reset_discrepancy_gram_identity': True,
        'reset_discrepancy_sigma_min_at_least_one': True,
        'attitude_estimate_enters_six_column_geometry': False,
        'attitude_error_enters_nominal_geometry_transfer': False,
        'same_cell_supplied_world_floor': str(floor),
        'injection_budget_supplied_NIS': str(budget['NIS']),
        'injection_budget_gain_action': True,
        'injection_budget_prior_attitude_covariance': True,
        'quiet_any_attitude_any_admitted_reference': {
            'scope': 'zero-residual quiet nominal record, any constant nominal attitude, committed reference |B|>=20 uT with horizontal fraction>=1/5, g>=9, acc/mag groups eight qualified predictions apart',
            'group_singular_floor_squared': str(quiet_c_squared),
            'group_singular_floor_used': str(quiet_c),
            'two_group_six_column_floor': str(quiet_s),
            'all_six_pivot_floor_12_rows': str(quiet_s/4),
            'covariance_action_generalized': False,
        },
        'nominal_attitude_column_transfer': {
            'window_s': '16',
            'normalized_gram_floor_zero_aw_error': str(column['normalized_gram_floor_without_injection']),
            'aw_error_threshold_16s_mps2': str(threshold16),
            'aw_error_threshold_long_window_mps2': str(long_window),
            'aw_variance_ceiling': str(ceiling(eps)),
            'retained_sqrt_V_upper_for_positive_16s_transfer': str(decimal_lower(storage_radius(threshold16, eps, 16))),
            'inherited_156_retained_sqrt_V_upper': str(threshold16/INHERITED_STD),
            'gyro_columns_certified': False,
        },
        'collinear_same_cell_witness': witness,
        'jerk_limited_collinear_cadence': {
            'horizontal_fraction': '1/5', 'window_s': '16',
            'length_weighted_mean_gap_floor_s': str(cadence),
            'regular_25Hz_all_collinear_excluded': cadence > F(1, 25),
        },
        'same_cell_floor_from_motion_and_bias_bounds_alone': False,
        'same_cell_route_refuted_under_three_assumptions': False,
        'aggregate_rows_retain_transverse_force_on_witness': True,
        'source_uniform_six_pivots': False,
        'uniform_historical_AG_readout_action': False,
        'rho0_certified': False,
        'theorem_closed': False,
    }
