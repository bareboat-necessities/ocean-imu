"""Signed nominal-mean attitude columns and the literal AW correction loop.

Proof: docs/ou3-world-frame-rows.md, sections 8--9. Exact rational checks
on supplied words; the carried falsification of the pointwise premise is in
aw_tracking_source_diagnostic.py. Enters V_next<=rho V+c_d|d|^2 only through
the attitude columns of the aggregate six-column floor (c, hence B_*, J_AG,
the full covariance upper bound and rho_0). Non-promoting: no source-uniform
bound on the nominal AW mean is claimed.
"""
from fractions import Fraction as F

from .aw_covariance_ceiling import ceiling
from .lin_path_certificate import inverse, small_x_source_defect
from .matrix_certificates import add, identity, is_psd, matmul, transpose
from .world_frame import (_column, decimal_lower, injection_budget,
                          nominal_attitude_column_floor, quaternion_rotation, skew)

G = F('9.80665')
FRACTION = F(1, 5)
SIGMA_FLOOR = F(1, 20)      # max(0.05, band floor) in apply_ou_tune_
STEP_MAX = F(6, 1000)
JERK = F(100)


def corollary_a_star(mean_transverse, mean_norm, fraction=FRACTION, gravity=G):
    """Normalized attitude-column Gram floor from the NOMINAL signed mean only.

    Rows [u_i]x, u_i=(a_hat_i-g e_z)/g, convex weights alpha_i, one magnetic row
    [b]x, |b|=1, horizontal fraction |e_z x b|=fraction. With mu=sum alpha_i
    a_hat_i, m=|mu| and m_perp=|mu x b|:
      sum alpha_i[u_i]x'[u_i]x >= [u_bar]x'[u_bar]x (convexity),
      lambda_min([u]x'[u]x+[b]x'[b]x) >= |u x b|^2/(|u|^2+1),
      |u_bar x b| >= fraction-m_perp/g, |u_bar| <= 1+m/g.
    No physical acceleration, attitude or attitude-error term appears.
    """
    mp, m, frac, g = map(F, (mean_transverse, mean_norm, fraction, gravity))
    if mp < 0 or m < mp or g <= 0 or not 0 < frac <= 1:
        raise ValueError('0<=m_perp<=m, g>0, fraction in (0,1] required')
    margin = frac-mp/g
    return {'transverse_margin': margin,
            'normalized_gram_floor': max(F(0), margin)**2/((1+m/g)**2+1),
            'nominal_transverse_mean_threshold': g*frac, 'positive': margin > 0}


def attitude_gram(weights, nominal_aw, field_unit, gravity=G):
    """Exact normalized Gram of weighted nominal accelerometer rows plus [b]x."""
    gram = matmul(transpose(skew(field_unit)), skew(field_unit))
    for w, a in zip(weights, nominal_aw):
        u = [F(a[0])/gravity, F(a[1])/gravity, F(a[2])/gravity-1]
        gram = add(gram, matmul(transpose(skew(u)), skew(u)), F(w))
    return gram


def nonlocal_tracking_example():
    """Pointwise |a_hat| up to 8 m/s^2 with zero nominal mean still gives gamma*.

    The pointwise Corollary A premise |a_hat-a|<=1.12383 is not needed: any
    physical history whose tracking error equals these nominal values fails it.
    """
    b = [F(7, 25), F(0), F(24, 25)]           # unit, horizontal fraction 7/25
    aw = [[F(8), F(0), F(0)], [F(-8), F(0), F(0)], [F(0), F(6), F(1)],
          [F(0), F(-6), F(-1)]]
    weights = [F(1, 4)]*4
    mean = [sum(w*a[c] for w, a in zip(weights, aw)) for c in range(3)]
    assert mean == [0, 0, 0]
    star = corollary_a_star(0, 0, F(7, 25))
    gram = attitude_gram(weights, aw, b)
    assert is_psd(add(gram, identity(3), -star['normalized_gram_floor']))
    return {'pointwise_nominal_aw_max_mps2': '8', 'nominal_mean': [str(x) for x in mean],
            'gamma_star': str(star['normalized_gram_floor']),
            'exact_gram_dominates_gamma_star': True,
            'pointwise_premise_satisfied_by_any_zero_error_history': False}


def physical_transfer(window_s, gravity=G, fraction=FRACTION):
    """|mu|<=2V/T+J h/4+|sum alpha_i e_i|: Corollary A is the pointwise case.

    The signed mean error threshold equals Corollary A's pointwise threshold,
    because only sum alpha_i e_i enters: e_i may be arbitrarily large.
    """
    pointwise = nominal_attitude_column_floor(window_s, 0, fraction, gravity)['aw_error_threshold']
    physical = 2*F('5.5')/F(window_s)+JERK*STEP_MAX/4
    signed = gravity*fraction-physical
    assert signed == pointwise
    return {'window_s': str(window_s), 'physical_mean_supply_mps2': str(physical),
            'signed_mean_error_threshold_mps2': str(signed),
            'equals_pointwise_corollary_A_threshold': True}


def nominal_mean_representation(operations, initial):
    """Exact: sum alpha_i a_hat_i = W_0 a_hat_0 + sum_c W_c Delta_c.

    operations: ('row', alpha) samples the AW mean before an accelerometer
    correction; ('corr', Delta) adds any applied AW correction (acc, S, mag);
    ('pred', phi) is the OU prediction a_hat<-phi a_hat. W_c=sum over later rows
    of alpha times the product of later phi, so 0<=W_c<=sum alpha<=1 and the
    nominal mean is a nonnegatively weighted SIGNED sum of AW corrections.
    """
    x = [F(v) for v in initial]
    direct = [F(0)]*3
    coeffs = []                     # (index or None, running W)
    ops = list(operations)
    weights = [F(0)]*(len(ops)+1)
    # backward accumulation of W for each correction and for the root
    acc = F(0)
    for i in range(len(ops)-1, -1, -1):
        kind, value = ops[i]
        if kind == 'row':
            acc += F(value)
        elif kind == 'pred':
            acc *= F(value)
        weights[i] = acc
    root_weight = acc
    for i, (kind, value) in enumerate(ops):
        if kind == 'row':
            direct = [d+F(value)*v for d, v in zip(direct, x)]
        elif kind == 'pred':
            if not 0 <= F(value) <= 1:
                raise ValueError('OU factor in [0,1] required')
            x = [F(value)*v for v in x]
        elif kind == 'corr':
            x = [v+F(d) for v, d in zip(x, value)]
            coeffs.append((weights[i], [F(d) for d in value]))
        else:
            raise ValueError('unknown operation')
    represented = [root_weight*v for v in map(F, initial)]
    for w, d in coeffs:
        represented = [r+w*dv for r, dv in zip(represented, d)]
    assert represented == direct
    total = sum(F(v) for k, v in ops if k == 'row')
    assert all(0 <= w <= total for w, _ in coeffs) and 0 <= root_weight <= total
    return {'weighted_mean': direct, 'root_weight': root_weight,
            'correction_weights': [w for w, _ in coeffs]}


def aw_loop_identity():
    """Exact signed AW error identity on a supplied chronological word.

    acc: e+=(I-Gamma)e+Gamma eta; other corrections: e+=e+xi;
    prediction: e+=phi e-(1-phi)a-Delta a. Summing e-e+ over the word gives
      sum_acc Gamma(e-eta) = e_0-e_N+sum xi-sum_pred[(1-phi)a_hat+Delta a],
    with sum Delta a=a_N-a_0 telescoping the physical acceleration.
    """
    gam = [[[F(1, 3), F(1, 50), F(0)], [F(-1, 40), F(2, 7), F(1, 90)], [F(0), F(1, 60), F(1, 4)]],
           [[F(1, 20), F(0), F(1, 30)], [F(0), F(1, 25), F(0)], [F(-1, 70), F(0), F(1, 16)]]]
    eta = [[F(1, 10), F(-3, 20), F(1, 5)], [F(-1, 8), F(0), F(1, 40)]]
    xi = [F(1, 200), F(-1, 300), F(1, 100)]
    a = [[F(1), F(0), F(-1, 2)], [F(6, 5), F(1, 10), F(-2, 5)], [F(3, 2), F(1, 5), F(0)]]
    phis = [F(199, 200), F(99, 100)]
    ahat = [F(9, 10), F(1, 5), F(-1, 3)]
    e = [x-y for x, y in zip(ahat, a[0])]
    e0 = list(e)
    lhs = [F(0)]*3
    rhs_extra = [F(0)]*3
    word = [('acc', 0), ('pred', 0), ('corr', None), ('acc', 1), ('pred', 1)]
    k = 0
    for kind, i in word:
        if kind == 'acc':
            diff = [x-y for x, y in zip(e, eta[i])]
            g_diff = [sum(gam[i][r][c]*diff[c] for c in range(3)) for r in range(3)]
            lhs = [x+y for x, y in zip(lhs, g_diff)]
            e = [x-y for x, y in zip(e, g_diff)]
        elif kind == 'corr':
            e = [x+y for x, y in zip(e, xi)]
            rhs_extra = [x+y for x, y in zip(rhs_extra, xi)]
        else:
            phi, now, nxt = phis[i], a[k], a[k+1]
            hat = [x+y for x, y in zip(e, now)]
            rhs_extra = [x-(1-phi)*h-(y-z) for x, h, y, z in zip(rhs_extra, hat, nxt, now)]
            e = [phi*x-(1-phi)*y-(z-y) for x, y, z in zip(e, now, nxt)]
            k += 1
    rhs = [x-y+z for x, y, z in zip(e0, e, rhs_extra)]
    assert lhs == rhs
    return {'exact_residual': '0', 'physical_increments_telescope': True,
            'gain_weighted_not_trapezoidal': True}


def world_innovation_factorization():
    """r=f-R_hat(a_hat-g e_z)-b_hat equals R_hat(y-a_hat), y=R_hat'(f-b_hat)+g e_z.

    Here f=R(a-g e_z)+b_a+n with the TRUE attitude R, so y-a=(R_hat'R-I)(a-g e_z)
    +R_hat'(b_a-b_hat+n): attitude, bias and residual enter only through eta.
    """
    rh, rt = quaternion_rotation((7, 1, -2, 1)), quaternion_rotation((13, 1, -2, 2))
    a, ahat = [F(1, 2), F(-3, 4), F(1, 5)], [F(2), F(1, 3), F(-1, 4)]
    ba, bh, n = [F(1, 10), F(0), F(-1, 20)], [F(1, 20), F(1, 50), F(0)], [F(1, 100), F(-1, 200), F(0)]
    gz = [F(0), F(0), G]
    f = [x+y+z for x, y, z in zip([r[0] for r in matmul(rt, _column([p-q for p, q in zip(a, gz)]))], ba, n)]
    r = [x-y-z for x, y, z in zip(f, [q[0] for q in matmul(rh, _column([p-s for p, s in zip(ahat, gz)]))], bh)]
    y = [x+z for x, z in zip([q[0] for q in matmul(transpose(rh), _column([p-s for p, s in zip(f, bh)]))], gz)]
    assert r == [q[0] for q in matmul(rh, _column([p-s for p, s in zip(y, ahat)]))]
    return {'exact_residual': '0', 'aw_row_linear_in_mean': True}


def aw_increment_cap():
    """Literal AW correction d_a=K_a r obeys d_a d_a'<=NIS K_a S K_a'<=NIS P_aw.

    Same Cauchy--Schwarz/Joseph argument as Lemma I, on the AW rows. With the
    commanded sigma at its .05 floor and the Lemma B excess decayed,
    |Delta a_hat|<=sqrt((1+eps)/400 NIS) per applied correction.
    """
    p = [[F(4), F(1), F(0), F(1, 2)], [F(1), F(3), F(1, 3), F(0)],
         [F(0), F(1, 3), F(2), F(1, 5)], [F(1, 2), F(0), F(1, 5), F(5)]]
    h = [[F(1), F(0), F(2), F(1)], [F(0), F(1), F(-1), F(1, 2)], [F(1, 3), F(1), F(0), F(2)]]
    rmeas = [[F(1, 2), F(1, 10), F(0)], [F(1, 10), F(1), F(0)], [F(0), F(0), F(3, 2)]]
    s = add(matmul(h, matmul(p, transpose(h))), rmeas)
    k = matmul(p, matmul(transpose(h), inverse(s)))
    # Treat the last three state coordinates as the AW block of this word.
    budget = injection_budget(k[1:], s, [F(3), F(-2), F(5, 2)], [row[1:] for row in p[1:]])
    assert budget['gain_action_dominates'] and budget['prior_dominates_gain_action']
    assert budget['prior_dominates_injection']
    eps = small_x_source_defect()[0]
    floor_variance = (1+eps)*SIGMA_FLOOR**2
    # Tracking a jerk-limited ramp needs |Delta a|=J h per step from AW moves.
    required_nis = (JERK*STEP_MAX)**2/floor_variance
    return {'supplied_NIS': str(budget['NIS']), 'gain_action_dominates': True,
            'prior_aw_dominates_increment': True,
            'aw_variance_at_sigma_floor': str(floor_variance),
            'per_step_NIS_needed_to_follow_jerk_limit_at_floor': str(decimal_lower(required_nis)),
            'global_aw_variance_ceiling': str(ceiling(eps))}


def certificate():
    star16 = corollary_a_star(0, 0)
    old16 = nominal_attitude_column_floor(16, 0, FRACTION)
    ops = [('row', F(1, 4)), ('corr', [F(1, 5), F(0), F(-1, 10)]), ('pred', F(99, 100)),
           ('row', F(1, 2)), ('corr', [F(-3, 10), F(1, 20), F(0)]), ('pred', F(49, 50)),
           ('corr', [F(1, 7), F(0), F(1, 9)]), ('row', F(1, 4))]
    rep = nominal_mean_representation(ops, [F(1), F(-1, 2), F(1, 3)])
    return {
        'qualification': 'OU3_AW_NOMINAL_MEAN_V1',
        'scope': 'real-arithmetic regular A21 world-frame nominal rows; supplied exact words audit algebra only',
        'corollary_A_star': {
            'premise': 'trapezoidal/convex nominal AW mean mu over the window; one applied magnetic row',
            'nominal_transverse_mean_threshold_mps2': str(star16['nominal_transverse_mean_threshold']),
            'normalized_gram_floor_zero_mean': str(star16['normalized_gram_floor']),
            'corollary_A_zero_error_16s_floor': str(old16['normalized_gram_floor_without_injection']),
            'physical_acceleration_enters': False,
            'attitude_error_enters': False,
        },
        'nonlocal_tracking_example': nonlocal_tracking_example(),
        'physical_transfer_16s': physical_transfer(16),
        'pointwise_tracking_premise_necessary': False,
        'nominal_mean_signed_correction_representation': {
            'weighted_mean': [str(x) for x in rep['weighted_mean']],
            'root_weight': str(rep['root_weight']),
            'correction_weights': [str(x) for x in rep['correction_weights']],
            'weights_nonnegative_at_most_row_mass': True,
        },
        'aw_loop_signed_identity': aw_loop_identity(),
        'world_innovation_factorization': world_innovation_factorization(),
        'aw_increment_cap': aw_increment_cap(),
        'source_uniform_nominal_mean_bound': False,
        'uniform_AW_tracking_bound': False,
        'theorem_closed': False,
    }


if __name__ == '__main__':
    import json
    print(json.dumps(certificate(), indent=2, sort_keys=True))
