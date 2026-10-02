"""Conditional six-column bridge retaining chronological transport defects.

Exact algebra subordinate to the historical-reader covariance comparison.
Supplied scalar bounds are premises; this module does not obtain them from
physical span or convert finite diagnostics into all-history certificates.
"""
from fractions import Fraction as F
from operator import index


def same_cell_geometry_squared_floor(force_floor, field_floor, cosine_upper):
    """Conditional raw row floor for [skew(f); skew(b) G].

    G is a chronological product of literal resets, sigma_min(G)>=1.
    The cosine is between f and G^-1 b, NOT between f and b. These
    actual-history bounds are premises, not new physical assumptions.
    """
    force, field, cosine = map(F, (force_floor, field_floor, cosine_upper))
    if min(force, field) <= 0 or not 0 <= cosine <= 1:
        raise ValueError('positive vector floors and absolute cosine in [0,1] required')
    return min(force**2, field**2) * (1-cosine)


def same_prediction_cell_groups(events):
    """Exact raw AG rows of actual acc/mag pairs with no intervening prediction.

    Events use the existing source observer format. All applied corrections
    stay in the optimal covariance comparison; their resets act literally on
    the auxiliary AG coordinates. The zero-lever-arm acc/mag BG rows must be
    zero. Sync and the read-only adaptive-state record have no AG transition.
    Unsupported hard events fail closed.
    """
    from .matrix_certificates import identity, matmul
    transport = identity(6)
    anchor = None
    groups = []
    for number, event in enumerate(events):
        kind = event['kind']
        if kind == 'prediction':
            f = [[F(x) for x in row] for row in event['F_AG']]
            if len(f) != 6 or any(len(r) != 6 for r in f):
                raise ValueError('six AG prediction coordinates required')
            if f[3:] != identity(6)[3:]:
                raise ValueError('literal constant gyro-error prediction required')
            transport = matmul(f, transport)
            anchor = None
        elif kind == 'reset':
            d = [F(row[0]) for row in event['d']]
            if len(d) != 3:
                raise ValueError('three reset coordinates required')
            x, y, z = d
            g = [[F(1), -z/2, y/2], [z/2, F(1), -x/2], [-y/2, x/2, F(1)]]
            f = identity(6)
            for i in range(3):
                f[i][:3] = g[i]
            transport = matmul(f, transport)
            if anchor is not None:
                anchor['local'] = matmul(g, anchor['local'])
        elif kind == 'correction':
            sensor = event['sensor']
            if sensor == 'S':
                continue
            if sensor not in ('acc', 'mag'):
                raise ValueError('unsupported applied observation')
            h = [[F(x) for x in row[:6]] for row in event['H']]
            if len(h) != 3 or any(len(r) != 6 or any(r[3:]) for r in h):
                raise ValueError('zero direct gyro columns required for same-cell factorization')
            if sensor == 'acc':
                anchor = {'C_acc': [r[:3] for r in h], 'local': identity(3),
                          'transport': transport, 'rows': matmul(h, transport),
                          'event': number}
            elif anchor is not None:
                c = anchor['C_acc'] + matmul([r[:3] for r in h], anchor['local'])
                raw = anchor['rows'] + matmul(h, transport)
                factored = matmul(c, anchor['transport'][:3])
                if raw != factored:
                    raise ArithmeticError('same-cell raw row factorization failed')
                groups.append({'C': c, 'anchor_transport': anchor['transport'],
                               'raw_rows': raw, 'acc_event': anchor['event'],
                               'mag_event': number})
        elif kind not in ('sync', 'sync_completion', 'adaptive_state'):
            raise ValueError('unretained hard event in regular AG group')
    return groups


def inverse_frame_gyro_prefixes(operations):
    """Literal real-source AG transport with C=A^-1 B, at every prefix.

    Events are ('predict', h, upper_norm_w_times_h) or ('reset', upper_norm_d).
    Both source rotation branches and G=I+[d]/2 have inverse norm <=1.
    Resets leave C unchanged. Their effect on later predictions is retained
    in beta>=||A^-1-I||, not charged against already accumulated gyro action.
    The supplied injection bounds are not obtained from physical excitation.
    """
    beta, defect, elapsed, attitude = F(0), F(0), F(0), F(1)
    out = []
    for op in operations:
        if len(op) == 3 and op[0] == 'predict':
            h, angle = map(F, op[1:])
            if h <= 0 or angle < 0:
                raise ValueError('positive step and nonnegative angle bound required')
            # ||R^-1 D-hI||: Rodrigues integral <=h*angle/2;
            # polynomial R^-1(D-hR) adds at most h*angle^2/3.
            defect += h * beta + h * (angle / 2 + angle**2 / 3)
            elapsed += h
            beta = min(F(2), beta + angle + angle**2 / 2)
            attitude *= 1 + angle**4 / 8
        elif len(op) == 2 and op[0] == 'reset':
            injection = F(op[1])
            if injection < 0:
                raise ValueError('nonnegative injection bound required')
            beta = min(F(2), beta + min(F(1), injection / 2))
            attitude *= 1 + injection**2 / 8
        else:
            raise ValueError('one literal prediction or reset per event required')
        out.append({'elapsed_time': elapsed, 'inverse_identity_defect_upper': beta,
                    'normalized_gyro_defect_upper': defect,
                    'gyro_singular_lower': max(F(0), elapsed - defect),
                    'attitude_norm_upper': attitude,
                    'inverse_attitude_norm_upper': F(1)})
    return out


def gyro_transport_prefixes(operations):
    """op=(h, ||R|| ceiling, ||R-I|| ceiling, ||D-hI|| ceiling).

    h=0 represents a reset and requires the fourth entry zero. Observations
    themselves are identity transitions in the latent historical array.
    h and D must be expressed in the same fixed gyro coordinate scaling.
    """
    a, b, t = F(1), F(0), F(0)
    rows = []
    for h, r, u, v in operations:
        h, r, u, v = map(F, (h, r, u, v))
        if min(h, r, u, v) < 0 or r == 0 or (h == 0 and v != 0):
            raise ValueError("invalid chronological transport bounds")
        a, b, t = r * a, r * b + u * t + v, t + h
        rows.append({"attitude_norm_upper": a, "gyro_identity_defect_upper": b,
                     "elapsed_scaled_time": t, "gyro_singular_lower": max(F(0), t - b)})
    return rows


def two_group_six_column_floor(attitude_row_floor, gyro_floor,
                                attitude_transport_upper, row_defect_upper):
    c, b, a, epsilon = map(F, (attitude_row_floor, gyro_floor,
                               attitude_transport_upper, row_defect_upper))
    if min(c, b) <= 0 or min(a, epsilon) < 0:
        raise ValueError("positive row/gyro floors and nonnegative ceilings required")
    inverse_upper = 1 + (a + 1) / b
    margin = c / inverse_upper - epsilon
    return {"block_inverse_upper": inverse_upper, "signed_margin": margin,
            "six_column_singular_lower": max(F(0), margin),
            "six_columns_certified_conditionally": margin > 0}


def six_pivot_floor_valid(singular_floor, row_count, proposed_pivot_floor):
    s, p = map(F, (singular_floor, proposed_pivot_floor))
    m = index(row_count)
    if s <= 0 or p <= 0 or m < 6:
        raise ValueError("positive floors and at least six actual rows required")
    return m * p * p <= s * s


def certificate():
    from .matrix_certificates import add, identity, ldlt, matmul, transpose
    # Supplied normalized exact example, NOT an observed shipping word.
    # C0=C1=I, A=I, B=I/2, E diagonal 1/100. Include a reset in the
    # separate scalar recurrence to audit chronological defect accounting.
    ops = [(F(1, 4), F(1), F(0), F(0)),
           (F(0), F(11, 10), F(1, 10), F(0)),
           (F(1, 4), F(1), F(0), F(0))]
    transport = gyro_transport_prefixes(ops)[-1]
    bound = two_group_six_column_floor(1, F(1, 2), 1, F(1, 100))
    s = bound["six_column_singular_lower"]
    o = [[F(0) for _ in range(6)] for _ in range(6)]
    for i in range(3):
        o[i][i] = F(101, 100)
        o[i + 3][i] = F(1)
        o[i + 3][i + 3] = F(51, 100)
    gram = matmul(transpose(o), o)
    _, pivots = ldlt(add(gram, identity(6), -s*s))
    assert all(x > 0 for x in pivots)
    p = s / 3  # 6/9 <= 1, no floating square root in the certificate.
    assert six_pivot_floor_valid(s, 6, p)
    failed = two_group_six_column_floor(1, F(1, 2), 1, F(1, 4))
    assert not failed["six_columns_certified_conditionally"]
    # Exact quiet nominal subcase: C=[skew(g ez);skew(B ex)], with
    # g>=9, B>=9; two coincident acc/mag groups eight qualified predictions apart, zero
    # corrections/resets. C'C >=81 I, so this is an analytic bound,
    # not an eigenvalue sampled from a Gramian or a physical bias theorem.
    quiet = two_group_six_column_floor(9, F(4, 125), 1, 0)
    quiet_s = quiet["six_column_singular_lower"]
    quiet_p = quiet_s / 4
    assert quiet_s == F(18, 127) and six_pivot_floor_valid(quiet_s, 12, quiet_p)
    inverse = inverse_frame_gyro_prefixes([
        ('predict', F(1, 200), 0), ('reset', 4)])
    assert inverse[-1]['gyro_singular_lower'] == F(1, 200)
    # Exact relaxed inter-anchor transport: not a reachable shipping word,
    # not a loss of all intermediate observation information.
    cancellation = [[F(3, 625) * x for x in row] for row in identity(3)]
    for d in (F(4), F(4), F(8, 3)):
        g = [[F(1), -d/2, F(0)], [d/2, F(1), F(0)], [F(0), F(0), F(1)]]
        cancellation = matmul(g, cancellation)
    cancellation = add(cancellation, identity(3), F(1, 25))
    assert cancellation == [[F(0), F(0), F(0)], [F(0), F(0), F(0)],
                             [F(0), F(0), F(28, 625)]]
    return {
        "qualification": "OU3_CHRONOLOGICAL_SIX_PIVOT_IMPLICATION_V1",
        "quiet_nominal_two_group_singular_floor": str(quiet_s),
        "quiet_nominal_all_six_pivot_floor": str(quiet_p),
        "quiet_nominal_scope": "zero-residual regular nominal source history, g>=9 and B>=9, actually applied acc/mag groups eight qualified predictions apart, raw AG scaling",
        "quiet_nominal_floor_identifies_physical_attitude_and_bias": False,
        "conditional_chronological_gyro_floor": True,
        "literal_reset_inverse_norm_at_most_one": True,
        "inverse_frame_reset_cancellation": True,
        "inverse_frame_every_prefix_gyro_bound": True,
        "terminal_reset_supplied_example_gyro_floor": str(inverse[-1]['gyro_singular_lower']),
        "inverse_frame_uniform_injection_budget_proved": False,
        "relaxed_reset_cancellation_gyro_block": [[str(x) for x in row] for row in cancellation],
        "reset_cancellation_shipping_reachability_proved": False,
        "reset_cancellation_full_historical_rank_loss_claimed": False,
        "conditional_two_group_six_column_floor": True,
        "same_prediction_cell_actual_row_factorization": True,
        "same_prediction_cell_row_defect_zero": True,
        "same_cell_reset_pulled_field_geometry_implication": True,
        "physical_span_supplies_reset_pulled_field_geometry": False,
        "singular_floor_to_all_six_greedy_pivots": True,
        "supplied_example_singular_floor": str(s),
        "supplied_example_all_six_pivot_floor": str(p),
        "supplied_example_exact_gram_residual_pivots": [str(x) for x in pivots],
        "supplied_reset_sequence_gyro_floor": str(transport["gyro_singular_lower"]),
        "supplied_failed_defect_margin": str(failed["signed_margin"]),
        "actual_chronological_resets_may_be_dropped": False,
        "physical_span_supplies_uniform_row_and_defect_budget": False,
        "two_temporal_margins_imply_six_pivots": False,
        "uniform_historical_AG_readout_action": False,
        "rho0_certified": False,
    }
