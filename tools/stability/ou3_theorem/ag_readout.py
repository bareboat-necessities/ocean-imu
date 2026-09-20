"""Historical AG readout action: exact conditional algebra, not a source certificate.

Enters the existing tail inequality by bounding carried AG covariance, then
certifying J from a matrix process comparison. Measurements have <=3 rows;
only a six-column raw observation array is solved. Backward action retains
all nuisance columns and every correlated process factor. No Riccati boxes,
independent-block information lifting, or gain substitutions are used.
"""
from fractions import Fraction as F

from .factor_certificates import _pivots
from .lin_path_certificate import inverse
from .matrix_certificates import add, encoded, identity, is_psd, ldlt, matmul, transpose


def _matrix(value, rows=None, cols=None):
    a = [[F(x) for x in row] for row in value]
    if (not a or not a[0] or any(len(row) != len(a[0]) for row in a)
            or (rows is not None and len(a) != rows)
            or (cols is not None and len(a[0]) != cols)):
        raise ValueError('nonempty matrix of the required shape needed')
    return a


def _zero(rows, cols):
    return [[F(0)]*cols for _ in range(rows)]


def _events(events, n):
    result = []
    for event in events:
        kind = event['kind']
        if kind == 'prediction':
            result.append((kind, _matrix(event['F'], n, n), _matrix(event['U'], n)))
        elif kind == 'reset':
            g = _matrix(event['G'], n, n)
            if len(_pivots(g)) != n:
                raise ValueError('nonsingular reset required')
            result.append((kind, g, None))
        elif kind == 'correction':
            if not event.get('accepted', True):
                continue
            h = _matrix(event['H'], cols=n)
            if not 1 <= len(h) <= 3:
                raise ValueError('one to three actual measurement rows required')
            v = _matrix(event['V'], len(h))
            ldlt(matmul(v, transpose(v)))
            result.append((kind, h, v))
        else:
            raise ValueError('unsupported operation; cannot omit a hard event')
    return result


def observation_array(events, n=21, ag=6):
    """Raw AG sensitivity including predictions/resets, before innovation filtering.

    Omitting the optimal correction from this auxiliary observation model
    does not omit it from the actual covariance. The readout bound follows
    from optimal linear estimation on exactly this model, with coefficients
    frozen from the actual execution. It is not a second nonlinear filter.
    """
    if not 0 < ag < n:
        raise ValueError('both AG and nuisance coordinates required')
    t = [row[:ag] for row in identity(n)]
    rows = []
    for kind, b, _ in _events(events, n):
        if kind == 'correction':
            rows.extend(matmul(b, t))
        else:
            t = matmul(b, t)
    return rows, t[:ag]


def factor_rows(rows):
    """Exact row-pivoted orthogonalization; no normal equations or tiny pivots.

    Select the largest remaining squared residual at each step. The old
    first-independent-row rule can invert a vanishing component even when
    the complete array remains well conditioned. This rule considers every
    available row. Its action still requires a separate uniform certificate.
    """
    residuals = _matrix(rows)
    selected = []
    for _ in range(len(residuals[0])):
        norms = [sum(x*x for x in row) for row in residuals]
        pivot = max(range(len(norms)), key=norms.__getitem__)
        if not norms[pivot]:
            break
        selected.append(pivot)
        v = residuals[pivot][:]
        for i, row in enumerate(residuals):
            scale = sum(x*y for x, y in zip(row, v))/norms[pivot]
            residuals[i] = [x-scale*y for x, y in zip(row, v)]
    return selected


def exact_readout(events, n=21, ag=6):
    """Select rows by factor pivoting and solve only the ag-by-ag system exactly.

    This is a witness constructor for a supplied word, not a uniform rank
    decision. No numerical rank threshold or normal-equation inversion.
    A uniform proof can use a coefficient-dependent reader with certified
    action; finite-word success is insufficient.
    """
    rows, endpoint = observation_array(events, n, ag)
    if not rows:
        raise ValueError('no applied observations')
    selected = factor_rows(rows)
    if len(selected) != ag:
        raise ValueError('raw AG observation array lacks full column rank')
    reader = _zero(ag, len(rows))
    reduced = matmul(endpoint, inverse([rows[i] for i in selected]))
    for i in range(ag):
        for j, index in enumerate(selected):
            reader[i][index] = reduced[i][j]
    return reader


def noise_action_lower(events, n=21, ag=6):
    """All-row necessary matrix test, independent of a selected minor.

    Every exact reader's action is >= T (sum O_i' R_i^-1 O_i)^-1 T'.
    Accumulate rank-at-most-three terms and solve only six coordinates. This
    is a lower bound on readout noise, NOT a covariance/information lifting.
    It neither upper-bounds the nuisance/process action nor certifies histories.
    """
    rows, terminal = observation_array(events, n, ag)
    gram, offset = _zero(ag, ag), 0
    for kind, h, v in _events(events, n):
        if kind != 'correction':
            continue
        block = rows[offset:offset+len(h)]
        offset += len(h)
        ri = inverse(matmul(v, transpose(v)))
        gram = add(gram, matmul(matmul(transpose(block), ri), block))
    ldlt(gram)
    return matmul(matmul(terminal, inverse(gram)), transpose(terminal))


def readout_action(events, reader, nuisance_upper, ag=6):
    """Verify exact AG root cancellation and return its full matrix action.

    If the actual historical root has P_nn <= nuisance_upper, its terminal
    AG covariance is <= action. No bound on root P_hh or P_hn is used.
    Q=U U', R_eff=V V' must be the actual factors (including safety/sync).
    Source coverage, factor enclosures and float32 are separate obligations.
    """
    upper = _matrix(nuisance_upper)
    ldlt(upper)
    n = ag+len(upper)
    sequence = _events(events, n)
    columns = sum(len(b) for kind, b, _ in sequence if kind == 'correction')
    reader = _matrix(reader, ag, columns)
    y = [row[:] for row in identity(n)[:ag]]
    action = _zero(ag, ag)
    offset = columns
    for kind, b, factor in reversed(sequence):
        if kind == 'correction':
            offset -= len(b)
            li = [row[offset:offset+len(b)] for row in reader]
            z = matmul(li, factor)
            action = add(action, matmul(z, transpose(z)))
            y = add(y, matmul(li, b), F(-1))
        else:
            if kind == 'prediction':
                z = matmul(y, factor)
                action = add(action, matmul(z, transpose(z)))
            y = matmul(y, b)
    if any(any(row[:ag]) for row in y):
        raise ValueError('uncancelled AG root: no AG prior ceiling is available')
    residual = [row[ag:] for row in y]
    action = add(action, matmul(matmul(residual, upper), transpose(residual)))
    if not is_psd(action):
        raise ArithmeticError('negative readout action')
    return {'action': action, 'nuisance_root_residual': residual,
            'AG_root_cancelled_exactly': True, 'source_uniform_verified': False}


def bootstrap(ag_upper, nuisance_upper, transition, process_factor, epsilon, eta=1):
    """Verify the matrix first-prediction implication, conditional on upper bounds.

    C_eta bounds the full covariance by block Cauchy--Schwarz. Check
    Q >= epsilon F C_eta F' with exact PSD elimination; retain both matrices.
    Then D >= delta P^-1 and J=delta (C_eta^-1)_hh. The resulting J can
    enter the retained nuisance/Schur implication. This check alone does NOT
    certify that any supplied upper comparison holds on shipping histories.
    """
    b, u = _matrix(ag_upper), _matrix(nuisance_upper)
    ldlt(b)
    ldlt(u)
    eps, eta = F(epsilon), F(eta)
    if min(eps, eta) <= 0:
        raise ValueError('strictly positive comparison constants required')
    ag, nu = len(b), len(u)
    c = [[(1+eta)*v for v in row]+[F(0)]*nu for row in b]
    c += [[F(0)]*ag+[(1+1/eta)*v for v in row] for row in u]
    f = _matrix(transition, ag+nu, ag+nu)
    qf = _matrix(process_factor, ag+nu)
    q = matmul(qf, transpose(qf))
    fcft = matmul(matmul(f, c), transpose(f))
    if not is_psd(add(q, fcft, -eps)):
        raise ValueError('Q - epsilon F C_eta F^T is not PSD')
    delta = eps/(1+eps)
    j = [[delta*v/(1+eta) for v in row] for row in inverse(b)]
    return {'full_upper': c, 'AG_loss_lower': j, 'delta': delta,
            'rho_upper': 1-delta, 'source_uniform_verified': False}


def supplied_fixture():
    """Varying, cross-coupled 21-state algebra audit, NOT a shipping history.

    Observation sparsity and literal reset form are retained. The LIN and
    process values are deliberately rational test inputs, not shipping OU
    discretization or an admitted physical trajectory.
    """
    def skew(v):
        x, y, z = map(F, v)
        return [[F(0), -z, y], [z, F(0), -x], [-y, x, F(0)]]

    def observation(v, acc=False):
        h = _zero(3, 21)
        for i, row in enumerate(skew(v)):
            h[i][:3] = [-x for x in row]
            if acc:
                h[i][15+i] = h[i][18+i] = F(1)
        return h

    events = []
    ha = observation(['0', '0', '-9.80665'], True)
    hm = observation(['35', '0', '20'])
    for step, d in [('.004', ['.01', '-.02', '.03']), ('.006', ['-.02', '.01', '.02'])]:
        f = identity(21)
        u = [[x/1000 for x in row] for row in identity(21)]
        for j in range(3):
            f[j][3+j] = F(step)
            f[6+j][15+j] = f[9+j][6+j] = f[12+j][9+j] = F(step)
            f[15+j][15+j], f[18+j][18+j] = F('.9'), F('.999')
            u[6+j][15+j] = F('.0002')
        events += [{'kind': 'prediction', 'F': f, 'U': u},
                   {'kind': 'correction', 'H': ha, 'V': [[x*F('.12') for x in row] for row in identity(3)]},
                   {'kind': 'correction', 'H': hm, 'V': [[x*F('.8') for x in row] for row in identity(3)]}]
        g = identity(21)
        for i, row in enumerate(skew(d)):
            for j, value in enumerate(row):
                g[i][j] += value/2
        events.append({'kind': 'reset', 'G': g})
    return events


def coefficient_relaxation_obstruction():
    """Source-shaped singular coefficient family, NOT an admitted history.

    Even bounding nominal AW by the physical acceleration ceiling would not
    repair this coefficient-only relaxation. The missing premise is its
    compatibility with the carried nominal mean/innovation/reset recursion.
    """
    g = F('9.80665')
    field = [F(45), F(0), F(45)]
    aw = [-g/2, F(0), g/2]
    force = [-g/2, F(0), -g/2]
    def observation(v, acc=False):
        x, y, z = v
        h = _zero(3, 21)
        h[0][:3], h[1][:3], h[2][:3] = [0, z, -y], [-z, 0, x], [y, -x, 0]
        if acc:
            for i in range(3):
                h[i][15+i] = h[i][18+i] = F(1)
        return h
    transition = identity(21)
    for i in range(3):
        transition[i][3+i] = F(1, 25)
    # Only the AG observation array is at issue: nuisance prediction/process
    # coefficients do not enter these root columns for the regular source
    # block structure. Identity on nuisance is not a claimed shipping step.
    events = []
    for _ in range(2):
        events.extend([
            {'kind': 'prediction', 'F': transition, 'U': [[0] for _ in range(21)]},
            {'kind': 'correction', 'H': observation(force, True), 'V': identity(3)},
            {'kind': 'correction', 'H': observation(field), 'V': identity(3)},
            {'kind': 'reset', 'G': identity(21)},
        ])
    rows, endpoint = observation_array(events)
    null = [[F(i in (0, 2)), F(i in (3, 5))] for i in range(6)]
    if any(any(row) for row in matmul(rows, null)):
        raise ArithmeticError('claimed annihilator is not exact')
    if len(factor_rows(rows)) != 4 or not any(any(row) for row in matmul(endpoint, null)):
        raise ArithmeticError('coefficient relaxation obstruction failed')
    return {'qualification': 'OU3_NOMINAL_COEFFICIENT_RELAXATION_OBSTRUCTION_V1',
            'nominal_aw': list(map(str, aw)), 'nominal_force': list(map(str, force)),
            'magnetic_field': list(map(str, field)), 'nominal_aw_norm_squared': str(g*g/2),
            'full_AG_array_rank': 4, 'AG_null_columns': encoded(null),
            'candidate_Gram_floor': 'I6', 'null_Rayleigh_margin': '-1',
            'LO_equals_terminal_AG_map_possible': False,
            'nominal_mean_recursion_satisfied': False,
            'same_history_magnetic_service_certified': False,
            'shipping_counterexample': False,
            'invalidated_method': 'independent coefficient ranges without nominal-history linkage'}


def gyro_alias_obstruction():
    """Mean-recursion-compatible regular-root relaxation, not a carried root.

    Quiet truth, zero residuals and the nominal bias -2*pi/h e_z give one
    complete nominal turn per sample. Literal real-arithmetic Rodrigues and
    B helpers give R=I and B=h e_z e_z'. Quaternion sign changes leave the
    observations unchanged. Nonzero process noise and arbitrary realized
    gains cannot reveal an exactly annihilated historical root column.
    Construction reachability and all-time magnetic service are NOT proved.
    """
    h = F(1, 200)
    transition = identity(21)
    transition[2][5] = h
    ha, hm = _zero(3, 21), _zero(3, 21)
    g, field = F('9.80665'), F(75)
    # -skew((0,0,-g)) and -skew((75,0,0)).
    ha[0][1], ha[1][0] = -g, g
    hm[1][2], hm[2][1] = field, -field
    for i in range(3):
        ha[i][15+i] = ha[i][18+i] = F(1)
    hs = _zero(3, 21)
    for i in range(3):
        hs[i][12+i] = F(1)
    events = []
    for _ in range(2):
        # Noise/nuisance placeholders are not source Q/FLIN certificates:
        # neither enters the six root columns in this block structure.
        events.append({'kind': 'prediction', 'F': transition, 'U': _zero(21, 1)})
        for observation in (hs, ha, hm):
            events.append({'kind': 'correction', 'H': observation, 'V': identity(3)})
            events.append({'kind': 'reset', 'G': identity(21)})
    rows, terminal = observation_array(events)
    null = [[F(i == 3), F(i == 4)] for i in range(6)]
    if any(any(row) for row in matmul(rows, null)) or len(factor_rows(rows)) != 4:
        raise ArithmeticError('complete-turn annihilator failed')
    if matmul(terminal, null) != null:
        raise ArithmeticError('endpoint must preserve both hidden gyro columns')
    return {'qualification': 'OU3_NOMINAL_GYRO_ALIAS_RELAXATION_V1',
            'step_s': str(h), 'nominal_gyro_bias_rad_s': ['0', '0', '-400*pi'],
            'true_gyro_and_bias': 'zero', 'nominal_force': ['0', '0', '-9.80665'],
            'magnetic_field': ['75', '0', '0'],
            'nominal_rotation_per_sample': '2*pi',
            'exact_bias_transport': encoded([[0, 0, 0], [0, 0, 0], [0, 0, h]]),
            'full_AG_array_rank': 4, 'AG_null_columns': encoded(null),
            'null_Rayleigh_margin_against_I6': '-1',
            'LO_equals_terminal_AG_map_possible': False,
            'quiet_physical_motion_and_true_bias_limits_satisfied': True,
            'regular_nominal_mean_recursion_satisfied_real_arithmetic': True,
            'all_nominal_innovations': 'zero',
            'arbitrary_realized_gains_cannot_remove_annihilator': True,
            'shipping_construction_reachability_verified': False,
            'all_time_magnetic_service_verified': False,
            'literal_float32_alias_claimed': False,
            'shipping_stability_refuted': False,
            'invalidated_method': 'innovation bounds alone without construction-linked nominal gyro control'}


def certificate():
    events, upper = supplied_fixture(), identity(15)
    reader = exact_readout(events)
    audit = readout_action(events, reader, upper)
    step = bootstrap(audit['action'], upper, events[0]['F'], events[0]['U'], F(1, 10**10))
    return {'qualification': 'OU3_AG_HISTORICAL_READOUT_ACTION_V1',
            'conditional_algebra_verified': True,
            'scope': 'supplied rational 21-state word; no shipping reachability or uniform enclosure',
            'AG_root_cancelled_exactly': audit['AG_root_cancelled_exactly'],
            'AG_action': encoded(audit['action']),
            'conditional_AG_loss_lower': encoded(step['AG_loss_lower']),
            'conditional_delta': str(step['delta']),
            'rank_three_measurement_rows': True,
            'reader_selection': 'exact largest-residual factor pivot; all observation rows considered',
            'coefficient_relaxation_obstruction': coefficient_relaxation_obstruction(),
            'nominal_gyro_alias_obstruction': gyro_alias_obstruction(),
            'prior_scale_obstruction': {
                'family': 'P0 = diag(t I6, I15), t > 0',
                'necessary_AG_loss_ceiling': '(D_word)_hh <= I6/t',
                'tested_t': '1000000000000',
                'candidate_J': 'I6/1000000',
                'Rayleigh_margin_ceiling': '-999999/1000000000000',
                'shipping_reachability_claimed': False,
                'magnetic_service_membership_claimed': False,
                'shipping_stability_refuted': False},
            'uniform_historical_readout_verified': False,
            'six_column_source_uniform_loss_verified': False,
            'full_21_covariance_upper_verified': False,
            'source_uniform_A21_linear_dissipativity': False,
            'theorem_closed': False}
