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
    t = identity(n)
    rows = []
    for kind, b, _ in _events(events, n):
        if kind == 'correction':
            rows.extend(row[:ag] for row in matmul(b, t))
        else:
            t = matmul(b, t)
    return rows, [row[:ag] for row in t[:ag]]


def exact_readout(events, n=21, ag=6):
    """Select independent rows and solve only the ag-by-ag system exactly.

    This is a witness constructor for a supplied word, not a uniform rank
    decision. No numerical rank threshold or normal-equation inversion.
    A uniform proof can use a coefficient-dependent reader with certified
    action; finite-word success is insufficient.
    """
    rows, endpoint = observation_array(events, n, ag)
    if not rows:
        raise ValueError('no applied observations')
    selected = _pivots(transpose(rows))
    if len(selected) != ag:
        raise ValueError('raw AG observation array lacks full column rank')
    reader = _zero(ag, len(rows))
    reduced = matmul(endpoint, inverse([rows[i] for i in selected]))
    for i in range(ag):
        for j, index in enumerate(selected):
            reader[i][index] = reduced[i][j]
    return reader


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
