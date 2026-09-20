"""Optional native carried-history probe; no uniform/float32 proof promotion.

python -m tools.stability.ou3_theorem.ag_readout_source_diagnostic --eigen DIR
Needs g++, Eigen and mpmath. Observer taps are inserted in a temporary header;
the repository shipping sources are never edited. An untapped control must
produce the identical terminal mean, quaternion, covariance and stage times.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path
import subprocess
import tempfile

REPO = Path(__file__).resolve().parents[3]
HEADER = Path('src/kalman_ou_iii/Kalman3D_Wave_OU_III.h')


def instrument(source):
    def once(old, new):
        nonlocal source
        if source.count(old) != 1:
            raise ValueError(f'shipping observer anchor changed: {old[:60]}')
        source = source.replace(old, new)

    once('    apply_pending_aw_covariance_inflation_();', '''    if (recording) {
        const T trace_phi = std::exp(-Ts / std::max(T(1e-3), tau_bacc_));
        const T trace_qscale = -T(0.5)*std::max(T(1e-3),tau_bacc_)
            *std::expm1(-T(2)*Ts/std::max(T(1e-3),tau_bacc_));
        readout_prediction(F_AA, F_LL, Q_AA, Q_LL, trace_phi, (Q_bacc_*trace_qscale).eval());
    }
    apply_pending_aw_covariance_inflation_();''')
    once('    apply_pending_aw_covariance_inflation_();\n    symmetrize_Pext_();   // Symmetry hygiene', '''    Eigen::Matrix<T,NX,NX> trace_before_sync;
    if (recording) trace_before_sync = Pext;
    apply_pending_aw_covariance_inflation_();
    symmetrize_Pext_();   // Symmetry hygiene
    if (recording) events.push_back(std::string("{\\"kind\\":\\"sync_completion\\",\\"before\\":")
        +matrix_json(trace_before_sync)+",\\"after\\":"+matrix_json(Pext)+'}');''')
    once('        Pext.template block<3,3>(OFF_AW, OFF_AW) += Delta;',
         '        readout_sync(Delta);\n        Pext.template block<3,3>(OFF_AW, OFF_AW) += Delta;')
    once('    ocean_imu::kalman::ou_detail::apply_left_error_reset<T, NX>(Pext, dtheta_injected);',
         '    readout_reset(dtheta_injected);\n'
         '    ocean_imu::kalman::ou_detail::apply_left_error_reset<T, NX>(Pext, dtheta_injected);')
    # Record only actually applied updates. Preserve any innovation safety bump
    # as effective R, without confusing it with floating covariance roundoff.
    for start, stop, sensor, hcode, noise in (
        ('::measurement_update_acc_only(', '::measurement_update_mag_only(', 'acc',
         'trace_h.template block<3,3>(0,0)=J_att; '
         'trace_h.template block<3,3>(0,3)=J_bg; '
         'trace_h.template block<3,3>(0,15)=J_aw; '
         'trace_h.template block<3,3>(0,18).setIdentity();', 'Racc'),
        ('::measurement_update_mag_only(', '::accelerometer_measurement_func(', 'mag',
         'trace_h.template block<3,3>(0,0)=J_att;', 'Rmag'),
        ('::applyIntegralZeroPseudoMeas()', '::measurement_update_position_pseudo(', 'S',
         'trace_h.template block<3,3>(0,12).setIdentity();', 'R_S'),
    ):
        a, b = source.index(start), source.index(stop)
        chunk = source[a:b]
        anchor = '    Eigen::LDLT<Matrix3> ldlt;'
        if chunk.count(anchor) != 1 or chunk.count('    xext.noalias() += K * r;') != 1:
            raise ValueError('shipping correction observer anchor changed')
        chunk = chunk.replace(anchor, '    const Matrix3 trace_s = S_mat;\n'+anchor)
        chunk = chunk.replace('    xext.noalias() += K * r;', f'''    if (recording) {{
        Eigen::Matrix<T,3,NX> trace_h = Eigen::Matrix<T,3,NX>::Zero();
        {hcode}
        readout_correction("{sensor}",trace_h,({noise}+(S_mat-trace_s)).eval(),K);
    }}
    xext.noalias() += K * r;''')
        source = source[:a]+chunk+source[b:]
    return source


def analyze(trace, dps=80):
    """Six-column QR-style reader and full backward action on exported operands.

    No source reachability enclosure follows from freezing float coefficients.
    Q blocks are retained as exported; tiny negative arithmetic eigenvalues
    are reported, never clipped or promoted to positive process factors.
    """
    import mpmath as mp
    if dps < 80:
        raise ValueError('at least 80 decimal digits required')
    with mp.workdps(dps):
        mat = lambda a: mp.matrix([[mp.mpf(float(x)) for x in row] for row in a])
        rows, sequence = [], []
        t = mp.eye(21)[:, :6]
        qmin, injection = mp.inf, mp.mpf(0)
        sync_skew = Fraction(0)
        sync_rayleigh = Fraction(0)
        for e in trace['events']:
            kind = e['kind']
            if kind == 'prediction':
                f, q = mp.eye(21), mp.zeros(21)
                f[:6, :6], f[6:18, 6:18] = mat(e['F_AG']), mat(e['F_LIN'])
                f[18:, 18:] = mp.mpf(e['phi_BA'])*mp.eye(3)
                q[:6, :6], q[6:18, 6:18], q[18:, 18:] = mat(e['Q_AG']), mat(e['Q_LIN']), mat(e['Q_BA'])
                qmin = min(qmin, min(mp.eigsy((q+q.T)/2, eigvals_only=True)))
            elif kind == 'sync':
                sync_skew = max(sync_skew, max(abs(Fraction(e['Q'][i][j])-Fraction(e['Q'][j][i]))
                                               for i in range(3) for j in range(3)))
                continue
            elif kind == 'sync_completion':
                exact_increment = completed_sync_increment(e)
                sync_rayleigh = min(sync_rayleigh, min(
                    (exact_increment[i][i]+exact_increment[j][j])/2-abs(exact_increment[i][j])
                    for i in range(21) for j in range(i+1, 21)))
                f, q = mp.eye(21), mat(exact_increment)
            elif kind == 'reset':
                d = mat(e['d'])
                x, y, z = d
                f, q = mp.eye(21), None
                f[:3, :3] += mp.matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]])/2
                injection = max(injection, mp.norm(d))
            else:
                h, r = mat(e['H']), mat(e['R'])
                rows.extend([list((h*t)[i, :]) for i in range(3)])
                sequence.append(('correction', h, r))
                continue
            t = f*t
            sequence.append((kind, f, q))
        o, endpoint = mp.matrix(rows), t[:6, :]
        # High-precision factor pivoting. This is diagnostic rank, never an
        # exact-zero test or a source-uniform singular-value certificate.
        residual = o.copy()
        selected = []
        for _ in range(6):
            norms = [mp.norm(residual[i, :])**2 for i in range(residual.rows)]
            index = max(range(len(norms)), key=norms.__getitem__)
            selected.append(index)
            v = residual[index, :].copy()
            for i in range(residual.rows):
                residual[i, :] -= (residual[i, :]*v.T)[0]/norms[index]*v
        reduced = endpoint*mp.matrix([rows[i] for i in selected])**-1
        reader = mp.zeros(6, len(rows))
        for j, index in enumerate(selected):
            reader[:, index] = reduced[:, j]
        y, action, offset = mp.eye(21)[:6, :], mp.zeros(6), len(rows)
        for kind, b, q in reversed(sequence):
            if kind == 'correction':
                offset -= 3
                li = reader[:, offset:offset+3]
                action += li*q*li.T
                y -= li*b
            else:
                if q is not None:
                    action += y*q*y.T
                y = y*b
        from .nuisance_upper_certificate import bounds
        u = mp.diag([mp.mpf(x.numerator)/x.denominator for x in bounds()[-1] for _ in range(3)])
        action += y[:, 6:]*u*y[:, 6:].T
        p = mat(trace['terminal_covariance'])[:6, :6]
        fmt = lambda x: mp.nstr(x, 24)
        return {'factor_selected_rows': selected,
                'AG_root_residual_diagnostic': fmt(mp.norm(y[:, :6])),
                'raw_AG_array_sigma_min': fmt(min(mp.svd(o, compute_uv=False))),
                'AG_action_lambda_max': fmt(max(mp.eigsy((action+action.T)/2, eigvals_only=True))),
                'AG_action_minus_literal_covariance_lambda_min': fmt(min(mp.eigsy((action+action.T-p-p.T)/2, eigvals_only=True))),
                'exported_Q_lambda_min': fmt(qmin), 'maximum_injection_rad': fmt(injection),
                'exported_sync_asymmetry_exact': str(sync_skew),
                'minimum_complete_sync_pair_Rayleigh_quotient_exact': str(sync_rayleigh),
                'applied_rows': len(rows), 'events': len(sequence),
                'source_uniform_verified': False, 'rigorous_enclosure': False}


def rational_upper_factor(matrix, bits=192):
    """Exact LDL factor with upward rational square roots of its pivots.

    Retain every off-diagonal correlation: if Q=L D L', this returns
    U=L sqrt(D_upper), hence U U'-Q is PSD by congruence. This is a supplied
    rational matrix enclosure, not an interval covariance recursion.
    """
    from .matrix_certificates import identity, transpose
    a = [[Fraction(x) for x in row] for row in matrix]
    if not a or any(len(row) != len(a) for row in a) or a != transpose(a):
        raise ValueError('symmetric square matrix required')
    l, d = identity(len(a)), []
    for j in range(len(a)):
        pivot = a[j][j]-sum(l[j][k]**2*d[k] for k in range(j))
        column = [a[i][j]-sum(l[i][k]*l[j][k]*d[k] for k in range(j))
                  for i in range(j+1, len(a))]
        if pivot < 0 or (not pivot and any(column)):
            raise ValueError('matrix is not positive semidefinite')
        d.append(pivot)
        for i, value in enumerate(column, j+1):
            l[i][j] = value/pivot if pivot else Fraction(0)
    scale = 1 << bits
    roots = []
    for pivot in d:
        value = pivot*scale*scale
        ceiling = -(-value.numerator//value.denominator)
        root = isqrt(ceiling)
        if root*root < ceiling:
            root += 1
        roots.append(Fraction(root, scale))
    return [[v*roots[j] for j, v in enumerate(row)] for row in l]


def completed_sync_increment(event):
    """Exact change from symmetric pre-state through literal add/symmetry.

    The pre-state may itself be asymmetric from prediction arithmetic.
    Its antisymmetric part is recorded, not interpreted as a covariance.
    No assertion about the prediction's rounding error follows here.
    """
    before = [[Fraction(x) for x in row] for row in event['before']]
    after = [[Fraction(x) for x in row] for row in event['after']]
    from .matrix_certificates import transpose
    if len(before) != 21 or any(len(row) != 21 for row in before):
        raise ValueError('full 21-state operation boundary required')
    if len(after) != 21 or any(len(row) != 21 for row in after) or after != transpose(after):
        raise ValueError('literal completed covariance must be symmetric')
    return [[after[i][j]-(before[i][j]+before[j][i])/2
             for j in range(21)] for i in range(21)]


def signed_upper_factor(matrix):
    """Congruence upper factor of an arbitrary exact symmetric defect.

    Signed rank-one/two Schur elimination retains all off-diagonal terms.
    Keep positive terms, drop negative terms, and round only positive square
    roots upward. This is an arithmetic-defect envelope, never a scalar
    contraction reduction or a repair of the shipping covariance.
    """
    from .matrix_certificates import add, is_psd, matmul, transpose
    original = [[Fraction(x) for x in row] for row in matrix]
    n = len(original)
    if not n or any(len(row) != n for row in original) or original != transpose(original):
        raise ValueError('symmetric square defect required')
    a = [row[:] for row in original]
    positive = []
    while any(any(row) for row in a):
        pivots = [i for i in range(n) if a[i][i]]
        if pivots:
            j = max(pivots, key=lambda i: abs(a[i][i]))
            terms = [(a[j][j], [a[i][j]/a[j][j] for i in range(n)])]
        else:
            i, j = next((i, j) for i in range(n) for j in range(i+1, n) if a[i][j])
            b = a[i][j]
            sign = 1 if b > 0 else -1
            terms = [(Fraction(1, 2)/abs(b), [a[k][i]+sign*a[k][j] for k in range(n)]),
                     (-Fraction(1, 2)/abs(b), [a[k][i]-sign*a[k][j] for k in range(n)])]
        for weight, column in terms:
            if weight > 0:
                root = rational_upper_factor([[weight]])[0][0]
                positive.append([root*x for x in column])
            a = [[a[i][j]-weight*column[i]*column[j] for j in range(n)] for i in range(n)]
    factor = transpose(positive) if positive else [[Fraction(0)] for _ in range(n)]
    if not is_psd(add(matmul(factor, transpose(factor)), original, -1)):
        raise ArithmeticError('signed factor did not enclose literal operation')
    return factor


def enclose_exported_word(trace):
    """Exact frozen-coefficient audit after the high-precision feasibility run.

    Every literal reset is retained. Q/R factors retain the correlated LIN
    block; complete literal sync/symmetry defects have signed upper factors.
    Float32 operation errors and coverage of real histories remain open.
    """
    from .ag_readout import exact_readout, readout_action
    from .matrix_certificates import add, encoded, identity, is_psd
    from .nuisance_upper_certificate import bounds
    f = Fraction
    events = []
    completed = 0
    pending_prediction = False
    for e in trace['events']:
        kind = e['kind']
        if kind == 'prediction':
            if pending_prediction:
                raise ValueError('missing sync/symmetry boundary before next prediction')
            pending_prediction = True
            transition, factor = identity(21), [[f(0)]*21 for _ in range(21)]
            for start, name in ((0, 'AG'), (6, 'LIN')):
                block = e['F_'+name]
                for i, row in enumerate(block):
                    transition[start+i][start:start+len(row)] = list(map(f, row))
            for i in range(18, 21):
                transition[i][i] = f(e['phi_BA'])
            for start, name in ((0, 'AG'), (6, 'LIN'), (18, 'BA')):
                block = rational_upper_factor(e['Q_'+name])
                for i, row in enumerate(block):
                    factor[start+i][start:start+len(row)] = row
            events.append({'kind': 'prediction', 'F': transition, 'U': factor})
        elif kind == 'correction':
            if pending_prediction:
                raise ValueError('missing sync/symmetry boundary before correction')
            events.append({'kind': kind, 'H': e['H'], 'V': rational_upper_factor(e['R'])})
        elif kind == 'reset':
            if pending_prediction:
                raise ValueError('missing sync/symmetry boundary before reset')
            x, y, z = [f(row[0]) for row in e['d']]
            reset = identity(21)
            cross = [[0, -z, y], [z, 0, -x], [-y, x, 0]]
            for i in range(3):
                for j in range(3):
                    reset[i][j] += cross[i][j]/2
            events.append({'kind': kind, 'G': reset})
        elif kind == 'sync':
            if not pending_prediction:
                raise ValueError('sync operand outside a prediction boundary')
            # Keep the raw operand in provenance, but charge the complete
            # literal operation exactly once at its following boundary.
            continue
        elif kind == 'sync_completion':
            if not pending_prediction:
                raise ValueError('sync completion without a prediction')
            pending_prediction = False
            factor = signed_upper_factor(completed_sync_increment(e))
            events.append({'kind': 'prediction', 'F': identity(21), 'U': factor})
            completed += 1
        else:
            raise ValueError('unrecorded source operation')
    if pending_prediction or completed != sum(e['kind'] == 'prediction' for e in trace['events']):
        raise ValueError('every prediction needs its complete sync/symmetry boundary')
    upper = [[f(0)]*15 for _ in range(15)]
    for i, value in enumerate(x for x in bounds()[-1] for _ in range(3)):
        upper[i][i] = value
    root_nuisance = [[f(x) for x in row[6:]] for row in trace['root_covariance'][6:]]
    if not is_psd(add(upper, root_nuisance, f(-1))):
        raise ArithmeticError('literal root nuisance covariance exceeds supplied comparison')
    result = readout_action(events, exact_readout(events), upper)
    action = result['action']
    # A readable full matrix certificate. Check its matrix slack exactly;
    # rounded entries are never treated as certified eigenvalue bounds.
    ceiling = [[f(round(1000*x), 1000) for x in row] for row in action]
    for i in range(6):
        ceiling[i][i] += f(1, 100)
    if not is_psd(add(ceiling, action, f(-1))):
        raise ArithmeticError('rounded full action ceiling failed exact PSD check')
    literal = [[f(x) for x in row[:6]] for row in trace['terminal_covariance'][:6]]
    if not is_psd(add(ceiling, literal, f(-1))):
        raise ArithmeticError('ceiling does not bound this literal terminal AG covariance')
    return {'qualification': 'OU3_EXPORTED_WORD_ACTION_ENCLOSURE_V1',
            'exported_coefficients_sha256': hashlib.sha256(json.dumps(trace, sort_keys=True).encode()).hexdigest(),
            'AG_root_cancelled_exactly': True, 'correlated_Q_R_enclosed_by_rational_factors': True,
            'literal_root_nuisance_comparison_verified': True,
            'complete_sync_operation_enclosed': True, 'sync_boundaries_verified': completed,
            'all_other_float_operations_enclosed': False,
            'full_matrix_action_ceiling_verified': True, 'action_ceiling': encoded(ceiling),
            'literal_terminal_AG_covariance_below_ceiling': True,
            'scope': 'one exported float coefficient sequence interpreted as exact rationals',
            'real_arithmetic_shipping_trajectory_enclosed': False,
            'source_uniform_verified': False, 'theorem_closed': False}


def run(eigen, headings=('0', '0.001', '0.000001', 'wave')):
    source = (REPO/HEADER).read_text()
    cases = []
    with tempfile.TemporaryDirectory(prefix='ou3-ag-source-') as directory:
        tmp = Path(directory)
        include = tmp/'kalman_ou_iii'
        include.mkdir()
        (include/HEADER.name).write_text(instrument(source))
        driver = REPO/'tools/stability/ag_readout_source.cpp'
        binaries = []
        for name, inc in [('observed', ['-I'+str(tmp)]), ('control', [])]:
            binary = tmp/name
            subprocess.run(['g++', '-O2', '-std=c++20', *inc, '-I'+str(REPO/'src'),
                            '-isystem', str(eigen), str(driver), '-o', str(binary)], check=True)
            binaries.append(binary)
        for heading in headings:
            observed, control = [json.loads(subprocess.check_output([str(p), heading], text=True)) for p in binaries]
            for key in control:
                if key != 'events' and control[key] != observed[key]:
                    raise ValueError(f'observer changed literal source output: {key}')
            cases.append({'input_profile': heading,
                          'live_step': observed['live_step'], 'refined_step': observed['refined_step'],
                          'active_step': observed['active_step'],
                          'literal_terminal_parity': True, **analyze(observed)})
            if heading in ('0', 'wave'):
                cases[-1]['exact_exported_word_enclosure'] = enclose_exported_word(observed)
    return {'qualification': 'OU3_CARRIED_SOURCE_READOUT_DIAGNOSTIC_V1',
            'shipping_header_sha256': hashlib.sha256(source.encode()).hexdigest(),
            'decimal_digits': 80, 'profile': 'construction through wrapper release; 225 to 225.32 s',
            'wave_truth': 'p_z=.4 sin(.6t), roll=.02 sin(.5t), B_world=(60,0,30); zero physical biases',
            'cases': cases, 'all_time_magnetic_service_certified': False,
            'source_uniform_verified': False, 'theorem_closed': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--eigen', type=Path, default=Path('/usr/include/eigen3'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(args.eigen), indent=2, sort_keys=True)+'\n')
