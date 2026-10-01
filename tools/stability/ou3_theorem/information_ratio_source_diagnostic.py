"""Carried feasibility of the information-ratio word contraction (non-promoting).

Runs the falsifiable experiment of `docs/ou3-proof-research-state.md` on carried
16- and 64-s words of the unchanged shipping estimator. For each word it
computes the root data information J, the data-plus-terminal information A,
the exact contraction rho_W, and the information-ratio bound of
`docs/ou3-corrected-word-proof.md` section 6 with the ideal C=P_0 and with the
joint-reader full diffuse limit over the preceding 64 s. It also records the
word Riccati diameter kappa_W=lambda_max(J^-1 A)=lambda_max(Pi^-1 P_diff) of
section 7, its slow/fast factorization, and the rank-one kernel bound whose
only covariance input is the scalar nu'P_0 nu along the physical tilt/BA
kernel nu=(theta_hat, -J_att theta_hat).

A temporary copy of the shipping header gets binary read-only taps; an
untapped control compiled from the unchanged header must end in the identical
state. The replay uses the literal exported F, Q, sync increments, applied H
and effective R, and literal resets, with optimal gains in float64. It is a
feasibility replay, not an enclosure; nothing here certifies a source-uniform
floor, magnetic service or float32 arithmetic.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from .ag_readout_source_diagnostic import HEADER, REPO

DRIVER = REPO/'tools/stability/information_ratio_source.cpp'
STEP = .005
SIZES = {'P': 36+144+36+144+1+9, 'Y': 9, 'A': 63+9, 'M': 63+9, 'S': 63+9, 'G': 3,
         'K': 1+441+21+4}
# (profile, first recorded step): 225 s is 45 s after BA release (transient
# BA/kernel covariance); 5000 s is near the stationary BA variance.
CASES = (('quiet', 45001), ('wave', 45001), ('collinear-service', 45001),
         ('sync-locked', 45001), ('quiet', 1000001), ('wave', 1000001))
MOVING = ('wave', 'collinear-service', 'sync-locked')
HISTORY_S, WORDS_S = 64, (16, 64)
KILL_FACTOR = 10
TILT_PREMISES = ('1e-4', '1e-3')           # rad^2 ceilings along theta_hat
KERNEL_PREMISE = '1e-3'
DIFFUSE = 1e8
FINITE_DIAMETER = 1e6
RELATIVE_DIFFERENCE_FLOOR = 1e-4         # replay-noise floor for the diameter identity residual
BA_MARGINAL = 1/1600
KAPPAS = np.unique(np.concatenate([np.geomspace(1.0001, 1e7, 481), [1.5, 2, 3, 4, 5, 8, 10]]))
NAMES = ['theta']*3+['bg']*3+['v']*3+['p']*3+['S']*3+['aw']*3+['ba']*3
SLOW = list(range(6))+list(range(18, 21))
FAST = list(range(6, 18))


def instrument(source):
    """Binary read-only taps at the anchors of the JSON readout observer."""
    def once(old, new):
        nonlocal source
        if source.count(old) != 1:
            raise ValueError(f'shipping observer anchor changed: {old[:60]}')
        source = source.replace(old, new)

    once('    apply_pending_aw_covariance_inflation_();', '''    if (recording) {
        const T trace_phi = std::exp(-Ts / std::max(T(1e-3), tau_bacc_));
        const T trace_qscale = -T(0.5)*std::max(T(1e-3),tau_bacc_)
            *std::expm1(-T(2)*Ts/std::max(T(1e-3),tau_bacc_));
        readout_prediction(F_AA, F_LL, Q_AA, Q_LL, static_cast<double>(trace_phi), (Q_bacc_*trace_qscale).eval());
    }
    apply_pending_aw_covariance_inflation_();''')
    once('        Pext.template block<3,3>(OFF_AW, OFF_AW) += Delta;',
         '        readout_sync(Delta);\n        Pext.template block<3,3>(OFF_AW, OFF_AW) += Delta;')
    once('    ocean_imu::kalman::ou_detail::apply_left_error_reset<T, NX>(Pext, dtheta_injected);',
         '    readout_reset(dtheta_injected);\n'
         '    ocean_imu::kalman::ou_detail::apply_left_error_reset<T, NX>(Pext, dtheta_injected);')
    for start, stop, sensor, hcode, noise in (
        ('::measurement_update_acc_only(', '::measurement_update_mag_only(', 'A',
         'trace_h.template block<3,3>(0,0)=J_att; '
         'trace_h.template block<3,3>(0,3)=J_bg; '
         'trace_h.template block<3,3>(0,15)=J_aw; '
         'trace_h.template block<3,3>(0,18).setIdentity();', 'Racc'),
        ('::measurement_update_mag_only(', '::accelerometer_measurement_func(', 'M',
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
        readout_correction('{sensor}',trace_h,({noise}+(S_mat-trace_s)).eval());
    }}
    xext.noalias() += K * r;''')
        source = source[:a]+chunk+source[b:]
    return source


def compile_driver(tmp, eigen):
    tapped = Path(tmp)/'observed/kalman_ou_iii'
    tapped.mkdir(parents=True)
    (tapped/HEADER.name).write_text(instrument((REPO/HEADER).read_text()))
    out = {}
    for name, extra in (('observed', ['-DINFO_RATIO_TAP', '-I'+str(Path(tmp)/'observed')]),
                        ('control', [])):
        exe = Path(tmp)/(name+'.bin')
        subprocess.run(['g++', '-O2', '-std=c++20', *extra, '-I'+str(REPO/'src'),
                        '-isystem', str(eigen), str(DRIVER), '-o', str(exe)], check=True)
        out[name] = exe
    return out


def parse(path):
    """Decode the tagged binary stream into chronological records."""
    data = Path(path).read_bytes()
    records, i = [], 0
    while i < len(data):
        tag = chr(data[i])
        n = SIZES[tag]
        records.append((tag, np.frombuffer(data, dtype='<f8', count=n, offset=i+1)))
        i += 1+8*n
    return records


def operations(records):
    """Map records to operations of the frozen 21-state word and root markers."""
    ops, markers = [], []
    for tag, v in records:
        if tag == 'P':
            f, q = np.eye(21), np.zeros((21, 21))
            f[:6, :6] = v[:36].reshape(6, 6)
            f[6:18, 6:18] = v[36:180].reshape(12, 12)
            q[:6, :6] = v[180:216].reshape(6, 6)
            q[6:18, 6:18] = v[216:360].reshape(12, 12)
            f[18:, 18:] = v[360]*np.eye(3)
            q[18:, 18:] = v[361:370].reshape(3, 3)
            ops.append(('P', f, sym(q)))
        elif tag == 'Y':
            d = np.zeros((21, 21))
            d[15:18, 15:18] = v.reshape(3, 3)
            ops.append(('Y', sym(d)))
        elif tag in 'AMS':
            ops.append(('C', v[:63].reshape(3, 21), sym(v[63:].reshape(3, 3)), tag))
        elif tag == 'G':
            x, y, z = v
            g = np.eye(21)
            g[:3, :3] += np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])/2
            ops.append(('G', g))
        else:
            markers.append({'op_index': len(ops), 'time': float(v[0]),
                            'P': sym(v[1:442].reshape(21, 21))})
    return ops, markers


def sym(a):
    return (a+a.T)/2


def replay(p0, ops, smoother=False):
    """Frozen-coefficient optimal-gain Riccati replay with fixed-point smoother.

    C=Cov(x_0,x_k|y) evolves by C<-C F', C<-C (I-KH)', C<-C G' and Sigma_00
    loses C H' S^-1 H C' at each correction (proof section 6).
    """
    p = sym(p0)
    c, s00 = (p.copy(), p.copy()) if smoother else (None, None)
    for op in ops:
        kind = op[0]
        if kind == 'P':
            f, q = op[1], op[2]
            p = sym(f@p@f.T+q)
            if smoother:
                c = c@f.T
        elif kind == 'Y':
            p = sym(p+op[1])
        elif kind == 'C':
            h, r = op[1], op[2]
            ph = p@h.T
            s = sym(h@ph+r)
            k = np.linalg.solve(s, ph.T).T
            a = np.eye(21)-k@h
            if smoother:
                ch = c@h.T
                s00 = sym(s00-ch@np.linalg.solve(s, ch.T))
                c = c@a.T
            p = sym(a@p@a.T+k@r@k.T)
        else:
            g = op[1]
            p = sym(g@p@g.T)
            if smoother:
                c = c@g.T
    return p, c, s00


def nuisance_upper():
    from .nuisance_upper_certificate import bounds
    return np.diag([float(x) for x in bounds()[-1] for _ in range(3)])


def joint_reader_full_upper(ops, scales=(1e8, 1e10)):
    """Full 21x21 diffuse limit lim_t Ric_W(diag(t I6, U_n)) over a window."""
    upper = nuisance_upper()
    out = []
    for t in scales:
        prior = np.zeros((21, 21))
        prior[:6, :6] = t*np.eye(6)
        prior[6:, 6:] = upper
        out.append(replay(prior, ops)[0])
    return out[-1], float(np.linalg.norm(out[-1]-out[0])/np.linalg.norm(out[-1]))


def info_ratio_bound(kappa, k):
    """Closed form of sup_y 1/(1+y)-1/((1+k)(1+kappa y)) (proof section 7)."""
    c = 1/(1+max(k, 0.0))
    if kappa <= 1 or c*kappa <= 1:
        return 1-c
    return (math.sqrt(kappa)-math.sqrt(c))**2/(kappa-1)


def gmax(big, small):
    """lambda_max(small^-1 big) through a Cholesky factor of small."""
    li = np.linalg.inv(np.linalg.cholesky(sym(small)))
    return float(np.max(np.linalg.eigvalsh(sym(li@big@li.T))))


def schur(m, keep, drop):
    return m[np.ix_(keep, keep)]-m[np.ix_(keep, drop)]@np.linalg.solve(m[np.ix_(drop, drop)], m[np.ix_(drop, keep)])


def lemma_margin(y, z, chat):
    """Best information-ratio bound over kappa for whitened C-hat>=I."""
    ch = np.linalg.cholesky(sym(chat))
    best = (1.0, None, None)
    for kappa in KAPPAS:
        k = max(float(np.max(np.linalg.eigvalsh(sym(ch.T@(z-kappa*y)@ch)))), 0.0)
        bound = info_ratio_bound(kappa, k)
        if bound < best[0]:
            best = (bound, float(kappa), k)
    return best


def kernel_vector(word, at_end=False):
    """Physical tilt/BA kernel from the first (or last) applied acc/mag rows.

    theta_hat spans the null space of the magnetic attitude Jacobian (the body
    field axis) and nu_ba=-J_att theta_hat cancels its accelerometer effect.
    """
    acc = mag = None
    for op in (reversed(word) if at_end else word):
        if op[0] == 'C' and op[3] == 'A' and acc is None:
            acc = op[1]
        if op[0] == 'C' and op[3] == 'M' and mag is None:
            mag = op[1]
        if acc is not None and mag is not None:
            break
    jm, ja = mag[:, :3], acc[:, :3]
    theta = np.linalg.eigh(jm.T@jm)[1][:, 0]
    nu = np.zeros(21)
    nu[:3] = theta
    nu[18:] = -ja@theta
    return nu


def rank_one_kernel(y, z, p0, nu):
    """Least k with Z-kappa Y <= k u u', u the whitened unit kernel, best kappa.

    k=lambda nu'P_0 nu for the raw-coordinate premise A-kappa J<=lambda nu nu'.
    Returns (bound, kappa, k) with the exact P_0 variance along nu.
    """
    l = np.linalg.cholesky(sym(p0))
    u = l.T@nu
    u /= np.linalg.norm(u)
    basis = np.linalg.qr(np.column_stack([u, np.eye(21)[:, :20]]))[0]
    if abs(abs(basis[:, 0]@u)-1) > 1e-12:
        raise ArithmeticError('kernel basis construction failed')
    pairs = []
    for kappa in KAPPAS:
        x = basis.T@(z-kappa*y)@basis
        rest = x[1:, 1:]
        if np.max(np.linalg.eigvalsh(sym(rest))) >= 0:
            continue
        pairs.append((float(kappa), max(float(x[0, 0]-x[0, 1:]@np.linalg.solve(rest, x[1:, 0])), 0.0)))
    return pairs


def best_bound(pairs, scale=1.0):
    """min over kappa of the bound with k scaled by (ceiling / exact variance)."""
    return min((info_ratio_bound(kappa, k*scale), kappa, k*scale) for kappa, k in pairs)


def shares(vec, scale):
    out = {}
    for i, name in enumerate(NAMES):
        out[name] = out.get(name, 0.0)+(vec[i]/scale[i])**2
    total = sum(out.values())
    return {k: round(v/total, 3) for k, v in out.items() if v/total >= .01}


def gyro_cap(word_s):
    from .word_diameter import gyro_persistence_cap
    return float(gyro_persistence_cap(word_s)['kappa_lower'])


def analyze_split(ops, markers, root_index, word_s):
    root, start, end = markers[root_index], markers[root_index-HISTORY_S], markers[root_index+word_s]
    if abs(end['time']-root['time']-word_s) > 1e-6 or abs(root['time']-start['time']-HISTORY_S) > 1e-6:
        raise ValueError('word markers are not at the requested times')
    hist = ops[start['op_index']:root['op_index']]
    word = ops[root['op_index']:end['op_index']]
    p0 = root['P']
    pend, c, s00 = replay(p0, word, smoother=True)
    forget = sym(s00-c@np.linalg.solve(pend, c.T))
    l = np.linalg.cholesky(sym(p0))
    li = np.linalg.inv(l)
    white = lambda a: sym(li@a@li.T)
    wy, wz = white(s00), white(forget)
    y, z = sym(np.linalg.inv(wy))-np.eye(21), sym(np.linalg.inv(wz))-np.eye(21)
    ev, vec = np.linalg.eigh(sym(wy-wz))
    rho = float(ev[-1])
    ideal = lemma_margin(y, z, np.eye(21))
    cfull, conv = joint_reader_full_upper(hist)
    chat = white(cfull)
    defect = float(np.min(np.linalg.eigvalsh(chat)))-1
    chat_fixed = chat+max(-defect, 0.0)*np.eye(21)
    joint = lemma_margin(y, z, chat_fixed)
    # Diameter: generalized eigenvalue of (Z,Y) and of (P_diff, Pi).
    ly = np.linalg.inv(np.linalg.cholesky(y))
    gev, gvec = np.linalg.eigh(sym(ly@z@ly.T))
    kappa_j = float(gev[-1])
    direction = l@(ly.T@gvec[:, -1])
    pi = replay(1e-14*np.eye(21), word)[0]
    pdiff = replay(DIFFUSE*np.eye(21), word)[0]
    kappa_d = gmax(pdiff, pi)
    # A finite float64 diffuse prior reproduces kappa only when J is well
    # conditioned; along an unobservable (quiet) kernel kappa_W is infinite.
    observable = kappa_j < FINITE_DIAMETER
    kappa_slow = gmax(pdiff[np.ix_(SLOW, SLOW)], pi[np.ix_(SLOW, SLOW)])
    kappa_fast_given_slow = gmax(schur(pdiff, FAST, SLOW), schur(pi, FAST, SLOW))
    # Rank-one kernel bound and its source-style scalar ceilings.
    nu = kernel_vector(word)
    nu_var = float(nu@p0@nu)
    pairs = rank_one_kernel(y, z, p0, nu)
    kern = best_bound(pairs)
    theta = nu[:3]
    tilt_var = float(theta@p0[:3, :3]@theta)
    nu_ba = float(np.linalg.norm(nu[18:]))
    ba_max = float(np.max(np.linalg.eigvalsh(p0[18:, 18:])))

    def ceiling_margin(tau):
        c_nu = (math.sqrt(tau)+nu_ba*math.sqrt(BA_MARGINAL))**2
        return 1-best_bound(pairs, c_nu/nu_var)[0], c_nu
    actual_tilt = ceiling_margin(tilt_var)
    premises = {t: ceiling_margin(float(t)) for t in TILT_PREMISES}
    # Kernel-bounded diameter at the 1e-3 rad^2 tilt premise: root information
    # (1/c_nu) nu nu', diffuse elsewhere; invariance of {P: nu'P nu <= c_nu}.
    c_nu = premises[KERNEL_PREMISE][1]
    unit = nu/np.linalg.norm(nu)
    prior = 1e8*(np.eye(21)-np.outer(unit, unit))+(c_nu/float(nu@nu))*np.outer(unit, unit)
    p_nu = replay(prior, word)[0]
    kappa_nu = gmax(p_nu, pi)
    # Invariance is checked for the kernel at the next root (end of the word).
    nu_end = kernel_vector(word, at_end=True)
    c_next = (math.sqrt(float(KERNEL_PREMISE))+float(np.linalg.norm(nu_end[18:]))*math.sqrt(BA_MARGINAL))**2
    propagated = float(nu_end@p_nu@nu_end)
    sig = lambda x: float(f'{x:.4g}')
    scale = np.sqrt(np.diag(p0))
    return {
        'word_s': word_s, 'root_time_s': round(root['time'], 3),
        'exact_margin': sig(1-rho),
        'slowest_direction_block_shares': shares(l@vec[:, -1], scale),
        'ideal_C': {'margin': sig(1-ideal[0]), 'kappa': sig(ideal[1]), 'k': sig(ideal[2])},
        'joint_reader_C': {'margin': sig(1-joint[0]), 'kappa': sig(joint[1]), 'k': sig(joint[2]),
                           'over_P0_lambda_max': sig(float(np.max(np.linalg.eigvalsh(chat)))),
                           'dominates_P0_within_1e-3': defect > -1e-3,
                           'diffuse_scale_converged_1e-6': conv < 1e-6,
                           'loss_factor': sig((1-rho)/(1-joint[0]))},
        'diameter': {'finite_kernel_observable': observable,
                     'kappa_information': sig(kappa_j) if observable else None,
                     'kappa_covariance': sig(kappa_d) if observable else None,
                     'relative_difference': sig(abs(kappa_j-kappa_d)/kappa_j) if observable else None,
                     'kappa_information_exceeds': None if observable else FINITE_DIAMETER,
                     'k0_margin': sig(2/(math.sqrt(kappa_j)+1)) if observable else 0.0,
                     'controlling_direction_block_shares': shares(direction, scale),
                     'slow_marginal': sig(kappa_slow) if observable else None,
                     'fast_given_slow': sig(kappa_fast_given_slow),
                     'gyro_persistence_cap': sig(gyro_cap(word_s))},
        'rank_one_kernel': {'margin_exact_nu_variance': sig(1-kern[0]), 'kappa': sig(kern[1]),
                            'k': sig(kern[2]), 'nu_P0_nu': sig(nu_var),
                            'tilt_variance_along_theta_hat': sig(tilt_var), 'nu_ba_norm': sig(nu_ba),
                            'BA_marginal_lambda_max': sig(ba_max),
                            'ceiling_actual_tilt': {'c_nu': sig(actual_tilt[1]), 'margin': sig(actual_tilt[0])},
                            'ceiling_tilt_premise': {t: {'c_nu': sig(m[1]), 'margin': sig(m[0])}
                                                     for t, m in premises.items()}},
        'kernel_bounded_diameter': {'tilt_premise_rad2': KERNEL_PREMISE, 'c_nu': sig(c_nu),
                                    'kappa_nu': sig(kappa_nu), 'margin_mu_equals_inverse_c': sig(1/kappa_nu),
                                    'next_root_kernel_ceiling': sig(c_next),
                                    'propagated_next_kernel_variance': sig(propagated),
                                    'kernel_rotation_rad': sig(math.acos(min(1.0, abs(float(nu[:3]@nu_end[:3]))))),
                                    'kernel_set_invariant': propagated <= c_next},
    }


def run_case(exes, profile, record_start, directory):
    stop = record_start+int(round((HISTORY_S+max(WORDS_S))/STEP))
    trace = Path(directory)/f'{profile}-{record_start}.bin'
    argv = [profile, str(record_start), str(stop), str(int(round(1/STEP))), str(trace)]
    observed = subprocess.run([str(exes['observed']), *argv], check=True,
                              capture_output=True, text=True).stdout.strip()
    control = subprocess.run([str(exes['control']), *argv[:-1], '/dev/null'], check=True,
                             capture_output=True, text=True).stdout.strip()
    if observed != control:
        raise ValueError(f'{profile}: observer changed the carried execution')
    ops, markers = operations(parse(trace))
    trace.unlink()
    end = observed.split()
    return {'profile': profile, 'record_start_s': round((record_start-1)*STEP, 3),
            'moving': profile in MOVING, 'live_step': int(end[1]), 'refined_step': int(end[2]),
            'active_step': int(end[3]), 'literal_terminal_parity': True,
            'words': [analyze_split(ops, markers, HISTORY_S, w) for w in WORDS_S]}


def diagnostic(eigen):
    with tempfile.TemporaryDirectory(prefix='ou3-info-ratio-') as tmp:
        exes = compile_driver(tmp, eigen)
        cases = [run_case(exes, p, s, tmp) for p, s in CASES]
    moving = [w['joint_reader_C']['loss_factor'] for c in cases if c['moving'] for w in c['words']]
    quiet = {c['record_start_s']: max(w['joint_reader_C']['loss_factor'] for w in c['words'])
             for c in cases if c['profile'] == 'quiet'}
    kernel_loss = max(w['exact_margin']/w['rank_one_kernel']['margin_exact_nu_variance']
                      for c in cases for w in c['words'])
    return {
        'qualification': 'OU3_INFORMATION_RATIO_SOURCE_FEASIBILITY_V1',
        'role': 'non-promoting carried feasibility of the information-ratio word contraction',
        'driver_sha256': hashlib.sha256(DRIVER.read_bytes()).hexdigest(),
        'header_sha256': hashlib.sha256((REPO/HEADER).read_bytes()).hexdigest(),
        'history_s': HISTORY_S, 'words_s': list(WORDS_S),
        'kill_criterion': {'factor': KILL_FACTOR, 'applies_to': 'MOVING words (wave, collinear 25 Hz, sync-locked)',
                           'worst_moving_joint_reader_loss': max(moving),
                           'passed': max(moving) <= KILL_FACTOR,
                           'quiet_joint_reader_loss_by_record_start_s': {str(k): v for k, v in quiet.items()}},
        'worst_rank_one_kernel_loss': round(kernel_loss, 3),
        'kernel_set_invariant_on_every_word': all(w['kernel_bounded_diameter']['kernel_set_invariant']
                                                  for c in cases for w in c['words']),
        'cases': cases,
        'replay': 'float64 optimal-gain replay of literal exported coefficients',
        'source_uniform_verified': False, 'rigorous_enclosure': False, 'theorem_closed': False,
    }


def differences(a, b, rel, path='$'):
    """Leaves where two records differ beyond a relative tolerance.

    Native float32 replays differ at about 1e-6 relative between toolchains
    (see the CI native replay binding in the research ledger); the record
    keeps four significant digits and only well-posed quantities. The diameter
    identity residual ``relative_difference`` is itself replay noise; its
    contract (<1e-3) is enforced by verify_diagnostic, so here both records
    only have to stay below a noise floor an order of magnitude tighter.
    """
    if path.endswith('.relative_difference') and a is not None and b is not None:
        return [] if max(a, b) <= RELATIVE_DIFFERENCE_FLOOR else [f'{path}: {a} != {b}']
    if isinstance(a, dict) and isinstance(b, dict):
        if a.keys() != b.keys():
            return [f'{path}: keys {sorted(set(a) ^ set(b))}']
        return [d for k in a for d in differences(a[k], b[k], rel, f'{path}.{k}')]
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f'{path}: length {len(a)} != {len(b)}']
        return [d for i, (x, y) in enumerate(zip(a, b)) for d in differences(x, y, rel, f'{path}[{i}]')]
    numeric = (int, float)
    if isinstance(a, numeric) and isinstance(b, numeric) and not isinstance(a, bool) and not isinstance(b, bool):
        return [] if abs(a-b) <= rel*max(abs(a), abs(b), 1e-12) else [f'{path}: {a} != {b}']
    return [] if a == b else [f'{path}: {a!r} != {b!r}']


def verify_diagnostic(record):
    """Internal consistency of a committed record; no replay required."""
    if record.get('qualification') != 'OU3_INFORMATION_RATIO_SOURCE_FEASIBILITY_V1':
        raise ValueError('wrong qualification')
    if record['source_uniform_verified'] or record['rigorous_enclosure'] or record['theorem_closed']:
        raise ValueError('diagnostic cannot promote')
    if record['driver_sha256'] != hashlib.sha256(DRIVER.read_bytes()).hexdigest():
        raise ValueError('information-ratio driver changed')
    if record['header_sha256'] != hashlib.sha256((REPO/HEADER).read_bytes()).hexdigest():
        raise ValueError('shipping header changed since the committed audit')
    if [(c['profile'], round(c['record_start_s']/STEP)+1) for c in record['cases']] != list(CASES):
        raise ValueError('case set changed')
    moving = []
    for case in record['cases']:
        if not case['literal_terminal_parity']:
            raise ValueError('observer parity missing')
        for w in case['words']:
            exact = w['exact_margin']
            for key in ('ideal_C', 'joint_reader_C'):
                if w[key]['margin'] > exact*(1+2e-3):
                    raise ValueError(f'{key} bound exceeds the exact margin')
            if w['joint_reader_C']['margin'] > w['ideal_C']['margin']*(1+2e-3):
                raise ValueError('joint-reader bound beats the ideal comparison')
            d = w['diameter']
            if d['finite_kernel_observable'] and d['kappa_information'] < d['gyro_persistence_cap']*(1-2e-3):
                raise ValueError('gyro-bias persistence cap violated')
            if not w['joint_reader_C']['dominates_P0_within_1e-3'] or not w['joint_reader_C']['diffuse_scale_converged_1e-6']:
                raise ValueError('joint-reader upper comparison did not converge or dominate')
            if d['finite_kernel_observable']:
                if d['relative_difference'] > 1e-3:
                    raise ValueError('diameter identity violated')
                if d['kappa_information'] > d['slow_marginal']*d['fast_given_slow']*(1+2e-3):
                    raise ValueError('slow/fast factorization violated')
            if d['k0_margin'] > exact*(1+2e-3):
                raise ValueError('k=0 bound exceeds the exact margin')
            r = w['rank_one_kernel']
            if r['margin_exact_nu_variance'] > exact*(1+2e-3) or r['BA_marginal_lambda_max'] > BA_MARGINAL*(1+1e-6):
                raise ValueError('kernel bound or BA marginal ceiling violated')
            for m in (r['ceiling_actual_tilt'], *r['ceiling_tilt_premise'].values()):
                if m['c_nu'] < 0 or m['margin'] > exact*(1+2e-3):
                    raise ValueError('scalar-ceiling kernel bound inconsistent')
            if r['ceiling_actual_tilt']['c_nu'] < r['nu_P0_nu']*(1-2e-3):
                raise ValueError('scalar ceiling below the actual kernel variance')
            kb = w['kernel_bounded_diameter']
            if kb['margin_mu_equals_inverse_c'] > exact*(1+2e-3) or kb['kappa_nu'] < 1:
                raise ValueError('kernel-bounded diameter bound inconsistent')
            if kb['kernel_set_invariant'] != (kb['propagated_next_kernel_variance'] <= kb['next_root_kernel_ceiling']) and \
                    abs(kb['propagated_next_kernel_variance']-kb['next_root_kernel_ceiling']) > 5e-4*kb['c_nu']:
                raise ValueError('kernel invariance flag inconsistent')
            if case['moving']:
                moving.append(w['joint_reader_C']['loss_factor'])
    kill = record['kill_criterion']
    if kill['worst_moving_joint_reader_loss'] != max(moving) or kill['passed'] != (max(moving) <= kill['factor']):
        raise ValueError('kill criterion summary inconsistent')
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--eigen', type=Path, default=Path('/usr/include/eigen3'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expect', type=Path)
    args = parser.parse_args()
    record = diagnostic(args.eigen)
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')
    verify_diagnostic(record)
    if args.expect:
        diff = differences(json.loads(args.expect.read_text()), record, 5e-3)
        if diff:
            print('\n'.join(diff))
            raise SystemExit('carried information-ratio audit differs from the committed record')


if __name__ == '__main__':
    main()
