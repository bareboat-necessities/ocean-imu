#!/usr/bin/env python3
"""Validated complete-word translation dissipation on the limiting P3 cell.

The old P4 scalar route kept only one covariance seed and discarded later
positive process/Joseph contributions.  This backend instead accumulates a
complete translation word in the P3-conditioned coordinates

    D = diag(sigma_min*h, sigma_min*h^2, sigma_min*h^3, sigma_min).

The source-dependent transition remains interval-valued, but the propagated
Loewner lower is kept as an exact rational matrix.  This is important for the
strongly anisotropic recurrent worst cell: repeatedly converting that lower to
a binary64 point matrix was losing SPD from ~1e-18 representation error even
when the exact lower remained positive.

For prediction, write F = Fc + E.  Since L >= 0,

    F L F' >= Fc L Fc' - gamma I,

where gamma outward-bounds the symmetric Fc L E' + E L Fc' cross term.  The
midpoint product and all rank-one Woodbury corrections are evaluated exactly
as Fractions.  Only genuine source uncertainty enters through outward-rounded
binary64 bounds used to form gamma.  Endpoint positivity of
L - delta*Sigma_upper is then checked by exact rational LDL/Schur pivots.

The S pseudo acts only in S and the translational accelerometer correction only
in a_w.  Both possible corrections are still applied at every IMU sample, which
is conservative relative to shipping scheduling/rejection.  This file proves
only the previously limiting translation source cell; it cannot by itself
promote the full P4 certificate.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from fractions import Fraction
from pathlib import Path

from ou3_interval import Interval, matrix_mul
import ou3_validated_transcendentals as VT
import ou3_p4_worst_translation_cell as WORST

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / 'tools' / 'ou3_proof_operating_domain.json'
WRAPPER = REPO / 'src' / 'kalman_ou_iii' / 'SeaStateFusionFilter_OU_III.h'
DEFAULT_HORIZON_S = 1.0
MAX_TAU_SPLIT_DEPTH = 14
N = 4


def down(x):
    return math.nextafter(float(x), -math.inf)


def up(x):
    return math.nextafter(float(x), math.inf)


def I(x):
    return Interval.outward_bounds(float(x), float(x))


def _q(x):
    return x if isinstance(x, Fraction) else Fraction.from_float(float(x))


def _qm(A):
    return [[_q(x) for x in row] for row in A]


def _down_fraction(q):
    f = float(q)
    if not math.isfinite(f):
        raise OverflowError('exact rational does not fit binary64')
    if Fraction.from_float(f) > q:
        f = math.nextafter(f, -math.inf)
    return f


def _up_fraction(q):
    f = float(q)
    if not math.isfinite(f):
        raise OverflowError('exact rational does not fit binary64')
    if Fraction.from_float(f) < q:
        f = math.nextafter(f, math.inf)
    return f


def _qi(q):
    """Outward interval enclosure of one exact rational scalar."""
    q = _q(q)
    return Interval.outward_bounds(_down_fraction(q), _up_fraction(q))


def _qpm(A):
    return [[_qi(x) for x in row] for row in A]


def _rational_spd(A):
    """Exact positive-definiteness test by unpivoted symmetric Schur pivots."""
    M = _qm(A)
    n = len(M)
    for k in range(n):
        p = M[k][k]
        if p <= 0:
            return False
        for i in range(k + 1, n):
            for j in range(i, n):
                v = M[i][j] - M[i][k] * M[j][k] / p
                M[i][j] = v
                M[j][i] = v
    return True


def _F(tau, h):
    x = I(h) / tau
    if x.hi >= 1e-2:
        raise RuntimeError('limiting cell left shipping small-x branch')
    a = VT.exp_interval(-x)
    em1 = VT.expm1_interval(-x)
    pva = -(tau * em1)
    x2 = x * x
    x3 = x2 * x
    x4 = x3 * x
    x5 = x4 * x
    ppa = tau * tau * (I(.5) * x2 - I(1 / 6) * x3 + I(1 / 24) * x4)
    psa = tau * tau * tau * (I(1 / 6) * x3 - I(1 / 24) * x4 + I(1 / 120) * x5)
    z = Interval.point(0)
    o = Interval.point(1)
    return [
        [o, z, z, pva / I(h)],
        [o, o, z, ppa / I(h * h)],
        [I(.5), o, o, psa / I(h * h * h)],
        [z, z, z, a],
    ]


def _midrad(A):
    C, R = [], []
    for row in A:
        cr, rr = [], []
        for a in row:
            c = min(max(.5 * a.lo + .5 * a.hi, a.lo), a.hi)
            r = up(max(abs(a.lo - c), abs(a.hi - c)))
            cr.append(c)
            rr.append(r)
        C.append(cr)
        R.append(rr)
    return C, R


def _abs_upper(a):
    return up(max(abs(a.lo), abs(a.hi)))


def _predict(L, F, rho):
    """Exact-rational Loewner lower for F L F' + rho I."""
    if not _rational_spd(L):
        raise RuntimeError('exact lower entered prediction non-SPD')
    Fc, R = _midrad(F)

    # Enclose Fc*L only for the genuine interval cross-term bound.  L itself
    # stays exact rational; no binary64 point lower is formed.
    B = matrix_mul([[I(x) for x in row] for row in Fc], _qpm(L))
    gamma = 0.0
    for i in range(N):
        rowsum = 0.0
        for j in range(N):
            cij = 0.0
            for k in range(N):
                cij = up(cij + up(_abs_upper(B[i][k]) * R[j][k]))
                cij = up(cij + up(R[i][k] * _abs_upper(B[j][k])))
            rowsum = up(rowsum + cij)
        gamma = max(gamma, rowsum)

    A = [[Fraction.from_float(float(Fc[i][j])) for j in range(N)] for i in range(N)]
    Q = _qm(L)
    AQ = [[sum((A[i][k] * Q[k][j] for k in range(N)), Fraction(0)) for j in range(N)] for i in range(N)]
    M = [[sum((AQ[i][k] * A[j][k] for k in range(N)), Fraction(0)) for j in range(N)] for i in range(N)]
    shift = Fraction.from_float(float(rho)) - Fraction.from_float(float(gamma))
    for i in range(N):
        M[i][i] += shift
    if not _rational_spd(M):
        raise RuntimeError(f'exact prediction lower lost SPD (source radius={gamma:.3e})')
    return M, gamma


def _rank1_information_update_lower(L, beta, qidx):
    """Exact lower for (L^-1 + beta e_q e_q')^-1."""
    if not (0 <= qidx < N and beta >= 0 and math.isfinite(beta)):
        raise RuntimeError('invalid rank-one information update')
    if beta == 0:
        return _qm(L), 0.0
    A = _qm(L)
    if not _rational_spd(A):
        raise RuntimeError('exact lower entered correction non-SPD')
    b = Fraction.from_float(float(beta))
    den = Fraction(1) + b * A[qidx][qidx]
    if den <= 0:
        raise RuntimeError('rank-one information denominator lost positivity')
    M = [[A[i][j] - b * A[i][qidx] * A[qidx][j] / den for j in range(N)] for i in range(N)]
    if not _rational_spd(M):
        raise RuntimeError('exact rank-one posterior lower lost SPD')
    return M, 0.0


def _spd_delta(L, upper, d):
    A = _qm(L)
    qd = Fraction.from_float(float(d))
    for i in range(N):
        # upper[i] is already an outward upper binary64 bound.
        A[i][i] -= qd * Fraction.from_float(float(upper[i]))
    return _rational_spd(A)


def _delta(L, upper):
    lo = 0.0
    trial = 1e-36
    while trial < 1 and _spd_delta(L, upper, trial):
        lo = trial
        trial *= 10
    if lo == 0:
        return 0.0
    hi = min(1.0, trial)
    for _ in range(64):
        mid = math.sqrt(lo * hi)
        if _spd_delta(L, upper, mid):
            lo = mid
        else:
            hi = mid
    return down(lo)


def _member(text, name):
    m = re.search(rf'float\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f\s*;', text)
    if not m:
        raise RuntimeError(f'cannot source-bind {name}')
    return float(m.group(1))


def _prop(tau, h, n, rho, betaS, betaA, depth=0):
    try:
        F = _F(tau, h)
        L = [[Fraction(0) for _ in range(N)] for _ in range(N)]
        mr = 0.0
        for k in range(n):
            if k == 0:
                r = Fraction.from_float(float(rho))
                L = [[r if i == j else Fraction(0) for j in range(N)] for i in range(N)]
            else:
                L, radius = _predict(L, F, rho)
                mr = max(mr, radius)
            L, radius = _rank1_information_update_lower(L, betaS, 2)
            mr = max(mr, radius)
            L, radius = _rank1_information_update_lower(L, betaA, 3)
            mr = max(mr, radius)
        return [(tau, L, mr, depth)]
    except RuntimeError:
        if depth >= MAX_TAU_SPLIT_DEPTH:
            raise
        mid = math.sqrt(tau.lo * tau.hi)
        a = Interval.outward_bounds(tau.lo, mid)
        b = Interval.outward_bounds(mid, tau.hi)
        return (_prop(a, h, n, rho, betaS, betaA, depth + 1) +
                _prop(b, h, n, rho, betaS, betaA, depth + 1))


def _mode(mode, domain_path, horizon_s):
    c = WORST.build_cell(mode, domain_path)
    s = WORST.serializable(c)
    row = c['row']
    h = float(c['sched']['dt_s'])
    x = c['x']
    tau = Interval.outward_bounds(h / x.hi, h / x.lo)
    sigma = float(c['sigma'].lo)
    rho = float(c['rho_translation_lower'])

    scale2 = [(sigma * h) ** 2, (sigma * h * h) ** 2,
              (sigma * h * h * h) ** 2, sigma * sigma]
    u = list(map(float, row['Sigma_diagonal_upper']))
    physical = [u[6], u[9], u[12], u[15]]
    upper = [(I(physical[i]) / I(scale2[i])).hi for i in range(N)]

    text = WRAPPER.read_text()
    rh = min(_member(text, 'R_S_x_factor_'), _member(text, 'R_S_y_factor_'), 1.0)
    rs = (I(rh) * I(float(c['rs'].lo))).lo
    rsvar = I(rs).square().lo
    acc = float(c['vector']['configured_measurement_bounds']['acc_measurement_std_mps2'])
    accvar = I(acc).square().lo
    rS = (I(rsvar) / I(scale2[2])).lo
    rA = (I(accvar) / I(scale2[3])).lo
    betaS = (I(1) / I(rS)).hi
    betaA = (I(1) / I(rA)).hi

    n = int(math.ceil(horizon_s / h))
    leaves = _prop(tau, h, n, rho, betaS, betaA)
    cert = []
    for t, L, radius, dep in leaves:
        d = _delta(L, upper)
        if d <= 0 or not _spd_delta(L, upper, d):
            raise RuntimeError(f'nonpositive endpoint margin on tau leaf {t.as_list()}')
        cert.append({
            'tau_s': t.as_list(),
            'delta_lower': d,
            'max_conditioned_radius_removed': radius,
            'split_depth': dep,
        })

    w = min(cert, key=lambda q: q['delta_lower'])
    old = float(row['direct_translation_generalized_margin_lower'])
    return {
        'source_cell': s,
        'conditioned_coordinates': 'D^-1[v,p,S,a_w]',
        'tau_interval_s': tau.as_list(),
        'tau_leaf_count': len(cert),
        'max_tau_split_depth_used': max(q['split_depth'] for q in cert),
        'steps': n,
        'horizon_s': horizon_s,
        'process_injection_lower_conditioned': rho,
        'S_measurement_information_beta_conditioned': betaS,
        'accelerometer_aw_information_beta_conditioned': betaA,
        'measurement_information_geometry': 'rank_one_S_and_aw_each_sample_exact_rational',
        'prediction_enclosure': 'midpoint_plus_symmetric_cross_term_loewner',
        'exact_rational_lower_retained_through_word': True,
        'corrections_allowed_every_sample_for_lower_bound': True,
        'artificial_S_variance_conditioned': rS,
        'artificial_acc_aw_variance_conditioned': rA,
        'translation_covariance_upper_conditioned': upper,
        'tau_leaf_certificates': cert,
        'complete_word_translation_margin_lower': w['delta_lower'],
        'limiting_tau_leaf': w,
        'old_single_seed_translation_margin_lower': old,
        'margin_widening_factor_lower': down(w['delta_lower'] / old),
        'interval_ldlt_endpoint_recertified': True,
    }


def build(domain_path=DEFAULT_DOMAIN, horizon_s=DEFAULT_HORIZON_S):
    horizon_s = float(horizon_s)
    if not math.isfinite(horizon_s) or horizon_s <= 0:
        raise ValueError('horizon_s must be finite positive')
    p = Path(domain_path).resolve()
    modes, failures = {}, []
    for mode in ('H', 'A'):
        try:
            modes[mode] = _mode(mode, p, horizon_s)
        except Exception as e:
            failures.append(f'{mode}: {e}')
    return {
        'qualification': 'OU3_P4_VALIDATED_WORST_CELL_COMPLETE_WORD_TRANSLATION_DISSIPATION',
        'source_only': True,
        'trajectory_replay_used': False,
        'outward_rounded': True,
        'horizon_s': horizon_s,
        'modes': modes,
        'P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS': 'PASS' if not failures and len(modes) == 2 else 'NOT_ESTABLISHED',
        'P4_USABLE_CERTIFICATE_STATUS': 'NOT_ESTABLISHED',
        'remaining_obligation': 'extend directional complete-word propagation to every reachable source cell/edge and attitude-bias blocks, then validate exact nonlinear return map',
        'failures': failures,
    }


def validate(d):
    f = list(d.get('failures', []))
    if d.get('source_only') is not True or d.get('trajectory_replay_used') is not False or d.get('outward_rounded') is not True:
        f.append('qualification flags invalid')
    if not float(d.get('horizon_s', 0)) > 0:
        f.append('invalid horizon')
    for mode in ('H', 'A'):
        m = d.get('modes', {}).get(mode, {})
        if not float(m.get('complete_word_translation_margin_lower', 0)) > 0:
            f.append(f'{mode}: no complete-word translation margin')
        if m.get('interval_ldlt_endpoint_recertified') is not True:
            f.append(f'{mode}: endpoint not recertified')
        if m.get('measurement_information_geometry') != 'rank_one_S_and_aw_each_sample_exact_rational':
            f.append(f'{mode}: directional measurement geometry missing')
        if m.get('prediction_enclosure') != 'midpoint_plus_symmetric_cross_term_loewner':
            f.append(f'{mode}: structured prediction enclosure missing')
        if m.get('exact_rational_lower_retained_through_word') is not True:
            f.append(f'{mode}: exact rational lower was not retained through word')
        if m.get('corrections_allowed_every_sample_for_lower_bound') is not True:
            f.append(f'{mode}: lower no longer covers maximum correction frequency')
        if not float(m.get('margin_widening_factor_lower', 0)) > 1:
            f.append(f'{mode}: complete word did not widen seed')
    if d.get('P4_USABLE_CERTIFICATE_STATUS') != 'NOT_ESTABLISHED':
        f.append('partial result prematurely promoted P4')
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--domain', type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument('--horizon-s', type=float, default=DEFAULT_HORIZON_S)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain, a.horizon_s)
    f = validate(d)
    d['validation_failures'] = f
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True))
    print(json.dumps({
        'translation_status': d['P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS'],
        'horizon_s': d['horizon_s'],
        'modes': {x: {
            'delta': d.get('modes', {}).get(x, {}).get('complete_word_translation_margin_lower'),
            'factor': d.get('modes', {}).get(x, {}).get('margin_widening_factor_lower'),
            'tau_leaves': d.get('modes', {}).get(x, {}).get('tau_leaf_count'),
            'geometry': d.get('modes', {}).get(x, {}).get('measurement_information_geometry'),
            'prediction': d.get('modes', {}).get(x, {}).get('prediction_enclosure'),
        } for x in ('H', 'A')},
        'failures': f,
    }, indent=2))
    return 0 if not f else 2


if __name__ == '__main__':
    raise SystemExit(main())
