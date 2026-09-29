"""Aggregate world-frame six-column floor: tube lemma and Theorem G0.

Proof: docs/ou3-world-frame-rows.md, sections 12--14. Enters
V_next<=rho V+c_d|d|^2 only through the six-column floor s of the historical
reader (then B_*, J_AG, P_max, rho_0). Theorem G0 is stated for the
INJECTION-FREE nominal world array (A~=I, B~=int R_hat'); its nominal window
premises (signed mean, L1 force) are measured, not proved, on actual
executions. All constants are exact rationals; cos is bounded below by its
alternating Taylor polynomial and square roots by integer square roots.
The downstream chain is a non-promoting feasibility computation.
"""
from fractions import Fraction as F
from math import isqrt

from .aw_tracking import G, FRACTION
from .signed_temporal import separated_reader_action_implication
from .world_frame import decimal_lower

OMEGA_INVARIANT = F('0.6108652381980153')+F('0.02')+F('0.02')+F('0.5')
OMEGA_PHYSICAL = F('0.6108652381980153')+F('0.02')+F('0.02')
B_MIN = F(20)                   # uT, MAGNETIC SERVICE field floor
H_MAX = F(6, 1000)


def cos_lower(x):
    """1-x^2/2+x^4/24-x^6/720 <= cos x for every real x (integrated chain)."""
    x = F(x)
    return 1-x**2/2+x**4/24-x**6/720


def sqrt_lower(v, digits=15):
    v = F(v)
    if v <= 0:
        return F(0)
    scale = 10**digits
    return F(isqrt(v.numerator*scale*scale//v.denominator), scale)


def sqrt_upper(v, digits=15):
    low = sqrt_lower(v, digits)
    return low if low*low == F(v) else low+F(1, 10**digits)


def tube_constant(omega, gap, eps, room):
    """Lemma T: |b.T(t*)| lower bound at a point with room l>=Delta each side.

    T=R_hat(t)'beta/|beta| is the unit tangent of c(t)=theta+int R_hat' beta,
    turning at most Omega per second; magnetic rows at gaps <=Delta see
    |P_b c(t_j)|<=eps|beta|. With c=arcsin|b.T(t*)| and e the transverse
    direction of T(t*), take samples straddling t* in the span where
    Omega|u-t*|+c<=pi/2 and integrate e.T along it:
    (A) if (pi/2-c)/Omega<=l: (2/Omega)(cos(Omega Delta)-sin c)<=2 eps;
    (B) otherwise (possible only if Omega l<pi/2), for every eps<r<=l-Delta:
        2r cos(c+Omega r)<=2 eps, so sin c>=sqrt(1-x^2)cos(Omega r)-x sin(Omega r),
        x=eps/r; the best r of a grid is taken.
    The bound is the smaller of both (rational: cos_lower, sin y<=y).
    """
    om, d, e, ell = map(F, (omega, gap, eps, room))
    case_a = cos_lower(om*d)-om*e
    if om*ell >= F(15708, 10000):          # 1.5708>pi/2: case B cannot occur
        return case_a
    best = None
    for k in range(1, 9):
        r = (ell-d)*k/8
        if r <= e:
            continue
        x = e/r
        bound = sqrt_lower(1-x*x)*cos_lower(om*r)-x*om*r
        best = bound if best is None else max(best, bound)
    return min(case_a, best if best is not None else F(0))


def theorem_g0(omega, gap, eps, window, separation, sigma_w, m_perp, force_l1,
               rows_per_window, start_offset, gravity=G, field_min=B_MIN, scale=F(1)):
    """Injection-free aggregate six-column floor (Theorem G0), exact rational.

    The word extends start_offset>=Delta before W1 and after W2 (the room of
    Lemma T); W1, W2 have length L and are separated by a gap G; every
    interval of length Delta contains an applied magnetic row; each window has
    nominal transverse mean
    |mean(a_hat) x b|<=m_perp and L1 normalized force mean(|f_hat|/g)<=u1.
    Part I (beta): Q>=q_I|beta|^2, q_I=min(B^2 eps^2, n g^2 (c0 G mu/2-rho u1)^2),
    mu=sigma_w-m_perp/g, rho=eps+(Delta/2)sqrt(1-c0^2).
    Part II (theta): Q>=K(sqrt(gamma*)|theta|-C|beta|)_+^2, K=min(n g^2,B^2),
    gamma*=mu^2/(u1^2+1), C=sqrt(gamma*) t_c+sqrt((u1 L/2)^2+(Delta/2)^2),
    t_c the W1 centre from the word start.
    Combination: s^2 >= q_I k^2/(1+lambda^2 k^2), k=sqrt(K gamma*)/(sqrt(q_I)+sqrt(K) C),
    for coordinates (theta, lambda beta).
    """
    om, d, e, L, gap_len, sw, mp, u1 = map(F, (omega, gap, eps, window, separation,
                                                sigma_w, m_perp, force_l1))
    g, bmin, lam = F(gravity), F(field_min), F(scale)
    n = F(rows_per_window)
    c0 = tube_constant(om, d, e, start_offset)
    if c0 <= 0:
        return {'tube_constant': c0, 'positive': False, 'reason': 'Omega Delta too large or no room'}
    mu = sw-mp/g
    if mu <= 0:
        return {'tube_constant': c0, 'positive': False, 'reason': 'nominal mean premise'}
    rho = e+d/2*sqrt_upper(1-c0*c0)
    acc_margin = c0*gap_len*mu/2-rho*u1
    if acc_margin <= 0:
        return {'tube_constant': c0, 'positive': False, 'reason': 'separation too short',
                'acc_margin': acc_margin}
    q_i = min(bmin**2*e**2, n*g**2*acc_margin**2)
    gamma = mu*mu/(u1*u1+1)
    k_ii = min(n*g**2, bmin**2)
    t_c = F(start_offset)+L/2
    c_ii = sqrt_upper(gamma)*t_c+sqrt_upper((u1*L/2)**2+(d/2)**2)
    k = sqrt_lower(k_ii*gamma)/(sqrt_upper(q_i)+sqrt_upper(k_ii)*c_ii)
    s2 = q_i*k*k/(1+lam*lam*k*k)
    return {'tube_constant': c0, 'mu': mu, 'rho': rho, 'acc_margin': acc_margin,
            'q_I': q_i, 'gamma_star': gamma, 'K_II': k_ii, 'C_II': c_ii, 'k': k,
            'six_column_floor_squared': s2, 'positive': s2 > 0}


def optimize_eps(omega, gap, **kw):
    """Best eps on a rational grid (the floor is valid for every eps>0)."""
    best = None
    for i in range(1, 60):
        eps = F(i, 200)
        r = theorem_g0(omega, gap, eps, **kw)
        if r.get('positive') and (best is None or r['six_column_floor_squared'] > best[1]['six_column_floor_squared']):
            best = (eps, r)
    return best


def tube_audit(dps=40):
    """80-digit audits of Lemma T on supplied nominal rotation histories.

    (a) rotation about the field axis: tangent constant along b, the lemma's
    hypothesis holds with eps=0 and |b.T|=1; (b) steady rolling at Omega with
    beta transverse: magnetic residuals are large, so the hypothesis fails at
    some sample (the magnetic rows see beta); (c) the planar circle tangent to
    the b-line at t=0 returns to it after 2 pi/Omega: with that sampling gap the
    tangent's b-component changes sign, so some condition like Omega Delta<pi/2
    is necessary for monotonicity.
    """
    import mpmath as mp
    with mp.workdps(dps):
        b = mp.matrix([mp.mpf(7)/25, 0, mp.mpf(24)/25])
        om = mp.mpf(OMEGA_INVARIANT.numerator)/OMEGA_INVARIANT.denominator

        def proj(v):
            return v-(b.T*v)[0]*b
        # (c) circle: c(t)=(1/om)(1-cos om t) e+(1/om) sin(om t) b, e=(24,0,-7)/25.
        e = mp.matrix([mp.mpf(24)/25, 0, -mp.mpf(7)/25])
        period = 2*mp.pi/om
        samples = [k*period for k in range(4)]
        residual = max(mp.norm(proj((1-mp.cos(om*t))/om*e+mp.sin(om*t)/om*b)) for t in samples)
        tangent_b = [mp.cos(om*t) for t in (0, period/2)]
        # (b) roll about e_x at om; beta=e_y: tangent R'beta rotates in y-z plane.
        worst = mp.mpf(0)
        for k in range(1, 21):
            t = mp.mpf(k)/2
            ct = mp.matrix([0, mp.sin(om*t)/om, (1-mp.cos(om*t))/om])
            worst = max(worst, mp.norm(proj(ct)))
        return {'circle_gap_s': mp.nstr(period, 12),
                'circle_max_sample_residual': mp.nstr(residual, 5),
                'circle_tangent_b_component_at_0_and_half_period': [mp.nstr(x, 6) for x in tangent_b],
                'monotonicity_needs_gap_below_2pi_over_Omega': True,
                'rolling_transverse_beta_max_magnetic_residual': mp.nstr(worst, 12),
                'field_axis_rotation_tangent_b_component': '1'}


def synthetic_audit(seed=7):
    """Falsification audit: actual sigma_min of injection-free arrays vs G0.

    Nominal rotations (field-axis spin, transverse roll, frozen, random
    piecewise rates) at |omega_hat|<=Omega, magnetic rows every Delta=1 s,
    accelerometer rows every .05 s inside W1=[2,18] and W2=[82,98] with
    nominal forces whose window mean is transverse-limited and whose L1 size
    is bounded. The G0 floor is recomputed from each array's own premises
    (with its own row count) and must not exceed the actual floor.
    Float64 linear algebra; non-promoting.
    """
    import math
    import numpy as np
    rng = np.random.default_rng(seed)
    g, bnorm = float(G), 20.0
    b = np.array([7.0, 0.0, 24.0])/25.0
    om = float(OMEGA_INVARIANT)
    dt = 0.01
    times = np.arange(0.0, 100.0+1e-9, dt)

    def skew(v):
        return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

    def expm(w):
        th = np.linalg.norm(w)
        if th < 1e-15:
            return np.eye(3)
        k = skew(w/th)
        return np.eye(3)+np.sin(th)*k+(1-np.cos(th))*k@k

    def run(rate_fn, force_fn):
        r = np.eye(3)
        gam = np.zeros((3, 3))
        rot, gams = [], []
        for t in times:
            rot.append(r.copy())
            gams.append(gam.copy())
            w = rate_fn(t)
            # body rate; R_hat world->body: R+=Exp(-w dt)R; Gamma+=int R'.
            step = expm(-w*dt)
            gam = gam+r.T*dt
            r = step@r
        rows, forces_w = [], {1: [], 2: []}
        for i, t in enumerate(times):
            tt = round(t, 6)
            if abs(tt-round(tt)) < 1e-9:
                rows.append(bnorm*skew(b)@np.hstack([np.eye(3), gams[i]]))
            win = 1 if 2 <= tt <= 18 else 2 if 82 <= tt <= 98 else 0
            if win and abs(tt*20-round(tt*20)) < 1e-9:
                u = force_fn(tt, win)
                forces_w[win].append(u)
                rows.append(g*skew(u)@np.hstack([np.eye(3), gams[i]]))
        o = np.vstack(rows)
        actual = np.linalg.svd(o, compute_uv=False).min()
        mperp = max(np.linalg.norm(np.cross(np.mean(forces_w[w], axis=0)+np.array([0, 0, 1.0]), b))*g
                    for w in (1, 2))
        l1 = max(np.mean([np.linalg.norm(u) for u in forces_w[w]]) for w in (1, 2))
        n = min(len(forces_w[1]), len(forces_w[2]))
        # Round the measured premises up on a 1e-6 grid (platform-stable).
        up = lambda x: F(math.ceil(float(x)*10**6)+1, 10**6)
        best = optimize_eps(OMEGA_INVARIANT, 1, window=16, separation=64, sigma_w=FRACTION,
                            m_perp=up(mperp), force_l1=up(l1), rows_per_window=n, start_offset=2)
        floor = float(best[1]['six_column_floor_squared']) if best else 0.0
        return float(actual), floor**.5

    down = np.array([0.0, 0.0, -1.0])
    axis_b = b.copy()
    wobble = lambda t, w: down+np.array([2.5*np.cos(9.52*t), 0.0, 0.0])/g
    cases = {
        'field_axis_spin': (lambda t: om*axis_b, wobble),
        'transverse_roll': (lambda t: om*np.array([1.0, 0.0, 0.0]), wobble),
        'frozen_attitude': (lambda t: np.zeros(3), wobble),
    }
    rates = rng.normal(size=(101, 3))
    rates = om*rates/np.linalg.norm(rates, axis=1)[:, None]
    forces = rng.normal(size=(2, 400, 3))*0.3
    cases['random_piecewise'] = (lambda t: rates[int(t)],
                                 lambda t, w: down+forces[w-1][int(round((t-(2 if w == 1 else 82))*20)) % 400]/1.0*0.3)
    out = {}
    for name, (rate_fn, force_fn) in cases.items():
        actual, floor = run(rate_fn, force_fn)
        assert floor <= actual
        out[name] = {'actual_sigma_min': float(f'{actual:.6g}'), 'G0_floor': float(f'{floor:.6g}')}
    return out


def feasibility(nominal_mean_perp, nominal_force_l1):
    """Theorem G0 at the shipping-relevant parameter points.

    Delta=1 s reads MAGNETIC SERVICE on every sliding 1-s window, Delta=2 s on
    disjoint windows only. Omega_invariant uses the .5 rad/s estimated-bias
    invariant; Omega_physical assumes the bias estimate within the physical
    envelope (retained region). 16-s windows, n>=16/.006 acc rows each.
    """
    rows = []
    n = int(F(16)/H_MAX)
    for label, omega in (('invariant', OMEGA_INVARIANT), ('physical', OMEGA_PHYSICAL)):
        for gap in (F(1), F(2)):
            for separation in (F(32), F(64)):
                best = optimize_eps(omega, gap, window=16, separation=separation,
                                    sigma_w=FRACTION, m_perp=nominal_mean_perp,
                                    force_l1=nominal_force_l1, rows_per_window=n,
                                    start_offset=2)
                rows.append({'Omega': label, 'magnetic_gap_s': str(gap),
                             'window_separation_s': str(separation),
                             'positive': best is not None,
                             'eps_s': str(best[0]) if best else None,
                             'six_column_floor_squared_lower':
                                 str(decimal_lower(best[1]['six_column_floor_squared'], 15)) if best else '0'})
    return rows


def downstream(six_column_floor_squared, gyro_schur_floor, word_s=100,
               bias_rw=F(1, 10**11), actual_floor=F(37)):
    """Push s to pivots, the scalar B_*, a least-squares reader and rho_0.

    A word of T s at h in [.004,.006] has at most 3(2T/.004) AG rows (one acc
    and at most one applied mag correction per prediction cell) and at least
    2T/.006 predictions plus acc corrections. Six greedy pivots >= s/sqrt(m).
    The existing scalar implication multiplies a coefficient ceiling >= 20
    (a magnetic row norm) over every operation. A least-squares reader
    L=T_h(O'O)^-1 O' has measurement action <= R_max|T_h|^2/s^2 (R_max=2^2:
    magnetic residual bound; the accelerometer .30104^2 is smaller); its gyro
    block is <= R_max/q with q the gyro Schur floor. Iterating the
    corrected-word prediction comparison over every prediction gives
    1-rho_0 <= T b0/P_bg,max: tiny unless the ceiling is near the true P_bg.
    """
    import math
    s2, q = F(six_column_floor_squared), F(gyro_schur_floor)
    rows = 3*int(2*word_s/F(4, 1000))
    operations = int(2*word_s/F(6, 1000))
    pivot = sqrt_lower(s2/rows)
    log_bstar = operations*math.log10(20)
    r_max = F(4)
    ls_action, bias_block = r_max/s2, r_max/q
    implication = separated_reader_action_implication(max(pivot, F(1, 10**12)), 20, 1, 1, 1, 1)
    margin = lambda ceiling: float(f'{float(F(word_s)*bias_rw/ceiling):.3e}')
    return {'word_s': word_s, 'rows_upper': rows, 'operations_lower': operations,
            'pivot_floor_lower': str(decimal_lower(pivot, 15)),
            'scalar_B_star_log10_lower': round(log_bstar, 1),
            'scalar_implication_finite': implication['B_star_finite'],
            'least_squares_action_upper': str(decimal_lower(ls_action, 3)),
            'least_squares_gyro_block_upper': str(decimal_lower(bias_block, 6)),
            'rho0_margin_upper_at_G0_action': margin(ls_action),
            'rho0_margin_upper_at_G0_gyro_block': margin(bias_block),
            'rho0_margin_upper_at_actual_floor_37': margin(r_max/actual_floor**2),
            'practical_linear_margin': False}


def certificate(nominal_mean_perp=F(2, 5), nominal_force_l1=F(6, 5)):
    """Premise values cover the carried worst (.348 m/s^2, 1.091) of
    aw-tracking-source-feasibility.json; they are supplied, not proved."""
    table = feasibility(nominal_mean_perp, nominal_force_l1)
    point = optimize_eps(OMEGA_INVARIANT, 1, window=16, separation=64, sigma_w=FRACTION,
                         m_perp=nominal_mean_perp, force_l1=nominal_force_l1,
                         rows_per_window=int(F(16)/H_MAX), start_offset=2)
    s2 = point[1]['six_column_floor_squared']
    return {
        'qualification': 'OU3_AGGREGATE_SIX_COLUMN_FLOOR_V1',
        'scope': 'injection-free nominal world array (A~=I, B~=int R_hat\'), real arithmetic; nominal window premises supplied, not proved',
        'magnetic_gap_from_service_s': '1',
        'magnetic_gap_reason': 'MAGNETIC SERVICE holds on every certified interval of length T_M=1 s, so each contains an applied informative row',
        'tube_lemma': {'constant': 'cos(Omega Delta)-Omega eps', 'requires': 'Omega Delta<pi/2',
                       'invariant_Omega_times_service_gap': str(OMEGA_INVARIANT),
                       'cos_lower_at_invariant_Omega': str(decimal_lower(cos_lower(OMEGA_INVARIANT), 12)),
                       'audit': tube_audit()},
        'theorem_G0_invariant_Omega_service_gap': {
            'Omega': str(OMEGA_INVARIANT), 'magnetic_gap_s': '1', 'window_s': '16',
            'window_separation_s': '64', 'eps_s': str(point[0]),
            'nominal_transverse_mean_premise_mps2': str(nominal_mean_perp),
            'nominal_force_L1_premise': str(nominal_force_l1),
            'tube_constant': str(decimal_lower(point[1]['tube_constant'], 12)),
            'q_I': str(point[1]['q_I']),
            'six_column_floor_squared_lower': str(decimal_lower(s2, 15)),
            'six_pivots_lower_for_m_rows': 's/sqrt(m)',
        },
        'feasibility': table,
        'synthetic_falsification_audit': synthetic_audit(),
        'downstream': downstream(s2, point[1]['q_I']),
        'injection_transport_included': False,
        'nominal_window_premises_proved': False,
        'source_uniform_six_column_floor': False,
        'rho0_certified': False,
        'theorem_closed': False,
    }


if __name__ == '__main__':
    import json
    print(json.dumps(certificate(), indent=2, sort_keys=True))
