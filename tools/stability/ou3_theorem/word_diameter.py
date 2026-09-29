"""Word Riccati diameter: exact form of the A21 information-ratio contraction.

Role in V_(j+1) <= rho V_j + c_d |d|^2: bounds the linear factor rho_0 of one
complete regular A21 word, with its actual frozen coefficients, prior, gains,
resets and correlated sources. Proofs are in docs/ou3-corrected-word-proof.md
section 7. Every check below is exact rational algebra on supplied words; the
gyro-bias persistence cap is a source-uniform inequality for the default
profile. Nothing here certifies a source-uniform upper bound on the diameter.

Notation for a frozen word with unknown root x_0 and unit sources s:
  y = O x_0 + A s,   x_N = Phi x_0 + T s,   Sigma = A A'.
  J = O' Sigma^-1 O (root information in the data),
  A_info = root information in (y, x_N),
  Pi = Cov(x_N | y, x_0) = Ric_W(0)     (known root),
  P_diff = lim_t Ric_W(t I)             (diffuse root).
"""
from fractions import Fraction as F
from math import isqrt

from .lin_path_certificate import inverse
from .matrix_certificates import add, identity, is_psd, ldlt, matmul, transpose


def _m(a):
    return [[F(v) for v in row] for row in a]


def _zero(r, c):
    return [[F(0)]*c for _ in range(r)]


def _scale(a, s):
    return [[s*v for v in row] for row in a]


def _sub(a, b):
    return add(a, b, F(-1))


def _block(a, rows, cols):
    return [[a[i][j] for j in cols] for i in rows]


def full_root_design(n, events):
    """Chronological design of a frozen word with the whole root unknown.

    events: {'kind':'prediction','F','U'} (Q=U U'), {'kind':'correction','H','V'}
    (R=V V'), {'kind':'reset','G'}. Every prediction/noise factor is one block
    of unit-covariance source columns; nothing is dropped or split.
    """
    sources = sum(len(e['U'][0]) if e['kind'] == 'prediction' else
                  len(e['V'][0]) if e['kind'] == 'correction' else 0 for e in events)
    phi, t = identity(n), _zero(n, sources)
    rows_o, rows_a, col = [], [], 0
    for e in events:
        kind = e['kind']
        if kind == 'prediction':
            f, u = _m(e['F']), _m(e['U'])
            phi, t = matmul(f, phi), matmul(f, t)
            for i in range(n):
                t[i][col:col+len(u[0])] = u[i]
            col += len(u[0])
        elif kind == 'reset':
            g = _m(e['G'])
            phi, t = matmul(g, phi), matmul(g, t)
        elif kind == 'correction':
            h, v = _m(e['H']), _m(e['V'])
            rows_o.extend(matmul(h, phi))
            block = matmul(h, t)
            for i, row in enumerate(v):
                block[i][col:col+len(row)] = row
            rows_a.extend(block)
            col += len(v[0])
        else:
            raise ValueError('unsupported operation')
    return {'O': rows_o, 'A': rows_a, 'Phi': phi, 'T': t}


def terminal_pair(design):
    """Exact J, A_info, Pi, P_diff and the closed-loop map Phi_tilde.

    A_info = J + Phi_tilde' Pi^-1 Phi_tilde and P_diff = Pi + Phi_tilde J^-1 Phi_tilde',
    where Phi_tilde = Phi - T A' Sigma^-1 O is the root coefficient of
    E[x_N | y, x_0] and Pi = T (I - A' Sigma^-1 A) T'. Both identities are the
    Schur/Woodbury forms of the joint Gaussian (y, x_N); they are checked here.
    """
    o, a, phi, t = (design[k] for k in ('O', 'A', 'Phi', 'T'))
    sigma = matmul(a, transpose(a))
    ldlt(sigma)
    si = inverse(sigma)
    j = matmul(matmul(transpose(o), si), o)
    ldlt(j)
    ta = matmul(t, transpose(a))
    tilde = _sub(phi, matmul(matmul(ta, si), o))
    pi = _sub(matmul(t, transpose(t)), matmul(matmul(ta, si), transpose(ta)))
    ldlt(pi)
    pdiff = add(pi, matmul(matmul(tilde, inverse(j)), transpose(tilde)))
    stacked = o+phi
    cov = [row_a+row_b for row_a, row_b in zip(sigma, matmul(a, transpose(t)))]
    cov += [row_a+row_b for row_a, row_b in zip(matmul(t, transpose(a)), matmul(t, transpose(t)))]
    a_info = matmul(matmul(transpose(stacked), inverse(cov)), stacked)
    if a_info != add(j, matmul(matmul(transpose(tilde), inverse(pi)), tilde)):
        raise ArithmeticError('A = J + Phi_tilde\' Pi^-1 Phi_tilde failed')
    return {'J': j, 'A': a_info, 'Pi': pi, 'P_diff': pdiff, 'Phi_tilde': tilde}


def riccati(p, events):
    """Literal frozen-coefficient optimal-gain Riccati recursion (exact)."""
    p = _m(p)
    for e in events:
        kind = e['kind']
        if kind == 'prediction':
            f, u = _m(e['F']), _m(e['U'])
            p = add(matmul(matmul(f, p), transpose(f)), matmul(u, transpose(u)))
        elif kind == 'reset':
            g = _m(e['G'])
            p = matmul(matmul(g, p), transpose(g))
        else:
            h, v = _m(e['H']), _m(e['V'])
            s = add(matmul(matmul(h, p), transpose(h)), matmul(v, transpose(v)))
            k = matmul(matmul(p, transpose(h)), inverse(s))
            p = _sub(p, matmul(matmul(k, h), p))
    return p


def closed_loop(p, events):
    """Closed-loop root-to-end map M and end covariance of the optimal word."""
    p, m = _m(p), identity(len(p))
    for e in events:
        kind = e['kind']
        if kind == 'prediction':
            f, u = _m(e['F']), _m(e['U'])
            p = add(matmul(matmul(f, p), transpose(f)), matmul(u, transpose(u)))
            m = matmul(f, m)
        elif kind == 'reset':
            g = _m(e['G'])
            p, m = matmul(matmul(g, p), transpose(g)), matmul(g, m)
        else:
            h, v = _m(e['H']), _m(e['V'])
            s = add(matmul(matmul(h, p), transpose(h)), matmul(v, transpose(v)))
            k = matmul(matmul(p, transpose(h)), inverse(s))
            a = _sub(identity(len(p)), matmul(k, h))
            p = _sub(p, matmul(matmul(k, h), p))
            m = matmul(a, m)
    return m, p


def generalized_max_bracket(big, small, width=F(1, 10**9)):
    """Exact rational bracket lo < lambda_max(small^-1 big) <= hi by bisection.

    Every test is exact semidefinite elimination of kappa*small - big.
    """
    lo, hi = F(0), F(1)
    while not is_psd(_sub(_scale(small, hi), big)):
        lo, hi = hi, 2*hi
    while hi-lo > width*max(F(1), lo):
        mid = (lo+hi)/2
        if is_psd(_sub(_scale(small, mid), big)):
            hi = mid
        else:
            lo = mid
    return lo, hi


def sqrt_upper(q, bits=96):
    q = F(q)
    scale = 1 << bits
    num = q.numerator*scale*scale
    root = isqrt(-(-num//q.denominator))
    while F(root, scale)**2 < q:
        root += 1
    return F(root, scale)


def k0_bound_upper(kappa):
    """Rational upper bound of (sqrt(kappa)-1)/(sqrt(kappa)+1) = tanh(log(kappa)/4).

    The map is increasing in sqrt(kappa), so an upward square root is safe.
    """
    r = sqrt_upper(kappa)
    return (r-1)/(r+1)


def information_ratio_bound(kappa, k):
    """Closed form of sup_y [1/(1+y) - 1/((1+k)(1+kappa y))], kappa>1, k>=0.

    With c=1/(1+k): (sqrt(kappa)-sqrt(c))^2/(kappa-1) if c kappa>=1, else 1-c.
    The stationary point solves (1+kappa y)=sqrt(c kappa)(1+y); evaluate there.
    """
    kappa, k = float(kappa), max(float(k), 0.0)
    c = 1/(1+k)
    if c*kappa <= 1:
        return 1-c
    return (kappa**.5-c**.5)**2/(kappa-1)


def closed_form_rational_check():
    """Exact check of the closed form at rational kappa=c^-1 squares.

    At kappa=p^2, c=q^2 (p>1>=q, pq>1) the maximizer y*=(pq-1)/(p^2-pq) is
    rational; verify f(y*) equals the formula and dominates a rational grid.
    """
    out = []
    for p, q in ((F(3), F(1)), (F(5), F(1, 2)), (F(7, 2), F(2, 3)), (F(10), F(1, 3))):
        kappa, c = p*p, q*q
        ystar = (p*q-1)/(p*p-p*q)
        def f(y, kappa=kappa, c=c):
            return 1/(1+y)-c/(1+kappa*y)
        value = f(ystar)
        formula = (p-q)**2/(kappa-1)
        if value != formula:
            raise ArithmeticError('closed form of the information-ratio supremum failed')
        grid = [F(i, 64) for i in range(0, 64*8)]
        if any(f(y) > value for y in grid):
            raise ArithmeticError('stationary value is not the supremum on the grid')
        out.append({'kappa': str(kappa), 'c': str(c), 'argmax': str(ystar), 'sup': str(value)})
    return out


def sharp_scalar_certificate():
    """Scalar word where the k=0 bound is attained: kappa=2, sup rho = 3-2 sqrt2.

    Word: correct y=x+v (R=1), then predict x+=x+w (Q=1). Then J=1, A=2,
    Pi=1, P_diff=2 and rho(p)=p/((p+1)(2p+1)) for prior p. The exact identity
      (3-2 sqrt2)(2p^2+3p+1)-p = ((2-sqrt2) p-(sqrt2-1))^2
    shows sup_p rho(p)=(sqrt2-1)/(sqrt2+1), attained at p=1/sqrt2.
    """
    events = [{'kind': 'correction', 'H': [[1]], 'V': [[1]]},
              {'kind': 'prediction', 'F': [[1]], 'U': [[1]]}]
    pair = terminal_pair(full_root_design(1, events))
    if (pair['J'], pair['A'], pair['Pi'], pair['P_diff']) != ([[1]], [[2]], [[1]], [[2]]):
        raise ArithmeticError('scalar diameter word changed')
    for p in (F(1, 10), F(7, 10), F(1), F(3)):
        m, pend = closed_loop([[p]], events)
        rho = m[0][0]**2*p/pend[0][0]
        if rho != p/((p+1)*(2*p+1)) or rho > k0_bound_upper(2):
            raise ArithmeticError('scalar contraction exceeded its diameter bound')
    # Coefficients of the perfect square with s=sqrt2: compare in Q(sqrt2).
    lhs = {'quad': (6, -4), 'lin': (8, -6), 'const': (3, -2)}      # a+b*s
    rhs = {'quad': ((2, -1), (2, -1)), 'lin': ((2, -1), (1, -1)), 'const': ((1, -1), (1, -1))}

    def mul(x, y):
        return (x[0]*y[0]+2*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
    square = {'quad': mul(*rhs['quad']), 'lin': tuple(2*v for v in mul(*rhs['lin'])),
              'const': mul(*rhs['const'])}
    if square != lhs:
        raise ArithmeticError('sharpness identity failed')
    return {'kappa': '2', 'sup_rho': '3-2*sqrt(2)', 'attained_at_prior': '1/sqrt(2)',
            'k0_bound_sharp': True}


def fixture_events(steps=6):
    """Supplied five-state slow/fast word, exact rationals, not a shipping history.

    x=(theta, b, a, s, c): theta+=theta+h b, b+=b, a+=a/2, s+=s+h a, c+=(9/10)c.
    Rows: 'acc' theta+a+c, 'S' s every step; 'mag' theta every other step,
    followed by a literal-form reset theta<-(11/10) theta. Slow=(theta,b,c).
    """
    h = F(1, 2)
    f = [[1, h, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, F(1, 2), 0, 0],
         [0, 0, h, 1, 0], [0, 0, 0, 0, F(9, 10)]]
    u = [[F(1, 10), 0, 0, 0, 0], [0, F(1, 100), 0, 0, 0], [0, 0, F(1, 2), 0, 0],
         [0, 0, F(1, 10), F(1, 100), 0], [0, 0, 0, 0, F(1, 50)]]
    acc = [[1, 0, 1, 0, 1]]
    s_row = [[0, 0, 0, 1, 0]]
    mag = [[1, 0, 0, 0, 0]]
    g = identity(5)
    g[0][0] = F(11, 10)
    events = []
    for k in range(steps):
        events.append({'kind': 'prediction', 'F': f, 'U': u})
        events.append({'kind': 'correction', 'H': acc, 'V': [[F(1, 5)]]})
        events.append({'kind': 'correction', 'H': s_row, 'V': [[F(1)]]})
        if k % 2:
            events.append({'kind': 'correction', 'H': mag, 'V': [[F(1, 5)]]})
            events.append({'kind': 'reset', 'G': g})
    return events


def schur(m, keep, drop):
    a = _block(m, keep, keep)
    b = _block(m, keep, drop)
    return _sub(a, matmul(matmul(b, inverse(_block(m, drop, drop))), transpose(b)))


def diameter_certificate():
    """Theorem D on the supplied word: identities, bracket, contraction bound."""
    events = fixture_events()
    n = 5
    pair = terminal_pair(full_root_design(n, events))
    j, a_info, pi, pdiff = pair['J'], pair['A'], pair['Pi'], pair['P_diff']
    # Pi and P_diff are the literal Riccati recursions from 0 and t I, t->infinity.
    if riccati(_zero(n, n), events) != pi:
        raise ArithmeticError('Pi differs from the zero-prior Riccati recursion')
    previous = None
    for t in (F(10)**3, F(10)**6):
        ric = riccati(_scale(identity(n), t), events)
        woodbury = add(pi, matmul(matmul(pair['Phi_tilde'], inverse(add(_scale(identity(n), 1/t), j))),
                                  transpose(pair['Phi_tilde'])))
        if ric != woodbury or not is_psd(_sub(pdiff, ric)):
            raise ArithmeticError('finite-prior Riccati differs from its Woodbury form')
        if previous is not None and not is_psd(_sub(ric, previous)):
            raise ArithmeticError('Riccati image is not increasing in the prior')
        previous = ric
    lo_j, hi_j = generalized_max_bracket(a_info, j)
    lo_p, hi_p = generalized_max_bracket(pdiff, pi)
    hi = max(hi_j, hi_p)
    lo = min(lo_j, lo_p)
    if not (is_psd(_sub(_scale(j, hi), a_info)) and is_psd(_sub(_scale(pi, hi), pdiff))
            and not is_psd(_sub(_scale(j, lo), a_info)) and not is_psd(_sub(_scale(pi, lo), pdiff))):
        raise ArithmeticError('information and covariance diameters differ')
    bound = k0_bound_upper(hi)
    priors = [identity(n), [[F(10), 0, 0, 0, 0], [0, F(1, 10), 0, 0, 0], [0, 0, F(3), 0, 0],
                            [0, 0, 0, F(1, 2), 0], [0, 0, 0, 0, F(1, 7)]],
              [[F(1, 100), 0, 0, 0, 0], [0, F(1, 1000), 0, 0, 0], [0, 0, F(1, 50), 0, 0],
               [0, 0, 0, F(1, 20), 0], [0, 0, 0, 0, F(1, 1000)]]]
    rhos = []
    for p0 in priors:
        m, pend = closed_loop(p0, events)
        energy = matmul(matmul(transpose(m), inverse(pend)), m)
        # rho <= bound  <=>  M' P_end^-1 M <= bound P_0^-1 (exact PSD test).
        if not is_psd(_sub(_scale(inverse(p0), bound), energy)):
            raise ArithmeticError('word contraction exceeds tanh(log(kappa)/4)')
        if not (is_psd(_sub(pend, pi)) and is_psd(_sub(pdiff, pend))):
            raise ArithmeticError('actual end covariance leaves [Pi, P_diff]')
        rlo, rhi = generalized_max_bracket(energy, inverse(p0), F(1, 10**6))
        rhos.append(str(rhi))
    return {'kappa_bracket': [str(lo), str(hi)], 'k0_bound_upper': str(bound),
            'exact_rho_upper_for_supplied_priors': rhos,
            'identities': ['A = J + Phi_tilde\' Pi^-1 Phi_tilde', 'P_diff = Pi + Phi_tilde J^-1 Phi_tilde\'',
                           'Pi = Ric_W(0)', 'Ric_W(t I) = Pi + Phi_tilde (I/t + J)^-1 Phi_tilde\''],
            'lambda_max(J^-1 A) = lambda_max(Pi^-1 P_diff)': True,
            'every_prior_end_covariance_in_[Pi,P_diff]': True,
            'scope': 'supplied rational five-state word; exact algebra'}


def composition_certificate():
    """kappa(W2 o W1) <= kappa(W2): [Pi_W, P_diff_W] lies inside [Pi_W2, P_diff_W2]."""
    events = fixture_events()
    split = next(i for i, e in enumerate(events) if e['kind'] == 'prediction' and i >= len(events)//2)
    whole = terminal_pair(full_root_design(5, events))
    tail = terminal_pair(full_root_design(5, events[split:]))
    if not (is_psd(_sub(whole['Pi'], tail['Pi'])) and is_psd(_sub(tail['P_diff'], whole['P_diff']))):
        raise ArithmeticError('word image is not nested in its tail image')
    _, hi_whole = generalized_max_bracket(whole['P_diff'], whole['Pi'])
    lo_tail, _ = generalized_max_bracket(tail['P_diff'], tail['Pi'])
    if not hi_whole <= lo_tail*(1+F(1, 10**6)):
        raise ArithmeticError('diameter increased under prefixing')
    return {'kappa_word_upper': str(hi_whole), 'kappa_tail_lower': str(lo_tail),
            'nested_images': True, 'rho_multiplicative': 'rho_W <= rho_W1 rho_W2 on every prior'}


def factorization_certificate():
    """kappa <= kappa(slow marginal) kappa(fast | slow) and the reverse order.

    With J <= A_s <= A (A_s adds the slow terminal state), lambda_max(J^-1 A)
    <= lambda_max(A_s^-1 A) lambda_max(J^-1 A_s); each factor is the diameter of
    a terminal marginal or Schur complement of [Pi, P_diff].
    """
    events = fixture_events()
    pair = terminal_pair(full_root_design(5, events))
    pi, pdiff = pair['Pi'], pair['P_diff']
    slow, fast = [0, 1, 4], [2, 3]
    _, whole = generalized_max_bracket(pdiff, pi)
    out = {}
    for name, first, second in (('slow_then_fast', slow, fast), ('fast_then_slow', fast, slow)):
        _, marg = generalized_max_bracket(_block(pdiff, first, first), _block(pi, first, first))
        _, cond = generalized_max_bracket(schur(pdiff, second, first), schur(pi, second, first))
        if not whole <= marg*cond*(1+F(1, 10**6)):
            raise ArithmeticError('diameter factorization failed')
        out[name] = {'marginal_upper': str(marg), 'conditional_upper': str(cond), 'product': str(marg*cond)}
    out['kappa_upper'] = str(whole)
    return out


def compression_reader_certificate():
    """Dual certificates: compression floors on J and reader ceilings on P_diff.

    For any L: J >= (L O)'(L Sigma L')^-1 (L O) (information of a function of
    the data). For any reader with L O = Phi: P_diff <= (T - L A)(T - L A)',
    with equality at the minimum-action reader.
    """
    events = fixture_events()
    design = full_root_design(5, events)
    pair = terminal_pair(design)
    o, a, phi, t = (design[k] for k in ('O', 'A', 'Phi', 'T'))
    sigma = matmul(a, transpose(a))
    # Compression: pairwise sums of consecutive rows (a coarsened data record).
    rows = len(o)
    comp = [[F(1) if j in (2*i, 2*i+1) else F(0) for j in range(rows)] for i in range(rows//2)]
    lo_, ls = matmul(comp, o), matmul(matmul(comp, sigma), transpose(comp))
    floor = matmul(matmul(transpose(lo_), inverse(ls)), lo_)
    if not is_psd(_sub(pair['J'], floor)):
        raise ArithmeticError('compression floor exceeds J')
    # Minimum-action full-root reader and a perturbed feasible reader.
    si = inverse(sigma)
    info_inv = inverse(pair['J'])
    ta = matmul(t, transpose(a))
    reader = add(matmul(ta, si), matmul(matmul(pair['Phi_tilde'], info_inv), matmul(transpose(o), si)))
    if matmul(reader, o) != phi:
        raise ArithmeticError('minimum-action reader does not cancel the root')
    action = lambda l: (lambda z: matmul(z, transpose(z)))(_sub(t, matmul(l, a)))
    if action(reader) != pair['P_diff']:
        raise ArithmeticError('minimum reader action differs from P_diff')
    # Perturb within the null space of O': add D (I - O J^-1 O' Sigma^-1).
    proj = _sub(identity(rows), matmul(matmul(o, info_inv), matmul(transpose(o), si)))
    d = [[F((i+2*j) % 3-1, 7) for j in range(rows)] for i in range(5)]
    other = add(reader, matmul(d, proj))
    if matmul(other, o) != phi or not is_psd(_sub(action(other), pair['P_diff'])):
        raise ArithmeticError('feasible reader action is not above P_diff')
    return {'compression_floor_below_J': True, 'minimum_reader_action_equals_P_diff': True,
            'feasible_reader_action_above_P_diff': True}


def gyro_persistence_cap(word_s, sigma_g=F('0.00135'), bias_density=F('1e-10'),
                         h_max=F('0.006'), omega_max=None):
    """Source-uniform lower bound on the diameter of every regular A21 word.

    For a root gyro-bias perturbation beta, the data-only mimic keeps b=0 and
    feeds the literal attitude noise w_theta=B_step beta at each prediction.
    It reproduces every applied row, reset and nuisance state exactly, so
      u'J u <= sum_k |B_step beta|^2 / lambda_min(Q_AA / Q_bb)
            <= T (1+e)^2 |beta|^2 / (sigma_g^2 (1-delta)).
    Here |B_step|<=h(1+e) (degree-18 series, e<=2/22! for |omega|h<1) and the
    literal Schur complement Q_tt-Q_tb Q_bb^-1 Q_bt >= sigma_g^2 h(1-delta),
    delta=h^2 b0(1+e)^2/(4 sigma_g^2), because Simpson's I_BB is PSD and
    |int B|<=h^2(1+e)/2. The terminal bias alone carries information
    1/(b0 T) (identity bias transition, Q_bb=b0 h I), so A >= E_bg E_bg'/(b0 T).
    Hence kappa_W >= sigma_g^2 (1-delta)/((1+e)^2 b0 T^2) for every admitted
    history, gain, reset and nuisance coefficient (zero lever arm).
    """
    t = F(word_s)
    e = F(2, 1124000727777607680000)          # 2/22!
    delta = h_max**2*bias_density*(1+e)**2/(4*sigma_g**2)
    kappa = sigma_g**2*(1-delta)/((1+e)**2*bias_density*t*t)
    margin_cap = 1-k0_bound_upper(kappa) if kappa > 1 else F(1)
    return {'word_s': str(t), 'kappa_lower': kappa, 'delta': delta,
            'k0_margin_cap_upper': margin_cap,
            'meaning': 'no k=0 information-ratio certificate on a T-second word can certify a margin above 2/(1+sqrt(kappa_lower))'}


def kernel_fixture_events(steps=6):
    """Supplied word with an exact tilt/bias-like kernel (no field row, no decay).

    x=(theta, b, a, s, c) as in fixture_events, but c+=c and no 'mag' row or
    reset: the 'acc' row sees theta+c, so nu=(1,0,0,0,-1) is invisible and J is
    singular along nu. Scaled-down analogue of the quiet A21 tilt/BA kernel.
    """
    events = []
    for e in fixture_events(steps):
        if e['kind'] == 'prediction':
            f = [row[:] for row in e['F']]
            f[4][4] = F(1)
            events.append({'kind': 'prediction', 'F': f, 'U': e['U']})
        elif e['kind'] == 'correction' and e['H'] != [[1, 0, 0, 0, 0]]:
            events.append(e)
    return events


def kernel_bounded_diameter(design, pair, nu, mu):
    """P_nu = Ric_W(root information mu nu nu') and kappa_nu = lambda_max(Pi^-1 P_nu).

    Adding the fictitious root observation nu'x_0 (precision mu) to the word
    leaves Pi unchanged and turns J, A into J+mu nu nu', A+mu nu nu'. The
    diameter identity for that word gives A - kappa_nu J <= (kappa_nu-1) mu nu nu'.
    """
    del design
    n = len(nu)
    col = [[F(v)] for v in nu]
    jn = add(pair['J_raw'] if 'J_raw' in pair else pair['J'], _scale(matmul(col, transpose(col)), F(mu)))
    ldlt(jn)
    return add(pair['Pi'], matmul(matmul(pair['Phi_tilde'], inverse(jn)), transpose(pair['Phi_tilde']))), jn, n


def terminal_pair_singular(design):
    """Pi, Phi_tilde, J, A_info without requiring J>0 (kernel words)."""
    o, a, phi, t = (design[k] for k in ('O', 'A', 'Phi', 'T'))
    sigma = matmul(a, transpose(a))
    ldlt(sigma)
    si = inverse(sigma)
    j = matmul(matmul(transpose(o), si), o)
    ta = matmul(t, transpose(a))
    tilde = _sub(phi, matmul(matmul(ta, si), o))
    pi = _sub(matmul(t, transpose(t)), matmul(matmul(ta, si), transpose(ta)))
    ldlt(pi)
    a_info = add(j, matmul(matmul(transpose(tilde), inverse(pi)), tilde))
    return {'J': j, 'A': a_info, 'Pi': pi, 'Phi_tilde': tilde}


def kernel_bounded_certificate():
    """Rank-structured information-ratio bound on a word with an exact kernel.

    Corollary: if A <= kappa J + lambda nu nu' then, for every root with
    nu'P_0 nu <= c, rho_W <= sup_y[1/(1+y)-1/((1+lambda c)(1+kappa y))]; the
    only covariance input is the scalar c. With lambda=(kappa_nu-1) mu the
    premise holds at kappa=kappa_nu(mu). The full diameter is infinite here.
    """
    events = kernel_fixture_events()
    n = 5
    design = full_root_design(n, events)
    pair = terminal_pair_singular(design)
    nu = [F(1), F(0), F(0), F(0), F(-1)]
    col = [[v] for v in nu]
    if any(any(v for v in row) for row in matmul(pair['J'], col)):
        raise ArithmeticError('supplied kernel is not invisible to the data')
    c = F(1, 4)
    mu = 1/c
    p_nu, jn, _ = kernel_bounded_diameter(design, pair, nu, mu)
    lo, hi = generalized_max_bracket(p_nu, pair['Pi'])
    an = add(pair['A'], _scale(matmul(col, transpose(col)), mu))
    if not is_psd(_sub(_scale(jn, hi), an)):
        raise ArithmeticError('kernel-bounded diameter does not cover A')
    lam = (hi-1)*mu
    premise = _sub(add(_scale(pair['J'], hi), _scale(matmul(col, transpose(col)), lam)), pair['A'])
    if not is_psd(premise):
        raise ArithmeticError('A <= kappa J + lambda nu nu\' failed')
    checked = []
    for scale in (F(1), F(1, 3), F(1, 50)):
        # Priors with nu'P_0 nu = scale*c and arbitrary other structure.
        # Large theta and c variances, strongly correlated so that the
        # invisible combination theta-c has variance exactly scale*c.
        p0 = [[F(0)]*n for _ in range(n)]
        for i, value in enumerate((F(3), F(1, 10), F(2), F(1, 2), F(3))):
            p0[i][i] = value
        p0[0][4] = p0[4][0] = F(3)-scale*c/2
        nu_var = matmul(matmul(transpose(col), p0), col)[0][0]
        if nu_var > c:
            raise ArithmeticError('supplied prior violates the scalar kernel ceiling')
        m, pend = closed_loop(p0, events)
        energy = matmul(matmul(transpose(m), inverse(pend)), m)
        k = lam*nu_var
        cc = 1/(1+k)
        # Rational upper bound of the closed form (sqrt(kappa)-sqrt(c))^2/(kappa-1).
        if cc*hi >= 1:
            bound = (sqrt_upper(hi)-(sqrt_upper(cc)-F(1, 1 << 90)))**2/(hi-1)
        else:
            bound = 1-cc
        if not is_psd(_sub(_scale(inverse(p0), bound), energy)):
            raise ArithmeticError('rank-structured information-ratio bound failed')
        checked.append({'nu_variance': str(nu_var), 'k': str(k), 'bound_upper': str(bound)})
    return {'kernel': ['1', '0', '0', '0', '-1'], 'J_nu_is_zero': True,
            'scalar_ceiling_c': str(c), 'kappa_nu_bracket': [str(lo), str(hi)],
            'lambda': str(lam), 'priors_checked': checked,
            'full_diameter_finite': False}


def scalar_kernel_ceiling(tilt_variance_ceiling, nu_ba_norm, ba_marginal=F(1, 1600)):
    """Block Cauchy--Schwarz ceiling on nu'P nu for nu=(theta_hat, nu_ba).

    With theta_hat a unit attitude vector, P_ba,ba <= ba_marginal I (proved in
    ou3-nuisance-upper-proof.md) and theta_hat'P_tt theta_hat <= tau:
      nu'P nu <= (sqrt(tau)+|nu_ba| sqrt(ba_marginal))^2,
    for arbitrary attitude/BA and all other cross covariance.
    """
    tau, nb, pb = F(tilt_variance_ceiling), F(nu_ba_norm), F(ba_marginal)
    return (sqrt_upper(tau)+nb*sqrt_upper(pb))**2


def schain_identity_certificate():
    """Exact S-chain cancellation of the neutral/AW root, AW noise and syncs.

    Literal per-axis transition (v,p,S,a): v+=v+phi_va a, p+=p+h v+phi_pa a,
    S+=S+h p+(h^2/2) v+phi_Sa a, a+=alpha a, plus noise (n_v,n_p,n_S,n_a) and
    PSD sync jumps on a. For c annihilating (1,t,t^2) on the S times and
      omega_k = sum_{t_j>=t_{k+1}} c_j (phi_Sa+phi_pa d+phi_va d^2/2), d=t_j-t_{k+1},
    one has exactly
      sum_k omega_k a_k - sum_j c_j S_j = -sum_k g_k'(n_v,n_p,n_S)_k,
      g_k = sum_{t_j>=t_{k+1}} c_j (d^2/2, d, 1),
    whatever the literal phi coefficients: the root (v,p,S,a), every AW noise
    and every sync jump cancel. Checked here as exact linear forms.
    """
    steps = [F(1, 200), F(1, 250), F(3, 500), F(1, 200), F(1, 250), F(1, 200), F(3, 500), F(1, 200)]
    tau = [F(1, 3), F(1, 3), F(2), F(2), F(1, 50), F(1, 50), F(12), F(12)]
    # Literal-form coefficients need not be exact integrals; use rational stand-ins.
    phis = [(h*(1-h/(2*t)), h*h/2*(1-h/(3*t)), h**3/6*(1-h/(4*t)), 1-h/t+h*h/(2*t*t))
            for h, t in zip(steps, tau)]
    # Variables: v0,p0,S0,a0, then (nv,np,nS,na) per step, then one sync per step.
    nvars = 4+4*len(steps)+len(steps)
    def var(i):
        row = [F(0)]*nvars
        row[i] = F(1)
        return row
    z = [var(0), var(1), var(2), var(3)]
    times, states = [F(0)], [z]
    for k, (h, (pva, ppa, psa, alpha)) in enumerate(zip(steps, phis)):
        v, p, s, a = z
        nb = 4+4*k
        v_next = [x+pva*y+w for x, y, w in zip(v, a, var(nb))]
        p_next = [x+h*y+ppa*w+u for x, y, w, u in zip(p, v, a, var(nb+1))]
        s_next = [x+h*y+h*h/2*w+psa*u+q for x, y, w, u, q in zip(s, p, v, a, var(nb+2))]
        a_next = [alpha*x+y+w for x, y, w in zip(a, var(nb+3), var(4+4*len(steps)+k))]
        z = [v_next, p_next, s_next, a_next]
        times.append(times[-1]+h)
        states.append(z)
    s_idx = [0, 2, 5, 8]
    tj = [times[j] for j in s_idx]
    # c annihilates 1,t,t^2 on the S times (third divided difference).
    c = []
    for j in range(4):
        prod = F(1)
        for i in range(4):
            if i != j:
                prod *= tj[j]-tj[i]
        c.append(6/prod)
    for power in range(3):
        if sum(cj*t**power for cj, t in zip(c, tj)):
            raise ArithmeticError('S weights do not annihilate quadratics')
    lhs = [F(0)]*nvars
    for k in range(len(steps)):
        omega = F(0)
        g = [F(0), F(0), F(0)]
        for cj, j in zip(c, s_idx):
            if j >= k+1:
                d = times[j]-times[k+1]
                pva, ppa, psa, _ = phis[k]
                omega += cj*(psa+ppa*d+pva*d*d/2)
                g = [g[0]+cj*d*d/2, g[1]+cj*d, g[2]+cj]
        lhs = [x+omega*y for x, y in zip(lhs, states[k][3])]
        nb = 4+4*k
        for comp in range(3):
            lhs[nb+comp] += g[comp]
    for cj, j in zip(c, s_idx):
        lhs = [x-cj*y for x, y in zip(lhs, states[j][2])]
    if any(lhs):
        raise ArithmeticError('S-chain cancellation identity failed')
    return {'steps': len(steps), 'S_times': [str(t) for t in tj], 'weights': [str(v) for v in c],
            'root_cancelled': True, 'aw_noise_and_sync_cancelled': True,
            'residual': '-sum_k g_k\'(n_v,n_p,n_S)_k exactly'}


def schain_residual_ceiling(g_norms_sq, steps, q_max=F(1600), defect=F(0)):
    """Variance ceiling of the S-chain residual for independent step sources.

    Var(n_v)<=q h^3/3, Var(n_p)<=q h^5/20, Var(n_S)<=q h^7/252 for the OU
    chain with driving density q<=2 sigma^2/tau (undamped kernels dominate);
    the trace bounds lambda_max, and the source defect multiplies by 1+eps.
    """
    total = F(0)
    for g2, h in zip(g_norms_sq, steps):
        total += F(g2)*q_max*(h**3/3+h**5/20+h**7/252)
    return (1+F(defect))*total


def certificate():
    cap = [gyro_persistence_cap(t) for t in (16, 32, 64, 128, 256)]
    return {
        'qualification': 'OU3_WORD_RICCATI_DIAMETER_V1',
        'theorem': ('For a frozen word with J>0: lambda_max(J^-1 A) = lambda_max(Pi^-1 P_diff) =: kappa_W, '
                    'and rho_W <= tanh(log(kappa_W)/4) = (sqrt(kappa_W)-1)/(sqrt(kappa_W)+1) for every root covariance.'),
        'closed_form_information_ratio_bound': '(sqrt(kappa)-sqrt(c))^2/(kappa-1), c=1/(1+k), if c kappa>=1; else 1-c',
        'closed_form_rational_checks': closed_form_rational_check(),
        'sharp_scalar_example': sharp_scalar_certificate(),
        'diameter_identity': diameter_certificate(),
        'composition': composition_certificate(),
        'slow_fast_factorization': factorization_certificate(),
        'dual_certificates': compression_reader_certificate(),
        'kernel_bounded_corollary': kernel_bounded_certificate(),
        'scalar_kernel_ceiling': {
            'inequality': "nu'P nu <= (sqrt(tau_theta)+|nu_ba|/40)^2 from P_ba,ba <= I/1600",
            'example_tau_1e-3_nu_ba_9.80665': str(scalar_kernel_ceiling(F(1, 1000), F('9.80665'))),
            'BA_part_alone_nu_ba_9.80665': str(F('9.80665')**2/1600)},
        's_chain_cancellation': schain_identity_certificate(),
        'gyro_bias_persistence_cap': [
            {'word_s': c['word_s'], 'kappa_lower': str(c['kappa_lower']),
             'kappa_lower_decimal': float(c['kappa_lower']),
             'k0_margin_cap_upper_decimal': float(c['k0_margin_cap_upper'])} for c in cap],
        'root_covariance_bound_used': False,
        'source_uniform_diameter_upper_bound': False,
        'rho0_certified': False,
        'theorem_closed': False,
    }


if __name__ == '__main__':
    import json
    print(json.dumps(certificate(), indent=2, sort_keys=True))
