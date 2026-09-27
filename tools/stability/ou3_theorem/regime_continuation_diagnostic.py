"""Non-promoting high-precision feasibility of the new algebraic bounds.

Supplied coefficients/root covariances audit the constructions before exact
certification. They are not a replay, shipping reachability, or coverage of
every admitted source history. No finite ratio is used as a universal bound.
"""
import mpmath as mp


def diagnostic():
    with mp.workdps(80):
        def skew(v):
            x, y, z = map(mp.mpf, v)
            return mp.matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]])

        def transport(ops):
            a, b, eye = mp.eye(3), mp.zeros(3), mp.eye(3)
            beta = defect = elapsed = mp.mpf(0)
            for op in ops:
                if op[0] == 'predict':
                    h, w = mp.mpf(op[1]), mp.matrix(list(map(mp.mpf, op[2])))
                    wh, angle = h * skew(w), h * mp.norm(w)
                    r, d, power = mp.expm(-wh), h * eye, eye.copy()
                    for k in range(1, 50):
                        power = power * (-wh)
                        d += h * power / mp.factorial(k + 1)
                    a, b = r * a, r * b + d
                    defect += h * beta + h * (angle/2 + angle**2/3)
                    elapsed += h
                    beta = min(2, beta + angle + angle**2/2)
                else:
                    d = mp.matrix(list(map(mp.mpf, op[1])))
                    g = eye + skew(d)/2
                    a, b = g * a, g * b
                    beta = min(2, beta + min(1, mp.norm(d)/2))
            return {'gyro_floor': mp.nstr(max(0, elapsed - defect), 65),
                    'gyro_singular_value': mp.nstr(min(mp.svd(b, compute_uv=False)), 65),
                    'inverse_A_norm': mp.nstr(max(mp.svd(a**-1, compute_uv=False)), 65)}

        ops = [('predict', '.005', ['.2', '-.1', '.05']),
               ('reset', ['.001', '-.002', '.0005'])] * 8
        ordinary = transport(ops)
        terminal = transport([('predict', '.005', ['0', '0', '0']),
                              ('reset', ['4', '0', '0'])])
        t, g, field = mp.mpf('.04'), mp.mpf('9.80665'), mp.mpf(75)
        qg, qb = mp.mpf('.00135')**2, mp.mpf('1e-10')
        f, h, q = mp.eye(12), mp.zeros(3, 12), mp.zeros(12)
        r = mp.diag([mp.mpf('.1')/g**2]*2 + [mp.mpf('.64')/field**2])
        for j in range(3):
            f[j, j+3] = t
            f[j+6, j+6], f[j+9, j+9] = mp.exp(-t/6), mp.exp(-t/5000)
            h[j, j] = 1
            q[j, j] = qg*t + qb*t**3/3
            q[j, j+3] = q[j+3, j] = qb*t*t/2
            q[j+3, j+3] = qb*t
            q[j+6, j+6] = mp.mpf('4.84')*(1-f[j+6, j+6]**2) + 16
            q[j+9, j+9] = (1-f[j+9, j+9]**2)/1600
        h[0, 7] = h[0, 10] = 1/g
        h[1, 6] = h[1, 9] = -1/g
        end, l0, l1 = mp.eye(12)[:6, :], mp.zeros(6, 3), mp.zeros(6, 3)
        for j in range(3):
            l0[j+3, j], l1[j, j], l1[j+3, j] = -1/t, 1, 1/t
        root, fresh = end*f-l0*h-l1*h*f, end-l1*h
        eta = [((156+mp.mpf(1)/40+mp.mpf(1)/3)/9)**2]*2 + [mp.mpf(1)/81]
        diagonal = [2*x for x in eta] + [2*(4*x/t**2+mp.mpf('2e-6')/t+mp.mpf('1e-9')*t/3) for x in eta]
        white = mp.diag([1/mp.sqrt(x) for x in diagonal])
        ratios = []
        for scale in (mp.mpf(1), mp.mpf('1e12')):
            factor = mp.diag([scale]*6+[156]*3+[mp.mpf(1)/40]*3)
            for j in range(3):
                factor[j, j+6], factor[j+3, j+9] = scale/2, scale/3
            p = factor * factor.T
            action = root*p*root.T + fresh*q*fresh.T + l0*r*l0.T + l1*r*l1.T
            ratios.append(mp.nstr(max(mp.eigsy(white*action*white, eigvals_only=True)), 65))
        force, field_vector = mp.matrix([0, 0, 9]), mp.matrix([10, 0, 0])
        group_reset = ((mp.eye(3)+skew([0, '.2', 0])/2) *
                       (mp.eye(3)+skew(['.1', 0, 0])/2))
        pulled = group_reset**-1 * field_vector
        cosine = abs((force.T*pulled)[0]) / (mp.norm(force)*mp.norm(pulled))
        ca, cm = skew(force), skew(field_vector)*group_reset
        group = mp.matrix([[ca[i,j] for j in range(3)] for i in range(3)] +
                          [[cm[i,j] for j in range(3)] for i in range(3)])
        geometry = {'conditional_floor': mp.sqrt(81*(1-cosine)),
                    'actual_singular_value': min(mp.svd(group, compute_uv=False)),
                    'reset_pulled_field_cosine': cosine}
        return {'qualification': 'OU3_REGIME_CONTINUATION_FEASIBILITY_V1',
                'decimal_digits': 80, 'role': 'supplied algebra feasibility only',
                'noncommuting_transport': ordinary, 'terminal_large_reset': terminal,
                'quiet_action_upper_comparison_ratios': ratios,
                'quiet_reader_AG_root_residual': mp.nstr(mp.norm(root[:, :6]), 15),
                'same_cell_geometry': {k: mp.nstr(v, 65) for k, v in geometry.items()},
                'shipping_reachability_certified': False,
                'source_uniform_verified': False, 'theorem_closed': False}
