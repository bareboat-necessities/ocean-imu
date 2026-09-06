"""Exact-rational operation regressions, not SEA3 words or theorem evidence.

The fixtures deliberately test identities with dense H18/A21 covariances.
They do not establish source reachability, recurrent-word contraction or P4.
"""
from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "stability"))
import ou3_p4_complete_sea3_accelerometer_operation_coordinate as ACC
import ou3_p4_complete_sea3_correction_information_bound as CORR
import ou3_p4_exact_reset_transport as RESET
from ou3_interval import Interval, matrix_point


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def tr(A):
    return list(map(list, zip(*A)))


def add(A, B):
    return [[a+b for a, b in zip(ar, br)] for ar, br in zip(A, B)]


def scale(A, x):
    return [[x*a for a in row] for row in A]


def sub(A, B):
    return add(A, scale(B, F(-1)))


def mul(A, B):
    return [[sum((a*b for a, b in zip(row, col)), F(0)) for col in tr(B)] for row in A]


def inv(A):
    n = len(A)
    a = [list(row)+e for row, e in zip(A, eye(n))]
    for j in range(n):
        k = next(i for i in range(j, n) if a[i][j])
        a[j], a[k] = a[k], a[j]
        pivot = a[j][j]
        a[j] = [x/pivot for x in a[j]]
        for i in range(n):
            if i != j:
                factor = a[i][j]
                a[i] = [x-factor*y for x, y in zip(a[i], a[j])]
    return [row[n:] for row in a]


def vec(*x):
    return [[F(v)] for v in x]


def skew(c):
    x, y, z = (row[0] for row in c)
    return [[F(0), -z, y], [z, F(0), -x], [-y, x, F(0)]]


def cayley(c):
    C = skew(c)
    den = F(1)+sum((row[0]**2 for row in c), F(0))/4
    return add(eye(3), scale(add(C, scale(mul(C, C), F(1, 2))), 1/den))


def shift(c, R, da, f):
    E = cayley(c)
    Q = mul(mul(tr(R), E), R)
    e0 = mul(tr(R), mul(sub(sub(E, eye(3)), skew(c)), f))
    mixed = mul(sub(Q, eye(3)), da)
    return add(mixed, e0), e0, mixed, E, Q


def energy(x, precision):
    return mul(mul(tr(x), precision), x)[0][0]


def maxabs(A):
    return max(abs(x) for row in A for x in row)


def iv(A):
    # Each rational is enclosed, rather than rounded and treated as exact.
    return [[Interval.outward_bounds(float(x), float(x)) for x in row] for row in A]


def fixture(n):
    c = vec(F(1, 4), F(-1, 5), F(1, 10))
    R = cayley(vec(F(1, 7), F(1, 9), F(-1, 11)))
    ahat, gravity = vec(F(1, 10), F(-1, 8), F(1, 6)), vec(0, 0, F(196133, 20000))
    da = vec(F(1, 20), F(-1, 30), F(1, 40))
    f = mul(R, sub(ahat, gravity))
    epsilon, e0, mixed, E, Q = shift(c, R, da, f)
    z = [[F(i+1, 10000)] for i in range(n)]
    z[:3], z[15:18] = c, da
    phi = [row[:] for row in z]
    phi[15:18] = add(da, epsilon)
    H = [[F(0) for _ in range(n)] for _ in range(3)]
    for i in range(3):
        H[i][:3] = scale(skew(f), F(-1))[i]
        H[i][15:18] = R[i]
        if n == 21:
            H[i][18+i] = F(1)
    diagonal = [F(1, 10**8)]*3 + [F(1, 10000)]*12 + [F(1, 100)]*3
    if n == 21:
        diagonal += [F(1, 10**8)]*3
    v = [[F(i+1, 10**6)] for i in range(n)]
    P = add([[diagonal[i] if i == j else F(0) for j in range(n)] for i in range(n)], mul(v, tr(v)))
    noise = scale(eye(3), F(1, 100))
    S = add(mul(mul(H, P), tr(H)), noise)
    K = mul(mul(P, tr(H)), inv(S))
    y = mul(H, phi)
    d = mul(K, y)
    dt = d[:3]
    theta2 = sum((row[0]**2 for row in dt), F(0))
    if not theta2 < F(1, 10000):
        raise AssertionError("fixture left shipping small-angle quaternion branch")
    qw = 1-theta2/8+theta2**2/384
    qk = F(1, 2)-theta2/48+theta2**2/3840
    a = scale(dt, 2*qk/qw)
    correction = cayley(a)
    numer = add(sub(c, a), scale(mul(skew(a), c), F(1, 2)))
    denom = F(1)+mul(tr(a), c)[0][0]/4
    cp = scale(numer, 1/denom)
    Rp = mul(correction, R)
    dap = sub(da, d[15:18])
    fp = mul(Rp, sub(add(ahat, d[15:18]), gravity))
    epsp, e0p, mixedp, _, _ = shift(cp, Rp, dap, fp)
    t = sub(z, d)
    zp = [row[:] for row in t]
    zp[:3] = cp
    phip = [row[:] for row in zp]
    phip[15:18] = add(dap, epsp)
    G = eye(n)
    G3 = add(eye(3), scale(skew(dt), F(1, 2)))
    Gi = eye(n)
    G3i = inv(G3)
    for i in range(3):
        G[i][:3], Gi[i][:3] = G3[i], G3i[i]
    Aphi = sub(phi, d)
    rho = sub(zp, mul(G, t))
    xi = [row[:] for row in rho]
    xi[15:18] = add(xi[15:18], sub(epsp, epsilon))
    old_xi = [row[:] for row in rho]
    old_xi[15:18] = add(old_xi[15:18], sub(e0p, e0))
    Pinv = inv(P)
    PJinv = add(Pinv, mul(mul(tr(H), inv(noise)), H))
    Prinv = mul(mul(tr(Gi), PJinv), Gi)
    J = energy(y, inv(S))
    b = mul(Gi, xi)
    change = energy(phip, Prinv)-energy(phi, Pinv)
    ledger = -J+2*mul(mul(tr(Aphi), PJinv), b)[0][0]+energy(b, PJinv)
    T = eye(n)
    Ti = eye(n)
    for i in range(3):
        T[15+i][15:18], Ti[15+i][15:18] = Q[i], tr(Q)[i]
    return locals()


class CoordinateTransportAlgebraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = {n: fixture(n) for n in (18, 21)}

    def test_original_shipping_H_times_Phi_is_exact_physical_residual(self):
        for n, a in self.cases.items():
            with self.subTest(n=n):
                physical = add(mul(sub(a['E'], eye(3)), a['f']), mul(mul(a['E'], a['R']), a['da']))
                if n == 21:
                    physical = add(physical, a['z'][18:21])
                self.assertEqual(physical, a['y'])
                self.assertGreater(maxabs(a['mixed']), F(1, 1000))

    def test_frozen_congruence_does_not_remove_shipping_mixed_remainder(self):
        for n, a in self.cases.items():
            with self.subTest(n=n):
                Hu, zu = mul(a['H'], a['Ti']), mul(a['T'], a['z'])
                self.assertEqual(mul(Hu, zu), mul(a['H'], a['z']))
                eta = sub(a['y'], mul(Hu, zu))
                pure = mul(a['R'], a['e0'])
                self.assertGreater(maxabs(sub(eta, pure)), F(1, 1000))
                self.assertNotEqual(energy(a['phi'], a['Pinv']), energy(a['z'], a['Pinv']))

    def test_pure_rotation_geometry_does_not_claim_full_eta_orthogonality(self):
        for a in self.cases.values():
            yr = mul(sub(a['E'], eye(3)), a['f'])
            h = mul(skew(a['c']), a['f'])
            eta = sub(yr, h)
            q2 = sum(row[0]**2 for row in a['c'])
            self.assertEqual(mul(tr(yr), eta)[0][0], 0)
            self.assertEqual(energy(eta, eye(3)), q2/(4+q2)*energy(h, eye(3)))
            self.assertEqual(energy(yr, eye(3)), 4/(4+q2)*energy(h, eye(3)))
            full_eta = sub(a['y'], mul(a['H'], a['z']))
            self.assertNotEqual(mul(tr(a['y']), full_eta)[0][0], 0)

    def test_gravity_only_normal_form_keeps_full_nonlinear_shift(self):
        for n, a in self.cases.items():
            B = scale(mul(tr(a['R']), skew(mul(a['R'], a['ahat']))), F(-1))
            TB, TBi = eye(n), eye(n)
            for i in range(3):
                TB[15+i][:3] = B[i]
                TBi[15+i][:3] = scale(B, F(-1))[i]
            HB = mul(a['H'], TBi)
            self.assertEqual([row[:3] for row in HB], skew(mul(a['R'], a['gravity'])))
            self.assertEqual(mul(HB, mul(TB, a['phi'])), a['y'])
            wlin = add(a['da'], mul(B, a['c']))
            wexact = add(wlin, a['epsilon'])
            self.assertEqual(wexact, mul(TB, a['phi'])[15:18])
            equivalent = sub(sub(mul(a['Q'], add(a['ahat'], a['da'])), a['ahat']),
                             mul(tr(a['R']), mul(sub(sub(a['E'], eye(3)), skew(a['c'])),
                                                   mul(a['R'], a['gravity']))))
            self.assertEqual(wexact, equivalent)
            PB = mul(mul(TB, a['P']), tr(TB))
            self.assertEqual(mul(mul(HB, PB), tr(HB)), mul(mul(a['H'], a['P']), tr(a['H'])))

    def test_interval_evaluator_encloses_exact_full_shift(self):
        for n, a in self.cases.items():
            with self.subTest(n=n):
                d = ACC.evaluate_operation_coordinate(iv(a['c']), iv(a['R']), iv(a['da']), iv(a['f']))
                for key, expected in (("epsilon_aw", a['epsilon']), ("mixed_aw_shift", a['mixed']),
                                      ("Phi_aw", a['phi'][15:18])):
                    for got, want in zip(d[key], expected):
                        self.assertLessEqual(F(got[0].lo), want[0])
                        self.assertGreaterEqual(F(got[0].hi), want[0])
                for physical, linear in zip(d['physical_residual_without_ba'], d['H_Phi_without_ba']):
                    residual = physical[0]-linear[0]
                    self.assertLessEqual(residual.lo, 0)
                    self.assertGreaterEqual(residual.hi, 0)

    def test_exact_full_shift_transport_and_signed_energy_ledger(self):
        for n, a in self.cases.items():
            with self.subTest(n=n):
                self.assertEqual(a['phip'], add(mul(a['G'], a['Aphi']), a['xi']))
                self.assertEqual(a['change'], a['ledger'])
                self.assertGreater(maxabs(sub(a['xi'], a['old_xi'])), F(1, 10000))
                expected = sub(mul(sub(mul(mul(tr(a['Rp']), cayley(a['cp'])), a['Rp']), a['Q']), a['da']),
                               mul(sub(mul(mul(tr(a['Rp']), cayley(a['cp'])), a['Rp']), eye(3)), a['d'][15:18]))
                self.assertEqual(sub(a['mixedp'], a['mixed']), expected)

    def test_shipping_left_injection_requires_right_inverse_error_composition(self):
        for a in self.cases.values():
            self.assertEqual(cayley(a['cp']), mul(a['E'], tr(a['correction'])))
            q = sum(float(row[0])**2 for row in a['c'])**0.5
            delta = sum(float(row[0])**2 for row in a['dt'])**0.5
            bound = RESET.reset_defect_bound(q*(1+1e-12), delta*(1+1e-12))
            self.assertLessEqual(sum(row[0]**2 for row in a['rho'][:3]), F(bound['reset_attitude_defect_norm_upper'])**2)

    def test_information_correlated_bounds_keep_full_H18_A21_P(self):
        for a in self.cases.values():
            for indices in ((0, 1, 2), (15, 16, 17)):
                upper = CORR.correction_norm_squared_upper(iv(a['P']), Interval.outward_bounds(float(a['J']), float(a['J'])), indices)
                actual = sum(a['d'][i][0]**2 for i in indices)
                self.assertLessEqual(actual, F(upper))

    def test_full_inverse_block_is_not_inverse_of_marginal(self):
        for n in (18, 21):
            P = eye(n)
            P[0][1] = P[1][0] = F(9, 10)
            upper = CORR.defect_precision_trace_upper(iv(P), (0,))
            self.assertGreaterEqual(F(upper), F(100, 19))
            self.assertAlmostEqual(100/19, upper, places=10)
            self.assertGreater(F(100, 19), 1/P[0][0])

    def test_covariance_reset_is_metric_isometry_not_covariance_reduction(self):
        G = add(eye(3), scale(skew(vec(F(1, 5), 0, 0)), F(1, 2)))
        reset = mul(G, tr(G))
        self.assertEqual(reset[1][1], F(101, 100))
        self.assertGreater(reset[1][1], 1)
        self.assertEqual(mul(mul(tr(G), inv(reset)), G), eye(3))

    def test_bad_covariance_and_vector_inputs_fail_closed(self):
        with self.assertRaises(ValueError):
            ACC.evaluate_operation_coordinate([], eye(3), [], [])
        with self.assertRaises(ValueError):
            CORR.defect_precision_trace_upper(matrix_point(eye(2)), (0,))
        with self.assertRaises(ValueError):
            CORR.defect_precision_trace_upper(matrix_point(eye(18)), (0, 0))
        P = eye(18)
        P[0][0] = F(-1)
        with self.assertRaises(RuntimeError):
            CORR.defect_precision_trace_upper(matrix_point(P), (0,))


if __name__ == '__main__':
    unittest.main()
