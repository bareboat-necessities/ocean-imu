"""Universal coefficient proofs and exact graph regressions, never source admission."""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.stability.ou3_alt_contraction import finite_measurement_graph as G


class Poly:
    """Tiny exact polynomial ring for all-variable, non-sampled identity checks."""
    N = 8

    def __init__(self, value=0):
        if isinstance(value, Poly):
            self.terms = value.terms.copy()
        elif isinstance(value, dict):
            self.terms = {m: F(c) for m, c in value.items() if c}
        else:
            self.terms = {(0,)*self.N: F(value)} if value else {}

    @classmethod
    def var(cls, i):
        e = [0]*cls.N
        e[i] = 1
        return cls({tuple(e): F(1)})

    def __add__(self, other):
        terms = defaultdict(F, self.terms)
        for m, c in Poly(other).terms.items():
            terms[m] += c
        return Poly(dict(terms))

    __radd__ = __add__

    def __neg__(self):
        return Poly({m: -c for m, c in self.terms.items()})

    def __sub__(self, other):
        return self+-Poly(other)

    def __rsub__(self, other):
        return Poly(other)+-self

    def __mul__(self, other):
        out = defaultdict(F)
        for m, c in self.terms.items():
            for n, d in Poly(other).terms.items():
                out[tuple(a+b for a, b in zip(m, n))] += c*d
        return Poly(dict(out))

    __rmul__ = __mul__

    def __eq__(self, other):
        return self.terms == Poly(other).terms


def solve3(S, r):
    a = [row[:]+[x] for row, x in zip(S, r)]
    for j in range(3):
        pivot = next(i for i in range(j, 3) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        v = a[j][j]
        a[j] = [x/v for x in a[j]]
        for i in range(3):
            if i != j:
                v = a[i][j]
                a[i] = [x-v*y for x, y in zip(a[i], a[j])]
    return [row[3] for row in a]


def dense_covariance():
    u = [F((i % 4)-2, 20) for i in range(21)]
    return G.plus(G.eye(21), [[x*y for y in u] for x in u])


class FiniteMeasurementTests(unittest.TestCase):
    def test_cayley_secant_identity_all_polynomial_coefficients(self):
        c = [Poly.var(i) for i in range(3)]
        U, D, En = G.cayley_polynomials(c)
        C = G.skew(c)
        self.assertEqual(G.mm(G.plus(G.scaled(G.eye(3), 2), C, -1), U), G.scaled(G.eye(3), 2*D))
        self.assertEqual(G.mm(U, C), G.plus(En, G.scaled(G.eye(3), D), -1))
        # E is a rotation on the real Cayley chart, not a linear approximation.
        self.assertEqual(G.mm(G.transpose(En), En), G.scaled(G.eye(3), D*D))

    def test_finite_reset_identity_all_eight_indeterminates(self):
        c, d = [Poly.var(i) for i in range(3)], [Poly.var(i) for i in range(3, 6)]
        w, k = Poly.var(6), Poly.var(7)
        W, twice_V, Lnum = G.right_reset_polynomials(c, d, w, k)
        gap = [twice_V[i]-W*c[i]+G.mv(Lnum, d)[i] for i in range(3)]
        self.assertEqual(gap, [Poly(0)]*3)
        # Includes BOTH quaternion branches because w,k were not specialized.

    def test_H18_full_latent_bias_covariance_enters_S_not_N(self):
        P = G.eye(21)
        latent = F(1, 62500)  # (0.004)^2 in exact real arithmetic.
        for i in range(18, 21):
            P[i][i] = latent
        H, _ = G.residual_factors('accelerometer', [0]*3, f_hat=[0, 0, -10], R_hat=G.eye(3))
        N, S = G.measurement_operands(P, G.eye(3), H, active_bias=False)
        H18 = [row[:18] for row in H]
        missing = G.plus(G.mm(G.mm(H18, G.eye(18)), G.transpose(H18)), G.eye(3))
        self.assertEqual(G.plus(S, missing, -1), G.scaled(G.eye(3), latent))
        self.assertEqual(N[18:21], G.zeros(3, 3))
        self.assertEqual(G.measurement_operands(P, G.eye(3), H, active_bias=True)[1], S)

    def test_masks_keep_all_nonzero_cross_covariance_terms(self):
        P = dense_covariance()
        H, _ = G.residual_factors('accelerometer', [0]*3, f_hat=[1, 2, -9], R_hat=G.eye(3))
        N, S = G.measurement_operands(P, G.eye(3), H, active_bias=False)
        expected = G.mm(P, G.transpose([row[:18]+[F(0)]*3 for row in H]))
        expected[18:21] = G.zeros(3, 3)
        self.assertEqual(N, expected)
        self.assertEqual(S, G.plus(G.mm(G.mm(H, P), G.transpose(H)), G.eye(3)))

    def test_descriptor_is_exact_at_finite_errors_all_three_measurements(self):
        z = [F((i % 3)-1, 10000) for i in range(24)]
        for active in (False, True):
            for kind in ('accelerometer', 'magnetometer', 'S_zero'):
                kw = dict(f_hat=[1, 2, -9], R_hat=G.eye(3), m_hat=[22, 0, 43])
                H, Hb = G.residual_factors(kind, z[:3], **kw)
                N, S = G.measurement_operands(dense_covariance(), G.eye(3), H, active_bias=active)
                # Source-coordinate values exercise algebra, not BRMM admission.
                nu = [F(1, 10000), F(-1, 12000), F(1, 15000)]
                r = [x+y for x, y in zip(G.mv(Hb, z), nu)]
                q = solve3(S, r)
                d = G.mv(N, q)
                w, k = G.small_quaternion_parameters(d[:3])
                g = G.descriptor(z[:3], d[:3], w, k, N, S, Hb, alpha=1)
                chi = z+q+nu+[F(1)]
                self.assertEqual(G.mv(g.equality, chi), [F(0)]*9)
                W, twice_V, _ = G.right_reset_polynomials(z[:3], d[:3], w, k)
                direct = [x/W for x in twice_V]+[z[i]-d[i] for i in range(3, 21)]+z[21:24]
                self.assertEqual(G.mv(g.after, chi), direct)
                self.assertEqual(G.mv(g.before, chi), z)
                # Exact graph is generally NOT obtained by applying its Jacobian.
                self.assertNotEqual(direct[:3], [z[i]-d[i] for i in range(3)])

    def test_wrong_c_or_theta_correction_cannot_satisfy_descriptor(self):
        H, Hb = G.residual_factors('S_zero', [0]*3)
        N, S = G.measurement_operands(dense_covariance(), G.eye(3), H, active_bias=True)
        g = G.descriptor([F(1, 10), 0, 0], [F(1, 20), 0, 0], 1, F(1, 2), N, S, Hb, alpha=1)
        self.assertNotEqual(G.mv(g.equality, [F(0)]*30+[F(1)]), [F(0)]*9)

    def test_physical_S_forcing_is_nonzero_at_zero_error(self):
        z = [F(0)]*24
        H, Hb = G.residual_factors('S_zero', z[:3])
        N, S = G.measurement_operands(dense_covariance(), G.eye(3), H, active_bias=True)
        physical_S = [F(1, 10000), 0, 0]
        nu = [-s for s in physical_S]
        q = solve3(S, nu)
        d = G.mv(N, q)
        w, k = G.small_quaternion_parameters(d[:3])
        g = G.descriptor(z[:3], d[:3], w, k, N, S, Hb, alpha=1)
        self.assertNotEqual(G.mv(g.after, z+q+nu+[F(1)]), z)
        self.assertEqual(G.mv(g.equality, z+q+nu+[F(1)]), [F(0)]*9)

    def test_projection_graph_and_same_beta_cross_terms(self):
        e, beta = [F(-4, 5), 0, 0], [F(0)]*3
        self.assertTrue(G.projection_graph_holds(e, beta, F(1, 2), F(2, 5)))
        self.assertFalse(G.projection_graph_holds(e, beta, F(1, 3), F(2, 5)))
        P = G.projection_matrix(F(1, 2))
        z = [F(0)]*18+e+beta
        self.assertEqual(G.mv(P, z)[18:21], [F(-2, 5), 0, 0])
        self.assertEqual(P[18][21], F(1, 2))
        self.assertTrue(G.projection_graph_holds([0]*3, [0]*3, 1, F(2, 5)))

    def test_rank3_full_storage_ledger_matches_dense_exactly(self):
        u = [F((i % 7)-3, 30) for i in range(24)]
        M = G.plus(G.eye(24), [[x*y for y in u] for x in u])
        B = [[F(((i+j) % 5)-2, 30) for j in range(3)] for i in range(24)]
        z, q = u, [F(1, 5), F(-1, 6), F(1, 7)]
        for a in (F(1), F(2, 3)):
            P = G.projection_matrix(a)
            self.assertEqual(G.projected_metric(M, a), G.mm(G.mm(G.transpose(P), M), P))
            out = G.mv(P, [x-y for x, y in zip(z, G.mv(B, q))])
            direct = G.dot(out, G.mv(M, out))-G.dot(z, G.mv(M, z))
            self.assertEqual(G.finite_storage_change(M, z, q, B, alpha=a), direct)

    def test_invalid_data_and_chart_poles_fail_closed(self):
        with self.assertRaises(TypeError):
            G.rational(.1)
        with self.assertRaises(ValueError):
            G.small_quaternion_parameters([F(1, 100), 0, 0])
        with self.assertRaises(ValueError):
            G.projection_matrix(0)
        H, Hb = G.residual_factors('S_zero', [0]*3)
        N, S = G.measurement_operands(G.eye(21), G.eye(3), H, active_bias=True)
        with self.assertRaises(ValueError):
            G.descriptor([0]*3, [1, 0, 0], 0, 1, N, S, Hb, alpha=1)
        with self.assertRaises(ValueError):
            G.measurement_operands(G.eye(18), G.eye(3), H, active_bias=False)


if __name__ == '__main__':
    unittest.main()
