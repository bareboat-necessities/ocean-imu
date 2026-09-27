"""Historical reader audits with large, correlated unknown AG priors."""
from fractions import Fraction as F
import unittest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.stability.ou3_theorem.matrix_certificates import add, identity, is_psd, matmul, transpose
from tools.stability.ou3_theorem.stationary_covariance import certificate, quiet_covariance_bounds


def zeros(n, m):
    return [[F(0) for _ in range(m)] for _ in range(n)]


def diagonal(values):
    return [[v if i == j else F(0) for j in range(len(values))] for i, v in enumerate(values)]


def action(term, covariance):
    return matmul(matmul(term, covariance), transpose(term))


class QuietCovarianceTest(unittest.TestCase):
    def test_exact_root_cancellation_and_full_action_with_correlations(self):
        # A supplied rational auxiliary word audits the inequality, not
        # real-arithmetic reachability. AW/BA factors obey the proven bounds.
        t, g, field = F(1, 25), F('9.80665'), F(75)
        f, h = identity(12), zeros(3, 12)
        q, r = zeros(12, 12), diagonal([F('.1') / g**2] * 2 + [F('.64') / field**2])
        for j in range(3):
            f[j][j + 3] = t
            f[j + 6][j + 6], f[j + 9][j + 9] = F(99, 100), F(999992, 1000000)
            h[j][j] = 1
            q[j][j] = F('.00135')**2 * t + F('1e-10') * t**3 / 3
            q[j][j + 3] = q[j + 3][j] = F('1e-10') * t**2 / 2
            q[j + 3][j + 3] = F('1e-10') * t
            q[j + 6][j + 6] = F(16) * (1 - f[j + 6][j + 6]**2) + 16
            q[j + 9][j + 9] = (1 - f[j + 9][j + 9]**2) / 1600
        h[0][7] = h[0][10] = 1 / g
        h[1][6] = h[1][9] = -1 / g
        end, l0, l1 = identity(12)[:6], zeros(6, 3), zeros(6, 3)
        for j in range(3):
            l0[j + 3][j], l1[j][j], l1[j + 3][j] = -1 / t, F(1), 1 / t
        root = add(add(matmul(end, f), matmul(l0, h), -1), matmul(matmul(l1, h), f), -1)
        fresh = add(end, matmul(l1, h), -1)
        self.assertTrue(all(x == 0 for row in root for x in row[:6]))
        upper = diagonal(quiet_covariance_bounds()['historical_AG_action_upper_diagonal'])
        previous = None
        for scale in (F(1), F(10**12)):
            factor = diagonal([scale] * 6 + [F(156)] * 3 + [F(1, 40)] * 3)
            for j in range(3):
                factor[j][j + 6], factor[j + 3][j + 9] = scale / 2, scale / 3
            prior = matmul(factor, transpose(factor))
            trial = add(add(action(root, prior), action(fresh, q)),
                        add(action(l0, r), action(l1, r)))
            self.assertTrue(is_psd(add(upper, trial, -1)))
            if previous is not None:
                self.assertEqual(trial, previous)
            previous = trial

    def test_every_prefix_matrix_comparison_is_uniform_in_delay(self):
        bounds = quiet_covariance_bounds()
        pair = diagonal(bounds['historical_AG_action_upper_diagonal'])
        upper = diagonal(bounds['every_prefix_AG_upper_diagonal'])
        # Exact PSD audits complement the analytic all-delay monotone proof.
        for t in (F(0), F(1, 1000), F(6, 125)):
            f, q = identity(6), zeros(6, 6)
            for j in range(3):
                f[j][j + 3] = t
                q[j][j] = F(1, 500000) * t + F(1, 10**9) * t**3 / 3
                q[j][j + 3] = q[j + 3][j] = F(1, 10**9) * t*t / 2
                q[j + 3][j + 3] = F(1, 10**9) * t
            self.assertTrue(is_psd(add(upper, add(action(f, pair), q), -1)))

    def test_nuisance_and_full_matrix_cross_terms_are_not_dropped(self):
        bounds = quiet_covariance_bounds()
        self.assertEqual(len(bounds['full_21_upper_diagonal']), 21)
        self.assertEqual(bounds['full_21_upper_diagonal'][:6],
                         [2 * x for x in bounds['every_prefix_AG_upper_diagonal']])

    def test_subcase_certificate_does_not_close_physical_stability(self):
        c = certificate()
        self.assertTrue(c['uniform_quiet_nominal_historical_action'])
        self.assertTrue(c['quiet_nominal_homogeneous_linear_loss_exists'])
        for key in ('stationary_compatible_class_practical_stability',
                    'source_uniform_moving_A21_covariance_upper',
                    'nonlinear_every_prefix_retention', 'theorem_closed'):
            self.assertFalse(c[key])

    def test_invalid_bounds_fail_closed(self):
        for kw in ({'gap_min': 0}, {'gap_max': F(1, 1000)}, {'prefix_max': -1},
                   {'gravity_min': 0}, {'gyro_density': 0}, {'aw_std': -1}):
            with self.assertRaises(ValueError):
                quiet_covariance_bounds(**kw)


if __name__ == '__main__':
    unittest.main()
