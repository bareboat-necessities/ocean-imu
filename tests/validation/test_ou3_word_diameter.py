"""Word Riccati diameter, kernel-bounded contraction and gyro persistence cap."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tools.stability.ou3_theorem.word_diameter import (  # noqa: E402
    certificate, closed_loop, full_root_design, generalized_max_bracket, gyro_persistence_cap,
    information_ratio_bound, k0_bound_upper, scalar_kernel_ceiling, terminal_pair)
from tools.stability.ou3_theorem.matrix_certificates import is_psd  # noqa: E402


class WordDiameterTest(unittest.TestCase):
    def test_committed_certificate_matches_exact_reproduction(self):
        path = ROOT/'reports/results/ou3_stability/word-diameter-certificate.json'
        record = json.loads(path.read_text())
        self.assertEqual(record, certificate())
        for key in ('root_covariance_bound_used', 'source_uniform_diameter_upper_bound',
                    'rho0_certified', 'theorem_closed'):
            self.assertFalse(record[key])
        self.assertTrue(record['diameter_identity']['lambda_max(J^-1 A) = lambda_max(Pi^-1 P_diff)'])
        self.assertTrue(record['sharp_scalar_example']['k0_bound_sharp'])
        self.assertFalse(record['kernel_bounded_corollary']['full_diameter_finite'])

    def test_closed_form_matches_numerical_supremum(self):
        for kappa, k in ((4.0, 0.0), (50.0, 3.0), (7.0, 20.0), (1e4, 0.5)):
            ys = [0.0]+[10**(-8+10*i/200000) for i in range(200001)]
            grid = max(1/(1+y)-1/((1+k)*(1+kappa*y)) for y in ys)
            self.assertGreaterEqual(information_ratio_bound(kappa, k), grid-1e-12)
            self.assertAlmostEqual(information_ratio_bound(kappa, k), grid, places=7)
            self.assertLess(information_ratio_bound(kappa, k), 1)

    def test_scalar_word_attains_the_k0_bound(self):
        events = [{'kind': 'correction', 'H': [[1]], 'V': [[1]]},
                  {'kind': 'prediction', 'F': [[1]], 'U': [[1]]}]
        pair = terminal_pair(full_root_design(1, events))
        self.assertEqual(pair['P_diff'][0][0]/pair['Pi'][0][0], 2)
        bound = k0_bound_upper(2)
        best = F(0)
        for i in range(1, 400):
            p = F(i, 400)*2
            m, pend = closed_loop([[p]], events)
            rho = m[0][0]**2*p/pend[0][0]
            self.assertLessEqual(rho, bound)
            best = max(best, rho)
        self.assertGreater(best, bound-F(1, 10**4))

    def test_diameter_bracket_is_consistent(self):
        pi = [[F(1), F(0)], [F(0), F(2)]]
        pdiff = [[F(3), F(1)], [F(1), F(5)]]
        lo, hi = generalized_max_bracket(pdiff, pi, F(1, 10**6))
        self.assertTrue(is_psd([[hi*pi[i][j]-pdiff[i][j] for j in range(2)] for i in range(2)]))
        self.assertFalse(is_psd([[lo*pi[i][j]-pdiff[i][j] for j in range(2)] for i in range(2)]))

    def test_gyro_persistence_cap_scales_with_inverse_square_word(self):
        caps = {t: gyro_persistence_cap(t)['kappa_lower'] for t in (16, 32, 64)}
        self.assertGreater(caps[16], F(71))
        self.assertGreater(caps[64], F(44, 10))
        self.assertEqual(caps[16]/caps[32], 4)
        self.assertLess(gyro_persistence_cap(64)['k0_margin_cap_upper'], F(65, 100))

    def test_scalar_kernel_ceiling_is_dominated_by_the_ba_marginal(self):
        ba_only = F('9.80665')**2/1600
        with_tilt = scalar_kernel_ceiling(F(1, 1000), F('9.80665'))
        self.assertGreaterEqual(with_tilt, ba_only)
        self.assertLess(with_tilt, F(13, 10)*ba_only)
        self.assertGreater(scalar_kernel_ceiling(F(1, 100), F(196, 100)), 4*F(196, 100)**2/1600)


if __name__ == '__main__':
    unittest.main()
