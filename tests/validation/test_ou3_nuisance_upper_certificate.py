"""Check exact root elimination and conservative source-bound constants."""
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.stability.ou3_theorem.nuisance_upper_certificate import bounds, certificate, interpolation_rows
from tools.stability.ou3_theorem.matrix_certificates import matmul


class NuisanceUpperTests(unittest.TestCase):
    def test_interpolation_cancels_arbitrary_neutral_root(self):
        # Unequal gaps and a delayed target exercise all cross terms.
        for a, b, d in [(F(8), F('8.156'), F('.156')),
                        (F('8.071'), F('8.023'), F('.037')),
                        (F('8.156'), F(8), F(0))]:
            t = a+b+d
            obs = [[s*s/2, s, F(1)] for s in [F(0), a, a+b]]
            expected = [[F(1), F(0), F(0)], [t, F(1), F(0)], [t*t/2, t, F(1)]]
            self.assertEqual(matmul(interpolation_rows(a,b,d), obs), expected)
            for row, ceiling in zip(interpolation_rows(a,b,d), bounds()[2]):
                self.assertLessEqual(sum(map(abs,row)), ceiling)

    def test_scalar_induction_and_fresh_noise_majorants(self):
        aw, fresh, _, _, diagonal = bounds()
        self.assertLess(aw, F(156)**2)
        self.assertLess(fresh, 1)
        self.assertEqual(diagonal[-1], F(5,1600))
        # BA initial/release variance is below the invariant stationary ceiling.
        self.assertLess(F('.004')**2, F(1,1600))
        self.assertFalse(certificate()['full_21_covariance_upper_verified'])
        self.assertFalse(certificate()['theorem_closed'])

    def test_reject_degenerate_observation_nodes(self):
        with self.assertRaises(ValueError):
            interpolation_rows(F(0), F(8))

if __name__ == '__main__':
    unittest.main()
