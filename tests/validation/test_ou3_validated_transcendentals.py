import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from ou3_interval import Interval
import ou3_validated_transcendentals as mod


class ValidatedTranscendentalTests(unittest.TestCase):
    def test_exp_neg_contains_libm_reference_across_supported_domain(self):
        points = [0.0, 1.0 / 2400.0, 0.005 / 0.02, 0.5, 1.0]
        points += [i / 64.0 for i in range(65)]
        for x in points:
            lo, hi = mod.exp_neg_scalar_bounds(x)
            self.assertLessEqual(lo, math.exp(-x))
            self.assertGreaterEqual(hi, math.exp(-x))

    def test_interval_exp_neg_is_monotone_and_contains_endpoints(self):
        for lo, hi in ((0.0, 0.25), (1e-4, 0.01), (0.2, 0.9)):
            I = mod.exp_neg(Interval(lo, hi))
            self.assertLessEqual(I.lo, math.exp(-hi))
            self.assertGreaterEqual(I.hi, math.exp(-lo))

    def test_ou_coefficients_enclose_direct_formulas(self):
        h = Interval.outward_bounds(0.005, 0.005)
        tau = Interval.outward_bounds(0.02, 12.0)
        c = mod.ou_discrete_coefficients(h, tau)
        for t in (0.02, 0.05, 0.5, 1.1, 4.0, 12.0):
            x = 0.005 / t
            alpha = math.exp(-x)
            phi_pa = t * t * (x + math.expm1(-x))
            phi_sa = t ** 3 * (0.5 * x * x - x - math.expm1(-x))
            self.assertLessEqual(c["alpha"].lo, alpha)
            self.assertGreaterEqual(c["alpha"].hi, alpha)
            self.assertLessEqual(c["phi_pa"].lo, phi_pa)
            self.assertGreaterEqual(c["phi_pa"].hi, phi_pa)
            self.assertLessEqual(c["phi_Sa"].lo, phi_sa)
            self.assertGreaterEqual(c["phi_Sa"].hi, phi_sa)

    def test_rejects_unvalidated_large_argument(self):
        with self.assertRaises(ValueError):
            mod.exp_neg(Interval(0.0, 1.01))
        with self.assertRaises(ValueError):
            mod.ou_discrete_coefficients(
                Interval.point(0.1), Interval.point(0.05)
            )


if __name__ == "__main__":
    unittest.main()
