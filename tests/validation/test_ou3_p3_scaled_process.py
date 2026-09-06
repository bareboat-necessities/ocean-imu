#!/usr/bin/env python3
import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))
    sys.path.insert(0, str(TOOLS / "stability"))

from ou3_interval import Interval
import ou3_p3_scaled_process as P


class ScaledProcessTests(unittest.TestCase):
    def test_regression_cell_crossing_old_005_cutoff_is_certified(self):
        x = Interval(0.049999999999999996, 0.050000498317846004)
        pieces = P.split_x_cell(x)
        self.assertTrue(pieces)
        self.assertGreater(min(r for _cell, r in pieces), 0.0)

    def test_correlated_exact_series_covers_outward_deployed_max(self):
        hi = math.nextafter(P.DEPLOYED_X_MAX, math.inf)
        self.assertEqual(P.NEAR_EXACT_SERIES_MAX_X, hi)
        x = Interval(0.249999, hi)
        rho = P.certified_cell_rho(x)
        self.assertGreater(rho, 0.0)

    def test_range_reduced_exp_majorant_is_rigorous_past_half(self):
        q = Fraction.from_float(math.nextafter(0.5, math.inf))
        upper = P._validated_exp_upper_fraction(q)
        self.assertGreater(upper, 1)
        self.assertLess(upper, 2)

    def test_tail_bound_stays_positive_at_deployed_rate_two_endpoint(self):
        xmax = Fraction.from_float(P.NEAR_EXACT_SERIES_MAX_X)
        tail = P._exp_tail_bound(2, xmax, P.NEAR_EXACT_SERIES_ORDER)
        self.assertGreater(tail, 0)


if __name__ == "__main__":
    unittest.main()
