#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from ou3_interval import Interval
import ou3_p4_prediction_directional_transport as T


class PredictionDirectionalTransportTests(unittest.TestCase):
    def I(self, x):
        return Interval.point(float(x))

    def test_pullback_fixture_preserves_rank_one_direction(self):
        F = [[self.I(1.0), self.I(2.0)], [self.I(0.0), self.I(1.0)]]
        Q = [[self.I(1.0), self.I(0.0)], [self.I(0.0), self.I(0.0)]]
        R = T.pullback(F, Q)
        self.assertTrue(R[0][0].contains(1.0))
        self.assertTrue(R[0][1].contains(2.0))
        self.assertTrue(R[1][0].contains(2.0))
        self.assertTrue(R[1][1].contains(4.0))

    def test_source_fixed_mode_prediction_is_invertible(self):
        d = T.build()
        self.assertEqual(T.validate(d), [])
        self.assertTrue(d["P4_PREDICTION_DIRECTIONAL_TRANSPORT_ESTABLISHED"])
        self.assertTrue(d["P3_process_noise_metric_comparison_retained"])
        self.assertFalse(d["per_prediction_scalarization_used"])
        self.assertFalse(d["condition_number_conversion_used"])
        self.assertFalse(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"])
        self.assertFalse(d["P5_FINITE_CAPTURE_ESTABLISHED"])
        for mode, n in (("H", 18), ("A", 21)):
            m = d["modes"][mode]
            self.assertEqual(m["dimension"], n)
            self.assertTrue(m["deterministic_transition_invertible"])
            self.assertGreater(m["deterministic_transition_determinant_abs_lower"], 0.0)
            self.assertTrue(m["directional_rank_preserved_exactly"])
            self.assertTrue(m["directional_nullity_preserved_exactly"])
            self.assertFalse(m["process_covariance_dropped_from_metric_argument"])
            self.assertFalse(m["P4_PROMOTED"])

    def test_positive_exp_factors(self):
        lo, hi = T._positive_exp_factor(0.0, 0.1)
        self.assertGreater(lo, 0.0)
        self.assertGreaterEqual(hi, lo)
        self.assertLessEqual(hi, 1.00000000000001)
        self.assertTrue(math.isfinite(lo))


if __name__ == "__main__":
    unittest.main()
