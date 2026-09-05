from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_shipping_prediction_primitives as mod  # noqa: E402
from ou3_interval import Interval  # noqa: E402


class Sea3ShippingPredictionPrimitivesTest(unittest.TestCase):
    def test_contract_is_source_neutral_and_shipping_bound(self):
        d = mod.build()
        self.assertEqual(mod.validate(d), [])
        self.assertTrue(d["validated_matrix_primitives_ready"])
        self.assertTrue(d["consumes_only_SEA3_derived_sample_coordinates"])
        self.assertFalse(d["source_domain_created_here"])
        self.assertFalse(d["arbitrary_bounded_input_source_created_here"])
        self.assertFalse(d["independent_tau_sigma_source_created_here"])
        self.assertFalse(d["P3_promoted"])

    def test_zero_rate_attitude_block_has_exact_small_rate_limit(self):
        h = Interval.point(0.005)
        qg = [Interval.point(1.0e-4)] * 3
        qb = Interval.point(1.0e-11)
        F, Q = mod.attitude_gyro_bias_F_Q(
            [Interval.point(0.0)] * 3, h, qg, qb
        )
        self.assertEqual((len(F), len(F[0])), (6, 6))
        self.assertEqual((len(Q), len(Q[0])), (6, 6))
        self.assertTrue(F[0][0].lo <= 1.0 <= F[0][0].hi)
        self.assertTrue(F[0][3].lo <= 0.005 <= F[0][3].hi)
        self.assertGreater(Q[0][0].lo, 0.0)
        self.assertGreater(Q[3][3].lo, 0.0)

    def test_integrated_ou_transition_covers_both_shipping_branches(self):
        h = Interval.point(0.005)
        # h/tau exactly at the source branch boundary.
        tau = Interval.point(0.5)
        F = mod.translation_axis_transition(tau, h)
        Q = mod.translation_axis_process(tau, h, Interval.point(0.25))
        self.assertEqual((len(F), len(F[0])), (4, 4))
        self.assertEqual((len(Q), len(Q[0])), (4, 4))
        self.assertTrue(F[0][0].lo <= 1.0 <= F[0][0].hi)
        self.assertGreater(F[0][3].hi, 0.0)
        for i in range(4):
            for j in range(4):
                self.assertEqual(Q[i][j], Q[j][i])

    def test_active_bias_matches_exact_GM_form(self):
        phi, Q = mod.active_accel_bias_F_Q(
            Interval.point(0.005),
            Interval.point(300.0),
            Interval.point(1.0e-6),
        )
        self.assertGreater(phi.lo, 0.0)
        self.assertLess(phi.hi, 1.0)
        for i in range(3):
            self.assertGreater(Q[i][i].lo, 0.0)


if __name__ == "__main__":
    unittest.main()
