from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_acceleration_covariance_coupling as coupling  # noqa: E402


class Sea3AccelerationCovarianceCouplingTests(unittest.TestCase):
    def test_twenty_minute_coupling_is_nonempty_but_rejects_cartesian_extreme(self):
        d = coupling.build(240000)
        self.assertEqual(coupling.validate(d), [])
        self.assertEqual(d["repository_total_Hs_upper_m"], 8.5)
        self.assertTrue(d["coupling_predicate_defined"])
        self.assertFalse(d["physical_vessel_pairing_qualified"])
        self.assertFalse(d["finite_horizon_good_event_promoted"])
        self.assertFalse(d["deterministic_left_inclusion_closed"])

        cartesian = d["independent_cartesian_extreme"]["evaluation"]
        self.assertFalse(cartesian["candidate_pass"])
        self.assertGreater(
            cartesian["validated_acceleration_trace_covariance_upper_m2_s4"],
            cartesian["required_acceleration_trace_covariance_upper_m2_s4"],
        )

        nonempty = d["high_sea_low_corner_witness"]["evaluation"]
        self.assertTrue(nonempty["candidate_pass"])
        self.assertLessEqual(
            nonempty["validated_acceleration_trace_covariance_upper_m2_s4"],
            nonempty["required_acceleration_trace_covariance_upper_m2_s4"],
        )

    def test_bound_matches_closed_form_and_is_monotone(self):
        hs = 8.5
        gain = 4.0
        fc = 0.03
        validated = coupling.acceleration_trace_covariance_upper_m2_s4(hs, gain, fc)
        nearest = math.pi**4 * hs**2 * gain**2 * fc**4
        self.assertGreater(validated, nearest)
        self.assertGreater(
            coupling.acceleration_trace_covariance_upper_m2_s4(hs, gain, 0.04),
            validated,
        )
        self.assertGreater(
            coupling.acceleration_trace_covariance_upper_m2_s4(hs, 4.0, fc),
            coupling.acceleration_trace_covariance_upper_m2_s4(hs, 2.0, fc),
        )

    def test_candidate_evaluation_is_fail_closed(self):
        d = coupling.build(200)
        bad = coupling.evaluate_tuple(d, hs_m=-1.0, gain=1.0, corner_hz=0.03)
        self.assertFalse(bad["candidate_pass"])
        self.assertTrue(bad["validation_failures"])


if __name__ == "__main__":
    unittest.main()
