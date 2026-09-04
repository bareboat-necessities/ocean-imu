import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_acceleration_covariance_coupling as ACC  # noqa: E402
import ou3_sea3_body_rate_covariance_coupling as RATE  # noqa: E402
import ou3_sea3_finite_horizon_concentration as CONC  # noqa: E402


class Sea3BodyRateCovarianceCouplingTests(unittest.TestCase):
    def test_twenty_minute_body_rate_coupling_is_validated_and_nonpromoting(self):
        d = RATE.build(240000)
        self.assertEqual(RATE.validate(d), [])
        self.assertTrue(d["validated_arithmetic"])
        self.assertTrue(d["outward_rounded"])
        self.assertTrue(d["finite_horizon_body_rate_candidate_producer_ready"])
        self.assertFalse(d["physical_vessel_rotational_RAO_qualified"])
        self.assertFalse(d["finite_horizon_good_event_promoted"])
        self.assertFalse(d["deterministic_left_inclusion_closed"])
        self.assertFalse(
            d["rotational_response_interface"]["universal_rotational_parameter_box_declared"]
        )
        self.assertTrue(
            d["rotational_response_interface"]["actual_vessel_or_response_family_must_supply_K_fc_q"]
        )

    def test_exact_degree_conversion_bound_matches_closed_form_outwardly(self):
        hs = 8.5
        gain = 0.125
        corner = 0.03
        exact_binary_expression = 8100.0 * hs * hs * gain * gain * corner * corner
        bound = RATE.body_rate_trace_covariance_upper_deg2_s2(
            hs,
            gain,
            corner,
            1.0,
        )
        self.assertGreaterEqual(bound, exact_binary_expression)
        self.assertLess(bound - exact_binary_expression, 1.0e-10)

    def test_nonzero_high_sea_low_corner_interface_witness_passes_twenty_minutes(self):
        d = RATE.build(240000)
        w = d["constructed_nonphysical_interface_witness"]
        self.assertTrue(w["not_a_measured_or_declared_vessel_RAO"])
        e = w["evaluation"]
        self.assertTrue(e["candidate_pass"])
        self.assertEqual(e["rotation_gain_rad_per_m"], 0.125)
        self.assertEqual(e["rotation_corner_hz"], 0.03)
        self.assertLessEqual(
            e["validated_body_rate_trace_covariance_upper_deg2_s2"],
            e["required_body_rate_trace_covariance_upper_deg2_s2"],
        )
        self.assertFalse(e["physical_vessel_rotational_RAO_qualified"])

    def test_same_gain_with_high_corner_fails_without_relaxing_p1_cap(self):
        d = RATE.build(240000)
        e = RATE.evaluate_tuple(
            d,
            hs_m=8.5,
            rotation_gain_rad_per_m=0.125,
            corner_hz=1.2,
            rolloff_power=1.0,
        )
        self.assertFalse(e["candidate_pass"])
        self.assertTrue(any("exceeds finite-horizon" in x for x in e["validation_failures"]))
        threshold = d["required_body_rate_trace_covariance_upper_deg2_s2"]
        self.assertLess(threshold, 900.0 / 108.0)
        self.assertGreater(math.nextafter(900.0 / 108.0, -math.inf), threshold)

    def test_rolloff_slower_than_first_order_is_rejected(self):
        d = RATE.build(200)
        e = RATE.evaluate_tuple(
            d,
            hs_m=1.0,
            rotation_gain_rad_per_m=0.01,
            corner_hz=0.1,
            rolloff_power=0.999,
        )
        self.assertFalse(e["candidate_pass"])
        self.assertTrue(any("q>=1" in x for x in e["validation_failures"]))

    def test_acceleration_and_body_rate_candidates_compose_in_existing_concentration_gate(self):
        samples = 240000
        concentration = CONC.build(samples)
        acc = ACC.build(samples)
        rate = RATE.build(samples)
        acc_eval = acc["high_sea_low_corner_witness"]["evaluation"]
        rate_eval = rate["constructed_nonphysical_interface_witness"]["evaluation"]
        self.assertTrue(acc_eval["candidate_pass"])
        self.assertTrue(rate_eval["candidate_pass"])
        self.assertEqual(acc["response_parameter_box_sha256"], concentration["response_parameter_box_sha256"])
        self.assertEqual(rate["response_parameter_box_sha256"], concentration["response_parameter_box_sha256"])

        combined = CONC.evaluate_covariance_candidate(
            concentration,
            acceleration_trace_covariance_upper_m2_s4=acc_eval[
                "validated_acceleration_trace_covariance_upper_m2_s4"
            ],
            body_rate_trace_covariance_upper_deg2_s2=rate_eval[
                "validated_body_rate_trace_covariance_upper_deg2_s2"
            ],
            validated_covariance_trace_enclosures=True,
            response_parameter_box_sha256=concentration["response_parameter_box_sha256"],
        )
        self.assertTrue(combined["finite_horizon_good_event_candidate_pass"])
        self.assertGreater(combined["finite_horizon_good_event_probability_lower"], 0.95)
        self.assertFalse(combined["deterministic_left_inclusion_promoted"])

    def test_hs_outside_declared_domain_fails_closed(self):
        d = RATE.build(200)
        e = RATE.evaluate_tuple(
            d,
            hs_m=math.nextafter(d["repository_total_Hs_upper_m"], math.inf),
            rotation_gain_rad_per_m=0.0,
            corner_hz=0.03,
            rolloff_power=1.0,
        )
        self.assertFalse(e["candidate_pass"])
        self.assertTrue(any("outside the declared SEA3" in x for x in e["validation_failures"]))


if __name__ == "__main__":
    unittest.main()
