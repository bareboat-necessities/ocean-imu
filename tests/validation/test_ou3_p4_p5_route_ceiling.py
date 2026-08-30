from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_p5_route_ceiling_certificate as CEILING


class Ou3P4P5RouteCeilingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = CEILING.build()

    def test_route_ceiling_validates_and_reports_the_obstruction(self):
        d = self.d
        self.assertEqual(CEILING.validate(d), [])
        self.assertEqual(d["P4_P5_UNIFORM_TRANSPORT_ROUTE_CEILING"], "BELOW_P1_HANDOFF")
        self.assertTrue(d["ceiling_is_about_the_proof_route_not_the_filter"])
        self.assertFalse(d["source_replay_used"])
        self.assertEqual(set(d["modes"]), {"H", "A"})

    def test_ceiling_is_below_every_bounded_P1_handoff_node(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            self.assertFalse(m["route_can_reach_P1_handoff"])
            self.assertLess(m["route_ceiling_absolute"], m["smallest_bounded_P1_handoff_cayley_norm"])
            self.assertLess(m["route_ceiling_at_shipping_prefix_factor"], m["route_ceiling_absolute"])
            self.assertGreater(m["shortfall_factor_ceiling_vs_largest_handoff"], 50.0)
            self.assertGreater(m["shortfall_factor_ceiling_vs_smallest_handoff"], 20.0)

    def test_ceiling_matches_delta_over_prefix_factor_times_injections(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            n = m["word_samples_upper"]
            self.assertAlmostEqual(m["route_ceiling_absolute"], 1.0 / n, delta=1.0e-12)
            self.assertAlmostEqual(
                m["route_ceiling_at_shipping_prefix_factor"],
                1.0 / (m["prefix_W_factor_upper"] * n),
                delta=1.0e-12,
            )

    def test_ceiling_does_not_depend_on_the_covariance_or_metric_scale(self):
        # theta_capture <= delta/(F N eps) has no a_t in it: the attitude chart
        # scale cancels, so no covariance retightening can move the ceiling.
        self.assertTrue(self.d["attitude_chart_scale_cancels_from_the_ceiling"])
        h, a = self.d["modes"]["H"], self.d["modes"]["A"]
        self.assertNotAlmostEqual(h["attitude_chart_scale"], 0.0)
        self.assertAlmostEqual(h["route_ceiling_absolute"], a["route_ceiling_absolute"], delta=1.0e-15)

    def test_current_radius_is_inside_its_own_ceiling_and_far_below_it(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            self.assertGreater(m["certified_attitude_capture_radius_now"], 0.0)
            self.assertLess(
                m["certified_attitude_capture_radius_now"],
                m["route_ceiling_at_shipping_prefix_factor"],
            )
            self.assertGreater(m["shortfall_factor_now_vs_largest_handoff"], 1.0e30)

    def test_breakeven_shows_the_word_structure_is_the_binding_hypothesis(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            # A word would have to inject fewer than two accepted attitude
            # corrections per second for the route to ever reach the handoff.
            self.assertLess(m["breakeven_injecting_operations_per_word"], 2.0)
            self.assertLess(m["word_samples_upper"] / 10.0,
                            m["state_operation_count_upper"])
            self.assertLess(m["breakeven_injection_fraction_at_source_word"], 0.01)

    def test_required_structural_change_retires_the_kappa_search(self):
        d = self.d
        joined = " ".join(d["required_structural_change"])
        self.assertIn("information decrease of the same", joined)
        self.assertIn("finite-angle sector", joined)
        self.assertIn("kappa", d["retired_search_direction"])


if __name__ == "__main__":
    unittest.main()
