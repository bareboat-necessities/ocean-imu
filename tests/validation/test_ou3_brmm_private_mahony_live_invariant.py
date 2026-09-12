from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_brmm_private_mahony_live_invariant as mod  # noqa: E402
import ou3_brmm_gravity_direction_forcing_qualification as DIR  # noqa: E402


class BrmmPrivateMahonyTwoPhaseInvariantTest(unittest.TestCase):
    """Two nested levels of one metric close at the padded 8.8 m/s^2 envelope.

    The retired single-level formulation could not: the level needed to contain
    the first-sample accelerometer seed exceeded the largest level fitting
    inside the 87 deg proof chart.  The replacement keeps one metric and splits
    the roles - an outer level for seed containment and chart retention, an
    inner level for the ultimate bound - which works because inward flow holds
    on every level above the minimal certifiable one, not only on a boundary.
    """

    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_certificate_closes(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertTrue(self.d["continuous_all_live_PI_invariant_closed"])
        self.assertTrue(self.d["invariant_strictly_inside_87deg_chart"])
        self.assertFalse(self.d["invariant_strictly_inside_60deg_chart"])

    def test_seed_angle_is_bound_to_the_padded_source_not_a_constant(self):
        g = self.d["BRMM_direction_geometry"]
        # The declared rational bound must dominate the source's own angle.
        self.assertTrue(g["seed_tilt_dominates_source_angle"])
        self.assertGreaterEqual(
            g["seed_tilt_rad_upper"], DIR.build()["instantaneous_direction_angle_upper_rad"]
        )
        # asin(8.8/9.80665) = 1.1137 rad, not the retired 0.955 rad constant.
        self.assertLess(g["seed_tilt_rad_upper"], 1.115)
        self.assertGreater(g["seed_tilt_deg_upper"], 63.0)
        self.assertTrue(g["same_history_decomposition_retained"])

    def test_three_level_constraints_are_jointly_satisfied(self):
        lv = self.d["two_phase_levels"]
        self.assertLessEqual(lv["seed_level_required"], lv["outer_level_C"])
        self.assertLess(lv["outer_level_C"], lv["chart_level_ceiling"])
        self.assertTrue(lv["seed_contained_in_outer_level"])
        self.assertTrue(lv["outer_level_inside_chart"])
        self.assertGreater(lv["chart_headroom_relative"], 0.0)
        self.assertLess(lv["inner_level_C"], lv["outer_level_C"])

    def test_both_level_boundaries_have_strict_inward_flow(self):
        for key in ("boundary_validation", "inner_boundary_validation"):
            b = self.d[key]
            self.assertTrue(b["closed"], key)
            self.assertGreater(b["strict_inward_margin_lower"], 0.0, key)
            self.assertEqual(b["cells_per_chart"], mod.CELLS_PER_CHART)
            self.assertEqual(b["endpoint_sector_checks"], 2)

    def test_capture_is_exponential_and_reaches_the_inner_level(self):
        cap = self.d["capture"]
        self.assertTrue(cap["inner_level_above_minimal"])
        self.assertGreater(cap["contraction_rate_lower_per_s"], 0.0)
        self.assertTrue(cap["monotone_decrease_above_inner_level"])
        # W_in = 2*support/q, and the rate is q/2.
        self.assertAlmostEqual(
            cap["minimal_certifiable_sqrt_level"],
            2.0 * cap["support_upper_ceiling"] / cap["q_lower_floor"],
            places=9,
        )
        self.assertAlmostEqual(
            cap["contraction_rate_lower_per_s"], cap["q_lower_floor"] / 2.0, places=12
        )
        self.assertLess(cap["sqrt_level_at_horizon_upper"], self.d["sqrt_C"])
        self.assertFalse(cap["capture_uses_deployed_timeout_as_hypothesis"])

    def test_tilt_bounds_are_ordered_and_inside_the_chart(self):
        self.assertLess(
            self.d["ultimate_tilt_deg_upper"], self.d["tilt_at_deployed_horizon_deg_upper"]
        )
        self.assertLess(
            self.d["tilt_at_deployed_horizon_deg_upper"], self.d["actual_tilt_deg_upper"]
        )
        self.assertLess(self.d["actual_tilt_deg_upper"], 87.0)
        self.assertLess(self.d["actual_tilt_rad_upper"], self.d["chart_radius_rad"])

    def test_metric_is_generic_and_consistent(self):
        p = self.d["metric_skew_p"]
        c = self.d["metric_scale_c"]
        self.assertEqual(self.d["metric_cholesky_R"], [[1.0, -p], [0.0, c]])
        self.assertEqual(self.d["metric_P"], [[1.0, -p], [-p, p * p + c * c]])
        self.assertAlmostEqual(self.d["metric_det"], c * c, places=12)
        self.assertFalse(self.d["metric_tuned_toward_magnetic_requirement"])

    def test_exp_bound_is_outward_and_composed(self):
        # exp(-x) for x beyond the validated half-step range.
        for x in (0.25, 1.0, 3.0, 4.0):
            self.assertGreaterEqual(mod._exp_neg_upper(x), math.exp(-x))
            self.assertLess(mod._exp_neg_upper(x), math.exp(-x) * 1.01)
        with self.assertRaises(ValueError):
            mod._exp_neg_upper(-1.0)

    def test_same_brmm_forcing_not_an_alternate_source(self):
        self.assertTrue(self.d["same_BRMM_specific_force_direction_required"])
        self.assertTrue(self.d["same_BRMM_gyro_bias_forcing_required"])
        self.assertTrue(self.d["same_history_direction_primitive_required"])
        self.assertFalse(self.d["source_generator"])
        self.assertFalse(self.d["trajectory_replay_used"])
        self.assertFalse(self.d["arbitrary_bounded_input_source_used"])

    def test_downstream_obligations_remain_fail_closed(self):
        self.assertFalse(self.d["shipping_binary32_discrete_invariant_closed"])
        self.assertFalse(self.d["complete_BRMM_family_materialized_here"])
        self.assertFalse(self.d["P3_promoted"])


if __name__ == "__main__":
    unittest.main()
