from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_brmm_private_mahony_live_invariant as mod  # noqa: E402
import ou3_brmm_gravity_direction_forcing_qualification as DIR  # noqa: E402


class BrmmPrivateMahonyLiveInvariantTest(unittest.TestCase):
    """The 87 deg-chart PI invariant does not close for the padded envelope.

    The certificate is retained as evidence, not as a passing lemma: the level
    needed to contain the first-sample accelerometer seed and the largest level
    that still fits inside the proof chart form an empty window.  Any future
    closure attempt must move the metric or the forcing qualification, and must
    consciously update these assertions.
    """

    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_seed_angle_is_bound_to_the_padded_source_not_a_constant(self):
        g = self.d["BRMM_direction_geometry"]
        q = DIR.build()
        self.assertAlmostEqual(
            g["initial_seed_tilt_rad_upper"],
            q["instantaneous_direction_angle_upper_rad"],
            places=15,
        )
        # asin(8.8/9.80665) = 1.1137 rad, not the retired 0.955 rad constant.
        self.assertAlmostEqual(g["initial_seed_tilt_rad_upper"], 1.1137282529726653, places=12)
        self.assertTrue(g["initial_seed_angle_closed"])
        self.assertTrue(g["same_history_decomposition_retained"])

    def test_admissible_level_window_is_empty_at_the_padded_envelope(self):
        w = self.d["admissible_level_window"]
        self.assertFalse(w["window_nonempty"])
        self.assertGreater(w["level_lower_required_to_contain_seed"], w["level_upper_allowed_by_proof_chart"])
        self.assertGreater(w["emptiness_factor"], 1.0)
        self.assertAlmostEqual(w["declared_level_C"], 1.4641, places=12)

    def test_invariant_stays_fail_closed_with_its_failures_recorded(self):
        failures = mod.validate(self.d)
        self.assertFalse(self.d["initial_set_inside_invariant"])
        self.assertFalse(self.d["continuous_all_live_PI_invariant_closed"])
        self.assertIn("initial_set_inside_invariant is not true", failures)
        self.assertTrue(any("no metric level contains the seed" in f for f in failures))

    def test_retained_geometry_is_still_certified(self):
        # The boundary flow and chart containment of the declared level are
        # unaffected by the seed obstruction and stay available to a successor.
        self.assertTrue(self.d["invariant_strictly_inside_87deg_chart"])
        self.assertFalse(self.d["invariant_strictly_inside_60deg_chart"])
        self.assertLess(self.d["actual_tilt_deg_upper"], 87.0)
        self.assertGreater(self.d["boundary_validation"]["strict_inward_margin_lower"], 0.0)

    def test_certified_tilt_is_a_level_radius_far_above_magnetic_needs(self):
        # The certified number is the ellipse radius, not an accuracy bound, and
        # it is what the magnetic capture composition has to consume.
        self.assertGreater(self.d["actual_tilt_deg_upper"], 80.0)
        self.assertAlmostEqual(
            self.d["actual_tilt_rad_upper"],
            math.sqrt(self.d["z_tilt_projection_sq_upper"])
            + 0.1 * self.d["BRMM_direction_geometry"]["direction_primitive_norm_upper_s"],
            places=12,
        )

    def test_same_brmm_forcing_not_an_alternate_source(self):
        self.assertTrue(self.d["same_BRMM_specific_force_direction_required"])
        self.assertTrue(self.d["same_BRMM_gyro_bias_forcing_required"])
        self.assertTrue(self.d["same_history_direction_primitive_required"])
        self.assertFalse(self.d["source_generator"])
        self.assertFalse(self.d["trajectory_replay_used"])
        self.assertFalse(self.d["arbitrary_bounded_input_source_used"])

    def test_discrete_shipping_composition_remains_fail_closed(self):
        self.assertFalse(self.d["shipping_binary32_discrete_invariant_closed"])
        self.assertFalse(self.d["complete_BRMM_family_materialized_here"])
        self.assertFalse(self.d["P3_promoted"])


if __name__ == "__main__":
    unittest.main()
