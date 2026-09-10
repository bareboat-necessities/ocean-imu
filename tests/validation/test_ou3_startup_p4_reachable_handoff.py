from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_startup_p4_reachable_handoff as mod  # noqa: E402


class StartupP4ReachableHandoffTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_shipping_parity_and_fail_closed_status(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(self.d["shipping_source_parity_failures"], [])
        self.assertFalse(self.d["P4_MOTION_PASS"])
        self.assertFalse(self.d["P4_PASS"])
        self.assertFalse(self.d["P5_PASS"])

    def test_translation_handoff_is_a_correlated_fiber_not_legacy_box(self):
        err = self.d["true_error_at_handoff"]
        self.assertEqual(err["position"], "e_p=0 exactly in handoff-local displacement gauge")
        self.assertEqual(err["integral_displacement"], "e_S=0 exactly in handoff-local integral gauge")
        self.assertIn("v_true", err["velocity"])
        legacy = self.d["legacy_entry_model_audit"]
        self.assertEqual(legacy["legacy_independent_S_radius_m_s"], 300.0)
        self.assertFalse(legacy["independent_S_ball_shipping_reachable_at_handoff"])
        self.assertFalse(legacy["independent_position_ball_shipping_reachable_at_handoff"])
        self.assertFalse(legacy["legacy_300_m_s_may_promote_P4"])

    def test_states_are_zero_but_errors_are_not_assumed_zero(self):
        x = self.d["estimated_state_at_handoff"]
        self.assertEqual(x["velocity_mps"], [0.0, 0.0, 0.0])
        self.assertEqual(x["position_m"], [0.0, 0.0, 0.0])
        self.assertEqual(x["integral_displacement_m_s"], [0.0, 0.0, 0.0])
        self.assertEqual(x["latent_acceleration_mps2"], [0.0, 0.0, 0.0])
        self.assertFalse(x["linear_state_reset_at_goLive"])
        self.assertTrue(x["states_inherited_from_constructor_without_preLive_propagation"])
        self.assertIn("a_w_true", self.d["true_error_at_handoff"]["latent_acceleration"])

    def test_timeout_requires_post_handoff_capture_case(self):
        modes = self.d["live_capture_modes"]
        self.assertIn("north/tuner may still be unready", modes["timeout_handoff"])
        self.assertFalse(modes["timeout_requires_full_P4_membership_at_goLive"])
        self.assertTrue(modes["capture_theorem_must_cover_early_live_transient"])
        self.assertEqual(self.d["deployed_live_handoff_upper_bound_s"], 150.0)

    def test_covariance_seed_is_not_true_error_membership(self):
        c = self.d["covariance_at_handoff"]
        self.assertFalse(c["translation_state_covariance_propagated_before_live"])
        self.assertTrue(c["translation_cross_covariances_zero_before_first_live_prediction"])
        self.assertTrue(c["aw_marginal_reseeded_to_current_committed_stationary_covariance"])
        self.assertTrue(c["aw_cross_covariances_cleared"])
        self.assertEqual(c["scheduler_elapsed_s"], 0.0)
        self.assertFalse(
            self.d["handoff_translation_coordinate_gauge"]["covariance_confidence_used_as_true_error_set"]
        )


if __name__ == "__main__":
    unittest.main()
