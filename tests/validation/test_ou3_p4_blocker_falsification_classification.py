#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_p4_blocker_falsification_classification as CLASS


class BlockerFalsificationClassificationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = CLASS.build()

    def test_validates_and_promotes_nothing(self):
        self.assertEqual(CLASS.validate(self.d), [])
        self.assertFalse(self.d["P4_promoted_here"])
        self.assertFalse(self.d["P4_MOTION_PASS"])
        self.assertFalse(self.d["P4_PASS"])
        self.assertFalse(self.d["P5_MAY_START"])

    def test_every_open_blocker_carries_exactly_one_declared_class(self):
        table = self.d["blocker_classification"]
        self.assertEqual(set(table), set(self.d["remaining_mathematical_P4_blockers"]))
        self.assertTrue(table)
        for name, row in table.items():
            self.assertIn(row["class"], CLASS.CLASSES, name)
            self.assertTrue(row["reason"], name)

    def test_no_counterexample_and_no_rigorous_infeasibility(self):
        self.assertFalse(self.d["counterexample_class_A_found"])
        self.assertFalse(self.d["rigorous_infeasibility_class_B_established"])
        self.assertFalse(self.d["P4_unprovable_on_declared_domain_supported"])
        self.assertNotIn("A", self.d["classes_present"])
        self.assertNotIn("B", self.d["classes_present"])

    def test_frozen_map_diagnostic_is_not_promoted_to_a_counterexample(self):
        self.assertFalse(self.d["frozen_map_diagnostic_promoted_to_lower_bound"])
        self.assertFalse(self.d["diagnostic_4p5788_treated_as_counterexample"])

    def test_integral_ball_is_not_the_correction_domain_limiter(self):
        self.assertFalse(self.d["integral_entry_ball_limits_correction_domain"])
        for mode, row in self.d["correction_domain_attribution"].items():
            with self.subTest(mode=mode):
                # The declared entry set makes S=0 the worst event.
                self.assertEqual(row["limiting_event_declared_entry_set"], "S_zero")
                # Feeding the correlated relation instead moves the limiter away
                # from S=0 without closing the ceiling, so the entry ball is not
                # what blocks this obligation.
                self.assertNotEqual(row["limiting_event_with_correlated_integral_relation"], "S_zero")
                self.assertFalse(row["closes_on_correlated_relation_alone"])
                self.assertGreater(row["attitude_trace_excess_factor"], 1e9)

    def test_transverse_attitude_cap_is_conditional_not_prior_independent(self):
        # The prior-independent form of this cap was asserted in an earlier
        # revision and is false for the deployed residual, which carries the
        # latent-acceleration and accelerometer-bias blocks alongside attitude.
        # This test pins the retraction so the claim cannot come back silently.
        cap = self.d["conditional_transverse_attitude_cap"]
        self.assertFalse(cap["prior_independent"])
        self.assertTrue(cap["conditional_on_joint_latent_bias_block"])
        self.assertGreater(cap["refutation_exceedance_factor"], 1.0)
        self.assertGreater(cap["refutation_achieved_marginal_rad2"],
                           cap["retracted_prior_independent_value_rad2"])
        self.assertTrue(cap["per_sample_accelerometer_update_declared"])
        self.assertFalse(cap["yaw_about_specific_force_capped_here"])
        self.assertTrue(cap["yaw_needs_asynchronous_magnetometer_and_gyro_bias_transport"])
        r = cap["accelerometer_std_lower_mps2"]
        f = cap["specific_force_norm_lower_mps2"]
        lam = cap["joint_latent_bias_lambda_max"]
        self.assertGreaterEqual(cap["transverse_variance_upper_per_axis_rad2"],
                                (r * r + 2.0 * lam) / (f * f))
        # And it is far above the retracted value, so nothing downstream may
        # keep treating the attitude block as capped at 1e-3 rad^2.
        self.assertGreater(cap["transverse_variance_upper_per_axis_rad2"], 1e-2)

    def test_one_shot_measurement_route_does_not_close_the_transverse_route(self):
        # With the retracted cap this route appeared to close under the same-cell
        # S^{-1} accounting. Under the conditional cap that conclusion does not
        # survive: the S/R discount grows with the cap, so the product saturates
        # and a looser attitude bound cannot be traded back for a tighter
        # ceiling. Recorded so the retraction is visible in the executable gate.
        route = self.d["magnitude_route_after_transverse_cap"]
        cap = route["reset_utility_norm_max"]
        self.assertGreater(route["accelerometer_ceiling_with_Rinverse_relaxation"], cap)
        self.assertFalse(route["Rinverse_relaxation_alone_closes_accelerometer"])
        self.assertGreaterEqual(route["same_cell_S_over_R_lower"], 1.0)
        self.assertTrue(route["third_attitude_direction_still_open"])
        self.assertFalse(route["one_shot_measurement_route_closes_the_correction_domain"])

    def test_one_shot_route_carries_a_sharp_threshold_not_just_a_failure(self):
        # The value of the retraction is that the route has an exact threshold:
        # ceiling(P)^2 = 2 P E r/(f^2 P + r) is increasing with a finite
        # supremum, and that supremum is above the reset utility cap, so there
        # is a P* at which the route would close. Recording P* -- and the
        # lambda_max(P_(a_w,b_a)) target it implies -- is what turns this open
        # blocker into a checkable objective.
        th = self.d["magnitude_route_after_transverse_cap"]["one_shot_route_threshold"]
        cap = self.d["magnitude_route_after_transverse_cap"]["reset_utility_norm_max"]
        self.assertTrue(th["supremum_exceeds_reset_utility_cap"])
        self.assertFalse(th["saturates_below_reset_cap_for_free"])
        self.assertGreater(th["ceiling_supremum_over_all_attitude_bounds"], cap)
        # The attained same-cell ceiling must sit below the supremum and above
        # the cap: that is exactly the near miss the threshold explains.
        attained = self.d["magnitude_route_after_transverse_cap"][
            "accelerometer_ceiling_retaining_same_cell_Sinverse"]
        self.assertGreater(attained, cap)
        self.assertLessEqual(attained, th["ceiling_supremum_over_all_attitude_bounds"])
        self.assertTrue(th["requirement_is_achievable_in_principle"])
        self.assertGreater(th["required_transverse_variance_per_axis_rad2"], 0.0)
        self.assertLess(th["required_transverse_variance_per_axis_rad2"],
                        th["attained_transverse_variance_per_axis_rad2"])
        self.assertGreater(th["transverse_variance_shortfall_factor"], 1.0)
        self.assertGreater(th["required_joint_latent_bias_lambda_max"], 0.0)
        self.assertLess(th["required_joint_latent_bias_lambda_max"],
                        th["attained_joint_latent_bias_lambda_max"])
        self.assertGreater(th["joint_latent_bias_shortfall_factor"], 1.0)

    def test_correlated_relation_summary_is_carried_through(self):
        rel = self.d["correlated_integral_relation"]
        self.assertEqual(rel["declared_independent_radius_m_s"], 300.0)
        self.assertLess(rel["correlated_one_cadence_radius_m_s"],
                        rel["chart_retaining_radius_upper_m_s"])
        self.assertTrue(rel["one_cadence_meets_chart_threshold"])
        self.assertFalse(rel["S_zero_gives_uniform_anchor_without_dwell"])
        self.assertTrue(rel["unconditional_anchor_requires_P5_or_reachable_P_SS_lower_bound"])

    def test_finite_precision_layers_stay_separate(self):
        fp = self.d["finite_precision"]
        self.assertTrue(fp["conditional_mathematical_binary32_additive_ISS_closed"])
        self.assertFalse(fp["target_toolchain_qualified"])
        self.assertFalse(fp["missing_toolchain_qualification_is_a_mathematical_P4_blocker"])
        self.assertFalse(fp["host_arithmetic_check_is_deployment_qualification"])
        self.assertTrue(fp["remaining_deployment_blockers"])

    def test_every_blocker_the_gate_can_emit_is_classifiable(self):
        # The gate emits three BIAS-family strings the moment a family flag
        # regresses. They must be classifiable so the producer reports the
        # regression instead of aborting before it writes anything.
        for family in ("BIAS0", "BIAS1", "BIAS2"):
            key = "source-uniform %s same-history physical-driver family" % family
            with self.subTest(family=family):
                self.assertIn(key, CLASS.BLOCKER_CLASSES)
                self.assertIn(key, CLASS.BLOCKER_REASONS)
                self.assertEqual(CLASS.BLOCKER_CLASSES[key], "E")
        self.assertEqual(set(CLASS.BLOCKER_CLASSES), set(CLASS.BLOCKER_REASONS))

    def test_coverage_helper_reports_a_gap_instead_of_raising(self):
        self.assertEqual(CLASS.every_gate_blocker_is_classifiable(
            {"remaining_mathematical_P4_blockers": list(CLASS.BLOCKER_CLASSES)}), [])
        self.assertEqual(CLASS.every_gate_blocker_is_classifiable(
            {"remaining_mathematical_P4_blockers": ["a blocker nobody classified"]}),
            ["a blocker nobody classified"])

    def test_class_and_reason_maps_must_agree(self):
        d = dict(self.d)
        original = dict(CLASS.BLOCKER_REASONS)
        try:
            CLASS.BLOCKER_REASONS.pop(next(iter(CLASS.BLOCKER_REASONS)))
            self.assertIn("blocker class and reason maps disagree", CLASS.validate(d))
        finally:
            CLASS.BLOCKER_REASONS.clear()
            CLASS.BLOCKER_REASONS.update(original)

    def test_unknown_blocker_text_is_rejected_rather_than_ignored(self):
        d = dict(self.d)
        d["remaining_mathematical_P4_blockers"] = ["a blocker nobody classified"]
        self.assertIn("classification table does not cover exactly the open blockers",
                      CLASS.validate(d))


if __name__ == "__main__":
    unittest.main()
