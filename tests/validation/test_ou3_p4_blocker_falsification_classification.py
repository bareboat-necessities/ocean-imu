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

    def test_transverse_attitude_cap_is_prior_independent_and_partial(self):
        cap = self.d["prior_independent_transverse_attitude_cap"]
        self.assertTrue(cap["prior_independent"])
        self.assertTrue(cap["per_sample_accelerometer_update_declared"])
        self.assertFalse(cap["yaw_about_specific_force_capped_here"])
        self.assertTrue(cap["yaw_needs_asynchronous_magnetometer_and_gyro_bias_transport"])
        r = cap["accelerometer_std_lower_mps2"]
        f = cap["specific_force_norm_lower_mps2"]
        self.assertGreaterEqual(cap["transverse_variance_upper_per_axis_rad2"], r * r / (f * f))
        self.assertLess(cap["transverse_trace_upper_two_axes_rad2"], 1e-2)

    def test_same_cell_Sinverse_is_what_closes_the_transverse_route(self):
        route = self.d["magnitude_route_after_transverse_cap"]
        cap = route["reset_utility_norm_max"]
        self.assertGreater(route["accelerometer_ceiling_with_Rinverse_relaxation"], cap)
        self.assertLess(route["accelerometer_ceiling_retaining_same_cell_Sinverse"], cap)
        self.assertFalse(route["Rinverse_relaxation_alone_closes_accelerometer"])
        self.assertTrue(route["same_cell_Sinverse_closes_accelerometer"])
        self.assertGreaterEqual(route["same_cell_S_over_R_lower"], 1.0)
        self.assertTrue(route["third_attitude_direction_still_open"])

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
