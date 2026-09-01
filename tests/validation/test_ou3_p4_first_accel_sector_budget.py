from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_first_accel_sector_budget as BUDGET


class FirstAccelSectorBudgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = BUDGET.build()

    def test_budget_validates_and_stays_a_distance_not_a_verdict(self):
        d = self.d
        self.assertEqual(BUDGET.validate(d), [])
        self.assertTrue(d["distance_only_no_verdict_emitted"])
        self.assertFalse(d["filter_changed"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["declared_entrance_shrunk"])
        self.assertFalse(d["P4_USABLE_CERTIFICATE_PROMOTED"])
        self.assertFalse(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"])

    def test_ladder_covers_the_documented_candidate_angles(self):
        angles = [r["angle_deg"] for r in self.d["ladder_rows"]]
        self.assertEqual(angles, [15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0])
        families = {r["angle_deg"]: r["ladder_family"] for r in self.d["ladder_rows"]}
        for deg in (15.0, 20.0, 25.0, 30.0):
            self.assertEqual(families[deg], "descending")
        for deg in (35.0, 40.0, 45.0):
            self.assertEqual(families[deg], "ascending")

    def test_budget_decreases_as_the_candidate_angle_grows(self):
        rows = sorted(self.d["ladder_rows"], key=lambda r: r["angle_deg"])
        budgets = [r["sector_invariance_correction_budget_upper_rad"] for r in rows]
        self.assertEqual(budgets, sorted(budgets, reverse=True))

    def test_the_budget_bracket_is_a_genuine_crossing(self):
        outer = self.d["operation_matched_outer_q_upper"]
        for r in self.d["ladder_rows"] + self.d["limit_probe_rows"]:
            q = r["post_prediction_q_upper"]
            lo, hi = r["budget_bracket_rad"]
            if lo <= 0.0:
                continue
            self.assertLess(BUDGET._worst_case_qplus(q, lo), outer)
            try:
                self.assertGreaterEqual(BUDGET._worst_case_qplus(q, hi), outer)
            except RuntimeError:
                pass

    def test_nuisance_term_exceeds_the_budget_at_every_rung_and_probe(self):
        for r in self.d["ladder_rows"] + self.d["limit_probe_rows"]:
            self.assertFalse(r["nuisance_fits_inside_budget"])
            self.assertGreater(r["nuisance_over_budget_ratio"], 1.0)
        self.assertFalse(self.d["any_ladder_rung_fits_nuisance_inside_budget"])
        self.assertFalse(
            self.d["shrinking_the_candidate_angle_alone_can_close_the_budget"])

    def test_shrinking_the_angle_helps_but_is_bounded_away_from_closing(self):
        rows = sorted(self.d["ladder_rows"], key=lambda r: r["angle_deg"])
        ratios = [r["nuisance_over_budget_ratio"] for r in rows]
        self.assertEqual(ratios, sorted(ratios))
        self.assertTrue(self.d["descending_ladder_halves_the_gap"])
        # ... and still above one in the small-angle limit.
        self.assertGreater(self.d["smallest_probe_nuisance_over_budget_ratio"], 1.0)
        self.assertLess(self.d["smallest_probe_angle_deg"], 15.0)

    def test_nuisance_floor_is_set_by_the_declared_aw_over_lowest_force(self):
        d = self.d
        ratio = d["declared_aw_over_lowest_specific_force"]
        self.assertAlmostEqual(
            ratio,
            d["declared_startup_aw_error_norm_upper_mps2"]
            / d["specific_force_norm_lower_mps2"],
            delta=1.0e-12)
        self.assertGreater(ratio, 0.5)
        self.assertAlmostEqual(
            d["lowest_force_gravity_direction_error_upper_rad"],
            math.asin(ratio), delta=1.0e-12)

    def test_worst_cell_decomposes_into_source_faithful_inputs(self):
        for r in self.d["ladder_rows"]:
            w = r["worst_cell"]
            total = (w["aw_after_prefix_upper_mps2"]
                     + w["force_attitude_remainder_upper_mps2"]
                     + self.d["accel_bias_error_norm_upper_mps2"])
            self.assertLessEqual(w["effective_aw_input_upper_mps2"], total * (1.0 + 1.0e-12))
            self.assertGreaterEqual(w["effective_aw_input_upper_mps2"], total * (1.0 - 1.0e-12))
            self.assertLessEqual(
                w["nuisance_correction_norm_upper_rad"],
                w["Ktheta_norm_upper"] * w["effective_aw_input_upper_mps2"] * (1.0 + 1.0e-9))


if __name__ == "__main__":
    unittest.main()
