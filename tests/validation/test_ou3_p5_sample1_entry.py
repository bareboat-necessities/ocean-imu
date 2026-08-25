import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_entry as G


class Ou3P5Sample1EntryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build(source_pieces=2)

    def test_sample1_entry_stage_validates(self):
        d = self.d
        self.assertEqual(G.validate(d), [])
        self.assertEqual(d["P5_SAMPLE1_ENTRY_CERTIFICATE"], "PASS")
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])

    def test_exact_sample0_family_is_carried_through_joseph_and_reset(self):
        d = self.d
        self.assertTrue(d["starts_from_exact_source_sample0_family_not_broad_v3_initial_hull"])
        self.assertTrue(d["full_18x18_covariance_propagated"])
        self.assertTrue(d["shipping_Joseph_update_used_for_first_accelerometer"])
        self.assertTrue(d["immediate_left_error_reset_congruence_used"])
        self.assertTrue(d["first_due_S_zero_mean_but_covariance_update_retained"])
        self.assertTrue(d["first_accel_canonical_Jatt_gravity_and_Jaw_identity_used"])
        self.assertTrue(d["source_complete_identity_branch_hull_retained"])

    def test_every_source_phase_cell_reaches_sample1_inside_q8(self):
        d = self.d
        self.assertEqual(d["evaluated_source_phase_cells"], d["expected_source_phase_cells"])
        self.assertGreater(d["evaluated_source_phase_cells"], 0)
        self.assertIsNone(d["first_failure"])
        self.assertTrue(d["sample1_entry_inside_q8"])
        self.assertTrue(math.isfinite(d["sample1_pre_measurement_cayley_norm_upper"]))
        self.assertLess(d["sample1_pre_measurement_cayley_norm_upper"], d["q_chart_target"])

    def test_stage_stops_before_sample1_measurements_and_does_not_promote_word(self):
        d = self.d
        self.assertTrue(d["sample1_entry_is_before_sample1_measurements"])
        self.assertFalse(d["sample1_S_accel_mag_prefix_evaluated_here"])
        self.assertEqual(d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(d["deployed_correction_limit_increased"])
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])
        self.assertEqual(
            d["next_obligation"],
            "CLASSIFY_SAMPLE1_PSEUDO_PHASE_FROM_SAMPLE0_PHASE_AND_EVALUATE_SAMPLE1_S_ACCEL_PREFIX",
        )


if __name__ == "__main__":
    unittest.main()
