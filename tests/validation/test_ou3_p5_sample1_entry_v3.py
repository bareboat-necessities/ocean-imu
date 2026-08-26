import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_entry_v3 as G


class Ou3P5Sample1EntryV3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build(source_pieces=2)

    def test_linear_gain_tightened_entry_validates(self):
        d = self.d
        self.assertEqual(G.validate(d), [])
        self.assertEqual(d["P5_SAMPLE1_ENTRY_CERTIFICATE"], "PASS")
        self.assertEqual(d["evaluated_source_phase_cells"], d["expected_source_phase_cells"])
        self.assertIsNone(d["first_failure"])

    def test_first_prefix_linear_structure_is_retained(self):
        d = self.d
        self.assertTrue(d["first_prefix_attitude_linear_cross_exact_zero"])
        self.assertTrue(d["first_prefix_linear_covariance_axis_isotropic"])
        self.assertTrue(d["first_due_S_preserves_linear_axis_isotropy"])
        self.assertEqual(
            d["first_accel_linear_gain_exact_form"],
            "K_g=p_gaw*S_a^{-1}, g in {v,p,S,a_w}",
        )
        self.assertEqual(d["first_accel_innovation_loewner_floor"], "S_a>=(p_aw+r_acc)I")
        self.assertEqual(d["first_accel_aw_gain_operator_norm_upper"], 1.0)
        self.assertFalse(d["raw_entrywise_linear_Kr_used_as_state_bound"])

    def test_sample1_state_bounds_remain_finite(self):
        d = self.d
        for value in d["max_sample1_state_group_norm_uppers"].values():
            self.assertTrue(math.isfinite(float(value)))
        self.assertTrue(d["sample1_entry_inside_q8"])
        self.assertLess(d["sample1_pre_measurement_cayley_norm_upper"], 8.0)

    def test_no_premature_word_promotion(self):
        d = self.d
        self.assertEqual(d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(d["deployed_correction_limit_increased"])
        self.assertFalse(d["sample1_S_accel_mag_prefix_evaluated_here"])
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])


if __name__ == "__main__":
    unittest.main()
