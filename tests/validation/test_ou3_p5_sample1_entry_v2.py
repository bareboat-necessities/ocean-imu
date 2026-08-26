import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_entry_v2 as G


class Ou3P5Sample1EntryV2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build(source_pieces=2)

    def test_tightened_sample1_entry_validates(self):
        d = self.d
        self.assertEqual(G.validate(d), [])
        self.assertEqual(d["P5_SAMPLE1_ENTRY_CERTIFICATE"], "PASS")
        self.assertEqual(d["evaluated_source_phase_cells"], d["expected_source_phase_cells"])
        self.assertIsNone(d["first_failure"])

    def test_exact_aw_gain_contraction_is_active(self):
        d = self.d
        self.assertEqual(d["first_accel_aw_gain_exact_formula"], "K_aw=p_aw*S^{-1}")
        self.assertEqual(d["first_accel_innovation_loewner_floor"], "S>=p_aw*I+R_acc")
        self.assertEqual(d["first_accel_aw_gain_operator_norm_upper"], 1.0)
        self.assertTrue(d["first_accel_aw_state_correction_norm_bounded_by_residual_norm"])
        self.assertFalse(d["raw_entrywise_aw_Kr_used_as_state_bound"])

    def test_aw_state_hull_is_finite_and_no_longer_catastrophic(self):
        d = self.d
        aw = float(d["max_sample1_state_group_norm_uppers"]["aw"])
        self.assertTrue(math.isfinite(aw))
        # The V1 entrywise K*r hull was 2610.858... m/s^2.  The exact
        # K_aw contraction should remove that dependency blow-up by a wide margin.
        self.assertLess(aw, 500.0)

    def test_no_theorem_or_deployment_gate_is_relaxed(self):
        d = self.d
        self.assertEqual(d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(d["deployed_correction_limit_increased"])
        self.assertTrue(d["sample1_entry_inside_q8"])
        self.assertLess(d["sample1_pre_measurement_cayley_norm_upper"], 8.0)
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])


if __name__ == "__main__":
    unittest.main()
