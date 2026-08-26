import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_prefix_v2 as G


class Ou3P5Sample1PrefixV2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build(source_pieces=2)

    def test_v2_validates_without_helper_wiring_failure(self):
        d = self.d
        self.assertEqual(G.validate(d), [])
        self.assertTrue(d["angle_conversion_wiring_fixed"])
        self.assertTrue(d["same_SO3_triangle_cayley_formula_as_first_post_reset"])
        self.assertGreater(d["evaluated_sample1_paths"], 0)
        if d["first_failure"] is not None:
            self.assertNotIn("AttributeError", d["first_failure"]["reason"])

    def test_status_remains_fail_closed(self):
        d = self.d
        self.assertIn(d["P5_SAMPLE1_S_ACCEL_PREFIX_CERTIFICATE"], ("PASS", "NOT_ESTABLISHED"))
        if d["P5_SAMPLE1_S_ACCEL_PREFIX_CERTIFICATE"] == "NOT_ESTABLISHED":
            self.assertIsNotNone(d["first_failure"])

    def test_no_gate_relaxation_or_word_promotion(self):
        d = self.d
        self.assertEqual(d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(d["deployed_correction_limit_increased"])
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])


if __name__ == "__main__":
    unittest.main()
