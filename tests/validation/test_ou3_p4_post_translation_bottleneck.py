from copy import deepcopy
from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_post_translation_bottleneck as B


class P4PostTranslationBottleneckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Unit-test the post-translation algebra without recomputing the rigorous
        # one-second complete-word enclosure.  The focused workflow runs the
        # actual rigorous producer and feeds its artifact to this diagnostic.
        cls.translation = {
            "P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS": "PASS",
            "modes": {
                "H": {"complete_word_translation_margin_lower": 9.676620313503055e-25},
                "A": {"complete_word_translation_margin_lower": 9.676620313503055e-25},
            },
        }
        cls.d = B.build(cls.translation)

    def test_validates(self):
        self.assertEqual(B.validate(self.d), [])

    def test_is_fail_closed(self):
        self.assertEqual(self.d["P4_USABLE_CERTIFICATE_STATUS"], "NOT_ESTABLISHED")
        self.assertFalse(self.d["blockwise_min_is_final_certificate"])
        self.assertFalse(self.d["cross_block_budget_is_final_certificate"])
        for mode in ("H", "A"):
            self.assertFalse(self.d["modes"][mode]["full_state_cross_block_bound_validated"])
            self.assertFalse(self.d["modes"][mode]["full_state_complete_word_cross_blocks_propagated"])
            self.assertFalse(self.d["modes"][mode]["usable_P4_promoted"])

    def test_identifies_positive_post_translation_margin_and_safe_cross_block_target(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            a = m["validated_complete_word_translation_margin_lower"]
            b = m["existing_direct_nontranslation_margin_lower"]
            self.assertGreater(a, 0.0)
            self.assertGreater(b, 0.0)
            self.assertEqual(m["diagnostic_blockwise_margin_lower"], min(a, b))
            self.assertGreater(m["diagnostic_widening_vs_old_full_margin_lower"], 0.0)
            self.assertIn(
                m["post_translation_limiting_block"],
                ("translation_complete_word", "nontranslation_existing_direct"),
            )
            budget = m["normalized_full_state_cross_block_spectral_norm_budget_lower_open"]
            self.assertTrue(m["cross_block_budget_outward_lower_enclosed"])
            self.assertEqual(
                budget,
                m["normalized_full_state_cross_block_spectral_norm_budget_open"],
            )
            self.assertGreater(budget, 0.0)
            # It is a conservative lower enclosure of the mathematical Schur
            # threshold, never a rounded-up surrogate.
            self.assertLessEqual(budget, math.sqrt(a * b))

    def test_derived_field_mutations_fail_validation(self):
        mutations = (
            ("diagnostic_blockwise_margin_lower", 1.0),
            ("diagnostic_widening_vs_old_full_margin_lower", 1.0),
            ("post_translation_limiting_block", "wrong"),
            ("translation_still_limits_after_full_word_widening", "wrong"),
            ("normalized_full_state_cross_block_spectral_norm_budget_lower_open", 1.0),
            ("normalized_full_state_cross_block_spectral_norm_budget_open", 1.0),
            ("cross_block_budget_outward_lower_enclosed", False),
        )
        for mode in ("H", "A"):
            for key, value in mutations:
                with self.subTest(mode=mode, key=key):
                    d = deepcopy(self.d)
                    d["modes"][mode][key] = value
                    self.assertNotEqual(B.validate(d), [])

    def test_nonpass_translation_input_is_rejected(self):
        d = B.build({"P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS": "NOT_ESTABLISHED", "modes": {}})
        self.assertNotEqual(B.validate(d), [])


if __name__ == "__main__":
    unittest.main()
