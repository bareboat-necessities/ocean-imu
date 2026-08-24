import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_implementation_word_language as IMPL_WORDS


class Ou3ImplementationWordLanguageTests(unittest.TestCase):
    def test_declared_domain_closes_conditional_word_language(self):
        d = IMPL_WORDS.build()
        self.assertEqual(IMPL_WORDS.validate(d), [])
        self.assertTrue(d["pass"])
        self.assertTrue(d["source_complete_relative_to_declared_theorem_hypotheses"])
        word = d["word_contract"]
        self.assertTrue(word["conditional_word_language"]["ready"])
        self.assertTrue(word["source_complete_relative_to_theorem_hypotheses"])
        self.assertEqual(word["vector_persistent_excitation"]["recurrence_window_s"], 1.0)
        self.assertGreater(word["conditional_word_language"]["word_samples_upper_at_configured_dt"], 0)
        tr = word["translation_recurrence"]
        self.assertEqual(tr["primary_route"], "FOUR_S_SPREAD_COMPLETE_V_P_S_AW_UCO")
        self.assertEqual(tr["aligned_firing_count"], 4)
        self.assertTrue(tr["three_firing_integrator_detectability_is_supporting_only"])
        self.assertFalse(word["conditional_word_language"]["one_sample_decrease_required"])
        self.assertTrue(word["source_branch_language"]["joint_source_reachability_required"])

    def test_language_stage_does_not_claim_enclosure(self):
        d = IMPL_WORDS.build()
        self.assertFalse(d["continuous_word_enclosed"])
        self.assertFalse(d["nonlinear_word_enclosed"])
        self.assertEqual(d["theorem_promotion"], "NOT_ESTABLISHED")


if __name__ == "__main__":
    unittest.main()
