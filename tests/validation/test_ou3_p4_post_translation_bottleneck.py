from pathlib import Path
import json, math, sys, unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools'))

import ou3_p4_post_translation_bottleneck as B
import ou3_p4_translation_full_word_interval as T


class P4PostTranslationBottleneckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.translation = T.build()
        cls.d = B.build(cls.translation)

    def test_validates(self):
        self.assertEqual(B.validate(self.d), [])

    def test_is_fail_closed(self):
        self.assertEqual(self.d['P4_USABLE_CERTIFICATE_STATUS'], 'NOT_ESTABLISHED')
        self.assertFalse(self.d['blockwise_min_is_final_certificate'])
        self.assertFalse(self.d['cross_block_budget_is_final_certificate'])
        for mode in ('H', 'A'):
            self.assertFalse(self.d['modes'][mode]['full_state_cross_block_bound_validated'])
            self.assertFalse(self.d['modes'][mode]['full_state_complete_word_cross_blocks_propagated'])
            self.assertFalse(self.d['modes'][mode]['usable_P4_promoted'])

    def test_identifies_a_positive_post_translation_margin_and_cross_block_target(self):
        for mode in ('H', 'A'):
            m = self.d['modes'][mode]
            self.assertGreater(m['validated_complete_word_translation_margin_lower'], 0.0)
            self.assertGreater(m['existing_direct_nontranslation_margin_lower'], 0.0)
            self.assertGreater(m['diagnostic_blockwise_margin_lower'], 0.0)
            self.assertGreater(m['diagnostic_widening_vs_old_full_margin_lower'], 1.0)
            self.assertIn(m['post_translation_limiting_block'], (
                'translation_complete_word', 'nontranslation_existing_direct'))
            expected = math.sqrt(
                m['validated_complete_word_translation_margin_lower'] *
                m['existing_direct_nontranslation_margin_lower'])
            self.assertEqual(expected, m['normalized_full_state_cross_block_spectral_norm_budget_open'])
            self.assertGreater(m['normalized_full_state_cross_block_spectral_norm_budget_open'], 0.0)


if __name__ == '__main__':
    unittest.main()
