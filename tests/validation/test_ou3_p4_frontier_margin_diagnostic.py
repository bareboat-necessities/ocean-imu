from pathlib import Path
import os, sys, unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools'))

import ou3_p4_frontier_margin_diagnostic as D


@unittest.skipUnless(
    os.environ.get('OU3_RUN_OBSOLETE_P4_FRONTIER') == '1',
    'retired microscopic P4 frontier regression runs only in its non-gating focused job',
)
class P4FrontierMarginDiagnosticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = D.build()

    def test_diagnostic_passes_without_promoting_tiny_funnel(self):
        self.assertEqual(D.validate(self.d), [])
        self.assertEqual(self.d['P4_USABLE_CERTIFICATE_STATUS'], 'NOT_ESTABLISHED')
        for mode in ('H', 'A'):
            m = self.d['modes'][mode]
            self.assertFalse(m['current_scalar_small_gain_route_usable'])
            self.assertGreater(m['p3_word_endpoint_delta_lower'], 0.0)
            self.assertGreater(m['p3_cell_count'], 0)
            self.assertGreaterEqual(m['p3_best_to_worst_delta_ratio'], 1.0)

    def test_worst_cell_is_explicit(self):
        for mode in ('H', 'A'):
            w = self.d['modes'][mode]['p3_worst_cell']
            self.assertEqual(w['mode'], mode)
            self.assertGreater(w['delta_full_lower'], 0.0)
            self.assertIn(w['limiting_block'], ('translation_RL_inverse_block', 'attitude_bias_or_active_ba_block'))


if __name__ == '__main__':
    unittest.main()
