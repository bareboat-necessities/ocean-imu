from pathlib import Path
import os, sys, unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools'))

import ou3_p4_frontier_combined_certificate as F
import ou3_validate_enclosure as ENC


@unittest.skipUnless(
    os.environ.get('OU3_RUN_OBSOLETE_P4_FRONTIER') == '1',
    'retired microscopic P4 frontier regression runs only in its non-gating focused job',
)
class P4FrontierCombinedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = F.build()

    def test_passes_and_strictly_widens_direct(self):
        self.assertEqual(F.validate(self.d), [])
        self.assertEqual(self.d['P4_FRONTIER_COMBINED_CERTIFICATE'], 'PASS')
        for mode in ('H', 'A'):
            m = self.d['modes'][mode]
            self.assertGreater(m['certified_level_W'], m['certified_level_W_before_frontier'])
            self.assertGreater(m['frontier_W_widening_factor_lower'], 1.0)

    def test_self_consistent_bootstrap_is_tighter_than_four(self):
        for mode in ('H', 'A'):
            m = self.d['modes'][mode]
            self.assertLess(m['frontier_bootstrap_gamma_upper'], 4.0)
            self.assertGreater(m['frontier_B_reduction_factor_lower'], 1.0)

    def test_selected_candidate_closes_prefix_and_strict_endpoint(self):
        for mode in ('H', 'A'):
            m = self.d['modes'][mode]
            best = max(m['frontier_candidates'], key=lambda r: r['W'])
            self.assertEqual(m['certified_level_W'], best['W'])
            self.assertLess(best['defect_ratio_upper'], best['strict_gap'])
            self.assertLessEqual(best['prefix_factor_upper'] ** 2, best['bootstrap_gamma_upper'] * (1.0 + 1e-15))
            self.assertLess(best['qprefix'], m['cayley_norm_limit'])

    def test_schema4_compatibility_and_A_projection(self):
        for mode in ('H', 'A'):
            out = ENC.validate_mode(mode, self.d['modes'][mode], {'required_path_metric': 'CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC'})
            self.assertTrue(out['linear_pass'], out['failures'])
            self.assertTrue(out['nonlinear_pass'], out['failures'])
        self.assertFalse(self.d['modes']['A']['active_bias_projection']['projection_surface_reached_in_certified_funnel'])

    def test_source_only(self):
        text = (ROOT / 'tools' / 'ou3_p4_frontier_combined_certificate.py').read_text().casefold()
        self.assertNotIn('monte carlo', text)
        self.assertNotIn('sampled trajectory', text)
        self.assertNotIn('numpy', text)


if __name__ == '__main__':
    unittest.main()
