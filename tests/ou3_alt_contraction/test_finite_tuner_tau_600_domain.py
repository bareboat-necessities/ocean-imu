"""Bounded-startup plus canonical 600-step deployed tau-domain regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_tuner_tau_600_domain as X


class Tests(unittest.TestCase):
    def test_full_startup_plus_word_envelope_stays_inside_local_roundoff_domain(self):
        self.assertTrue(X.assert_induction_closes())
        self.assertGreater(X.lower_after(X.TOTAL_UPDATES),F(0))
        self.assertLess(X.upper_after(X.TOTAL_UPDATES),X.DOMAIN_ABS_MAX)
        self.assertEqual(X.lower_after(0),F(2,5))
        self.assertEqual(X.upper_after(0),F(12))
        self.assertEqual(X.TOTAL_UPDATES,30600)

    def test_envelope_is_monotone_conservative_over_startup_and_word(self):
        self.assertGreater(X.lower_after(1),X.lower_after(X.TOTAL_UPDATES))
        self.assertLess(X.upper_after(1),X.upper_after(X.TOTAL_UPDATES))
        with self.assertRaises(ValueError): X.lower_after(X.TOTAL_UPDATES+1)
        with self.assertRaises(ValueError): X.upper_after(-1)

    def test_readiness_closes_tau_domain_only_for_bounded_startup_plus_one_word(self):
        r=X.readiness()
        self.assertTrue(r['shipping_tau_initial_timeout_rate_and_target_domain_source_shape_matches'])
        self.assertTrue(r['shipping_tau_predecessor_domain_inductively_closed_through_startup_and_600_word'])
        self.assertTrue(r['source_uniform_tau_roundoff_supply_bound_closed_through_startup_and_600_word'])
        self.assertFalse(r['successive_words_tau_domain_tiled_indefinitely'])
        self.assertFalse(r['upstream_WPE_binary32_frequency_production_closed'])
        self.assertFalse(r['complete_word_finite_identity'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
