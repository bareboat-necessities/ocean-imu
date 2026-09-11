import unittest
from tools.stability.ou3_alt_contraction import hybrid_word as H

class HybridWordTests(unittest.TestCase):
    def test_configured_hybrid_factorization_closes_without_claiming_release_progress(self):
        d=H.build();self.assertEqual(H.validate(d),[])
        self.assertTrue(d['H18_same_mode_word_closed']);self.assertTrue(d['A21_same_mode_word_closed'])
        self.assertTrue(d['both_configured_guard_outcomes_covered']);self.assertTrue(d['H18_A21_edge_attached']);self.assertTrue(d['hybrid_H18_A21_word_closed'])
        self.assertTrue(d['actual_A21_release_covariance_comes_from_H18_held_covariance']);self.assertFalse(d['parallel_hypothetical_A21_covariance_used_at_release'])
        self.assertFalse(d['release_guard_eventual_reachability_closed_here']);self.assertFalse(d['startup_capture_consumed']);self.assertFalse(d['storage_search_allowed_here'])

if __name__=='__main__':unittest.main()
