from fractions import Fraction as F
import unittest
from tools.stability.ou3_alt_contraction import finite_aw_sync_clock_binary64 as X

class Tests(unittest.TestCase):
    def test_source_and_default_partition_close_on_bounded_prefix(self):
        r=X.readiness()
        self.assertTrue(r['shipping_inner_time_and_aw_sync_source_shape_matches'])
        self.assertTrue(r['binary64_clock_strictly_advances_through_startup_plus_word'])
        self.assertTrue(r['default_due_partition_is_exactly_gap_ge_21_samples'])
        self.assertTrue(r['canonical_aw_sync_binary64_predicate_closed'])
        self.assertLessEqual(r['max_deployed_elapsed_over_20_samples'],r['default_adapt_every_float'])
        self.assertGreater(r['min_deployed_elapsed_over_21_samples'],r['default_adapt_every_float'])
        self.assertFalse(r['storage_search_allowed'])

    def test_exact_tick_qualifier_matches_20_and_21_sample_boundaries(self):
        for start in (0,1600,30000,30579):
            a=X.qualify_exact_tick(F(start+20,200),F(start,200))
            b=X.qualify_exact_tick(F(start+21,200),F(start,200))
            self.assertFalse(a['exact_due']); self.assertFalse(a['deployed_due'])
            self.assertTrue(b['exact_due']); self.assertTrue(b['deployed_due'])

    def test_noncanonical_or_nondefault_inputs_fail_closed(self):
        with self.assertRaises(ValueError): X.qualify_exact_tick(F(1,3),0)
        with self.assertRaises(ValueError): X.qualify_exact_tick(F(1,200),0,F(1,5))

if __name__=='__main__': unittest.main()
