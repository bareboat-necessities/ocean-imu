from fractions import Fraction as F
import unittest
from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as X

class Tests(unittest.TestCase):
    def test_exact_binary32_5ms_literal(self):
        self.assertEqual(X.DT_FLOAT,F(5368709,1073741824))
        self.assertLess(X.DT_FLOAT,F(1,200))

    def test_timeout_plus_one_word_is_exhaustively_strictly_advancing(self):
        r=X.build()
        self.assertTrue(r.shipping_source_shape_matches)
        self.assertTrue(r.all_updates_strictly_advance)
        self.assertEqual(r.startup_timeout_clock,F(9829793,65536))
        self.assertEqual(r.latest_word_end_clock,F(10026593,65536))

    def test_all_prefix_grid_and_three_second_elapsed_errors_are_bounded_exactly(self):
        r=X.build()
        self.assertEqual(r.max_absolute_grid_error,F(50319,1638400))
        self.assertEqual(r.max_absolute_grid_error_step,25607)
        self.assertEqual(r.max_three_second_elapsed_error,F(3,1024))
        self.assertEqual(r.max_three_second_elapsed_error_start_step,25607)
        self.assertLess(r.max_three_second_elapsed_error,F(1,250))  # < 4 ms

    def test_exact_late_time_stall_witness_blocks_indefinite_clock_claim(self):
        self.assertEqual(X.STALL_WITNESS_TIME,F(1<<17))
        self.assertEqual(X.STALL_WITNESS_NEXT,X.STALL_WITNESS_TIME)
        r=X.readiness()
        self.assertTrue(r['late_time_2pow17_stall_witness_present'])
        self.assertEqual(r['late_time_stall_witness_s'],F(1<<17))
        self.assertFalse(r['indefinite_wrapper_clock_lifetime_closed'])

    def test_readiness_closes_only_canonical_prefix_not_indefinite_lifetime(self):
        r=X.readiness()
        self.assertTrue(r['shipping_outer_clock_is_float_and_incremented_by_dt'])
        self.assertTrue(r['all_binary32_clock_updates_strictly_advance_through_timeout_plus_one_word'])
        self.assertTrue(r['canonical_5ms_wrapper_clock_prefix_binary32_closed'])
        self.assertFalse(r['arbitrary_dt_wrapper_clock_closed'])
        self.assertFalse(r['indefinite_wrapper_clock_lifetime_closed'])
        self.assertFalse(r['magnetic_counter_lifetime_closed'])
        self.assertFalse(r['ALT_LIVE_PASS'])
        self.assertFalse(r['ALT_END_TO_END_PASS'])

if __name__=='__main__': unittest.main()
