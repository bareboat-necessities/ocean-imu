from fractions import Fraction as F
import unittest
from tools.stability.ou3_alt_contraction import finite_magnetic_wrapper_clock as X

class Tests(unittest.TestCase):
    def test_timestamp_is_exact_descendant_of_canonical_physical_endpoint(self):
        t=X.at_physical_time(F(7))
        self.assertEqual(t.source_step,1400)
        self.assertEqual(t.wrapper_time,F(3670075,524288))
        self.assertTrue(X.outer_mag_delay_passed(t))
        with self.assertRaisesRegex(ValueError,'5 ms physical grid'):
            X.at_physical_time(F(7001,1000))

    def test_refinement_threshold_differs_from_ideal_physical_clock(self):
        at90=X.at_physical_time(F(90))
        self.assertEqual(at90.wrapper_time,F(11795193,131072))
        self.assertLess(at90.wrapper_time,F(90))
        self.assertFalse(X.refinement_due(at90))

        at90005=X.at_physical_time(F(18001,200))
        self.assertFalse(X.refinement_due(at90005))
        at90010=X.at_physical_time(F(9001,100))
        self.assertEqual(at90010.source_step,18002)
        self.assertEqual(at90010.wrapper_time,F(11796503,131072))
        self.assertTrue(X.refinement_due(at90010))

    def test_shipping_elapsed_uses_wrapper_difference_or_literal_fallback(self):
        a=X.at_physical_time(F(90)); b=X.at_physical_time(F(18001,200))
        self.assertEqual(X.shipping_elapsed(b,a.wrapper_time),F(655,131072))
        self.assertEqual(X.shipping_elapsed(a,None),F(1,200))
        self.assertEqual(X.shipping_elapsed(a,a.wrapper_time),F(1,200))
        self.assertEqual(X.shipping_elapsed(a,a.wrapper_time+1),F(1,200))

    def test_readiness_refuses_to_promote_old_single_clock_magnetic_word(self):
        r=X.readiness()
        self.assertTrue(r['physical_and_outer_wrapper_clocks_are_distinct_coordinates'])
        self.assertTrue(r['exact_binary32_wrapper_timestamp_derived_from_canonical_physical_endpoint'])
        self.assertTrue(r['shipping_nonadvancing_clock_fallback_dt_materialized'])
        self.assertTrue(r['canonical_prefix_wrapper_clock_arithmetic_closed'])
        self.assertFalse(r['dual_clock_magnetic_word_composed'])
        self.assertFalse(r['indefinite_wrapper_clock_lifetime_closed'])
        self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
