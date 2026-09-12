"""Deployment-lifetime obstruction regressions; not instability evidence."""
import unittest

from tools.stability.ou3_alt_contraction import finite_indefinite_runtime_obstruction as X


class Tests(unittest.TestCase):
    def test_current_shipping_exposes_both_indefinite_lifetime_blockers(self):
        r=X.build()
        self.assertTrue(r.wrapper_clock_finite_prefix_closed)
        self.assertTrue(r.exact_clock_stall_witness_present)
        self.assertFalse(r.wrapper_clock_indefinite_closed)
        self.assertTrue(r.signed_mag_counter_increment_present)
        self.assertFalse(r.signed_mag_counter_saturation_present)
        self.assertFalse(r.mag_schedule_uniform_call_upper_present)
        self.assertFalse(r.signed_mag_counter_lifetime_closed)
        self.assertFalse(r.indefinite_current_shipping_execution_closed)

    def test_indefinite_promotion_guard_fails_closed(self):
        with self.assertRaisesRegex(RuntimeError,'indefinite current-shipping execution'):
            X.assert_indefinite_shipping_execution_closed()

    def test_obstruction_does_not_claim_finite_word_or_dynamical_failure(self):
        r=X.readiness()
        self.assertFalse(r['finite_word_contraction_invalidated_by_this_obstruction'])
        self.assertFalse(r['dynamical_instability_claimed'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_END_TO_END_PASS'])


if __name__=='__main__': unittest.main()
