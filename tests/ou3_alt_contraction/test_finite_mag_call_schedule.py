from fractions import Fraction as F
import unittest
from tools.stability.ou3_alt_contraction import finite_mag_call_schedule as X

class Tests(unittest.TestCase):
    def test_default_25hz_schedule_forces_250_call_guard_within_10s(self):
        r=X.release_reachability()
        self.assertEqual(r.unlock_count,250)
        self.assertEqual(r.first_call_after_live_max,F(1,25))
        self.assertEqual(r.elapsed_first_to_unlock_max,F(249,25))
        self.assertEqual(r.live_to_unlock_max,F(10))
        self.assertTrue(r.strict_one_second_guard_satisfied)
        self.assertTrue(X.require_release_reachable(r))

    def test_schedule_is_explicit_assumption_and_fail_closed(self):
        with self.assertRaisesRegex(ValueError,'wrong magnetometer'): X.Schedule(F(1,25),F(1,25),'OTHER')
        slow=X.Schedule(F(1,2),F(1,1000))
        r=X.release_reachability(slow,2)
        self.assertFalse(r.strict_one_second_guard_satisfied)
        with self.assertRaisesRegex(ValueError,'one-second'): X.require_release_reachable(r)

    def test_readiness_does_not_claim_complete_edge(self):
        r=X.readiness()
        self.assertTrue(r['named_deterministic_mag_call_schedule_declared'])
        self.assertTrue(r['strict_one_second_guard_forced'])
        self.assertFalse(r['external_acc_bias_hold_excluded_here'])
        self.assertFalse(r['H18_A21_complete_word_edge_attached'])
        self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
