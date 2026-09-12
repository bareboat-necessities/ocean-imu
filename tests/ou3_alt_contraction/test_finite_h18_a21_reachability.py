import unittest
from fractions import Fraction as F
from tools.stability.ou3_alt_contraction import finite_h18_a21_reachability as X

class Tests(unittest.TestCase):
    def test_unlock_is_forced_within_10s(self):
        r=X.reach(external_hold=False)
        self.assertEqual(r.mag_updates_applied,250)
        self.assertEqual(r.live_to_unlock_upper,F(10))
        self.assertTrue(r.accel_bias_lock_forced_clear)
        self.assertTrue(r.A21_enable_forced)

    def test_external_hold_does_not_block_lock_clear_but_blocks_immediate_A21(self):
        r=X.reach(external_hold=True)
        self.assertTrue(r.accel_bias_lock_forced_clear)
        self.assertFalse(r.A21_enable_forced)

    def test_readiness_keeps_same_history_edge_open(self):
        r=X.readiness()
        self.assertTrue(r['qualified_mag_schedule_to_literal_unlock_guard_attached'])
        self.assertTrue(r['accel_bias_lock_forced_clear_within_10s_of_gauged_Live'])
        self.assertFalse(r['H18_A21_same_history_covariance_edge_composed_with_reachability'])
        self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
