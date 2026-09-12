"""Qualified H18/A21 reachability + hold/release composition regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_h18_a21_transition as X
from tools.stability.ou3_alt_contraction import finite_mag_bias_gate as G
import test_finite_core as FC


class Tests(unittest.TestCase):
    def test_no_hold_reaches_A21_within_10_seconds_and_applies_floor(self):
        cfg=G.Config(); h=FC.root('H')
        out=X.after_qualified_unlock(h,cfg,external_hold=False)
        self.assertEqual(out.mode_after_unlock,'A'); self.assertEqual(out.state.mode,'A')
        self.assertLessEqual(out.live_to_unlock_upper,F(10))
        floor=F(4,1000)**2
        for i in range(18,21): self.assertGreaterEqual(out.state.covariance[i][i],floor)

    def test_hold_may_keep_H18_indefinitely_then_release_exactly_to_A21(self):
        cfg=G.Config(); h=FC.root('H')
        held=X.after_qualified_unlock(h,cfg,external_hold=True)
        self.assertEqual(held.state.mode,'H'); self.assertFalse(held.control.locked)
        self.assertTrue(held.control.hold); self.assertTrue(held.held_continuation_allowed)
        active=X.release_external_hold(held,cfg)
        self.assertEqual(active.state.mode,'A'); self.assertFalse(active.control.hold)

    def test_active_branch_can_be_held_again_without_forcing_eventual_A21(self):
        cfg=G.Config(); h=FC.root('H')
        active=X.after_qualified_unlock(h,cfg,external_hold=False)
        held=X.assert_external_hold(active,cfg)
        self.assertEqual(held.state.mode,'H'); self.assertTrue(held.control.hold)
        for i in range(18,21):
            for j in range(18): self.assertEqual(held.state.covariance[i][j],0)

    def test_wrong_start_mode_fails_closed(self):
        with self.assertRaisesRegex(ValueError,'H18'):
            X.after_qualified_unlock(FC.root('A'),G.Config(),external_hold=False)

    def test_readiness_keeps_clock_and_complete_word_open(self):
        r=X.readiness()
        self.assertTrue(r['no_hold_unlock_to_A21_covariance_edge_composed'])
        self.assertTrue(r['indefinite_external_hold_H18_branch_retained'])
        self.assertTrue(r['eventual_A21_not_assumed_for_arbitrary_hold_history'])
        self.assertFalse(r['deployment_clock_roundoff_closed'])
        self.assertFalse(r['complete_same_history_Live_word']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
