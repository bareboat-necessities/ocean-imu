"""Deployed periodic a_w covariance synchronization regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_periodic_aw_sync as X

class Tests(unittest.TestCase):
    def test_strict_cadence_predicate_matches_shipping(self):
        s=X.State(False,0)
        at=X.tick(s,time=F(1,10),adapt_every=F(1,10),live=True)
        self.assertFalse(at.requested_now); self.assertIs(at.state,s)
        after=X.tick(s,time=F(101,1000),adapt_every=F(1,10),live=True)
        self.assertTrue(after.requested_now); self.assertTrue(after.state.pending)
        self.assertEqual(after.state.last_sync_time,F(101,1000))

    def test_nonlive_or_disabled_never_queue(self):
        s=X.State(False,0)
        self.assertFalse(X.tick(s,time=1,adapt_every=F(1,10),live=False).state.pending)
        self.assertFalse(X.tick(s,time=1,adapt_every=F(1,10),live=True,enabled=False).state.pending)

    def test_existing_pending_survives_until_prediction_and_clock_updates_only_on_new_tick(self):
        s=X.State(True,F(1,5))
        quiet=X.tick(s,time=F(1,4),adapt_every=F(1,10),live=True)
        self.assertTrue(quiet.state.pending); self.assertEqual(quiet.state.last_sync_time,F(1,5))
        consumed=X.consume_at_prediction(quiet.state)
        self.assertFalse(consumed.pending); self.assertEqual(consumed.last_sync_time,F(1,5))

    def test_nondefault_immediate_policies_fail_closed(self):
        for kw in ({'legacy':True},{'congruent':True}):
            with self.assertRaises(NotImplementedError):
                X.tick(X.State(),time=1,adapt_every=F(1,10),live=True,**kw)

    def test_readiness_is_explicit_about_unrepresented_policy_branches(self):
        r=X.readiness()
        self.assertTrue(r['deployed_periodic_queue_predicate_materialized'])
        self.assertFalse(r['legacy_immediate_replacement_branch_attached'])
        self.assertFalse(r['congruent_immediate_sync_branch_attached'])
        self.assertFalse(r['clock_binary64_roundoff_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
