"""Deployed periodic a_w covariance synchronization regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_periodic_aw_sync as X

SIG=((F(1),0,0),(0,F(2),0),(0,0,F(3)))

class Tests(unittest.TestCase):
    def test_strict_cadence_predicate_matches_shipping_and_snapshots_target(self):
        s=X.State(False,0)
        at=X.tick(s,time=F(1,10),adapt_every=F(1,10),live=True)
        self.assertFalse(at.requested_now); self.assertIs(at.state,s)
        after=X.tick(s,time=F(101,1000),adapt_every=F(1,10),live=True,active_sigma=SIG)
        self.assertTrue(after.requested_now); self.assertTrue(after.state.pending)
        self.assertEqual(after.state.last_sync_time,F(101,1000)); self.assertEqual(after.state.target,SIG)

    def test_nonlive_or_disabled_never_queue_or_consume_sigma(self):
        s=X.State(False,0)
        self.assertFalse(X.tick(s,time=1,adapt_every=F(1,10),live=False).state.pending)
        self.assertFalse(X.tick(s,time=1,adapt_every=F(1,10),live=True,enabled=False).state.pending)
        with self.assertRaises(ValueError):
            X.tick(s,time=1,adapt_every=F(1,10),live=False,active_sigma=SIG)

    def test_existing_pending_target_survives_until_prediction(self):
        s=X.State(True,F(1,5),SIG)
        quiet=X.tick(s,time=F(1,4),adapt_every=F(1,10),live=True)
        self.assertTrue(quiet.state.pending); self.assertEqual(quiet.state.target,SIG)
        consumed=X.consume_at_prediction(quiet.state)
        self.assertFalse(consumed.pending); self.assertIsNone(consumed.target); self.assertEqual(consumed.last_sync_time,F(1,5))

    def test_later_active_sigma_cannot_rewrite_already_queued_target(self):
        s=X.State(True,F(1,5),SIG)
        self.assertEqual(X.floor_target(s),SIG)
        # Until another due tick is actually executed, no new active Sigma is accepted.
        with self.assertRaises(ValueError):
            X.tick(s,time=F(1,4),adapt_every=F(1,10),live=True,
                   active_sigma=((9,0,0),(0,9,0),(0,0,9)))

    def test_due_branch_requires_current_active_sigma(self):
        with self.assertRaisesRegex(ValueError,'current active'):
            X.tick(X.State(),time=1,adapt_every=F(1,10),live=True)

    def test_nondefault_immediate_policies_fail_closed(self):
        for kw in ({'legacy':True},{'congruent':True}):
            with self.assertRaises(NotImplementedError):
                X.tick(X.State(),time=1,adapt_every=F(1,10),live=True,**kw)

    def test_readiness_is_explicit_about_unrepresented_policy_branches(self):
        r=X.readiness()
        self.assertTrue(r['queued_floor_target_snapshotted_at_request'])
        self.assertTrue(r['next_boundary_tuner_commit_cannot_rewrite_queued_target'])
        self.assertFalse(r['legacy_immediate_replacement_branch_attached'])
        self.assertFalse(r['congruent_immediate_sync_branch_attached'])
        self.assertFalse(r['clock_binary64_roundoff_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
