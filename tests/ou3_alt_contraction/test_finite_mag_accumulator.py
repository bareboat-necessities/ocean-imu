"""Default unweighted MagAutoTuner accumulator regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_accumulator as X


class Tests(unittest.TestCase):
    def test_accepted_sample_updates_same_statistics_with_unit_weight(self):
        s=X.State(); cfg=X.Config(min_samples=2,min_window=F(1,100),sample_dt=F(1,200))
        a=X.accepted_sample(s,cfg,world_sample=(2,3,4),world_norm=5,dt=F(1,200))
        self.assertEqual(a.state.world_sum,(2,3,4)); self.assertEqual(a.state.weight_sum,1)
        self.assertEqual(a.state.accepted_count,1); self.assertFalse(a.finalize_due)
        b=X.accepted_sample(a.state,cfg,world_sample=(4,5,6),world_norm=7,dt=F(1,200))
        self.assertEqual(b.state.world_sum,(6,8,10)); self.assertEqual(b.state.weight_sum,2)
        self.assertEqual(b.state.accepted_window,F(1,100)); self.assertTrue(b.finalize_due)
        self.assertEqual(X.mean(b.state),(3,4,5))

    def test_nonpositive_dt_uses_configured_sample_period(self):
        cfg=X.Config(min_samples=1,min_window=0,sample_dt=F(1,200))
        for dt in (0,F(-1,10)):
            out=X.accepted_sample(X.State(),cfg,world_sample=(1,0,0),world_norm=1,dt=dt)
            self.assertEqual(out.state.accepted_window,F(1,200)); self.assertTrue(out.finalize_due)

    def test_timeout_can_finalize_before_min_window_but_not_before_count(self):
        cfg=X.Config(min_samples=2,min_window=10,max_window=1,sample_dt=F(1,2))
        a=X.accepted_sample(X.State(),cfg,world_sample=(1,0,0),world_norm=1,dt=F(1,2))
        self.assertFalse(a.finalize_due)
        b=X.accepted_sample(a.state,cfg,world_sample=(1,0,0),world_norm=1,dt=F(1,2))
        self.assertTrue(b.finalize_due)

    def test_weighted_branch_and_empty_mean_fail_closed(self):
        with self.assertRaises(NotImplementedError):
            X.accepted_sample(X.State(),X.Config(quality_weighting=True),world_sample=(1,0,0),world_norm=1,dt=F(1,200))
        with self.assertRaises(ValueError): X.mean(X.State())

    def test_readiness_does_not_claim_acceptance_or_binary32(self):
        r=X.readiness()
        self.assertTrue(r['accepted_world_sum_count_time_weight_recurrence_materialized'])
        self.assertFalse(r['raw_mag_sample_acceptance_gates_attached'])
        self.assertFalse(r['world_sample_tilt_quaternion_binary32_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
