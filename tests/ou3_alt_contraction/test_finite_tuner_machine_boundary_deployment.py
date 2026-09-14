"""Common machine TuneState pending-boundary regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as T
from tools.stability.ou3_alt_contraction import finite_tuner_machine_boundary_deployment as X


class Tests(unittest.TestCase):
    def pending(self):
        s=T.initial()
        return T.State(s.tau,s.sigma,s.rs,True)

    def test_no_pending_is_literal_identity_and_consumes_nothing(self):
        s=T.initial(); out=X.imu_boundary(s,live=True,min_RS=B.rn32(F(15,100)),max_RS=B.rn32(100))
        self.assertFalse(out.consumed); self.assertEqual(out.state,s)
        self.assertIsNone(out.separate); self.assertIsNone(out.fma)

    def test_live_pending_commits_whole_carried_snapshot(self):
        s=self.pending(); out=X.imu_boundary(s,live=True,min_RS=B.rn32(F(15,100)),max_RS=B.rn32(100))
        self.assertTrue(out.consumed); self.assertFalse(out.state.pending)
        self.assertEqual(out.separate.tau,s.tau.separate)
        self.assertEqual(out.separate.sigma,s.sigma.separate)
        self.assertEqual(out.separate.RS,s.rs.separate)
        self.assertEqual(out.separate.rs_commit.stored_RS,s.rs.separate)
        self.assertEqual(out.fma.tau,s.tau.fma)
        self.assertEqual(out.fma.sigma,s.sigma.fma)
        self.assertEqual(out.fma.RS,s.rs.fma)
        self.assertEqual(out.fma.rs_commit.stored_RS,s.rs.fma)

    def test_prelive_pending_clears_common_transaction_without_RS_write(self):
        s=self.pending(); out=X.imu_boundary(s,live=False,min_RS=B.rn32(F(15,100)),max_RS=B.rn32(100))
        self.assertTrue(out.consumed); self.assertFalse(out.state.pending)
        self.assertIsNone(out.separate.rs_commit); self.assertIsNone(out.fma.rs_commit)
        self.assertEqual((out.separate.tau,out.separate.sigma,out.separate.RS),
                         (s.tau.separate,s.sigma.separate,s.rs.separate))

    def test_mag_and_hold_preserve_pending_whole_state(self):
        s=self.pending(); self.assertIs(X.mag_or_hold(s),s)

    def test_crossed_or_detached_snapshot_fails_closed(self):
        s=self.pending(); out=X.imu_boundary(s,live=True,min_RS=B.rn32(F(15,100)),max_RS=B.rn32(100))
        with self.assertRaisesRegex(ValueError,'tau commit detached'):
            X.Boundary(out.state,
                X.Snapshot('separate',B.rn32(2),out.separate.sigma,out.separate.RS,out.separate.rs_commit),
                out.fma,True,True)

    def test_readiness_closes_common_boundary_not_master_word(self):
        r=X.readiness()
        for k in ('full_tau_sigma_RS_machine_TuneState_pending_transaction_composed',
                  'pending_consumed_only_at_following_IMU_boundary',
                  'tau_and_sigma_commit_from_same_carried_machine_snapshot',
                  'Live_RS_commit_uses_same_carried_machine_snapshot',
                  'preLive_pending_clears_without_RS_covariance_write',
                  'separate_and_FMA_histories_remain_global_and_coherent',
                  'MAG_and_HOLD_are_literal_machine_TuneState_identity'):
            self.assertTrue(r[k])
        for k in ('Eigen_set_RS_noise_execution_correspondence_closed',
                  'startup_frontend_machine_TuneState_product_attached',
                  'Live_600_step_machine_TuneState_product_attached',
                  'source_uniform_complete_600_step_word_qualified',
                  'storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
