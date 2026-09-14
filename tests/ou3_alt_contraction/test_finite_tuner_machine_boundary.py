"""Common machine TuneState boundary regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as TAU
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_deployment_ledger as SIG
from tools.stability.ou3_alt_contraction import finite_tuner_rs_deployment_ledger as RS
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as M
from tools.stability.ou3_alt_contraction import finite_tuner_machine_boundary as X
from tools.stability.ou3_alt_contraction import finite_tuner_commit as C


def cfg():
    # Non-cubic branch isolates the common pending transaction from the still-open
    # sqrt/libm correspondence while retaining tau-derived cadence and Live R_S.
    return C.CommitConfig(F(3,2),F(1,100),2,F(3,200),True,False,
                          F(2,5),5,F(4,5),F(6,5),F(7,5))


def state(*,pending=True):
    # Distinct global compiler values ensure the boundary cannot accidentally
    # cross or collapse the two coherent machine histories.
    return M.State(
      TAU.State(B.rn32(1),B.rn32(F(5,4)),7),
      SIG.State(B.rn32(F(1,10)),B.rn32(F(1,8)),7),
      RS.State(B.rn32(2),B.rn32(F(9,4)),7),
      pending)


class Tests(unittest.TestCase):
    def test_pending_live_boundary_commits_each_whole_track_atomically(self):
        s=state(); out=X.imu_boundary(s,cfg(),live=True,separate_band_noise_floor_sigma=B.rn32(F(1,20)),fma_band_noise_floor_sigma=B.rn32(F(3,50)))
        self.assertTrue(out.consumed); self.assertFalse(out.state.pending)
        self.assertEqual(out.separate_commit.tau,s.tau.separate)
        self.assertEqual(out.fma_commit.tau,s.tau.fma)
        self.assertEqual(out.separate_commit.Sigma_aw[2][2],B.mul(s.sigma.separate,s.sigma.separate))
        self.assertEqual(out.fma_commit.Sigma_aw[2][2],B.mul(s.sigma.fma,s.sigma.fma))
        self.assertNotEqual(out.separate_commit.R_S,out.fma_commit.R_S)

    def test_preLive_pending_consumes_common_transaction_but_writes_no_RS(self):
        out=X.imu_boundary(state(),cfg(),live=False,separate_band_noise_floor_sigma=B.rn32(F(1,20)),fma_band_noise_floor_sigma=B.rn32(F(3,50)))
        self.assertTrue(out.consumed); self.assertFalse(out.state.pending)
        self.assertIsNone(out.separate_commit.R_S); self.assertIsNone(out.fma_commit.R_S)
        self.assertEqual(out.separate_commit.tau,B.rn32(1))
        self.assertEqual(out.fma_commit.tau,B.rn32(F(5,4)))

    def test_nonpending_boundary_is_identity_and_consumes_no_witness(self):
        s=state(pending=False); out=X.imu_boundary(s,cfg(),live=True)
        self.assertFalse(out.consumed); self.assertIs(out.state,s)
        self.assertIsNone(out.separate_commit); self.assertIsNone(out.fma_commit)
        with self.assertRaisesRegex(ValueError,'no-pending boundary consumes no'):
            X.imu_boundary(s,cfg(),live=True,separate_band_noise_floor_sigma=B.rn32(F(1,20)))

    def test_missing_machine_floor_does_not_fall_back_to_exact_shadow(self):
        with self.assertRaisesRegex(ValueError,'each compiler track'):
            X.imu_boundary(state(),cfg(),live=True)
        with self.assertRaisesRegex(ValueError,'each compiler track'):
            X.imu_boundary(state(),cfg(),live=True,separate_band_noise_floor_sigma=B.rn32(F(1,20)))

    def test_machine_floor_witnesses_and_rounded_squares_are_retained(self):
        s=state(); sb=B.rn32(F(3,10)); fb=B.rn32(F(2,5))
        out=X.imu_boundary(s,cfg(),live=True,separate_band_noise_floor_sigma=sb,fma_band_noise_floor_sigma=fb)
        self.assertEqual(out.arithmetic.separate.band_noise_floor_sigma,sb)
        self.assertEqual(out.arithmetic.fma.band_noise_floor_sigma,fb)
        self.assertEqual(out.separate_commit.Sigma_aw[2][2],B.mul(sb,sb))
        self.assertNotEqual(out.separate_commit.Sigma_aw[2][2],sb*sb)
        self.assertEqual(out.fma_commit.Sigma_aw[2][2],B.mul(fb,fb))
        old=C.commit(C.TuneState(s.tau.separate,s.sigma.separate,s.rs.separate),cfg(),
                     pending=True,live=True,band_noise_floor_sigma=sb)
        self.assertNotEqual(out.separate_commit.Sigma_aw,old.Sigma_aw)
        with self.assertRaisesRegex(ValueError,'rounded binary32 outputs'):
            replace(out,separate_commit=old)

    def test_tau_floor_and_sync_use_the_same_rounded_commit(self):
        s=state(); s=replace(s,tau=replace(s.tau,separate=B.rn32(F(1,10000))))
        out=X.imu_boundary(s,cfg(),live=False,separate_band_noise_floor_sigma=B.rn32(F(1,20)),
                           fma_band_noise_floor_sigma=B.rn32(F(1,20)),sync_covariance=True)
        self.assertEqual(out.separate_commit.tau,B.rn32(F(1,1000)))
        self.assertEqual(out.separate_commit.aw_floor_target,out.separate_commit.Sigma_aw)
        self.assertIsNone(out.separate_commit.R_S)

    def test_MAG_and_HOLD_are_whole_product_identity(self):
        s=state(); self.assertIs(X.mag_or_hold(s),s)

    def test_mismatched_scalar_update_counts_cannot_reach_boundary(self):
        with self.assertRaisesRegex(ValueError,'different update counts'):
            M.State(TAU.State(updates=1),SIG.State(updates=2),RS.State(updates=1),True)

    def test_RS_ledger_covers_same_bounded_startup_plus_word_horizon(self):
        self.assertEqual(RS.MAX_UPDATES,TAU.MAX_UPDATES)
        self.assertEqual(RS.MAX_UPDATES,SIG.MAX_UPDATES)
        M.State(TAU.State(updates=30600),SIG.State(updates=30600),RS.State(updates=30600),False)
        with self.assertRaises(ValueError): RS.State(updates=30601)

    def test_readiness_closes_common_boundary_only(self):
        r=X.readiness()
        self.assertTrue(r['next_boundary_common_machine_commit_attached'])
        self.assertTrue(r['applied_parameters_packaged_from_binary32_operation_graph'])
        self.assertTrue(r['one_common_pending_bit_controls_tau_sigma_RS_transaction'])
        self.assertTrue(r['MAG_and_HOLD_preserve_whole_machine_TuneState_and_pending_by_identity'])
        for k in ('compiler_FP_contraction_mode_qualified','exact_commit_to_binary32_shipping_correspondence_closed',
                  'target_libm_correspondence_closed','startup_frontend_machine_TuneState_product_attached',
                  'Live_600_step_machine_TuneState_product_attached','source_uniform_complete_600_step_word_qualified',
                  'storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
