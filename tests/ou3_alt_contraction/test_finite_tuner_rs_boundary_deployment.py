"""Dual-compiler R_S pending-boundary staging regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as CAND
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_binary32 as AB
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_machine_real_join as AJ
from tools.stability.ou3_alt_contraction import finite_tuner_rs_boundary_deployment as X
from tools.stability.ou3_alt_contraction import finite_tuner_rs_deployment_ledger as L
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_machine_real_join as T


def pair():
    c=D.shipping_defaults(qeff_pow_result=B.rn32(1))
    def one(tau,sigma,f):
        target=CAND.TargetState(F(f),F(1),B.rn32(F(tau)),B.rn32(F(sigma)))
        tj=T.join(c,target,pow_result=B.rn32(16),sqrt_result=B.rn32(F(1,5)))
        dt=B.rn32(F(1,200)); safe=min(max(target.tau_target,AB.TIME_MIN),AB.TIME_MAX)
        requested=B.mul(c.adapt_RS_mult,safe); lo=min(max(dt,AB.HORIZON_MIN),AB.HORIZON_MAX)
        rssec=min(max(requested,lo),AB.HORIZON_MAX); x=B.div(dt,rssec)
        elo,ehi,_,_=EXP.enclosure(x); e=B.rn32((elo+ehi)/2)
        return tj,AJ.join(c,tau_target=target.tau_target,dt=dt,exp_decay=e)
    return one(F(5,2),F(9,10),F(1,5)),one(F(2),F(4,5),F(1,4))


def advanced():
    (st,sa),(ft,fa)=pair(); s=L.initial()
    r=L.step(s,separate_target=st,separate_alpha=sa,fma_target=ft,fma_alpha=fa)
    return s,r


class Tests(unittest.TestCase):
    def test_sample_candidate_is_not_committed_on_same_sample(self):
        initial,result=advanced(); state=X.State(initial,False)
        nxt=X.after_sample(state,result,pending_after=True)
        self.assertEqual(nxt.ledger,result.state)
        self.assertTrue(nxt.pending)
        self.assertNotEqual(nxt.ledger.separate,initial.separate)
        self.assertNotEqual(nxt.ledger.fma,initial.fma)

    def test_following_live_boundary_commits_carried_snapshot(self):
        initial,result=advanced(); carried=X.after_sample(X.State(initial),result,pending_after=True)
        out=X.imu_boundary(carried,live=True,min_RS=B.rn32(F(15,100)),max_RS=B.rn32(100))
        self.assertTrue(out.consumed); self.assertFalse(out.state.pending)
        self.assertEqual(out.separate_commit.stored_RS,carried.ledger.separate)
        self.assertEqual(out.fma_commit.stored_RS,carried.ledger.fma)

    def test_no_pending_boundary_is_no_RS_write(self):
        s=X.begin(); out=X.imu_boundary(s,live=True,min_RS=B.rn32(F(15,100)),max_RS=B.rn32(100))
        self.assertFalse(out.consumed); self.assertIsNone(out.separate_commit); self.assertIsNone(out.fma_commit)
        self.assertEqual(out.state.ledger,s.ledger)

    def test_prelive_pending_clears_here_without_RS_covariance_write(self):
        s=X.State(L.initial(),True)
        out=X.imu_boundary(s,live=False,min_RS=B.rn32(F(15,100)),max_RS=B.rn32(100))
        self.assertFalse(out.consumed); self.assertFalse(out.state.pending)
        self.assertIsNone(out.separate_commit); self.assertIsNone(out.fma_commit)

    def test_mag_hold_preserves_candidate_and_pending_by_identity(self):
        s=X.State(L.State(B.rn32(F(3,5)),B.rn32(F(7,10)),3),True)
        self.assertIs(X.mag_or_hold(s),s)

    def test_detached_sample_result_is_rejected(self):
        initial,result=advanced(); other=L.State(B.rn32(F(3,5)),B.rn32(F(3,5)),0)
        with self.assertRaisesRegex(ValueError,'detached from boundary predecessor'):
            X.after_sample(X.State(other),result,pending_after=True)

    def test_readiness_closes_staging_not_full_common_boundary_or_master(self):
        r=X.readiness()
        for k in ('sample_k_RS_candidate_is_carried_without_same_sample_commit',
                  'pending_RS_is_consumed_only_at_following_IMU_boundary',
                  'following_Live_boundary_commits_carried_separate_and_FMA_snapshots',
                  'MAG_and_HOLD_are_literal_RS_ledger_and_pending_identity',
                  'preLive_RS_boundary_performs_no_RS_covariance_write'):
            self.assertTrue(r[k])
        for k in ('full_common_tuner_pending_transaction_attached','active_MEKF_RS_Eigen_correspondence_closed',
                  'startup_frontend_RS_machine_history_attached','Live_600_step_RS_machine_history_attached',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
