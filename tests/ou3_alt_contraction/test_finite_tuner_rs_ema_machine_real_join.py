"""One-step machine versus exact-real R_S EMA residual regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_candidate_interval as I
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_binary32 as A
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_machine_real_join as AJ
from tools.stability.ou3_alt_contraction import finite_tuner_rs_ema_machine_real_join as X
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_machine_real_join as T


def target_join():
    cfg=D.shipping_defaults(qeff_pow_result=B.rn32(1))
    target=C.TargetState(F(1,5),F(1),B.rn32(F(5,2)),B.rn32(F(9,10)))
    return T.join(cfg,target,pow_result=B.rn32(16),sqrt_result=B.rn32(F(1,5)))


def exp_for(joined,*,tau=None):
    mult=joined.cfg.adapt_RS_mult; tau=joined.target.tau_target if tau is None else B.rn32(tau)
    dt=B.rn32(F(1,200)); safe=min(max(tau,A.TIME_MIN),A.TIME_MAX)
    requested=B.mul(mult,safe); lo=min(max(dt,A.HORIZON_MIN),A.HORIZON_MAX)
    rssec=min(max(requested,lo),A.HORIZON_MAX); x=B.div(dt,rssec)
    elo,ehi,_,_=EXP.enclosure(x)
    return dt,B.rn32((elo+ehi)/2)


def alpha_for(joined,*,tau=None):
    tau=joined.target.tau_target if tau is None else B.rn32(tau); dt,e=exp_for(joined,tau=tau)
    return A.step(mult=joined.cfg.adapt_RS_mult,tau_target=tau,dt=dt,exp_decay=e)


def alpha_join_for(joined):
    dt,e=exp_for(joined)
    return AJ.join(joined.cfg,tau_target=joined.target.tau_target,dt=dt,exp_decay=e)


class Tests(unittest.TestCase):
    def test_residual_recurrence_keeps_alpha_supply_separate(self):
        prev=I.IntervalTuneState(F(11,10),F(1,2),F(49,100),F(51,100))
        pm=B.rn32(F(1,2)); ae=F(1,100); ma=B.rn32(F(11,1000))
        out=X.join(prev,pm,target_join(),exact_alpha=ae,machine_alpha=ma)
        self.assertEqual(out.previous_residual_lo,pm-prev.RS_hi)
        self.assertEqual(out.previous_residual_hi,pm-prev.RS_lo)
        self.assertEqual(out.alpha_supply,ma-ae)
        self.assertEqual(out.next_residual_lo,out.machine.next-out.exact_next_hi)
        self.assertEqual(out.next_residual_hi,out.machine.next-out.exact_next_lo)

    def test_source_owned_alpha_feeds_same_candidate_recurrence(self):
        joined=target_join(); alpha=alpha_for(joined)
        prev=I.IntervalTuneState(F(11,10),F(1,2),F(49,100),F(51,100)); pm=B.rn32(F(1,2)); ae=F(1,100)
        out=X.join_from_alpha_step(prev,pm,joined,exact_alpha=ae,alpha_step=alpha)
        self.assertIs(out.alpha_source,alpha)
        self.assertEqual(out.machine.alpha,alpha.alpha)
        self.assertEqual(out.alpha_supply,alpha.alpha-ae)

    def test_point_alpha_dual_history_component_join(self):
        joined=target_join(); alpha=alpha_for(joined)
        prev=I.IntervalTuneState(F(11,10),F(1,2),F(49,100),F(51,100)); pm=B.rn32(F(1,2)); ae=F(1,100)
        out=X.join_compiler_modes(prev,pm,joined,exact_alpha=ae,alpha_step=alpha)
        self.assertEqual(out.separate_residual_lo,out.machine.next_separate-out.exact_next_hi)
        self.assertEqual(out.fma_residual_hi,out.machine.next_fma-out.exact_next_lo)

    def test_strong_join_uses_same_source_exact_alpha_interval_and_both_compiler_modes(self):
        joined=target_join(); aj=alpha_join_for(joined)
        prev=I.IntervalTuneState(F(11,10),F(1,2),F(49,100),F(51,100)); pm=B.rn32(F(1,2))
        out=X.join_from_alpha_interval(prev,pm,joined,alpha_join=aj)
        self.assertIs(out.alpha,aj)
        self.assertIs(out.machine.alpha_source,aj.machine)
        self.assertEqual(out.separate_residual_lo,out.machine.next_separate-out.exact_next_hi)
        self.assertEqual(out.separate_residual_hi,out.machine.next_separate-out.exact_next_lo)
        self.assertEqual(out.fma_residual_lo,out.machine.next_fma-out.exact_next_hi)
        self.assertEqual(out.fma_residual_hi,out.machine.next_fma-out.exact_next_lo)
        self.assertLessEqual(out.exact_next_lo,out.exact_next_hi)

    def test_alpha_from_different_tau_candidate_is_rejected(self):
        joined=target_join(); detached=alpha_for(joined,tau=1)
        prev=I.IntervalTuneState(F(1),F(1),F(1,2),F(1,2))
        with self.assertRaisesRegex(ValueError,'SAME candidate tau'):
            X.join_from_alpha_step(prev,B.rn32(F(1,2)),joined,exact_alpha=F(1,100),alpha_step=detached)
        with self.assertRaisesRegex(ValueError,'SAME candidate tau'):
            X.join_compiler_modes(prev,B.rn32(F(1,2)),joined,exact_alpha=F(1,100),alpha_step=detached)

    def test_equal_alpha_is_allowed_but_not_assumed(self):
        prev=I.IntervalTuneState(F(1),F(1),F(1,2),F(1,2)); a=B.rn32(F(1,100))
        out=X.join(prev,B.rn32(F(1,2)),target_join(),exact_alpha=a,machine_alpha=a)
        self.assertEqual(out.alpha_supply,0)

    def test_nonbinary_machine_predecessor_fails_closed(self):
        prev=I.IntervalTuneState(F(1),F(1),F(1,2),F(1,2))
        with self.assertRaisesRegex(ValueError,'machine predecessor'):
            X.join(prev,F(1,3),target_join(),exact_alpha=F(1,100),machine_alpha=B.rn32(F(1,100)))

    def test_readiness_keeps_source_uniform_alpha_and_commit_open(self):
        r=X.readiness()
        for k in ('exact_RS_interval_EMA_image_materialized',
                  'exact_RS_EMA_consumes_same_source_alpha_interval_not_free_point',
                  'multi_affine_predecessor_target_alpha_box_image_closed_by_exact_vertices',
                  'binary32_RS_EMA_step_composed_with_same_target_join',
                  'source_owned_binary32_RS_alpha_can_feed_persistent_recurrence','same_candidate_tau_owns_RS_alpha_horizon',
                  'separate_and_FMA_RS_histories_both_joined_to_same_exact_interval','compiler_mode_not_guessed_by_theorem_join',
                  'machine_minus_exact_predecessor_RS_residual_carried','machine_minus_exact_next_RS_residual_interval_exposed'):
            self.assertTrue(r[k])
        for k in ('alpha_RS_target_libm_correspondence_closed','shipping_compiler_FP_contraction_mode_qualified',
                  'source_uniform_target_and_alpha_supply_bounds_closed','binary32_RS_commit_correspondence_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
