"""One-step machine versus exact-real R_S EMA residual regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_candidate_interval as I
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_rs_ema_machine_real_join as X
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_machine_real_join as T


def target_join():
    cfg=D.shipping_defaults(qeff_pow_result=B.rn32(1))
    target=C.TargetState(F(1,5),F(1),B.rn32(F(5,2)),B.rn32(F(9,10)))
    return T.join(cfg,target,pow_result=B.rn32(16),sqrt_result=B.rn32(F(1,5)))


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

    def test_equal_alpha_is_allowed_but_not_assumed(self):
        prev=I.IntervalTuneState(F(1),F(1),F(1,2),F(1,2))
        a=B.rn32(F(1,100))
        out=X.join(prev,B.rn32(F(1,2)),target_join(),exact_alpha=a,machine_alpha=a)
        self.assertEqual(out.alpha_supply,0)
        self.assertLessEqual(out.next_residual_lo,out.next_residual_hi)

    def test_nonbinary_machine_predecessor_fails_closed(self):
        prev=I.IntervalTuneState(F(1),F(1),F(1,2),F(1,2))
        with self.assertRaisesRegex(ValueError,'machine predecessor'):
            X.join(prev,F(1,3),target_join(),exact_alpha=F(1,100),machine_alpha=B.rn32(F(1,100)))

    def test_readiness_keeps_source_uniform_alpha_and_commit_open(self):
        r=X.readiness()
        for k in ('exact_RS_interval_EMA_image_materialized','binary32_RS_EMA_step_composed_with_same_target_join',
                  'machine_minus_exact_predecessor_RS_residual_carried','machine_alpha_minus_exact_alpha_supply_carried',
                  'machine_minus_exact_next_RS_residual_interval_exposed'):
            self.assertTrue(r[k])
        for k in ('alpha_RS_target_libm_correspondence_closed','source_uniform_target_and_alpha_supply_bounds_closed',
                  'binary32_RS_commit_correspondence_closed','source_uniform_complete_600_step_word_qualified',
                  'storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
