"""Separate nested/final Qaxis covariance exp-root regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_qaxis_binary32_branch as QB
from tools.stability.ou3_alt_contraction import finite_qaxis_exp_binary32 as Q


class Tests(unittest.TestCase):
    def test_general_branch_retains_two_independent_binary32_exp_results(self):
        b=QB.branch(F(1,4),F(1,200))
        self.assertFalse(b.small)
        m=B.rn32(F(9801,10000)); f=B.rn32(F(19603,20000))
        p=Q.pair(b,marginal=m,final=f)
        self.assertEqual(p.marginal,m); self.assertEqual(p.final,f)
        self.assertNotEqual(p.marginal,p.final)

    def test_small_branch_rejects_ghost_exp_witnesses(self):
        b=QB.branch(1,F(1,200)); self.assertTrue(b.small)
        with self.assertRaisesRegex(ValueError,'executes no covariance'):
            Q.pair(b,marginal=B.rn32(F(995,1000)),final=B.rn32(F(995,1000)))

    def test_detached_general_exp_is_rejected(self):
        b=QB.branch(F(1,4),F(1,200))
        with self.assertRaisesRegex(ValueError,'nested Qaxis exp detached'):
            Q.pair(b,marginal=B.rn32(F(9,10)),final=B.rn32(F(9801,10000)))

    def test_readiness_keeps_libm_fail_closed(self):
        r=Q.readiness()
        self.assertTrue(r['shipping_nested_and_final_Qaxis_covariance_exp_calls_present'])
        self.assertTrue(r['Qaxis_nested_exp_retained_separately_from_OU_transition_exp'])
        self.assertTrue(r['Qaxis_final_exp_retained_separately_from_OU_transition_exp'])
        self.assertFalse(r['Qaxis_nested_and_final_exp_bit_identity_assumed'])
        self.assertTrue(r['both_Qaxis_exp_results_bound_to_same_binary32_x_real_enclosure'])
        self.assertFalse(r['Qaxis_exp_libm_binary32_correspondence_closed'])
        self.assertFalse(r['ALT_LIVE_PASS'])

    def test_general_branch_roundoff_budget_uses_compiled_threshold(self):
        r=Q.general_branch_exp_error_budget()
        self.assertEqual(r['argument_interval'][0],QB.SMALL_THRESHOLD)
        self.assertLess(QB.SMALL_THRESHOLD,F(1,100))
        self.assertEqual(r['sufficient_absolute_exp_error'],F(1,2**24))
        self.assertGreater(r['lower_boundary_slack_after_error'],0)
        self.assertGreater(r['upper_boundary_slack_after_error'],0)
        self.assertFalse(r['target_libm_satisfies_budget_proved'])


if __name__=='__main__': unittest.main()
