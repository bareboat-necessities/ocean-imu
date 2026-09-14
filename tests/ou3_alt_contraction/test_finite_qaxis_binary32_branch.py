"""Exact binary32 Qaxis coefficient-branch regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_qaxis_binary32_branch as Q


class Tests(unittest.TestCase):
    def test_small_and_general_branches_are_both_reachable(self):
        self.assertTrue(Q.branch(1,F(1,200)).small)
        self.assertFalse(Q.branch(F(1,4),F(1,200)).small)

    def test_branch_carries_same_rounded_tau_h_operands(self):
        b=Q.branch(F(11,10),F(1,200))
        self.assertEqual(b.tau_float,Q.rn32_pos(F(11,10)))
        self.assertEqual(b.h_float,Q.rn32_pos(F(1,200)))
        self.assertEqual(b.x,Q.mul32_pos(b.h_float,Q.div32_pos(Q.ONE,b.tau_eff)))
        self.assertEqual(b.small,b.x<Q.SMALL_THRESHOLD)

    def test_detached_branch_claim_is_rejected(self):
        b=Q.branch(1,F(1,200))
        with self.assertRaisesRegex(ValueError,'comparison'):
            Q.Branch(b.tau_input,b.h_input,b.tau_float,b.h_float,b.tau_eff,b.inv,b.x,not b.small)

    def test_shipping_source_shape_and_fail_closed_readiness(self):
        r=Q.readiness()
        self.assertTrue(r['shipping_nested_and_final_Qaxis_share_literal_tau_h_branch_shape'])
        self.assertTrue(r['small_general_comparison_binary32_attached'])
        self.assertFalse(r['source_owned_tau_deployment_commit_correspondence_closed'])
        self.assertFalse(r['coefficient_formula_binary32_roundoff_closed'])
        self.assertFalse(r['Eigen_PSD_hygiene_deployment_closed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
