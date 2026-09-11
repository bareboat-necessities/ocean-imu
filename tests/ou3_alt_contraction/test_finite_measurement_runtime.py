"""safe-LDLT control-flow regressions; not finite-precision qualification."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_measurement_runtime as R
import test_finite_core as FC


class Tests(unittest.TestCase):
    def test_first_success_uses_no_shift(self):
        b=R.SafeLDLT(True,None,F(3,2),F(1,10**7)); self.assertEqual(b.innovation_shift,0); self.assertTrue(b.accepted)
        s=FC.root('H'); out=R.measurement(s,'S_zero',ldlt=b,R=M.eye(3))
        self.assertTrue(out.accepted); self.assertIsNotNone(out.accepted_graph)
        self.assertEqual(out.accepted_graph.innovation,tuple(map(tuple,out.accepted_graph.innovation)))

    def test_retry_bump_is_literal_max_formula(self):
        b=R.SafeLDLT(False,True,F(2),F(1,10**9)); self.assertEqual(b.bump,F(3,10**6)); self.assertEqual(b.innovation_shift,b.bump)
        s=FC.root('H'); base=R.measurement(s,'S_zero',ldlt=R.SafeLDLT(True,None,2,F(1,10**9)),R=M.eye(3)); retried=R.measurement(s,'S_zero',ldlt=b,R=M.eye(3))
        self.assertTrue(retried.accepted); self.assertNotEqual(base.state.covariance,retried.state.covariance)

    def test_machine_epsilon_can_dominate_bump(self):
        b=R.SafeLDLT(False,True,0,F(1,100)); self.assertEqual(b.bump,F(1,100))

    def test_double_failure_rejects_without_any_update(self):
        s=FC.root('A'); b=R.SafeLDLT(False,False,1,F(1,10**7)); out=R.measurement(s,'S_zero',ldlt=b,R=M.eye(3))
        self.assertFalse(out.accepted); self.assertEqual(out.state,s); self.assertIsNone(out.accepted_graph)

    def test_invalid_branch_declarations_fail_closed(self):
        with self.assertRaises(ValueError): R.SafeLDLT(True,False,1,F(1,10**7))
        with self.assertRaises(TypeError): R.SafeLDLT(False,None,1,F(1,10**7))
        with self.assertRaises(ValueError): R.SafeLDLT(True,None,-1,F(1,10**7))
        with self.assertRaises(ValueError): R.SafeLDLT(True,None,1,0)

    def test_readiness_stays_fail_closed(self):
        r=R.readiness(); self.assertTrue(r['double_LDLT_failure_rejection_branch']); self.assertTrue(r['same_retry_shift_used_by_gain_and_Joseph']); self.assertFalse(r['Eigen_LDLT_outcomes_finite_precision_attached']); self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
