"""Exact post-prediction prefix regressions; not source/stability evidence."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_post_prediction as P
import test_finite_core as FC


def aw_block(state):
    c=list(map(list,state.covariance)); return [r[15:18] for r in c[15:18]]


class Tests(unittest.TestCase):
    def test_not_pending_is_identity_and_consumes_no_solver_witness(self):
        s=FC.root('A'); out=P.aw_floor(s,pending=False,target=M.eye(3))
        self.assertEqual(out.state,s); self.assertFalse(out.pending_after); self.assertIsNone(out.solver_success); self.assertFalse(out.inflation_applied)
        with self.assertRaises(ValueError): P.aw_floor(s,pending=False,target=M.eye(3),solver_success=True,eigenvectors=M.eye(3),eigenvalues=[1,1,1])

    def test_pending_solver_failure_clears_request_without_covariance_change(self):
        s=FC.root('A'); out=P.aw_floor(s,pending=True,target=M.eye(3),solver_success=False)
        self.assertEqual(out.state,s); self.assertFalse(out.pending_after); self.assertFalse(out.solver_success); self.assertFalse(out.inflation_applied)

    def test_success_adds_only_positive_spectral_part_to_aw_block(self):
        s=FC.root('A'); old=aw_block(s)
        delta=[[F(1,5),0,0],[0,F(-1,7),0],[0,0,F(2,9)]]
        target=M.plus(old,delta)
        out=P.aw_floor(s,pending=True,target=target,solver_success=True,eigenvectors=M.eye(3),eigenvalues=[F(1,5),F(-1,7),F(2,9)])
        new=list(map(list,out.state.covariance)); expected=M.plus(old,[[F(1,5),0,0],[0,0,0],[0,0,F(2,9)]])
        self.assertEqual([r[15:18] for r in new[15:18]],expected)
        before=list(map(list,s.covariance))
        for i in range(21):
            for j in range(21):
                if not (15<=i<18 and 15<=j<18): self.assertEqual(new[i][j],before[i][j])
        self.assertTrue(out.inflation_applied); self.assertFalse(out.pending_after)

    def test_eigensystem_must_belong_to_same_delta(self):
        s=FC.root('A'); old=aw_block(s); target=M.plus(old,M.eye(3))
        with self.assertRaises(ValueError): P.aw_floor(s,pending=True,target=target,solver_success=True,eigenvectors=M.eye(3),eigenvalues=[1,1,2])
        bad=[[1,1,0],[0,1,0],[0,0,1]]
        with self.assertRaises(ValueError): P.aw_floor(s,pending=True,target=target,solver_success=True,eigenvectors=bad,eigenvalues=[1,1,1])

    def test_scheduler_preserves_credit_until_due(self):
        q=P.Scheduler(F(1,10),F(3,100))
        due,q=q.step(F(1,100)); self.assertFalse(due); self.assertEqual(q.elapsed,F(1,25))
        due,q=q.step(F(3,50)); self.assertTrue(due); self.assertEqual(q.elapsed,0)

    def test_scheduler_keeps_remainder_after_overshoot(self):
        q=P.Scheduler(F(1,10),F(9,100)); due,q=q.step(F(3,100))
        self.assertTrue(due); self.assertEqual(q.elapsed,F(1,50))

    def test_scheduler_tolerance_can_service_just_before_deadline(self):
        q=P.Scheduler(F(1,10),F(9,100),F(1,1000)); due,q=q.step(F(9,1000))
        self.assertTrue(due); self.assertEqual(q.elapsed,0)

    def test_composed_prefix_orders_floor_then_scheduler(self):
        s=FC.root('H'); old=aw_block(s); target=M.plus(old,M.eye(3)); sched=P.Scheduler(F(1,100),0)
        out=P.post_prediction_prefix(s,h=F(1,100),pending_aw_floor=True,aw_floor_target=target,scheduler=sched,floor_solver_success=True,floor_eigenvectors=M.eye(3),floor_eigenvalues=[1,1,1])
        self.assertTrue(out.S_service_due); self.assertEqual(out.scheduler.elapsed,0); self.assertEqual(out.floor.state,out.state)
        self.assertEqual(aw_block(out.state),target)

    def test_readiness_remains_fail_closed(self):
        r=P.readiness(); self.assertTrue(r['pending_aw_floor_positive_part_success_branch']); self.assertTrue(r['scheduler_due_not_due_branches']); self.assertFalse(r['S_service_ldlt_branch_attached']); self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
