"""Committed ActiveParameters -> prediction/S-service composition regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_active_runtime_word as X
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_runtime as PR
from tools.stability.ou3_alt_contraction import finite_qaxis_runtime as Q
from tools.stability.ou3_alt_contraction import finite_ou_runtime_primitives as O
from tools.stability.ou3_alt_contraction import finite_post_prediction as PP
from tools.stability.ou3_alt_contraction import finite_measurement_runtime as MR
import test_finite_prediction_runtime as FP
import test_finite_runtime_parameters as RP

PASS=Q.PSDWitness(True)


def setup_prediction():
    helper=FP.Tests(); s,seg,g,a,_,bias,_=helper.make('A',False)
    raw=helper.raw_sample(s,g); active=RP.committed(True)
    ou=O.OUDecay(seg.h,active.tau,F(199,200))
    q=PR.QAxisBranch(False,active.Sigma_aw,(PASS,PASS,PASS),(PASS,PASS,PASS),F(1,10**7))
    return active,s,seg,raw,a,ou,bias,q


class Tests(unittest.TestCase):
    def test_active_commit_is_mandatory_ancestor_of_raw_prediction(self):
        active,s,seg,raw,a,ou,bias,q=setup_prediction()
        out=X.prediction_from_active(active,s,seg,raw,angular=a,Qbase=M.zeros(6,6),ou=ou,bias=bias,qaxis=q)
        self.assertIs(out.active,active); self.assertEqual(out.state.reference,seg.after)

    def test_prediction_rejects_correlated_or_detached_tuner_parameters_before_runtime(self):
        active,s,seg,raw,a,ou,bias,q=setup_prediction()
        bad=PR.QAxisBranch(True,active.Sigma_aw,(PASS,),(PASS,),F(1,10**7))
        with self.assertRaises(ValueError):
            X.prediction_from_active(active,s,seg,raw,angular=a,Qbase=M.zeros(6,6),ou=ou,bias=bias,qaxis=bad)
        with self.assertRaises(ValueError):
            X.prediction_from_active(active,s,seg,raw,angular=a,Qbase=M.zeros(6,6),ou=O.OUDecay(seg.h,2,F(199,200)),bias=bias,qaxis=q)

    def test_active_period_roots_scheduler_and_due_live_RS(self):
        active,s,seg,raw,a,ou,bias,q=setup_prediction()
        pred=X.prediction_from_active(active,s,seg,raw,angular=a,Qbase=M.zeros(6,6),ou=ou,bias=bias,qaxis=q)
        sched=PP.Scheduler(active.pseudo_period,active.pseudo_period-seg.h)
        prefix=X.post_prediction_from_active(pred,h=seg.h,pending_aw_floor=False,
            aw_floor_target=M.zeros(3,3),scheduler=sched)
        self.assertTrue(prefix.S_service_due)
        serviced=X.service_S_from_active(active,prefix,ldlt=MR.SafeLDLT(True,None,1,F(1,10**7)))
        self.assertIsNotNone(serviced.measurement); self.assertTrue(serviced.measurement.accepted)

    def test_detached_scheduler_period_is_rejected(self):
        active,s,seg,raw,a,ou,bias,q=setup_prediction()
        pred=X.prediction_from_active(active,s,seg,raw,angular=a,Qbase=M.zeros(6,6),ou=ou,bias=bias,qaxis=q)
        with self.assertRaises(ValueError):
            X.post_prediction_from_active(pred,h=seg.h,pending_aw_floor=False,
                aw_floor_target=M.zeros(3,3),scheduler=PP.Scheduler(active.pseudo_period+1,0))

    def test_not_due_branch_consumes_no_ldlt(self):
        active,s,seg,raw,a,ou,bias,q=setup_prediction()
        pred=X.prediction_from_active(active,s,seg,raw,angular=a,Qbase=M.zeros(6,6),ou=ou,bias=bias,qaxis=q)
        prefix=X.post_prediction_from_active(pred,h=seg.h,pending_aw_floor=False,
            aw_floor_target=M.zeros(3,3),scheduler=PP.Scheduler(active.pseudo_period,0))
        self.assertFalse(prefix.S_service_due)
        out=X.service_S_from_active(active,prefix); self.assertIsNone(out.measurement)
        with self.assertRaises(ValueError): X.service_S_from_active(active,prefix,ldlt=MR.SafeLDLT(True,None,1,F(1,10**7)))

    def test_readiness_fail_closed(self):
        r=X.readiness()
        self.assertTrue(r['same_history_tuner_parameters_attached_at_event_interface'])
        self.assertTrue(r['independent_Qaxis_branch_enforced_after_tuner_commit'])
        self.assertFalse(r['runtime_transcendental_and_solver_finite_precision_attached'])
        self.assertFalse(r['sensor_residual_source_bounds_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
