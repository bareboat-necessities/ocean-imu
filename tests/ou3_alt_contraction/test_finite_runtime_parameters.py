"""Tuner commit to runtime-parameter ancestry regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_tuner_commit as T
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as R
from tools.stability.ou3_alt_contraction import finite_prediction_runtime as PR
from tools.stability.ou3_alt_contraction import finite_qaxis_runtime as Q
from tools.stability.ou3_alt_contraction import finite_ou_runtime_primitives as O
from tools.stability.ou3_alt_contraction import finite_post_prediction as PP

PASS=Q.PSDWitness(True)

def committed(live=True):
    cfg=T.CommitConfig(F(3,2),F(1,100),2,F(3,200),True,True,F(2,5),5,F(4,5),F(6,5),F(7,5))
    tune=T.TuneState(1,F(1,10),2)
    c=T.commit(tune,cfg,pending=True,live=live,band_noise_floor_sigma=F(1,20),rs_sqrt_scale=F(1,10) if live else None)
    return R.ActiveParameters.from_commit(c)

class Tests(unittest.TestCase):
    def test_commit_forces_same_tau_diagonal_sigma_and_independent_qaxis(self):
        active=committed(); ou=O.OUDecay(F(1,200),active.tau,F(199,200))
        q=PR.QAxisBranch(False,active.Sigma_aw,(PASS,PASS,PASS),(PASS,PASS,PASS),F(1,10**7))
        self.assertTrue(active.require_prediction(ou=ou,qaxis=q))
        badcorr=PR.QAxisBranch(True,active.Sigma_aw,(PASS,),(PASS,),F(1,10**7))
        with self.assertRaises(ValueError): active.require_prediction(ou=ou,qaxis=badcorr)
        with self.assertRaises(ValueError): active.require_prediction(ou=O.OUDecay(F(1,200),2,F(199,200)),qaxis=q)
        badsig=PR.QAxisBranch(False,M.eye(3),(PASS,PASS,PASS),(PASS,PASS,PASS),F(1,10**7))
        with self.assertRaises(ValueError): active.require_prediction(ou=ou,qaxis=badsig)

    def test_commit_forces_scheduler_period(self):
        active=committed(); self.assertTrue(active.require_scheduler(PP.Scheduler(active.pseudo_period,0)))
        with self.assertRaises(ValueError): active.require_scheduler(PP.Scheduler(active.pseudo_period+1,0))

    def test_period_retarget_preserves_elapsed_credit_below_new_deadline(self):
        active=committed(); old=PP.Scheduler(2,F(1,3),F(1,1000))
        new=R.retarget_scheduler(active,old)
        self.assertEqual(new.period,active.pseudo_period); self.assertEqual(new.elapsed,F(1,3)); self.assertEqual(new.tolerance,F(1,1000))

    def test_overdue_retarget_requires_explicit_nextafter_witness(self):
        active=committed(); old=PP.Scheduler(2,F(3,2))
        # committed() has new period 3/2, so equality is already due/overdue.
        with self.assertRaisesRegex(ValueError,'nextafter witness'):
            R.retarget_scheduler(active,old)
        parked=R.retarget_scheduler(active,old,park=R.NextafterParkWitness(F(149,100)))
        self.assertEqual(parked.period,F(3,2)); self.assertEqual(parked.elapsed,F(149,100))
        with self.assertRaises(ValueError):
            R.retarget_scheduler(active,old,park=R.NextafterParkWitness(F(3,2)))

    def test_live_commit_forces_same_RS(self):
        active=committed(True); self.assertTrue(active.require_RS(active.R_S))
        with self.assertRaises(ValueError): active.require_RS(M.eye(3))
        cold=committed(False)
        with self.assertRaises(ValueError): cold.require_RS(M.eye(3))

    def test_readiness_fail_closed(self):
        r=R.readiness(); self.assertTrue(r['tuner_commit_forces_independent_Qaxis_branch']); self.assertTrue(r['committed_live_RS_to_S_measurement']); self.assertTrue(r['period_change_preserves_subdeadline_elapsed_credit']); self.assertFalse(r['nextafter_binary32_ancestry_attached']); self.assertFalse(r['TuneState_frontend_history_attached']); self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
