"""Highest finite prediction runtime composition regressions; no source promotion."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_runtime as R
from tools.stability.ou3_alt_contraction import finite_qaxis_runtime as Q
from tools.stability.ou3_alt_contraction import finite_attitude_runtime as A
from tools.stability.ou3_alt_contraction import finite_ou_runtime_primitives as O
import test_finite_core as FC

PASS=Q.PSDWitness(True)

class Tests(unittest.TestCase):
    def make(self,mode='A',correlated=True):
        s=FC.root(mode); segment,kw=FC.physical_successor(s); h=segment.h
        gyro=[s.reference.gyro_bias[i]-s.z[3+i] for i in range(3)]
        angular=A.AngularRuntime((0,0,0),h); ou=O.OUDecay(h,1,F(199,200))
        bias=O.BiasDecay(mode=='A',1,F(199,200) if mode=='A' else 1,M.scaled(M.eye(3),F(1,100000)))
        qaxis=R.QAxisBranch(correlated,M.zeros(3,3),(PASS,) if correlated else (PASS,PASS,PASS),(PASS,) if correlated else (PASS,PASS,PASS),F(1,10**7))
        return s,segment,gyro,angular,ou,bias,qaxis

    def test_active_correlated_prediction_has_no_free_matrices(self):
        s,seg,g,a,ou,b,q=self.make('A',True)
        out=R.prediction(s,seg,gyro_body=g,angular=a,Qbase=M.scaled(M.eye(6),F(1,100000)),ou=ou,bias=b,qaxis=q)
        self.assertEqual(out.reference,seg.after); self.assertEqual(out.z[21:24],seg.after.beta)

    def test_held_independent_prediction_keeps_ba_estimate_fixed(self):
        s,seg,g,a,ou,b,q=self.make('H',False)
        out=R.prediction(s,seg,gyro_body=g,angular=a,Qbase=M.zeros(6,6),ou=ou,bias=b,qaxis=q)
        self.assertEqual(out.reference,seg.after)
        # H18 holds b_a_hat. Since e_ba=beta_true-b_a_hat, truth decay/driver
        # moves e_ba by exactly the truth increment; do not incorrectly add
        # the physical driver alone.
        before_hat=tuple(s.reference.beta[i]-s.z[18+i] for i in range(3))
        after_hat=tuple(seg.after.beta[i]-out.z[18+i] for i in range(3))
        self.assertEqual(after_hat,before_hat)
        expected=tuple(s.z[18+i]+(seg.phi_true-1)*s.reference.beta[i]+seg.bias_driver[i] for i in range(3))
        self.assertEqual(out.z[18:21],expected)

    def test_covariance_angular_rate_must_match_nominal_gyro(self):
        s,seg,g,a,ou,b,q=self.make('A',True)
        wrong=A.AngularRuntime((F(1,10**8),0,0),seg.h)
        with self.assertRaises(ValueError): R.prediction(s,seg,gyro_body=g,angular=wrong,Qbase=M.zeros(6,6),ou=ou,bias=b,qaxis=q)

    def test_step_must_be_shared_by_physical_attitude_and_ou_paths(self):
        s,seg,g,a,ou,b,q=self.make('A',True); bad=O.OUDecay(seg.h*2,1,F(99,100))
        with self.assertRaises(ValueError): R.prediction(s,seg,gyro_body=g,angular=a,Qbase=M.zeros(6,6),ou=bad,bias=b,qaxis=q)

    def test_qaxis_witness_count_matches_correlated_branch(self):
        with self.assertRaises(ValueError): R.QAxisBranch(True,M.eye(3),(PASS,PASS),(PASS,),F(1,10**7))
        with self.assertRaises(ValueError): R.QAxisBranch(False,M.eye(3),(PASS,),(PASS,),F(1,10**7))

    def test_readiness_fail_closed(self):
        r=R.readiness(); self.assertFalse(r['free_prediction_transition_matrices_at_entry']); self.assertTrue(r['same_bias_corrected_gyro_for_nominal_and_covariance']); self.assertTrue(r['Qaxis_runtime_graph_composed']); self.assertFalse(r['same_history_tuner_parameters_attached']); self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
