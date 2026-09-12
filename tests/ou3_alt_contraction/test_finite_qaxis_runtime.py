"""Integrated-OU Qaxis formula regressions; not deployment/source qualification."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_qaxis_runtime as Q

PASS=Q.PSDWitness(True)

class Tests(unittest.TestCase):
    def test_zero_sigma_is_zero_through_both_hygiene_layers(self):
        out=Q.qaxis4(1,F(1,200),0,F(199,200),marginal_psd=PASS,final_psd=PASS,machine_epsilon=F(1,10**7))
        self.assertEqual(out,M.zeros(4,4))

    def test_small_branch_marginal_matches_literal_leading_formula(self):
        tau=1; h=F(1,1000); s=F(3,2); a=F(999,1000); out=Q.marginal_raw(tau,h,s,a)
        inv=F(1); h2=h*h; h3=h2*h; h4=h3*h; h5=h4*h; h6=h5*h; h7=h6*h; h8=h7*h; h9=h8*h
        expected=s*(F(2,3)*h3-F(1,2)*h4+F(7,30)*h5-F(1,12)*h6+F(31,1260)*h7-F(1,160)*h8+F(127,90720)*h9)
        self.assertEqual(out[0][0],expected); self.assertEqual(out,M.transpose(out))

    def test_general_branch_retains_independent_nested_and_final_alpha_roots(self):
        tau=F(1,2); h=F(1,100); s=F(2); a=F(49,50)
        out=Q.marginal_raw(tau,h,s,a)
        inv=2; qc=2*s*inv; K22=tau*(1-a*a)/2
        self.assertEqual(out[2][2],qc*K22)
        changed=Q.marginal_raw(tau,h,s,F(48,50)); self.assertNotEqual(out,changed)
        base=Q.qaxis4(tau,h,1,a,marginal_psd=PASS,final_psd=PASS,machine_epsilon=F(1,10**7),small_branch=False,
                      marginal_alpha=a,final_alpha=a)
        split=Q.qaxis4(tau,h,1,a,marginal_psd=PASS,final_psd=PASS,machine_epsilon=F(1,10**7),small_branch=False,
                       marginal_alpha=F(979,1000),final_alpha=F(981,1000))
        self.assertNotEqual(base,split)

    def test_qaxis_copies_regularized_marginal_into_v_p_a_indices(self):
        tau=1; h=F(1,200); a=F(199,200); eps=F(1,10**7); marg=Q.regularize_psd(Q.marginal_raw(tau,h,1,a),PASS,machine_epsilon=eps)
        out=Q.qaxis4(tau,h,1,a,marginal_psd=PASS,final_psd=PASS,machine_epsilon=eps); idx=(0,1,3)
        for i in range(3):
            for j in range(3): self.assertEqual(out[idx[i]][idx[j]],marg[i][j])
        self.assertEqual(out,M.transpose(out))

    def test_ldlt_accept_cannot_bypass_same_matrix_inertia(self):
        S=[[1,0,0],[0,-2,0],[0,0,3]]
        with self.assertRaisesRegex(ValueError,'LDLT-accept witness detached'):
            Q.regularize_psd(S,Q.PSDWitness(True),machine_epsilon=F(1,10**7))

    def test_psd_eigen_branch_clips_negative_eigenvalues(self):
        S=[[1,0,0],[0,-2,0],[0,0,3]]; w=Q.PSDWitness(False,True,tuple(map(tuple,M.eye(3))),(1,-2,3))
        out=Q.regularize_psd(S,w,machine_epsilon=F(1,10**7)); self.assertEqual(out,[[1,0,0],[0,0,0],[0,0,3]])

    def test_psd_eigen_failure_adds_shipping_tolerance(self):
        S=M.eye(4); eps=F(1,10**7); w=Q.PSDWitness(False,False); out=Q.regularize_psd(S,w,machine_epsilon=eps)
        tol=64*eps
        for i in range(4): self.assertEqual(out[i][i],1+tol)

    def test_eigensystem_must_match_same_matrix(self):
        S=M.eye(3); w=Q.PSDWitness(False,True,tuple(map(tuple,M.eye(3))),(1,1,2))
        with self.assertRaises(ValueError): Q.regularize_psd(S,w,machine_epsilon=F(1,10**7))

    def test_readiness_stays_fail_closed(self):
        r=Q.readiness()
        self.assertTrue(r['free_Qaxis_matrix_removed_by_this_lemma'])
        self.assertTrue(r['LDLT_accept_branch_has_same_matrix_exact_inertia_guard'])
        self.assertTrue(r['nested_and_final_Qaxis_general_exp_roots_can_be_retained_separately'])
        self.assertFalse(r['arbitrary_LDLT_accept_boolean_can_bypass_matrix_relation'])
        self.assertFalse(r['Qaxis_covariance_exp_runtime_source_attached'])
        self.assertFalse(r['Eigen_LDLT_eigensolver_outcomes_attached'])
        self.assertFalse(r['complete_word_finite_identity'])
        self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
