"""Exact attitude runtime algebra regressions; not source/stability evidence."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_attitude_runtime as A
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_ou_runtime_primitives as O
import test_finite_core as FC


class Tests(unittest.TestCase):
    def test_zero_rate_small_branch_is_literal_polynomial_limit(self):
        r=A.AngularRuntime((0,0,0),F(1,200)); R,B=r.rot_B(r.h)
        self.assertEqual(R,M.eye(3)); self.assertEqual(B,M.scaled(M.eye(3),r.h))
        Faa=r.F_AA(); self.assertEqual([row[:3] for row in Faa[:3]],M.eye(3)); self.assertEqual([row[3:] for row in Faa[:3]],B)

    def test_general_branch_requires_same_rate_and_times(self):
        full=A.TrigWitness(1,0,1,1); half=A.TrigWitness(F(1,2),0,1,1)
        r=A.AngularRuntime((1,0,0),1,full,half); self.assertFalse(r.small_rate)
        with self.assertRaises(ValueError): A.AngularRuntime((2,0,0),1,full,half)
        with self.assertRaises(ValueError): A.AngularRuntime((1,0,0),1,A.TrigWitness(F(1,2),0,1,1),half)

    def test_isotropic_branch_predicate_matches_shipping_tolerance(self):
        self.assertTrue(A.is_isotropic_shipping(M.eye(3)))
        q=M.eye(3); q[0][1]=q[1][0]=F(1,10**8)
        self.assertFalse(A.is_isotropic_shipping(q))

    def test_structured_Q_with_zero_bias_noise_reduces_isotropic_gyro_integral(self):
        r=A.AngularRuntime((0,0,0),F(1,100)); qb=M.zeros(6,6)
        for i in range(3): qb[i][i]=F(2,1000)
        q=A.structured_Q_pre_hygiene(r,qb)
        expected=M.zeros(6,6)
        for i in range(3): expected[i][i]=F(2,1000)*r.h
        self.assertEqual(q,expected)

    def test_fast_Q_branch_is_Qbase_times_h(self):
        r=A.AngularRuntime((0,0,0),F(1,200)); qb=M.scaled(M.eye(6),F(1,1000)); Faa,Qaa=A.attitude_blocks(r,qb,use_exact_Q=False)
        self.assertEqual(Qaa,M.scaled(qb,r.h)); self.assertEqual(Faa,r.F_AA())

    def test_psd_hygiene_success_and_retry_branches(self):
        S=M.eye(6); self.assertEqual(A.psd_hygiene_6(S,first_ldlt_success=True),S)
        bad=M.scaled(M.eye(6),-1)
        out=A.psd_hygiene_6(bad,first_ldlt_success=False,second_ldlt_success=True)
        self.assertTrue(all(out[i][i]>=A.EPS_PSD for i in range(6)))
        out2=A.psd_hygiene_6(bad,first_ldlt_success=False,second_ldlt_success=False)
        self.assertTrue(all(out2[i][i]>out[i][i] for i in range(6)))

    def test_highest_prediction_entry_has_no_free_attitude_or_linear_matrices(self):
        s=FC.root('A'); segment,kw=FC.physical_successor(s)
        angular=A.AngularRuntime((0,0,0),segment.h)
        ou=O.OUDecay(segment.h,1,F(199,200)); bias=O.BiasDecay(True,1,F(199,200),M.scaled(M.eye(3),F(1,100000)))
        qbase=M.scaled(M.eye(6),F(1,100000)); qaxis=[M.zeros(4,4)]*3
        out=A.runtime_paired_prediction(s,segment,gyro_body=kw['gyro_body'],angular=angular,Qbase=qbase,ou=ou,bias=bias,independent_qaxis=qaxis,use_exact_Q=True,first_ldlt_success=True)
        self.assertEqual(out.reference,segment.after); self.assertEqual(out.z[21:24],segment.after.beta)

    def test_readiness_stays_fail_closed_on_trig_and_ldlt_ancestry(self):
        r=A.readiness(); self.assertTrue(r['free_F_AA_Q_AA_removed_at_runtime_entry']); self.assertFalse(r['trig_values_same_runtime_source_attached']); self.assertFalse(r['Eigen_LDLT_outcomes_finite_precision_attached']); self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
