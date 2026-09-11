"""Exact shared-runtime-root regressions; not BRMM/source admission."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_ou_runtime_primitives as R
from tools.stability.ou3_alt_contraction import finite_prediction_covariance as C
import test_finite_core as FC


class Tests(unittest.TestCase):
    def test_general_branch_matches_shipping_closed_form(self):
        r=R.OUDecay(F(1,200),F(1,2),F(99,100))
        self.assertFalse(r.small_coeff_branch)
        coeff=r.axis_coefficients()[0]; h,tau,a=r.h,r.tau,r.alpha; x=h/tau; em1=a-1
        self.assertEqual(coeff[0],tau*(1-a))
        self.assertEqual(coeff[1],tau*tau*(x+em1))
        self.assertEqual(coeff[2],tau**3*(F(1,2)*x*x-x-em1))
        self.assertEqual(r.axis_coefficients(),(coeff,coeff,coeff))

    def test_small_branch_matches_safe_phi_A_polynomial(self):
        r=R.OUDecay(F(1,1000),F(1),F(999,1000))
        self.assertTrue(r.small_coeff_branch)
        va,pa,Sa,a,h=r.axis_coefficients()[0]; x=h/r.tau
        self.assertEqual(va,r.tau*(1-a))
        self.assertEqual(pa,r.tau**2*(F(1,2)*x**2-F(1,6)*x**3+F(1,24)*x**4))
        self.assertEqual(Sa,r.tau**3*(F(1,6)*x**3-F(1,24)*x**4+F(1,120)*x**5))

    def test_branch_boundary_is_literal_strict_less_than(self):
        self.assertTrue(R.OUDecay(F(99,10000),1,F(99,100)).small_coeff_branch)
        self.assertFalse(R.OUDecay(F(1,100),1,F(99,100)).small_coeff_branch)

    def test_active_BA_QBB_uses_same_phi_squared(self):
        q=M.scaled(M.eye(3),F(1,1000)); b=R.BiasDecay(True,F(3,2),F(9,10),q)
        expected=M.scaled(q,F(3,4)*(1-F(81,100)))
        self.assertEqual(b.covariance_increment(),expected)
        blocks=R.runtime_blocks(F_AA=M.eye(6),Q_AA=M.zeros(6,6),ou=R.OUDecay(F(1,200),1,F(199,200)),bias=b,independent_qaxis=[M.zeros(4,4)]*3)
        self.assertEqual(blocks.phi_b,b.phi_b)
        self.assertEqual(list(map(list,blocks.Q_BB)),expected)

    def test_held_branch_has_identity_factor_and_zero_QBB(self):
        b=R.BiasDecay(False,1,1,M.eye(3))
        self.assertEqual(b.covariance_increment(),M.zeros(3,3))
        with self.assertRaises(ValueError): R.BiasDecay(False,1,F(99,100),M.eye(3))

    def test_runtime_blocks_remove_free_FLL_and_QBB(self):
        ou=R.OUDecay(F(1,200),1,F(199,200)); b=R.BiasDecay(True,1,F(199,200),M.eye(3))
        qaxis=M.scaled(M.eye(4),F(1,10000)); sig=M.eye(3)
        blocks=R.runtime_blocks(F_AA=M.eye(6),Q_AA=M.zeros(6,6),ou=ou,bias=b,qaxis_unit=qaxis,sigma_aw=sig)
        self.assertEqual(list(map(list,blocks.F_LL)),C.linear_transition(ou.axis_coefficients()))
        self.assertEqual(list(map(list,blocks.Q_LL)),C.correlated_linear_process(qaxis,sig))
        self.assertEqual(list(map(list,blocks.Q_BB)),b.covariance_increment())

    def test_runtime_paired_prediction_uses_no_free_axis_or_phi_arguments(self):
        s=FC.root('A'); segment,kw=FC.physical_successor(s)
        ou=R.OUDecay(segment.h,1,F(199,200)); b=R.BiasDecay(True,1,F(199,200),M.scaled(M.eye(3),F(1,100000)))
        out=R.runtime_paired_prediction(s,segment,gyro_body=kw['gyro_body'],ou=ou,bias=b,F_AA=M.eye(6),Q_AA=M.zeros(6,6),independent_qaxis=[M.zeros(4,4)]*3)
        self.assertEqual(out.reference,segment.after)
        self.assertEqual(out.z[21:24],segment.after.beta)

    def test_detached_duration_and_mode_are_rejected(self):
        s=FC.root('A'); segment,kw=FC.physical_successor(s); b=R.BiasDecay(True,1,F(199,200),M.eye(3))
        with self.assertRaises(ValueError):
            R.runtime_paired_prediction(s,segment,gyro_body=kw['gyro_body'],ou=R.OUDecay(segment.h*2,1,F(99,100)),bias=b,F_AA=M.eye(6),Q_AA=M.zeros(6,6),independent_qaxis=[M.zeros(4,4)]*3)
        with self.assertRaises(ValueError):
            R.runtime_paired_prediction(FC.root('H'),segment,gyro_body=kw['gyro_body'],ou=R.OUDecay(segment.h,1,F(199,200)),bias=b,F_AA=M.eye(6),Q_AA=M.zeros(6,6),independent_qaxis=[M.zeros(4,4)]*3)

    def test_invalid_branch_inputs_fail_closed(self):
        for args in ((0,1,1),(1,0,1),(1,1,0),(1,1,F(3,2))):
            with self.assertRaises(ValueError): R.OUDecay(*args)
        with self.assertRaises(ValueError): R.BiasDecay(True,F(1,10000),F(99,100),M.eye(3))
        ou=R.OUDecay(F(1,200),1,F(199,200)); b=R.BiasDecay(False,1,1,M.eye(3))
        with self.assertRaises(ValueError): R.runtime_blocks(F_AA=M.eye(6),Q_AA=M.zeros(6,6),ou=ou,bias=b)
        with self.assertRaises(ValueError): R.runtime_blocks(F_AA=M.eye(6),Q_AA=M.zeros(6,6),ou=ou,bias=b,qaxis_unit=M.eye(4),sigma_aw=M.eye(3),independent_qaxis=[M.eye(4)]*3)

    def test_readiness_stays_fail_closed(self):
        r=R.readiness()
        self.assertTrue(r['free_axis_coefficients_removed_at_runtime_entry'])
        self.assertTrue(r['free_phi_hat_removed_at_runtime_entry'])
        self.assertTrue(r['free_Q_BB_removed_at_runtime_entry'])
        self.assertFalse(r['alpha_exp_runtime_source_attached'])
        self.assertFalse(r['attitude_F_Q_primitives_source_attached'])
        self.assertFalse(r['complete_word_finite_identity'])
        self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
