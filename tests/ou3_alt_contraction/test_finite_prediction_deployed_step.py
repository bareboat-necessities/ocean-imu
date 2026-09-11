import unittest
from fractions import Fraction as F
from tools.stability.ou3_alt_contraction import finite_prediction_deployed_step as D
from tools.stability.ou3_alt_contraction import finite_prediction_graph as G


class FinitePredictionDeployedStepTests(unittest.TestCase):
    def test_zero_increment_step_is_identity_quaternion(self):
        self.assertEqual(D.step_quaternion_polynomial([0,0,0]),[F(1),F(0),F(0),F(0)])

    def test_threshold_is_strict(self):
        with self.assertRaises(ValueError):D.step_quaternion_polynomial([F(1,100),0,0])

    def test_same_omega_bg_h_drive_shadow_and_nominal(self):
        z=[F(0)]*24;z[:3]=[F(1,4),F(-1,7),F(1,9)];z[3:6]=[F(1,100),F(-1,120),F(1,140)]
        omega=[F(3,10),F(-1,5),F(1,4)];h=F(1,200)
        coeff=G.translation_axis_coefficients(0,0,0,1,h)
        d=D.prediction('H',z,omega_hat=omega,h=h,q15=[F(0)]*15,coeff=coeff,w_bias=[0,0,0],phi_true=1)
        dn=[-x*h for x in omega];ds=[(-omega[i]+z[3+i])*h for i in range(3)]
        expected=G.relative_attitude_prediction(z[:3],D.step_quaternion_polynomial(ds),D.step_quaternion_polynomial(dn))
        self.assertEqual(d.state_after[:3],expected)

    def test_readiness_closes_branch_binding_but_not_source_word(self):
        d=D.readiness()
        self.assertTrue(d['finite_prediction_graph_consumed'])
        self.assertTrue(d['source_uniform_polynomial_quaternion_branch_consumed'])
        self.assertTrue(d['same_omega_bg_h_generate_nominal_and_shadow_steps'])
        self.assertFalse(d['free_step_quaternion_coefficients_used'])
        self.assertFalse(d['q15_same_history_source_relation_attached_to_finite_word'])
        self.assertFalse(d['complete_word_finite_identity'])
        self.assertFalse(d['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
