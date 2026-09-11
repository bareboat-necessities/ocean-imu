import unittest
from fractions import Fraction as F
from tools.stability.ou3_alt_contraction import finite_prediction_graph as G


class FinitePredictionGraphTests(unittest.TestCase):
    def test_identity_steps_preserve_finite_cayley_error(self):
        c=[F(3,10),F(-1,5),F(1,7)]
        self.assertEqual(G.relative_attitude_prediction(c,[1,0,0,0],[1,0,0,0]),c)

    def test_projective_step_quaternion_scaling_cancels(self):
        c=[F(1,3),F(-2,7),F(1,5)]
        qs=[F(9,10),F(1,10),F(-1,20),F(1,30)]
        qn=[F(19,20),F(-1,12),F(1,15),F(1,18)]
        a=G.relative_attitude_prediction(c,qs,qn)
        b=G.relative_attitude_prediction(c,[3*x for x in qs],[5*x for x in qn])
        self.assertEqual(a,b)

    def test_output_cayley_reconstructs_same_projective_quaternion(self):
        c=[F(2,9),F(-1,8),F(1,11)]
        qs=[F(7,8),F(1,7),F(-1,13),F(1,17)]
        qn=[F(11,12),F(-1,9),F(1,14),F(1,19)]
        qp=G.quat_mul(G.quat_mul(qs,[F(2),*c]),G.quat_conj(qn))
        cp=G.relative_attitude_prediction(c,qs,qn)
        # (2,c_plus) is exactly a scalar multiple of q_plus.
        for i in range(3):
            self.assertEqual(qp[0]*cp[i],2*qp[1+i])

    def test_translation_is_exact_true_minus_estimate_recurrence(self):
        h=F(1,200); pva=F(1,250); ppa=F(1,100000); pSa=F(1,50000000); alpha=F(99,100)
        coeff=G.translation_axis_coefficients(pva,ppa,pSa,alpha,h)
        e=[F(i-5,20) for i in range(12)]
        q=[F(i-7,30) for i in range(15)]
        out=G.translation_prediction(e,q,coeff)
        for a in range(3):
            ev,ep,eS,ea=e[a],e[3+a],e[6+a],e[9+a]
            a0,a1,J0,J1,J2=q[a],q[3+a],q[6+a],q[9+a],q[12+a]
            self.assertEqual(out[a],ev+pva*ea+J0-pva*a0)
            self.assertEqual(out[3+a],ep+h*ev+ppa*ea+J1-ppa*a0)
            self.assertEqual(out[6+a],eS+h*ep+h*h*ev/2+pSa*ea+J2-pSa*a0)
            self.assertEqual(out[9+a],alpha*ea+a1-alpha*a0)

    def test_bias_driver_is_shared_not_duplicated(self):
        e=[F(1,5),F(-1,7),F(1,9)]; beta=[F(1,11),F(1,13),F(-1,15)]; w=[F(1,100),F(-1,120),F(1,140)]
        ph=F(9,10);pt=F(19,20)
        ep,bp=G.bias_prediction(e,beta,w,ph,pt)
        for i in range(3):
            self.assertEqual(ep[i]-ph*e[i]-(pt-ph)*beta[i],w[i])
            self.assertEqual(bp[i]-pt*beta[i],w[i])

    def test_H_mode_holds_estimated_bias_but_not_physical_truth(self):
        z=[F(0)]*24
        z[18:21]=[F(1,5),F(-1,6),F(1,7)]
        z[21:24]=[F(1,8),F(1,9),F(-1,10)]
        w=[F(1,100),F(-1,110),F(1,120)]
        pt=F(49,50)
        q15=[F(0)]*15
        coeff=G.translation_axis_coefficients(0,0,0,1,F(1,200))
        d=G.prediction('H',z,q_shadow=[1,0,0,0],q_nominal=[1,0,0,0],q15=q15,coeff=coeff,w_bias=w,phi_true=pt)
        for i in range(3):
            self.assertEqual(d.state_after[18+i],z[18+i]+(pt-1)*z[21+i]+w[i])
            self.assertEqual(d.state_after[21+i],pt*z[21+i]+w[i])

    def test_A_mode_uses_configured_estimator_factor(self):
        z=[F(0)]*24;z[18]=F(2,5);z[21]=F(1,4)
        coeff=G.translation_axis_coefficients(0,0,0,1,F(1,200))
        d=G.prediction('A',z,q_shadow=[1,0,0,0],q_nominal=[1,0,0,0],q15=[F(0)]*15,coeff=coeff,w_bias=[F(1,100),0,0],phi_true=F(19,20),phi_hat=F(9,10))
        self.assertEqual(d.state_after[18],F(9,10)*F(2,5)+F(1,20)*F(1,4)+F(1,100))

    def test_cayley_pole_and_missing_A_factor_fail_closed(self):
        with self.assertRaises(ValueError):
            G.relative_attitude_prediction([0,0,0],[0,1,0,0],[1,0,0,0])
        coeff=G.translation_axis_coefficients(0,0,0,1,F(1,200))
        with self.assertRaises(TypeError):
            G.prediction('A',[F(0)]*24,q_shadow=[1,0,0,0],q_nominal=[1,0,0,0],q15=[F(0)]*15,coeff=coeff,w_bias=[0,0,0],phi_true=1)

    def test_readiness_closes_only_algebra(self):
        d=G.readiness()
        self.assertTrue(d['finite_attitude_prediction_identity'])
        self.assertTrue(d['finite_translation_q15_forcing_identity'])
        self.assertTrue(d['finite_shared_bias_driver_identity'])
        self.assertFalse(d['source_uniform_step_quaternion_graph_attached'])
        self.assertFalse(d['complete_word_finite_identity'])
        self.assertFalse(d['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
