import unittest
import numpy as np
from tools.stability.ou3_theorem.adaptive_shaped_storage import (
    scaling, normalized_prediction_exact, physical_prediction,
    normalized_factor, tuner_coboundary, prediction_parity, normalized_S_correction21)

class AdaptiveShapedStorageTests(unittest.TestCase):
    def test_exact_real_prediction_normalization(self):
        for h,tau,sigma in [(.005,.02,.05),(.005,1.1,.4),(.006,12.,4.)]:
            self.assertLess(prediction_parity(h,tau,sigma),2e-12)

    def test_tuner_coboundary_telescopes(self):
        vals=[(.8,.2),(1.1,.4),(3.2,.7),(2.4,.5)]
        C=np.eye(4)
        for (ta,sa),(tb,sb) in zip(vals,vals[1:]):
            C=tuner_coboundary(ta,sa,tb,sb)@C
        direct=scaling(*vals[-1])@np.linalg.inv(scaling(*vals[0]))
        self.assertLess(np.linalg.norm(C-direct,ord=np.inf),2e-12)

    def test_generated_coordinate_factor_identity(self):
        h=.005; a=(1.1,.4); b=(1.12,.41)
        A=physical_prediction(h,a[0])
        left=normalized_factor(A,*a,*b)
        right=tuner_coboundary(*a,*b)@normalized_prediction_exact(h,a[0])
        self.assertLess(np.linalg.norm(left-right,ord=np.inf),2e-12)

    def test_S_factor_is_similarity_of_literal_gain(self):
        H=np.zeros((3,21));H[:,12:15]=np.eye(3)
        K=np.zeros((21,3));K[12:15,:]=.2*np.eye(3);K[15:18,:]=.03*np.eye(3)
        A=normalized_S_correction21(K,H,1.7,.42)
        self.assertTrue(np.all(np.isfinite(A)))
        self.assertLess(np.linalg.norm(A[12:15,12:15]-.8*np.eye(3)),1e-12)

if __name__=="__main__":
    unittest.main()
