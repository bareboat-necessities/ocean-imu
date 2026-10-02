import unittest
import numpy as np
from tools.stability.ou3_theorem.adaptive_shaped_storage import (
    scaling, normalized_prediction_exact, physical_prediction,
    normalized_factor, tuner_coboundary, prediction_parity)

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

if __name__=="__main__":
    unittest.main()
