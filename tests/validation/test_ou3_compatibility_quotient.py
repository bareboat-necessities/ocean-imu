import unittest
import numpy as np
from tools.stability.ou3_theorem.compatibility_quotient import decompose, quotient_terminal_map

class CompatibilityQuotientTests(unittest.TestCase):
    def test_covariance_metric_decomposition(self):
        P=np.diag([2.,3.,5.,7.]); r=np.array([1.,0.,2.,0.]); e=np.array([.3,-.2,.4,.7])
        c=decompose(P,r,e)
        self.assertAlmostEqual(c["V"],c["alpha"]**2+c["V_perp"],places=12)
        self.assertTrue(all(c["projector_checks"].values()))
    def test_rotating_line_is_explicit_forcing(self):
        P=np.eye(3); r0=np.array([1.,0.,0.]); r1=np.array([0.,1.,0.])
        q=quotient_terminal_map(P,r0,P,r1,np.eye(3),np.zeros(3))
        self.assertEqual(q["quotient_dimension"],2)
        self.assertGreater(q["gauge_injection_norm"],.99)
    def test_preserved_line_has_no_transverse_injection(self):
        P=np.eye(3); r=np.array([1.,0.,0.])
        q=quotient_terminal_map(P,r,P,r,np.eye(3),np.zeros(3))
        self.assertLess(q["gauge_injection_norm"],1e-12)

if __name__=="__main__": unittest.main()
