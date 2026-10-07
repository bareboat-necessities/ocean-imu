import unittest
import numpy as np
from tools.stability.ou3_theorem.compatibility_quotient import decompose, quotient_terminal_map, quotient
from tools.stability.ou3_theorem.planar_compatibility_quotient_mean import compatibility_line, legacy_line

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

    def test_conditioned_metric_is_scale_invariant(self):
        rng=np.random.default_rng(653)
        U,_=np.linalg.qr(rng.normal(size=(21,21)))
        P=U@np.diag(np.geomspace(2e-8,.03,21))@U.T
        r=legacy_line(); M=np.eye(21)*.8
        for scale in (1e-12,1.,1e12):
            _,_,_,_,checks=quotient(scale*P,r)
            self.assertTrue(all(checks.values()))
            q=quotient_terminal_map(scale*P,r,scale*P,r,M,np.zeros(21))
            self.assertAlmostEqual(np.linalg.norm(q['M_Q'],2),.8,places=10)
            self.assertLess(q['gauge_injection_norm'],1e-10)

    def test_invalid_metrics_fail_closed(self):
        for P in (np.diag([1.,-1.]),np.zeros((2,2)),np.array([[1.,.5],[0.,1.]]),np.diag([1.,np.nan])):
            with self.assertRaises((ValueError,np.linalg.LinAlgError)):
                quotient(P,[1.,1.])
        with self.assertRaises(ValueError): quotient(np.eye(2),[0.,0.])

    def test_physical_chart_line_annihilates_actual_measurement_models(self):
        def skew(v):
            x,y,z=v; return np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
        for k in (40000,40500,41000,41500,44000):
            phase=.02*np.sin(np.pi*(k/200)/10)
            c,s=np.cos(phase),np.sin(phase)
            R=np.array([[c,0,-s],[0,1,0],[s,0,c]])
            ax=-.02*(np.pi/10)**2*np.sin(np.pi*(k/200)/10)
            H=np.zeros((3,21)); H[:,:3]=-skew(R@np.array([ax,0,-9.80665])); H[:,18:]=np.eye(3)
            r=compatibility_line(k)
            np.testing.assert_allclose(H@r,0,atol=3e-15)
            np.testing.assert_allclose(-skew(75*R[:,0])@r[:3],0,atol=3e-15)
            if k==40000: self.assertGreater(np.linalg.norm(H@legacy_line()),19.)

if __name__=="__main__": unittest.main()
