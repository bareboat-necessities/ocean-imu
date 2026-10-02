import unittest,numpy as np
from tools.stability.ou3_theorem.closed_loop_release_stream import *
class ClosedLoopReleaseTests(unittest.TestCase):
 def test_ou_mean_literal_coefficients(self):
  x=np.zeros(21);x[15]=1.
  y,a=predict_mean(x,.005,1.1)
  pva,ppa,psa,aa=ou_coeffs(.005,1.1)
  self.assertAlmostEqual(y[6],pva);self.assertAlmostEqual(y[9],ppa)
  self.assertAlmostEqual(y[12],psa);self.assertAlmostEqual(y[15],aa)
 def test_world_acc_row_has_aw_ba(self):
  x=np.zeros(21);H=exact_world_acc_H(x)
  self.assertTrue(np.allclose(H[:,15:18],np.eye(3)))
  self.assertTrue(np.allclose(H[:,18:21],np.eye(3)))
 def test_reset_zeros_local_theta(self):
  x=np.ones(21);P=covariance_interval();dx=np.zeros(21);dx[0]=.01
  y,_,e=reset_after_correction(x,P,dx)
  self.assertTrue(np.allclose(y[:3],0));self.assertAlmostEqual(e["dtheta"][0],.01)
 def test_interval_H_uses_same_aw_radius(self):
  xm=np.zeros(21);xr=np.zeros(21);xr[15]=.2
  H=interval_world_acc_H(xm,xr)
  self.assertGreater(H.rad[1][2],0.)
 def test_interval_reset_uses_same_dtheta_radius(self):
  G=interval_reset_G(np.array([.1,0,0]),np.array([.02,0,0]))
  self.assertGreater(G.rad[1][2],0.)
 def test_interval_prediction_zero_radius_point_tau(self):
  xm=np.zeros(21);xr=np.zeros(21);xm[15]=1.
  m,r,_=interval_predict_mean(xm,xr,.005,1.1,1.1)
  p,_=predict_mean(xm,.005,1.1)
  self.assertTrue(np.allclose(m,p));self.assertLess(np.max(r),1e-12)
if __name__=="__main__":unittest.main()
