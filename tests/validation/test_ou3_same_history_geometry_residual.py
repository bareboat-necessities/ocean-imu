import unittest
import numpy as np
from tools.stability.ou3_theorem.same_history_geometry_residual import *
class GeometryResidualTests(unittest.TestCase):
 def test_zero_error_exact_specific_force(self):
  R=np.eye(3);a=np.array([.2,-.1,.3]);aw=a.copy();z=np.zeros(3)
  q=acc_geometry_and_residual(R_bw_phys=R,theta_error=z,physical_a=a,
      accel_error=z,aw_hat=aw,ba_hat=z)
  self.assertLess(np.linalg.norm(q["residual"]),1e-12)
 def test_s_identity(self):
  q=s_geometry_and_residual(np.array([1.,2.,3.]),np.array([.2,.3,.4]))
  self.assertLess(q["identity_error"],1e-12)
 def test_injection_is_same_Kr(self):
  K=np.zeros((21,3));K[0,0]=2.
  q=correction_injection(K,np.array([.3,0,0]))
  self.assertAlmostEqual(q["dtheta"][0],.6)
 def test_gravity_alone_refuses_yaw(self):
  with self.assertRaises(ArithmeticError):
   physical_rotation_from_gravity_heading(np.array([0.,0.,1.]))
 def test_world_geometry_eliminates_absolute_yaw(self):
  q=world_acc_row(aw_hat=np.array([.2,.1,.3]))
  self.assertFalse(q["absolute_attitude_required"])
  self.assertTrue(np.allclose(q["force_world"],np.array([.2,.1,.3-G])))
 def test_world_residual_signed_source(self):
  r=world_acc_residual_source(physical_a=np.array([1.,0,0]),aw_hat=np.array([.7,0,0]),
      world_sensor_error=np.array([-.1,0,0]),ba_world=np.array([.05,0,0]))
  self.assertAlmostEqual(r[0],.15)
if __name__=="__main__":unittest.main()
