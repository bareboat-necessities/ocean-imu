import unittest,numpy as np
from tools.stability.ou3_theorem.interval_riccati_21 import diagonal_interval,selector_h
from tools.stability.ou3_theorem.rank_loss_interval_factor import exact
from tools.stability.ou3_theorem.reachable_history_enclosure import HistoryCell,HistoryInterval
from tools.stability.ou3_theorem.causal_covariance_supply import new_state
class CausalCovarianceSupplyTests(unittest.TestCase):
 def seed(self):
  c=HistoryCell((("accel_history",HistoryInterval(-1,1)),),"h")
  P=diagonal_interval((1.,)*21,(1.,)*21)
  return new_state(c,P,2)
 def test_S_gain_generated_and_joint_updated(self):
  s=self.seed();H=selector_h(12);R=diagonal_interval((.25,)*3,(.25,)*3);D=np.zeros((21,2));D[12,0]=.1
  z=s.correct(H,R,D,"S")
  self.assertTrue(z["innovation_inverse"]["verified"]);self.assertEqual(len(s.prefix),1)
  self.assertIsNotNone(s.interval_joint);self.assertEqual(s.interval_joint.dependency_token,"h")
  self.assertNotEqual(float(s.joint.B[12,0]),0.)
 def test_verified_full_precision_used(self):
  s=self.seed();Q,cert=s._Qinterval();self.assertTrue(cert["verified"]);self.assertEqual(Q.shape,(21,21))
 def test_prediction_and_reset_preserve_token(self):
  s=self.seed();F=exact(np.eye(21).tolist());Q=diagonal_interval((1e-6,)*21,(1e-6,)*21)
  s.predict(F,Q,np.zeros((21,2)));s.reset(F,np.zeros((21,2)))
  self.assertEqual(s.joint.dependency_token,"h");self.assertEqual(len(s.prefix),2)
 def test_covariance_only_sync_does_not_add_mean_supply(self):
  s=self.seed();n=len(s.prefix);s.covariance_sync_aw(2.);self.assertEqual(len(s.prefix),n)
if __name__=="__main__":unittest.main()
