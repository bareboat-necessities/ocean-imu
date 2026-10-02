# ruff: noqa: F401, F811
import unittest,numpy as np
from tools.stability.ou3_theorem.verified_linked_qcqp import Box,verified_supply_bnb,verified_interval_supply_bnb\nfrom tools.stability.ou3_theorem.rank_loss_interval_factor import exact
class VerifiedLinkedQCQPTests(unittest.TestCase):
 def test_scalar_supply_upper(self):
  B=np.array([[.2]]);C=np.array([[-.1]]);e=Box(np.array([-1.]),np.array([1.]));u=Box(np.array([-.5]),np.array([.5]))
  A=np.array([[1.],[-1.]]);b=np.array([.5,.5])
  z=verified_supply_bnb(B,C,e,u,A,b,lambda box:True,tol=.02,max_leaves=10000)
  self.assertTrue(z.verified);self.assertGreaterEqual(z.upper,.225)
 def test_interval_BC_radius_is_included(self):
  e=Box(np.array([-1.]),np.array([1.]));u=Box(np.array([-.5]),np.array([.5]))
  B=exact([[.2]]);C=exact([[-.1]])
  # Add explicit coefficient radius.
  from tools.stability.ou3_theorem.interval_riccati_21 import IMat
  B=IMat(B.mid,((.01,),));C=IMat(C.mid,((.02,),))
  z=verified_interval_supply_bnb(B,C,e,u,np.zeros((0,1)),np.zeros(0),lambda box:True,tol=.05,max_leaves=10000)
  self.assertTrue(z.verified);self.assertGreater(z.upper,.225)
 def test_unknown_side_fails_closed(self):
  B=np.zeros((1,1));C=np.zeros((1,1));e=Box(np.array([0.]),np.array([0.]));u=Box(np.array([0.]),np.array([0.]))
  z=verified_supply_bnb(B,C,e,u,np.zeros((0,1)),np.zeros(0),lambda box:None,tol=1e-6)
  self.assertFalse(z.verified)
if __name__=="__main__":unittest.main()
