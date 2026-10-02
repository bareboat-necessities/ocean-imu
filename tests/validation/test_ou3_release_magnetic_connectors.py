import unittest,numpy as np
from tools.stability.ou3_theorem.aggregate_magnetic_service import *
from tools.stability.ou3_theorem.explicit_release_set import certificate
class ReleaseMagneticConnectorTests(unittest.TestCase):
 def test_aggregate_service_action(self):
  s=AggregateMagneticService(1,1);E=np.zeros((2,21));E[0,0]=1;E[1,2]=1
  A,c=service_action_lower(s,E);self.assertTrue(c["aggregate_service_used"]);self.assertFalse(c["event_schedule_enumerated"])
  self.assertAlmostEqual(A.mid[0][0],1.)
 def test_release_fails_closed_on_missing_AG_numeric(self):
  c=certificate();self.assertTrue(c["constructive_leaf_seed_available"]);self.assertTrue(c["release_image_must_be_computed_by_causal_propagator"])
if __name__=="__main__":unittest.main()
