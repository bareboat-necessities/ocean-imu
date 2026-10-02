import math,unittest,numpy as np
from tools.stability.ou3_theorem.force_direction_cone import *
class ForceConeTests(unittest.TestCase):
 def test_skew_norm_identity_bound(self):
  c=ForceCone(10,12,np.array([1.,0,0]),.1)
  self.assertGreaterEqual(c.skew_operator_radius(),1+12*2*math.sin(.05))
 def test_six_cones_cover_cube_directions(self):
  cs=octahedral_direction_cones()
  for v in [np.array(x,float) for x in ((1,1,1),(-1,1,1),(1,-1,1),(1,1,-1))]:
   v/=np.linalg.norm(v)
   self.assertTrue(any(math.acos(min(1,max(-1,float(v@c.axis))))<=c.half_angle+1e-12 for c in cs))
if __name__=="__main__":unittest.main()
