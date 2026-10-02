import unittest,math,numpy as np
from tools.stability.ou3_theorem.inj_block_algebra import INJ
from tools.stability.ou3_theorem.changing_axis_inj import rodrigues,common_rotation_transport
from tools.stability.ou3_theorem.same_history_axis_mismatch import *
class ChangingAxisTests(unittest.TestCase):
 def test_common_rotation_exact(self):
  x=INJ(2,.3,-.2);n=np.array([0.,0.,1.]);U=rodrigues([1,2,3],.2)
  y,m,e=common_rotation_transport(x,U,n);self.assertLess(e,1e-12);self.assertEqual(y,x)
 def test_exact_difference_norms(self):
  a=.3;self.assertAlmostEqual(projector_difference_norm(a),math.sin(a))
  self.assertAlmostEqual(skew_difference_norm(a),2*math.sin(a/2))
if __name__=="__main__":unittest.main()
