import unittest
import numpy as np
from tools.stability.ou3_theorem.planar_parity import EVEN,ODD,permutation,service_probe_transform,split_information
class PlanarParityTests(unittest.TestCase):
 def test_partition(self):
  self.assertEqual((len(EVEN),len(ODD)),(12,9));self.assertEqual(sorted(EVEN+ODD),list(range(21)))
  p=permutation();self.assertTrue(np.allclose(p@p.T,np.eye(21)))
 def test_probe_transform_is_orthogonal(self):
  t=service_probe_transform();self.assertTrue(np.allclose(t@t.T,np.eye(4)))
 def test_block_information_split(self):
  # +/- representation of two independent scalar parity blocks.
  t=service_probe_transform();d=np.diag([2.,3.,5.,7.]);i=t.T@d@t
  a,b,c=split_information(i);self.assertLess(c,1e-12)
  self.assertTrue(np.allclose(a,np.diag([2.,3.])));self.assertTrue(np.allclose(b,np.diag([5.,7.])))
if __name__=="__main__":unittest.main()
