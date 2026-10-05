import struct,tempfile,unittest
from pathlib import Path
import numpy as np
from tools.stability.ou3_theorem.planar_service_stream import MAGIC,SIZES
from tools.stability.ou3_theorem.planar_every_root_service import initial_phi
class EveryRootServiceTest(unittest.TestCase):
 def test_initial_physical_columns_have_expected_norms(self):
  for k in (1,200,4000):
   p=initial_phi(k)
   self.assertAlmostEqual(np.linalg.norm(p[:3,0]),1.0,places=12)
   self.assertAlmostEqual(np.linalg.norm(p[:3,2]),1.0,places=12)
   self.assertAlmostEqual(np.linalg.norm(p[3:6,1]),.02,places=12)
   self.assertAlmostEqual(np.linalg.norm(p[3:6,3]),.02,places=12)
if __name__=="__main__":unittest.main()
