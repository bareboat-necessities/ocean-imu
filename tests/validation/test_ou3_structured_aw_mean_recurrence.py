import unittest,math
import numpy as np
from tools.stability.ou3_theorem.structured_aw_mean_recurrence import *
class AWMeanTests(unittest.TestCase):
 def test_zero_source_decays(self):
  z=step(2,.9,1,.2,0,0);self.assertAlmostEqual(z['A_next'],1.8)
 def test_low_aw_axis_small(self):
  z=weighted_inj_defect(1,0,.1,.9);self.assertLess(z['D'],.01)
 def test_source_terms_linked(self):
  z=aw_mean_step(1,.8,.2,.3,.1,.4);self.assertAlmostEqual(z['A_next'],.9)
 def test_vector_order(self):
  aw=np.array([1.,0,0]);I=np.eye(3)
  z=vector_chronology_step(aw,.8,[('acc',I,np.array([.2,0,0])),('S',I,np.array([-.2,0,0]))])
  self.assertAlmostEqual(z['A_next'],.8)
 def test_mag_order(self):
  I=np.eye(3);z=vector_chronology_step(np.zeros(3),1.,[('acc',I,np.array([1.,0,0])),('mag',I,np.array([0.,1.,0]))])
  self.assertAlmostEqual(z['A_next'],math.sqrt(2))
if __name__=='__main__':unittest.main()
