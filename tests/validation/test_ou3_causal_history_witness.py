import unittest,numpy as np
from tools.stability.ou3_theorem.reachable_history_enclosure import HistoryCell,HistoryInterval
from tools.stability.ou3_theorem.causal_history_witness import CausalWitnessState
from tools.stability.ou3_theorem.marine_magnetic_qcqp import VectorBox,MagneticEventBox
class CausalHistoryWitnessTests(unittest.TestCase):
 def vb(self,x,r=0):x=np.asarray(x,float);return VectorBox(x-r,x+r)
 def test_physical_history_and_kernel_cone(self):
  c=HistoryCell((("physical",HistoryInterval(0,1)),),"h");w=CausalWitnessState(c)
  w.append_physical(0,self.vb([0,0,0]),self.vb([0,0,0]),self.vb([0,0,0]),self.vb([0,0,0]),self.vb([0,0,1]))
  w.append_physical(30,self.vb([0,0,0]),self.vb([.04,0,0]),self.vb([0,0,0]),self.vb([0,0,0]),self.vb([.04,0,.9992]))
  w.append_kernel_line([1,0,0],np.eye(3),np.ones(6)*1e-4)
  r,rad,z=w.compatibility_cone();self.assertTrue(z["verified"]);self.assertEqual(len(r),21)
 def test_magnetic_event_retained(self):
  c=HistoryCell((("physical",HistoryInterval(0,1)),),"h");w=CausalWitnessState(c)
  w.append_magnetic(MagneticEventBox(.5,.5,self.vb([30,0,0]),self.vb([0,0,0]),True))
  self.assertEqual(len(w.magnetic_events),1)
if __name__=="__main__":unittest.main()
