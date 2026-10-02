import unittest,numpy as np
from tools.stability.ou3_theorem.history_witness_builder import *
from tools.stability.ou3_theorem.marine_magnetic_qcqp import VectorBox,MagneticEventBox
class HistoryWitnessBuilderTests(unittest.TestCase):
 def vb(self,x,r=.001):x=np.array(x,float);return VectorBox(x-r,x+r)
 def test_builds_shared_marine_kernel_witness(self):
  ts=[0.,30.,60.]
  v=[TimedVectorBox(t,t,self.vb([0,0,0])) for t in ts]
  p=[TimedVectorBox(t,t,self.vb([.04*(i%2),0,0])) for i,t in enumerate(ts)]
  a=[TimedVectorBox(t,t,self.vb([0,0,0])) for t in ts];j=a
  g=[AttitudeGravityBox(t,t,self.vb([np.sin(.04*i),0,np.cos(.04*i)])) for i,t in enumerate(ts)]
  me=[MagneticEventBox(t,t,self.vb([30,0,0]),self.vb([0,0,0]),True) for t in np.arange(.5,60,1)]
  w=build_history_witness(dependency_token="h",velocity=v,position=p,acceleration=a,jerk=j,gravity=g,magnetic_events=me,start=0,end=60,attitude_line_mid=[1,0,0],attitude_line_rad=[1e-4]*3,ba_graph_mid=np.eye(3),ba_graph_rad=np.zeros((3,3)))
  self.assertEqual(w.dependency_token,"h");self.assertTrue(w.kernel_line.verified);self.assertTrue(w.marine_magnetic.gravity_windows);self.assertTrue(all(w.marine_magnetic.gravity_windows))
 def test_missing_long_window_witness_fails(self):
  x=[TimedVectorBox(0,0,self.vb([0,0,0]))];g=[AttitudeGravityBox(0,0,self.vb([0,0,1]))]
  with self.assertRaises(ArithmeticError):build_history_witness(dependency_token="h",velocity=x,position=x,acceleration=x,jerk=x,gravity=g,magnetic_events=[],start=0,end=60,attitude_line_mid=[1,0,0],attitude_line_rad=[0,0,0],ba_graph_mid=np.eye(3),ba_graph_rad=np.zeros((3,3)))
if __name__=="__main__":unittest.main()
