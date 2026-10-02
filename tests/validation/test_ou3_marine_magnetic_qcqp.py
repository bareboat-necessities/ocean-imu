import unittest,numpy as np,json
from pathlib import Path
from tools.stability.ou3_theorem.marine_magnetic_qcqp import *
C=json.loads((Path(__file__).resolve().parents[2]/"tools/stability/ou3_theorem/constants.json").read_text())
class MarineMagneticTests(unittest.TestCase):
 def vb(self,x,r=0):x=np.array(x,float);return VectorBox(x-r,x+r)
 def test_norm_cap(self):
  self.assertIs(norm_cap(self.vb([1,0,0],.01),2),True)
  self.assertIs(norm_cap(self.vb([3,0,0],.01),2),False)
 def test_displacement_witness(self):
  self.assertIs(difference_span_witness(self.vb([0,0,0],.001),self.vb([.04,0,0],.001),.03),True)
 def test_gravity_span_witness(self):
  th=.04;self.assertIs(gravity_direction_span_witness(self.vb([0,0,1]),self.vb([math.sin(th),0,math.cos(th)]),.034),True)
 def test_magnetic_recurring_service(self):
  es=[MagneticEventBox(t,t,self.vb([30,0,0]),self.vb([0,0,0]),True) for t in (.5,1.5,2.5,3.5)]
  self.assertIs(recurring_magnetic_service(es,0,4,1,20,75,2),True)
 def test_existential_window_group(self):
  good=(self.vb([0,0,1]),self.vb([.04,0,.9992]))
  bad=(self.vb([0,0,1]),self.vb([0,0,1]))
  w=MarineMagneticWitness([self.vb([0,0,0])],[self.vb([0,0,0])],[self.vb([0,0,0])],[self.vb([0,0,0])],[[bad,good]],[[ (self.vb([0,0,0]),self.vb([.04,0,0])) ]],[],0,0)
  # magnetic interval vacuous at zero horizon; one good pair is sufficient.
  self.assertIs(callback(w,C),True)
if __name__=="__main__":unittest.main()
