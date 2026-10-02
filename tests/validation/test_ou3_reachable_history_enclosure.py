import math,unittest
from tools.stability.ou3_theorem.reachable_history_enclosure import *

class ReachableHistoryEnclosureTests(unittest.TestCase):
 def test_generated_output_cannot_be_split(self):
  with self.assertRaises(ValueError): HistoryCell((("tau",HistoryInterval(.1,1)),),"h")
 def test_dependency_must_remain_shared(self):
  c=HistoryCell((("fast_accel_primitive",HistoryInterval(-.05,.05)),),"h")
  def bad(x):
   p=JointPairBound(1,2,0,.1,"different",True)
   return PropagatedCell(x,None,p,())
  with self.assertRaises(ArithmeticError): propagate(c,bad)
 def test_joint_quotient(self):
  c=HistoryCell((("slow_accel_endpoint",HistoryInterval(-.1,.1)),),"h")
  p=JointPairBound(.4,.8,-.1,.02,"h",True)
  x=PropagatedCell(c,{},p,(p,))
  z=certify((x,),.1)
  self.assertAlmostEqual(z["linked_supply_to_dissipation_upper"],.05)
  self.assertTrue(z["entry_budget_passed"])
 def test_zero_set_not_silently_divided(self):
  p=JointPairBound(0,1,0,.1,"h",False)
  self.assertTrue(math.isinf(p.quotient_upper()))
if __name__=="__main__": unittest.main()
