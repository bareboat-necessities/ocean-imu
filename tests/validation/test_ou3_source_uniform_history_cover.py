import unittest
from tools.stability.ou3_theorem.source_uniform_history_cover import generate,max_ratio
from tools.stability.ou3_theorem.reachable_history_enclosure import HistoryCell,HistoryInterval
class SourceUniformCoverTests(unittest.TestCase):
 def test_splits_only_input_coordinate(self):
  r=HistoryCell((("physical_x",HistoryInterval(-1,1)),),"h")
  def prop(c):return c
  def cert(c):
   x=dict(c.coordinates)["physical_x"]
   if x.width<.3:
    class R: ratio_upper=.1
    return {"verified":True,"ratio":R()}
   return {"verified":False,"input_sensitivity":{"physical_x":1}}
  z=generate((r,),prop,cert,max_cells=100,max_depth=5)
  self.assertTrue(z.verified);self.assertTrue(max_ratio(z)["verified"])
 def test_unresolved_fails_closed(self):
  r=HistoryCell((("physical_x",HistoryInterval(-1,1)),),"h")
  z=generate((r,),lambda c:c,lambda c:{"verified":False},max_cells=2,max_depth=0)
  self.assertFalse(z.verified)
if __name__=="__main__":unittest.main()
