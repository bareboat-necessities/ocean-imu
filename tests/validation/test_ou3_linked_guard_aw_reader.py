import unittest,numpy as np
from tools.stability.ou3_theorem.interval_riccati_21 import IMat
from tools.stability.ou3_theorem.linked_guard_aw_reader import *
def P(a):
 a=np.asarray(a,float);z=np.zeros_like(a)
 return IMat(tuple(map(tuple,a)),tuple(map(tuple,z)))
class ReaderTests(unittest.TestCase):
 def test_single_acc(self):
  K=np.zeros((21,3));K[15,0]=.4;H=np.zeros((3,21))
  r=terminal_aw_reader_interval([{"kind":"acc","sample":0,"K":P(K),"H":P(H)}],1,(1,0,0))
  self.assertAlmostEqual(r["mid"][0,0],.4);self.assertEqual(r["rad"][0,0],0.)
 def test_future_prediction(self):
  K=np.zeros((21,3));K[15,0]=.4;H=np.zeros((3,21));F=np.eye(21);F[15,15]=.5
  ev=[{"kind":"acc","sample":0,"K":P(K),"H":P(H)},{"kind":"prediction","sample":1,"F":P(F)}]
  r=terminal_aw_reader_interval(ev,2,(1,0,0));self.assertAlmostEqual(r["mid"][0,0],.2)
 def test_paired_product(self):
  m,r=product(np.array([2.]),np.array([.1]),np.array([.3]),np.array([.02]))
  self.assertGreaterEqual(m[0]+r[0],2.1*.32)
if __name__=="__main__":unittest.main()
