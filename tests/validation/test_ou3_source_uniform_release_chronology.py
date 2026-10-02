import unittest
from tools.stability.ou3_theorem.source_uniform_release_chronology import *
class ReleaseChronologyTests(unittest.TestCase):
 def test_seed_is_numeric(self):
  z=release_seed();self.assertEqual(len(z["mean_state"]),21)
  self.assertEqual(len(z["covariance_diagonal_upper"]),21)
 def test_scheduler_fixed_tau(self):
  class X:
   def __init__(self,x):self.lo=x;self.hi=x
  tr=[{"tau":X(1.1)} for _ in range(100)]
  ev=schedule_from_adaptation(tr,.005)
  self.assertTrue(all(x["resolved"] for x in ev))
  self.assertGreater(len(ev),0)
 def test_scheduler_interval_fails_closed(self):
  class X:
   lo=1.;hi=2.
  ev=schedule_from_adaptation([{"tau":X()}],.005)
  self.assertFalse(ev[0]["resolved"])
if __name__=="__main__":unittest.main()
