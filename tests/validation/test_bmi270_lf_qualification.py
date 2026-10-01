import math, unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.stability.ou3_theorem.bmi270_lf_qualification import design_fir, response, G, THETA_MARGIN
class LFQualificationTests(unittest.TestCase):
 def test_joint_intercepts(self):
  self.assertAlmostEqual(G*math.sin(THETA_MARGIN/2),.055578868,places=8)
  self.assertAlmostEqual(THETA_MARGIN/60,1.88916587367928e-4,places=12)
 def test_filter_passband(self):
  fs=200.; h=design_fir(fs)
  self.assertGreaterEqual(response(h,fs,.15),.99)
  self.assertAlmostEqual(sum(h),1.,places=12)
if __name__=="__main__": unittest.main()
