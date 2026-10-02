import unittest,math
from tools.stability.ou3_theorem.structured_covariance_defect import *
class StructuredDefectTests(unittest.TestCase):
 def test_zero_stays_zero(self):
  self.assertEqual(prediction_defect(0.,2.),0.)
  self.assertEqual(correction_defect(0.,10.,2.,.1)["D_next"],0.)
 def test_prediction_congruence(self):
  self.assertAlmostEqual(prediction_defect(.1,2.,.03),.43)
 def test_resolvent_fails_closed(self):
  z=correction_defect(1.,10.,2.,.3)
  self.assertFalse(z["closed"]);self.assertTrue(math.isinf(z["D_next"]))
 def test_literal_order(self):
  z=chronology_step(.001,[{"kind":"prediction","F_norm":1.1},
    {"kind":"correction","P_struct_norm":2.,"H_norm":1.,"S_struct_inv_norm":.1},
    {"kind":"reset","G_norm":1.01,"new_structured_reset_defect":.0002}])
  self.assertTrue(z["closed"]);self.assertGreater(z["D_next"],0.)
if __name__=="__main__":unittest.main()
