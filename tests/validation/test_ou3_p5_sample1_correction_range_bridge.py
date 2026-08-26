import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_correction_range_bridge as B

class T(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=B.build()
 def test_pass(self):
  self.assertEqual(B.validate(self.d),[])
  self.assertEqual(self.d["P5_SAMPLE1_SCALAR_CORRECTION_RANGE_BRIDGE"],"PASS")
  self.assertTrue(self.d["scalar_core_inside_validated_quaternion_range"])
  self.assertGreater(self.d["scalar_range_headroom_rad_lower"],0.0)
 def test_no_overclaim(self):
  self.assertFalse(self.d["reset_process_tangent_force_perturbations_included"])
  self.assertFalse(self.d["sample1_S_due_not_due_family_closed_here"])
  self.assertFalse(self.d["complete_sample1_branch_closed_here"])
  self.assertFalse(self.d["whole_word_promoted_here"])
  self.assertFalse(self.d["N_H_words_set_here"])
 def test_numbers(self):
  self.assertLess(self.d["scalar_core_correction_norm_upper_rad"],self.d["deployed_quaternion_v2_range_upper_rad"])
  self.assertTrue(math.isfinite(self.d["scalar_core_correction_norm_upper_rad"]))

if __name__=="__main__": unittest.main()
