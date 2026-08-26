import math
import sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_rotation_gauge_refinement_v2 as G

class T(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=G.build(delta_pieces=4,force_pieces=4)
 def test_validates_fail_closed(self):
  d=self.d; self.assertEqual(G.validate(d),[]); self.assertGreater(d["evaluated_joint_cells"],0)
  self.assertIn(d["P5_SAMPLE1_DELTA_FORCE_SUBDIVIDED_WITNESS_REFINEMENT"],("PASS","NOT_ESTABLISHED"))
  if d["P5_SAMPLE1_DELTA_FORCE_SUBDIVIDED_WITNESS_REFINEMENT"]=="PASS": self.assertIsNone(d["first_unclosed_joint_cell"])
  else: self.assertIsNotNone(d["first_unclosed_joint_cell"])
 def test_structure(self):
  d=self.d
  self.assertTrue(d["sample0_canonical_correction_magnitude_subdivided"])
  self.assertTrue(d["sample1_force_component_cube_subdivided"])
  self.assertTrue(d["sample1_H_theta_exact_skew_structure_retained"])
  self.assertTrue(d["sample1_J_aw_exact_identity_in_transported_body_gauge"])
 def test_numerics(self):
  d=self.d; self.assertTrue(math.isfinite(d["max_sample1_residual_norm_upper_mps2"])); self.assertTrue(math.isfinite(d["max_sample1_correction_norm_upper_rad"])); self.assertGreaterEqual(d["fixed_pivot_inverse_count"],0); self.assertGreaterEqual(d["spectral_fallback_inverse_count"],0)
 def test_no_promotion(self):
  d=self.d; self.assertEqual(d["deployed_correction_limit_rad"],6.0); self.assertFalse(d["deployed_correction_limit_increased"]); self.assertFalse(d["complete_source_cell_refined_here"]); self.assertFalse(d["whole_word_promoted_here"]); self.assertFalse(d["N_H_words_set_here"])
if __name__=="__main__": unittest.main()
