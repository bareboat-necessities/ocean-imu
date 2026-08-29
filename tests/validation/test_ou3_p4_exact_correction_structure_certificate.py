from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import ou3_p4_exact_correction_structure_certificate as C
import ou3_validate_enclosure as ENC
class P4ExactCorrectionStructureTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=C.build()
 def test_passes(self):
  self.assertEqual(C.validate(self.d),[]); self.assertEqual(self.d['P4_EXACT_CORRECTION_STRUCTURE_WORD_CERTIFICATE'],'PASS')
 def test_monotone_vs_thirdgen(self):
  for mode in ('H','A'):
   m=self.d['modes'][mode]; self.assertGreaterEqual(m['certified_level_W'],m['certified_level_W_before_exact_structure']); self.assertGreaterEqual(m['exact_structure_W_widening_factor_vs_thirdgen_lower'],1.0)
 def test_p5_exact_identities_are_bound(self):
  self.assertEqual(self.d['p5_effective_vector_input_certificate'],'PASS'); self.assertEqual(self.d['p5_exact_correction_transport_certificate'],'PASS')
  for mode in ('H','A'):
   m=self.d['modes'][mode]; self.assertTrue(m['p5_S_zero_eta_exact_zero_bound']); self.assertTrue(m['p5_magnetometer_radial_gain_action_exact_zero_bound']); self.assertTrue(m['p5_accelerometer_eta_effective_aw_input_bound']); self.assertFalse(m['reset_condition_number_multiplier_used'])
 def test_safety_and_schema4(self):
  for mode in ('H','A'):
   m=self.d['modes'][mode]; self.assertLess(m['accepted_correction_norm_prefix_upper'],1e-2); o=ENC.validate_mode(mode,m,{'required_path_metric':'CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC'}); self.assertTrue(o['nonlinear_pass'],o['failures'])
  self.assertFalse(self.d['modes']['A']['active_bias_projection']['projection_surface_reached_in_certified_funnel'])
if __name__=='__main__': unittest.main()
