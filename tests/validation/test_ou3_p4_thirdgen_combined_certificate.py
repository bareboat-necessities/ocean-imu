from pathlib import Path
import os,sys,unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import ou3_p4_thirdgen_combined_certificate as C
import ou3_validate_enclosure as ENC
@unittest.skipUnless(os.environ.get('OU3_RUN_OBSOLETE_P4_FRONTIER')=='1','retired microscopic P4 frontier regression')
class ThirdGenCombinedP4Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=C.build()
 def test_combined_passes(self):
  self.assertEqual(C.validate(self.d),[]); self.assertEqual(self.d['P4_THIRDGEN_COMBINED_WORD_CERTIFICATE'],'PASS'); self.assertTrue(self.d['thirdgen_stacks_431_and_432_refinements'])
 def test_combined_never_regresses_431(self):
  for mode in ('H','A'):
   m=self.d['modes'][mode]; self.assertGreaterEqual(m['certified_level_W'],m['certified_level_W_before_combined']); self.assertGreaterEqual(m['thirdgen_W_widening_factor_vs_431_lower'],1.0-2e-15); self.assertGreater(m['thirdgen_total_W_widening_factor_vs_legacy_lower'],1.0)
 def test_radius_continuation_is_real_and_maximal_on_grid(self):
  for mode in ('H','A'):
   m=self.d['modes'][mode]; self.assertTrue(m['thirdgen_combined_radius_continuation']); self.assertGreater(m['thirdgen_candidate_count_certified'],0); self.assertEqual(m['certified_level_W'],max(x['W'] for x in m['thirdgen_candidates']))
 def test_exact_prefix_and_safety(self):
  for mode in ('H','A'):
   m=self.d['modes'][mode]; self.assertGreaterEqual(m['thirdgen_exact_prefix_factor_upper'],1.0); self.assertLess(m['prefix_canonical_error_norm_upper'],m['cayley_norm_limit']); self.assertLess(m['accepted_correction_norm_prefix_upper'],1e-2)
 def test_schema4_compatibility(self):
  for mode in ('H','A'):
   o=ENC.validate_mode(mode,self.d['modes'][mode],{'required_path_metric':'CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC'}); self.assertTrue(o['linear_pass'],o['failures']); self.assertTrue(o['nonlinear_pass'],o['failures'])
 def test_A_projection_interior(self):
  self.assertFalse(self.d['modes']['A']['active_bias_projection']['projection_surface_reached_in_certified_funnel'])
 def test_source_only(self):
  text=(ROOT/'tools'/'ou3_p4_thirdgen_combined_certificate.py').read_text().casefold(); self.assertNotIn('monte carlo',text); self.assertNotIn('sampled trajectory',text); self.assertNotIn('numpy',text)
if __name__=='__main__':unittest.main()
