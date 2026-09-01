from pathlib import Path
import os, sys, unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import ou3_p4_full_word_dissipation_route as R

@unittest.skipUnless(
 os.environ.get('OU3_RUN_OBSOLETE_P4_FRONTIER') == '1',
 'retired microscopic P4 frontier regression runs only in its non-gating focused job',
)
class FullWordDissipationRouteTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=R.build()
 def test_route_validates(self):
  self.assertEqual(R.validate(self.d),[])
 def test_old_tiny_frontier_is_not_promoted(self):
  self.assertEqual(self.d['P4_USABLE_CERTIFICATE_STATUS'],'NOT_ESTABLISHED')
  self.assertTrue(self.d['old_scalar_frontier_is_regression_baseline_only'])
 def test_current_delta_is_seed_only(self):
  for mode in ('H','A'):
   m=self.d['modes'][mode]
   self.assertTrue(m['current_numeric_margin_is_single_seed_then_preserved'])
   self.assertFalse(m['later_additive_PSD_terms_quantified_in_current_delta'])
 def test_replacement_forbids_old_bottlenecks(self):
  self.assertTrue(self.d['path_metric_required'])
  self.assertTrue(self.d['source_reachability_required'])
  self.assertTrue(self.d['scalar_uniform_min_delta_for_all_cells_forbidden_as_final_route'])
  self.assertTrue(self.d['scalar_BW_small_gain_forbidden_as_final_route'])

if __name__=='__main__': unittest.main()
