from pathlib import Path
import math, sys, unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import ou3_p4_translation_full_word_design as D

class P4TranslationFullWordDesignTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=D.build()
 def test_probe_is_design_only_and_valid(self):
  self.assertEqual(D.validate(self.d),[])
  self.assertTrue(self.d['ordinary_floating_point_design_only'])
  self.assertFalse(self.d['validated_for_theorem_promotion'])
  self.assertEqual(self.d['P4_USABLE_CERTIFICATE_STATUS'],'NOT_ESTABLISHED')
 def test_complete_word_margin_is_reported(self):
  for mode in ('H','A'):
   for h,row in self.d['modes'][mode]['horizons'].items():
    q=row['worst_grid_point']['translation_complete_word_generalized_margin_design']
    self.assertTrue(math.isfinite(q) and q>0.0,(mode,h,q))
    self.assertGreater(row['design_worst_to_old_margin_ratio'],0.0)
 def test_no_replay(self): self.assertFalse(self.d['trajectory_replay_used'])

if __name__=='__main__': unittest.main()
