from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import ou3_p4_translation_full_word_interval as C

class P4TranslationFullWordIntervalTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=C.build()
 def test_validated_worst_cell_translation_passes(self):
  self.assertEqual(C.validate(self.d),[])
  self.assertEqual(self.d['P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS'],'PASS')
  self.assertEqual(self.d['P4_USABLE_CERTIFICATE_STATUS'],'NOT_ESTABLISHED')
 def test_margin_widens_seed(self):
  for mode in ('H','A'):
   m=self.d['modes'][mode]
   self.assertTrue(m['interval_ldlt_endpoint_recertified'])
   self.assertGreater(m['complete_word_translation_margin_lower'],m['old_single_seed_translation_margin_lower'])
   self.assertGreater(m['margin_widening_factor_lower'],1.0)
 def test_source_only(self):
  self.assertTrue(self.d['source_only']); self.assertFalse(self.d['trajectory_replay_used']); self.assertTrue(self.d['outward_rounded'])

if __name__=='__main__': unittest.main()
