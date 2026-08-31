from pathlib import Path
import os, sys, unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import ou3_p4_translation_full_word_interval_fast as C


@unittest.skipUnless(
 os.environ.get('OU3_RUN_EXPENSIVE_P4_FULL_WORD') == '1',
 'rigorous full-word P4 enclosure runs only in focused ou3-p4-frontier-combined CI',
)
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
 def test_measurement_information_keeps_translation_directions(self):
  for mode in ('H','A'):
   m=self.d['modes'][mode]
   self.assertEqual(m['measurement_information_geometry'],'rank_one_S_and_aw_each_sample_exact_rational')
   self.assertTrue(m['corrections_allowed_every_sample_for_lower_bound'])
   self.assertGreater(m['S_measurement_information_beta_conditioned'],0.0)
   self.assertGreater(m['accelerometer_aw_information_beta_conditioned'],0.0)
 def test_exact_transition_and_lower_are_retained_through_word(self):
  for mode in ('H','A'):
   m=self.d['modes'][mode]
   self.assertEqual(m['prediction_enclosure'],'exact_rational_transition_interval_rowwise_loewner_dyadic192')
   self.assertTrue(m['exact_rational_transition_enclosure'])
   self.assertTrue(m['exact_rational_lower_retained_through_word'])
   self.assertTrue(m['dyadic_loewner_compression'])
   self.assertEqual(m['dyadic_loewner_bits'],192)
   self.assertLess(m['dyadic_loewner_max_added_diagonal_loss'],1e-50)
 def test_source_only(self):
  self.assertTrue(self.d['source_only']); self.assertFalse(self.d['trajectory_replay_used']); self.assertTrue(self.d['outward_rounded'])

if __name__=='__main__': unittest.main()
