import unittest
from tools.stability.ou3_alt_contraction import finite_master_guard as G

class Tests(unittest.TestCase):
    def test_closed_progress_is_consumed_but_storage_stays_blocked(self):
        x=G.build()
        self.assertEqual(G.validate(x),[])
        for v in x['closed_subobligations'].values(): self.assertTrue(v)
        self.assertIn('WPE_log_and_exp_libm_correspondence',x['open_obligations'])
        self.assertIn('complete_source_uniform_600_step_word',x['open_obligations'])
        self.assertNotIn('all_admitted_ungauged_entries_have_a_finite_attitude_representation',x['open_obligations'])
        self.assertTrue(x['finite_storage_status']['all_admitted_startup_entries_represented'])
        self.assertTrue(x['closed_subobligations']['all_nonzero_fresh_attitudes_represented_by_joint24_atlas'])
        self.assertNotIn('signed_magnetic_counter_safety_on_the_finite_word',x['open_obligations'])
        self.assertEqual(x['falsified_prerequisites'],{})
        self.assertTrue(x['finite_storage_status']['finite_word_counter_safety_closed'])
        self.assertEqual(x['research_outcome'],'finite_master_qualification_incomplete')
        self.assertTrue(x['open_obligations']['startup_accumulation_to_handoff_full_frame_bound'])
        self.assertEqual(x['conditional_timeout_plus_word_last_sample'],30602)
        self.assertEqual(x['startup_entry_obstruction']['Cayley_denominator_at_south'],0)
        self.assertFalse(x['finite_master_guard_closed'])
        self.assertFalse(x['storage_search_allowed'])
        self.assertIn('finite-state storage blocked',x['finite_storage_guard_error'])

    def test_closed_counter_prerequisite_must_reach_master_status(self):
        x=G.build()
        x['finite_storage_status']['finite_word_counter_safety_closed']=False
        self.assertIn('proved counter safety not consumed by storage guard',G.validate(x))

if __name__=='__main__': unittest.main()
