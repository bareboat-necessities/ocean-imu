import unittest
from tools.stability.ou3_alt_contraction import finite_master_guard as G

class Tests(unittest.TestCase):
    def test_closed_progress_is_consumed_but_storage_stays_blocked(self):
        x=G.build()
        self.assertEqual(G.validate(x),[])
        for v in x['closed_subobligations'].values(): self.assertTrue(v)
        self.assertIn('WPE_log_and_exp_libm_correspondence',x['open_obligations'])
        self.assertIn('complete_source_uniform_600_step_word',x['open_obligations'])
        self.assertFalse(x['finite_master_guard_closed'])
        self.assertFalse(x['storage_search_allowed'])
        self.assertIn('finite-state storage blocked',x['finite_storage_guard_error'])

if __name__=='__main__': unittest.main()
