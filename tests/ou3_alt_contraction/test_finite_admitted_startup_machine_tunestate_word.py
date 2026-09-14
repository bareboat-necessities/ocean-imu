import unittest
from tools.stability.ou3_alt_contraction import finite_admitted_startup_machine_tunestate_word as X

class Tests(unittest.TestCase):
    def test_readiness_closes_boundary_attachment_but_not_universal_startup(self):
        r=X.readiness()
        self.assertTrue(r['startup_to_Live_machine_history_attachment_closed'])
        self.assertTrue(r['guard_private_Mahony_WPE_band_stats_stillness_machine_memory_crosses_goLive_without_reseed'])
        self.assertTrue(r['whole_tau_sigma_RS_machine_TuneState_and_applied_parameters_cross_goLive'])
        self.assertTrue(r['same_goLive_state_seeds_admitted_600_edge_machine_product'])
        self.assertFalse(r['every_admitted_startup_history_reaches_this_boundary_product'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])

if __name__=='__main__': unittest.main()
