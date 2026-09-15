import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_startup_wpe_machine_word as X


class Tests(unittest.TestCase):
    def test_readiness_closes_attachment_not_universal_arithmetic(self):
        r=X.readiness()
        for k in (
            'literal_reset_full_WPE_machine_history_attached_through_startup',
            'same_private_Mahony_vertical_drives_startup_WPE_machine',
            'full_WPE_machine_history_preserved_by_same_goLive_edge',
            'admitted_Live_constructor_accepts_only_startup_produced_WPE_state',
            'admitted_source_private_Mahony_startup_to_Live_invariant_closed',
            'startup_guard_Mahony_WPE_machine_to_admitted_Live_attachment_closed'):
            self.assertTrue(r[k],k)
        for k in (
            'every_admitted_startup_history_reaches_this_product',
            'source_uniform_startup_WPE_supply_bounds_closed',
            'target_WPE_libm_and_compiler_profile_correspondence_closed',
            'source_uniform_complete_600_step_word_qualified',
            'storage_search_allowed','ALT_STARTUP_PASS','ALT_LIVE_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[k],k)


if __name__=='__main__': unittest.main()
