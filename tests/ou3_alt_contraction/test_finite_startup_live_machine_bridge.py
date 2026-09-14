"""Whole-machine goLive bridge guard regressions."""
import unittest

from tools.stability.ou3_alt_contraction import finite_startup_live_machine_bridge as X


class Tests(unittest.TestCase):
    def test_bridge_rejects_non_startup_product(self):
        with self.assertRaisesRegex(TypeError,'whole-machine product'):
            X.bridge(object(),object(),object())

    def test_readiness_closes_identity_not_machine_commit_correspondence(self):
        r=X.readiness()
        for k in ('coherent_startup_WPE_tau_sigma_RS_product_consumed',
                  'goLive_executes_no_WPE_or_machine_TuneState_adaptation_recurrence',
                  'goLive_preserves_WPE_and_whole_machine_TuneState_by_identity',
                  'goLive_preserves_common_pending_bit_into_first_Live_boundary',
                  'construction_to_TunerReady_whole_machine_history_shape_available'):
            self.assertTrue(r[k])
        for k in ('machine_TuneState_to_goLive_active_parameters_binary32_correspondence_closed',
                  'every_admitted_startup_history_reaches_this_goLive_product',
                  'all_target_libm_correspondence_closed','storage_search_allowed',
                  'ALT_STARTUP_PASS','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
