import unittest
from tools.stability.ou3_alt_contraction import live_mahony_regional as M

class LiveMahonyRegionalTests(unittest.TestCase):
    def test_regional_invariance_closes_without_claiming_startup_entry(self):
        d=M.build();self.assertEqual(M.validate(d),[])
        self.assertTrue(d['regional_continuous_invariance_closed'])
        self.assertTrue(d['regional_source_order_binary32_discrete_invariance_closed'])
        self.assertGreater(d['binary32_robust_boundary_margin_lower'],0)
        self.assertGreater(d['discrete_V_margin_lower'],0)
        self.assertGreater(d['one_step_chart_round_margin_rad'],0)
        self.assertFalse(d['startup_entry_into_regional_set_closed_here'])
        self.assertFalse(d['startup_initial_seed_membership_required_for_Live_invariance'])
        self.assertFalse(d['ALT_LIVE_PASS'])

if __name__=='__main__':unittest.main()
