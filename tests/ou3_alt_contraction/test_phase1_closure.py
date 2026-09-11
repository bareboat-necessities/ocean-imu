import unittest
from tools.stability.ou3_alt_contraction import phase1_closure as P

class Phase1ClosureTests(unittest.TestCase):
    def test_phase1_closes_master_graph_but_not_stability(self):
        d=P.build();self.assertEqual(P.validate(d),[])
        self.assertTrue(all(d['phase1_status'].values()))
        self.assertTrue(d['storage_search_allowed'])
        self.assertTrue(d['common_joint24_storage_must_be_searched_first'])
        self.assertFalse(d['parameter_dependent_storage_search_allowed_initially'])
        self.assertFalse(d['release_guard_eventual_reachability_closed'])
        self.assertFalse(d['startup_capture_closed'])
        self.assertFalse(d['every_prefix_chart_retention_closed'])
        self.assertFalse(d['deployment_finite_precision_closed'])
        self.assertFalse(d['ALT_LIVE_PASS'])
        self.assertFalse(d['ALT_END_TO_END_PASS'])

if __name__=='__main__':unittest.main()
