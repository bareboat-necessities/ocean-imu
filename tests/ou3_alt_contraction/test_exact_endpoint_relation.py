import unittest
from tools.stability.ou3_alt_contraction import exact_endpoint_relation as E


class ExactEndpointRelationTest(unittest.TestCase):
    def test_exact_relation_is_universal_not_enumerated(self):
        d=E.build(); self.assertEqual(E.validate(d),[])
        self.assertTrue(d['exact_universal_endpoint_map_relation_defined'])
        self.assertTrue(d['all_branch_successors_retained'])
        self.assertTrue(d['all_bias_lineages_retained'])
        self.assertEqual(d['endpoint_map_set_symbol'],'W_ALT_joint24(O^601_BRMM)')
        self.assertFalse(d['finite_source_enumeration_used'])
        self.assertFalse(d['trajectory_replay_used'])
        self.assertFalse(d['independent_sample_boxes_used'])
        self.assertFalse(d['numeric_interval_endpoint_representation_closed'])
        self.assertFalse(d['common_M_source_uniform_projected_LDLT_closed'])
        self.assertFalse(d['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
