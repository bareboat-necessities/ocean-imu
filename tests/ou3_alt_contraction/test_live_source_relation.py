import unittest
from tools.stability.ou3_alt_contraction import live_source_relation as R

class LiveSourceRelationTests(unittest.TestCase):
    def test_regional_universal_relation_closes_without_startup(self):
        d=R.build();self.assertEqual(R.validate(d),[])
        self.assertEqual(d['physical_outer_relation'],'O^601_BRMM')
        self.assertEqual(d['sample_count'],601)
        self.assertTrue(d['COMPLETE_BRMM_left_inclusion_consumed'])
        self.assertTrue(d['regional_frontend_predecessor_invariant_consumed'])
        self.assertTrue(d['same_history_joint_frontend_transition_operator_bound'])
        self.assertTrue(d['trusted_Riccati_event_attachment_operator_bound'])
        self.assertTrue(d['all_branch_successors_retained_by_operator'])
        self.assertTrue(d['regional_universal_Normal_Live_source_relation_closed'])
        self.assertFalse(d['startup_capture_consumed'])
        self.assertFalse(d['startup_entry_membership_closed_here'])
        self.assertFalse(d['replay_or_seeded_realization_used'])
        self.assertFalse(d['complete_600_step_physical_joint24_cocycle_closed_here'])
        self.assertFalse(d['storage_search_allowed'])
        self.assertFalse(d['ALT_LIVE_PASS'])

if __name__=='__main__':unittest.main()
