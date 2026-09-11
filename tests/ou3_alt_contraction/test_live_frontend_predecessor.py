import unittest
from tools.stability.ou3_alt_contraction import live_frontend_predecessor as F

class LiveFrontendPredecessorTests(unittest.TestCase):
    def test_live_predecessor_closes_without_startup_capture(self):
        d=F.build();self.assertEqual(F.validate(d),[])
        self.assertTrue(d['regional_normal_live_frontend_predecessor_set_closed'])
        self.assertTrue(d['regional_Live_entry_membership_is_theorem_premise'])
        self.assertTrue(d['regional_Mahony_invariant_consumed'])
        self.assertFalse(d['startup_capture_consumed'])
        self.assertFalse(d['startup_entry_closed_here'])
        self.assertTrue(all(d['strict_invariance_checks'].values()))
        self.assertFalse(d['ALT_LIVE_PASS'])

if __name__=='__main__':unittest.main()
