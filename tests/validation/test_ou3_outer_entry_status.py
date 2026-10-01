import unittest
from tools.stability.ou3_theorem.outer_entry_status import certificate

class T(unittest.TestCase):
    def test_outer_absorption_does_not_claim_inner_entry(self):
        r=certificate()
        self.assertTrue(r["compact_outer_release_storage_bound_exists"])
        self.assertTrue(r["finite_practical_absorbing_storage_radius_exists"])
        self.assertEqual(r["local_target_storage"],"9/400")
        self.assertFalse(r["strict_entry_inequality_certified"])
        self.assertFalse(r["entry_into_local_ball_proved"])

if __name__=="__main__":
    unittest.main()
