import unittest
from tools.stability.ou3_theorem.jump_twist_remainder import certificate

class T(unittest.TestCase):
    def test_source_vs_finite_error_classification(self):
        r=certificate()
        self.assertFalse(r["injection_squared_times_error_is_source_only"])
        self.assertEqual(r["normalized_reader_norm_ceiling"],"4")
        self.assertTrue(r["signed_twist_source_uniform_bound_certified"])
        self.assertFalse(r["theorem_closed"])

if __name__=="__main__":
    unittest.main()
