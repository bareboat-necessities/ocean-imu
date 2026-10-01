import unittest
from fractions import Fraction as F
from tools.stability.ou3_theorem.aw_signed_mean_reader_diagnostic import certificate

class T(unittest.TestCase):
    def test_carried_stress_suite_fits_new_reader_budget(self):
        r=certificate()
        self.assertEqual(F(r["worst_carried_signed_mean_error_mps2"]),F("0.371143"))
        self.assertTrue(r["carried_feasible"])
        self.assertGreater(F(r["carried_margin_mps2"]),F("0.524"))
        self.assertFalse(r["source_uniform_verified"])
        self.assertFalse(r["theorem_closed"])

if __name__=="__main__":
    unittest.main()
