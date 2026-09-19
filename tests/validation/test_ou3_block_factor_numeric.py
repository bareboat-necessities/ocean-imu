from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.block_factor_numeric import certificate, sharpened_factor_coordinate_certificate

class NumericBlockFactorTests(unittest.TestCase):
    def test_constructive_real_arithmetic_mu_rho(self):
        r=certificate()
        self.assertTrue(r["verified"])
        self.assertGreater(r["ell_LIN"],1.19e-6)
        self.assertGreater(r["mu_cov"],0.0)
        self.assertLess(r["rho0"],1.0)
        self.assertEqual(r["root_metric_gamma"],1.0)
        self.assertFalse(r["float32_whole_word_supply_closed"])
        self.assertFalse(r["theorem_closed"])\n\n    def test_sharpened_coordinates_improve_margin_without_fitting_rho(self):\n        old=certificate(); r=sharpened_factor_coordinate_certificate()\n        self.assertTrue(r["verified"])\n        self.assertGreater(r["mu_cov"],old["mu_cov"]*50)\n        self.assertGreater(r["mu_translation"],2.1e-3)\n        self.assertEqual(r["selection_role"],"fixed proof-coordinate congruence, not sampled/fitted rho")\n        self.assertFalse(r["float32_whole_word_supply_closed"])

if __name__=="__main__":unittest.main()
