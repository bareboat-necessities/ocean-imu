from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.block_factor_numeric import certificate


class NumericBlockFactorTests(unittest.TestCase):
    def test_partial_lin_certificate_cannot_promote_full_state_contraction(self):
        r=certificate()
        self.assertTrue(r['lin_matrix_verified'])
        self.assertTrue(r['full_word_energy_algebra_verified'])
        self.assertFalse(r['verified'])
        self.assertIsNone(r['mu_cov'])
        self.assertIsNone(r['rho0'])
        self.assertFalse(r['restricted_service_lifting_valid'])
        self.assertFalse(r['constructive_full_A21_mu_rho_enclosure'])
        self.assertFalse(r['float32_whole_word_supply_closed'])
        self.assertFalse(r['theorem_closed'])


if __name__=='__main__':
    unittest.main()
