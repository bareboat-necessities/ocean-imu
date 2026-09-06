#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_sea3_full_word_reset_congruence as RESET


class Sea3FullWordResetCongruenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = RESET.self_test()
        cls.failures = RESET.validate()

    def test_shipping_formula_and_nonsingularity(self):
        self.assertEqual([], self.failures)
        self.assertTrue(self.d["source_parity_pass"])
        self.assertTrue(self.d["reset_determinant_lower_at_least_one"])
        self.assertFalse(self.d["small_angle_needed_for_nonsingularity"])

    def test_joint_P_Psi_Omega_and_margin_congruence(self):
        self.assertTrue(self.d["joint_decomposition_identity_enclosed"])
        self.assertTrue(self.d["M_delta_reset_congruence_identity_enclosed"])

    def test_prior_free_factorization_survives_reset(self):
        self.assertTrue(self.d["prior_free_D_unchanged_by_reset"])
        self.assertTrue(self.d["prior_free_reconstruction_P_identity_enclosed"])
        self.assertTrue(self.d["prior_free_reconstruction_Psi_identity_enclosed"])
        self.assertTrue(self.d["prior_free_reconstruction_Omega_identity_enclosed"])

    def test_reset_primitive_cannot_promote_P3(self):
        self.assertFalse(self.d["P3_promoted"])


if __name__ == "__main__":
    unittest.main()
