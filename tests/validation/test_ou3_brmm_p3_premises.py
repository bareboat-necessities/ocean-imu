from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/stability"))
import ou3_brmm_p3_premises as premises


class BrmmP3PremisesTest(unittest.TestCase):
    def test_execution_regime_is_not_runtime_flag_or_physical_recurrence(self):
        d = premises.build()
        self.assertEqual([], premises.validate(d))
        for key in ("runtime_Live_flag_implies_all_premises", "BRMM_alone_implies_vector_PE",
                    "P3_conclusion_is_an_assumption", "bias_error_decay_required",
                    "spectral_membership_required", "physical_execution_admission_proved_here"):
            bad = deepcopy(d)
            bad[key] = True
            self.assertTrue(premises.validate(bad), key)

    def test_projection_covariance_write_cannot_inherit_identity_proof(self):
        text = premises.MEKF.read_text()
        self.assertTrue(premises.projection_covariance_parity(text))
        mutated = text.replace("b *= (acc_bias_limit_ / n);", "b *= (acc_bias_limit_ / n); P.setIdentity();")
        self.assertFalse(premises.projection_covariance_parity(mutated))

    def test_all_branches_and_modes_are_needed_without_physical_promotion(self):
        d = premises.conditional_coverage(True, True)
        self.assertTrue(premises.coverage_closed(d))
        for branch in premises.BRANCHES:
            missing = deepcopy(d)
            del missing[branch]
            self.assertFalse(premises.coverage_closed(missing))
            for key in ("H18", "A21", "physical_admission_certified"):
                bad = deepcopy(d)
                bad[branch][key] = not bad[branch][key]
                self.assertFalse(premises.coverage_closed(bad))


if __name__ == "__main__":
    unittest.main()
