import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
MOD = ROOT / "tools" / "stability" / "ou3_innovation_psd_plus_R_inverse.py"
spec = importlib.util.spec_from_file_location("ou3_innovation_psd_plus_R_inverse", MOD)
M = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(M)


class InnovationPsdPlusRInverseTests(unittest.TestCase):
    def test_pathological_rectangular_hull_is_not_the_physical_family(self):
        d = M.build()
        self.assertTrue(d["generic_rectangular_inverse_rejected"])
        self.assertTrue(d["provenance_aware_inverse_finite"])
        self.assertTrue(d["point_regression_inverses_contained"])
        self.assertFalse(d["P4_PASS"])
        self.assertFalse(d["P5_MAY_START"])

    def test_proof_records_structural_source_of_invertibility(self):
        d = M.build()
        p = d["proof"]
        self.assertEqual(p["retained_identity"], "S=H P H^T+R")
        self.assertEqual(p["derived_Q_property"], "H P H^T>=0")
        self.assertTrue(p["uniform_R_SPD_validated"])
        self.assertTrue(p["det_S_lower_from_Loewner_monotonicity"] > 0.0)
        self.assertTrue(p["singular_entrywise_hull_members_are_not_admitted"])
        self.assertFalse(p["same_history_K_dependency_closed_here"])

    def test_validation_is_fail_closed(self):
        d = M.build()
        self.assertEqual(M.validate(d), [])
        bad = dict(d)
        bad["P4_PASS"] = True
        self.assertTrue(M.validate(bad))


if __name__ == "__main__":
    unittest.main()
