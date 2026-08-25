import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_large_angle_sector_certificate as SECTOR


class Ou3P5LargeAngleSectorCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = SECTOR.build()

    def test_actual_source_correlated_tuple_is_used(self):
        d = self.d
        self.assertEqual(SECTOR.validate(d), [])
        self.assertTrue(d["source_correlated_covariance_gain_tuple"])
        self.assertTrue(d["validated_interval_Kalman_gain_used"])
        self.assertTrue(d["exact_deployed_quaternion_backend_used"])
        self.assertTrue(d["exact_Rodrigues_backend_used"])
        self.assertFalse(d["independent_gain_extrema_used"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])

    def test_counterexample_lies_inside_both_gauged_nodes(self):
        w = self.d["validated_counterexample"]
        self.assertTrue(w["source_admissible"])
        self.assertEqual(w["xi_norm"], 0.0)
        q = w["cayley_norm_interval"][1]
        self.assertLess(q, w["normal_gauged_q_upper"])
        self.assertLess(q, w["timeout_gauged_q_upper"])
        self.assertGreaterEqual(w["magnetic_sine_separation"], 0.1)
        self.assertGreaterEqual(w["specific_force_norm_mps2"], 5.0)
        self.assertLessEqual(w["specific_force_norm_mps2"], 30.0)

    def test_raw_VR_positive_alpha_sector_is_strictly_disproved(self):
        d = self.d
        w = d["validated_counterexample"]
        self.assertGreater(w["V_R_before_interval"][0], 0.0)
        self.assertLess(w["D_R_deployed_interval"][1], 0.0)
        self.assertLess(w["D_R_Rodrigues_interval"][1], 0.0)
        self.assertTrue(w["deployed_energy_increases_strictly"])
        self.assertTrue(w["Rodrigues_energy_increases_strictly"])
        self.assertTrue(w["requested_positive_alpha_sector_disproved"])
        self.assertTrue(d["beta_cannot_repair_xi_zero_counterexample"])
        self.assertEqual(
            d["P5_RAW_VR_LARGE_ANGLE_SECTOR"],
            "DISPROVED_ON_DECLARED_SOURCE_FAMILY",
        )

    def test_a_false_sector_is_not_promoted_by_increasing_beta(self):
        d = self.d
        self.assertIn("source-shaped Cayley/information", d["required_theorem_correction"])
        self.assertIn("do not tune beta", d["next_obligation"])
        self.assertTrue(d["single_vector_positive_alpha_sector_also_impossible"])


if __name__ == "__main__":
    unittest.main()
