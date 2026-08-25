import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_outer_information_geometry as G


class Ou3P5OuterInformationGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build()

    def test_two_gauged_nodes_have_strict_finite_angle_information(self):
        d = self.d
        self.assertEqual(G.validate(d), [])
        self.assertEqual(d["P5_FINITE_ANGLE_INFORMATION_GEOMETRY_CERTIFICATE"], "PASS")
        n = d["nodes"]["normal_gauged"]
        t = d["nodes"]["timeout_gauged"]
        self.assertGreater(n["exact_pair_residual_information_vs_goLive_attitude_metric_lower"], 0.0)
        self.assertGreater(t["exact_pair_residual_information_vs_goLive_attitude_metric_lower"], 0.0)
        self.assertLess(n["cayley_norm_upper"], t["cayley_norm_upper"])
        self.assertGreater(n["exact_cayley_residual_factor_lower"], t["exact_cayley_residual_factor_lower"])

    def test_geometry_is_source_correlated_and_uses_actual_handoff_covariance(self):
        d = self.d
        self.assertTrue(d["actual_goLive_covariance_tuple_used"])
        self.assertTrue(d["joint_source_tuple_required"])
        self.assertFalse(d["independent_gain_extrema_product_used"])
        self.assertTrue(d["identity_requires_source_correlated_H_P_R_K_S"])
        info = d["goLive_attitude_information"]
        self.assertGreater(info["lambda_min_lower"], 0.0)
        self.assertGreaterEqual(info["lambda_max_upper"], info["lambda_min_lower"])

    def test_exact_finite_angle_factor_not_small_angle_taylor(self):
        d = self.d
        self.assertIn("4||c||^2/(4+||c||^2)", d["exact_cayley_vector_residual_identity"])
        p = d["packet_geometry"]
        self.assertGreater(p["linear_pair_information_mu_lower"], 0.0)
        self.assertGreater(p["angular_factor_lower"], 0.0)
        self.assertGreater(p["a_f_lower"], 0.0)
        self.assertGreater(p["a_m_lower"], 0.0)

    def test_exact_joseph_information_identity_is_retained(self):
        ident = self.d["exact_joseph_tangent_information_identity"]
        self.assertIn("S^-1", ident)
        self.assertIn("eta^T R^-1 eta", ident)

    def test_complete_outer_word_is_fail_closed_until_exact_transport(self):
        d = self.d
        self.assertEqual(d["P5_GAUGED_OUTER_CAYLEY_INFORMATION_WORD_SECTOR"], "NOT_ESTABLISHED")
        self.assertIn("deployed Cayley/quaternion", d["remaining_word_term"])
        self.assertIn("S-to-attitude", d["remaining_word_term"])
        text = (ROOT / "tools" / "ou3_p5_outer_information_geometry.py").read_text(encoding="utf-8")
        self.assertNotIn("ou3_exact_replay", text)
        self.assertNotIn("numpy", text)


if __name__ == "__main__":
    unittest.main()
