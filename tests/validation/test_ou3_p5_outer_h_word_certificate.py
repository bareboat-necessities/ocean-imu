import json
import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_outer_h_word_certificate as OUTER


class Ou3P5OuterHWordCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = OUTER.build()
        print("P5_OUTER_H_WORD=" + json.dumps({
            "status": cls.d.get("P5_OUTER_H_WORD_CERTIFICATE"),
            "first_failure": cls.d.get("first_failure"),
            "normal": cls.d.get("node_word_tests", {}).get("normal"),
            "timeout": cls.d.get("node_word_tests", {}).get("timeout"),
        }, sort_keys=True), flush=True)

    def test_outer_backend_is_source_bound_and_not_local_P4_radius_extension(self):
        d = self.d
        self.assertEqual(OUTER.validate(d), [])
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertFalse(d["local_P4_BW_recurrence_reused_as_outer_certificate"])
        self.assertEqual(d["metric_route"], "SAME_NORMALIZED_CAYLEY_SOURCE_INFORMATION_GEOMETRY_AS_P4")

    def test_outer_backend_is_anisotropic_and_retains_S_to_attitude(self):
        d = self.d
        self.assertEqual(d["exact_linear_coordinates_not_charged_as_nonlinearity"], ["v", "p"])
        self.assertTrue(d["S_innovation_treated_as_exact_linear_selector"])
        self.assertTrue(d["full_S_to_attitude_gain_retained"])
        for name in ("normal", "timeout"):
            t = d["node_word_tests"][name]
            self.assertFalse(t["vector_defect"]["v_or_p_charged_as_vector_nonlinearity"])
            self.assertFalse(t["vector_defect"]["S_charged_as_vector_measurement_nonlinearity"])
            self.assertTrue(t["S_to_attitude_prefix"]["full_S_to_attitude_gain_retained"])
            self.assertTrue(t["S_to_attitude_prefix"]["S_innovation_is_exactly_linear"])

    def test_normal_and_timeout_nodes_use_actual_finite_angle_handoff_bounds(self):
        nodes = self.d["handoff_nodes"]
        qn = nodes["normal"]["rotation"]["cayley_norm_upper"]
        qt = nodes["timeout"]["rotation"]["cayley_norm_upper"]
        self.assertGreater(qn, 0.0)
        self.assertGreater(qt, qn)
        self.assertLess(qt, 1.0)
        for name in ("normal", "timeout"):
            r = nodes[name]["rotation"]
            self.assertGreaterEqual(r["R_minus_I_minus_skew_norm_upper"], 0.0)
            self.assertGreater(r["R_minus_I_norm_upper"], 0.0)
            self.assertTrue(math.isfinite(r["R_minus_I_minus_skew_linearized_ratio_upper"]))

    def test_word_test_uses_source_safe_branch_counts_and_direct_gap(self):
        for name in ("normal", "timeout"):
            t = self.d["node_word_tests"][name]
            c = t["word_counts"]
            self.assertGreaterEqual(c["accepted_accel_corrections_upper"], 200)
            self.assertGreater(c["accepted_mag_corrections_upper"], 0)
            self.assertGreater(c["S_zero_corrections_upper"], 0)
            self.assertGreater(t["P3_homogeneous_sqrt_decrease_lower"], 0.0)
            self.assertTrue(math.isfinite(t["outer_vector_nonlinear_information_ratio_upper"]))
            self.assertTrue(math.isfinite(t["outer_vector_word_decrease_margin_lower"]))

    def test_certificate_never_promotes_without_both_prefix_and_word_decrease(self):
        tests = self.d["node_word_tests"]
        expected = all(
            t["prefix_group_correction_pass"] and t["vector_word_decrease_pass"]
            for t in tests.values()
        )
        self.assertEqual(self.d["outer_word_decrease_all_nodes"], expected)
        self.assertEqual(
            self.d["P5_OUTER_H_WORD_CERTIFICATE"],
            "PASS" if expected else "NOT_ESTABLISHED",
        )
        if not expected:
            self.assertNotEqual(self.d["first_failure"], "NONE")
            self.assertIn("do not enlarge the local P4 radius", self.d["next_widening_if_not_established"])


if __name__ == "__main__":
    unittest.main()
