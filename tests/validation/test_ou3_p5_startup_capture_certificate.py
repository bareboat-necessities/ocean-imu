import json
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_startup_capture_certificate as P5


class Ou3P5StartupCaptureCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = P5.build()
        print("P5_CAPTURE_IDENTIFICATION=" + json.dumps({
            "P5_FINITE_CAPTURE_CERTIFICATE": cls.d.get("P5_FINITE_CAPTURE_CERTIFICATE"),
            "P5_OBSTRUCTION_IDENTIFIED": cls.d.get("P5_OBSTRUCTION_IDENTIFIED"),
            "P5_OUTER_BRIDGE_REQUIREMENTS_IDENTIFIED": cls.d.get("P5_OUTER_BRIDGE_REQUIREMENTS_IDENTIFIED"),
            "first_obstruction": cls.d.get("first_obstruction"),
            "W_inner": cls.d.get("P4_H_inner_level_W"),
            "W_capture": cls.d.get("P4_H_strict_decrease_W_threshold_lower"),
            "q_design": cls.d.get("P4_H_nonlinear_design_canonical_norm_radius"),
            "weakest_axis_witness": cls.d.get("weakest_axis_witness"),
            "largest_axis_witness": cls.d.get("largest_axis_witness"),
            "outer_bridge_requirements": cls.d.get("outer_bridge_requirements"),
        }, sort_keys=True), flush=True)

    def test_current_P5_first_obstruction_is_source_certified(self):
        d = self.d
        self.assertEqual(P5.validate(d), [])
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertEqual(d["P1_STARTUP_CERTIFICATE"], "PASS")
        self.assertEqual(d["P4_EXACT_NONLINEAR_WORD_CERTIFICATE"], "PASS")
        self.assertEqual(d["P5_OBSTRUCTION_IDENTIFIED"], "PASS")
        self.assertEqual(d["P5_OUTER_BRIDGE_REQUIREMENTS_IDENTIFIED"], "PASS")
        self.assertEqual(d["P5_FINITE_CAPTURE_CERTIFICATE"], "NOT_ESTABLISHED")
        self.assertEqual(d["first_obstruction"], "P1_HANDOFF_OUTSIDE_P4_CERTIFIED_CAPTURE_DOMAIN")
        self.assertFalse(d["first_required_P5_inequality_holds"])
        self.assertFalse(d["finite_capture_iteration_permitted"])
        self.assertIsNone(d["N_H_words"])

    def test_P4_inner_seed_is_inside_but_startup_witnesses_are_outside_capture_domain(self):
        d = self.d
        self.assertGreater(d["P4_H_inner_level_W"], 0.0)
        self.assertGreater(d["P4_H_strict_decrease_W_threshold_lower"], d["P4_H_inner_level_W"])
        self.assertGreater(d["P4_inner_seed_to_decrease_threshold_W_factor"], 15.0)
        self.assertLess(d["P4_inner_seed_to_decrease_threshold_W_factor"], 17.0)
        self.assertAlmostEqual(d["P4_H_strict_decrease_W_threshold_lower"], 5.26676120278833e-140, delta=1e-153)
        self.assertAlmostEqual(d["P4_H_nonlinear_design_canonical_norm_radius"], 8.650578521054014e-13, delta=1e-25)
        rows = d["axis_witnesses"]
        self.assertEqual({r["group"] for r in rows}, {"b_g", "v", "p", "S", "a_w"})
        for row in rows:
            self.assertTrue(row["outside_P4_inner_seed"])
            self.assertTrue(row["outside_P4_strict_decrease_domain"])
            self.assertTrue(row["outside_P4_nonlinear_design_radius"])
            self.assertGreater(row["W_lower_over_strict_decrease_threshold"], 1.0)

    def test_even_weakest_declared_axis_witness_misses_capture_by_many_orders(self):
        d = self.d
        weak = d["weakest_axis_witness"]
        self.assertEqual(weak["group"], "b_g")
        self.assertAlmostEqual(weak["axis_witness_W_lower"], 1.0e-4, delta=1.0e-18)
        self.assertGreater(weak["W_lower_over_strict_decrease_threshold"], 1.0e135)
        strong = d["largest_axis_witness"]
        self.assertEqual(strong["group"], "S")
        self.assertAlmostEqual(strong["axis_witness_W_lower"], 9.0e4, delta=1.0e-8)
        self.assertGreater(strong["W_lower_over_strict_decrease_threshold"], 1.0e144)

    def test_normal_and_timeout_handoffs_are_finite_angle_but_far_outside_q_design(self):
        b = self.d["outer_bridge_requirements"]
        normal = b["normal_handoff_cayley_norm_upper"]
        timeout = b["timeout_handoff_cayley_norm_upper"]
        self.assertGreater(normal, 0.09)
        self.assertLess(normal, 0.10)
        self.assertGreater(timeout, 0.40)
        self.assertLess(timeout, 0.42)
        self.assertGreater(timeout, normal)
        self.assertTrue(b["normal_inside_current_promoted_cayley_norm_limit"])
        self.assertTrue(b["timeout_inside_current_promoted_cayley_norm_limit"])
        self.assertGreater(b["normal_over_current_P4_design_radius_factor"], 1.0e11)
        self.assertGreater(b["timeout_over_current_P4_design_radius_factor"], 4.0e11)

    def test_unchanged_isotropic_P4_recurrence_is_not_a_plausible_outer_bridge(self):
        b = self.d["outer_bridge_requirements"]
        self.assertGreater(b["uniform_B_reduction_factor_needed_at_weakest_witness"], 1.0e67)
        self.assertGreater(b["uniform_B_reduction_factor_needed_at_largest_witness"], 1.0e72)
        self.assertGreater(
            b["uniform_B_reduction_factor_needed_at_largest_witness"],
            b["uniform_B_reduction_factor_needed_at_weakest_witness"],
        )
        self.assertIn("Simply enlarging q_design", b["interpretation"])
        text = " ".join(b["required_proof_structure"])
        self.assertIn("exact SO(3)", text)
        self.assertIn("anisotropic", text)
        self.assertIn("source-node subdivision", text)
        self.assertIn("overlaps the existing P4 inner seed", text)

    def test_P5_does_not_fake_a_word_count_outside_P4_domain(self):
        d = self.d
        self.assertIn("would extrapolate a local certificate outside its proof domain", d["reason_iteration_is_not_permitted"])
        self.assertIn("outer H capture bridge", d["required_next_certificate"])
        self.assertIn("do not proceed", d["next_obligation"])


if __name__ == "__main__":
    unittest.main()
