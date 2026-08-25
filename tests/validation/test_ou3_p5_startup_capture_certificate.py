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
            "first_obstruction": cls.d.get("first_obstruction"),
            "outer_bridge_requirements": cls.d.get("outer_bridge_requirements"),
        }, sort_keys=True), flush=True)

    def test_current_local_P5_obstruction_is_source_certified(self):
        d = self.d
        self.assertEqual(P5.validate(d), [])
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertEqual(d["P1_STARTUP_CERTIFICATE"], "PASS")
        self.assertEqual(d["P4_EXACT_NONLINEAR_WORD_CERTIFICATE"], "PASS")
        self.assertEqual(d["heading_branch_contract"], "PASS")
        self.assertEqual(d["P5_OBSTRUCTION_IDENTIFIED"], "PASS")
        self.assertEqual(d["P5_OUTER_BRIDGE_REQUIREMENTS_IDENTIFIED"], "PASS")
        self.assertEqual(d["P5_FINITE_CAPTURE_CERTIFICATE"], "NOT_ESTABLISHED")
        self.assertEqual(d["first_obstruction"], "P1_HANDOFF_OUTSIDE_P4_CERTIFIED_CAPTURE_DOMAIN")
        self.assertFalse(d["first_required_P5_inequality_holds"])
        self.assertFalse(d["finite_capture_iteration_permitted"])
        self.assertIsNone(d["N_H_words"])
        self.assertEqual(d["completion_object"], "ou3_p5_outer_h_bridge_certificate.py")

    def test_P4_inner_seed_is_inside_but_startup_witnesses_are_outside_capture_domain(self):
        d = self.d
        self.assertGreater(d["P4_H_inner_level_W"], 0.0)
        self.assertGreater(d["P4_H_strict_decrease_W_threshold_lower"], d["P4_H_inner_level_W"])
        self.assertGreater(d["P4_inner_seed_to_decrease_threshold_W_factor"], 15.0)
        self.assertLess(d["P4_inner_seed_to_decrease_threshold_W_factor"], 17.0)
        rows = d["axis_witnesses"]
        self.assertEqual({r["group"] for r in rows}, {"b_g", "v", "p", "S", "a_w"})
        for row in rows:
            self.assertTrue(row["outside_P4_inner_seed"])
            self.assertTrue(row["outside_P4_strict_decrease_domain"])
            self.assertTrue(row["outside_P4_nonlinear_design_radius"])

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

    def test_gauged_handoffs_use_composed_full_attitude_not_tilt_only(self):
        b = self.d["outer_bridge_requirements"]
        self.assertTrue(b["P1_gravity_cosines_are_tilt_only"])
        normal = b["normal_gauged_full_attitude_cayley_norm_upper"]
        timeout = b["timeout_gauged_full_attitude_cayley_norm_upper"]
        self.assertGreater(normal, 0.20)
        self.assertLess(normal, 0.35)
        self.assertGreater(timeout, 0.50)
        self.assertLess(timeout, 0.70)
        self.assertGreater(timeout, normal)
        self.assertTrue(b["normal_gauged_inside_current_promoted_cayley_norm_limit"])
        self.assertTrue(b["timeout_gauged_inside_current_promoted_cayley_norm_limit"])
        self.assertGreater(b["normal_gauged_over_current_P4_design_radius_factor"], 1.0e11)
        self.assertGreater(b["timeout_gauged_over_current_P4_design_radius_factor"], 5.0e11)

    def test_ungauged_timeout_has_no_fake_full_heading_radius(self):
        b = self.d["outer_bridge_requirements"]
        self.assertFalse(b["timeout_ungauged_full_heading_cayley_bound_available"])
        self.assertIn("YAW_QUOTIENT", b["timeout_ungauged_required_route"])
        self.assertIn("yaw-quotient", " ".join(b["required_proof_structure"]))

    def test_unchanged_isotropic_P4_recurrence_is_not_outer_bridge(self):
        b = self.d["outer_bridge_requirements"]
        self.assertGreater(b["uniform_B_reduction_factor_needed_at_weakest_witness"], 1.0e67)
        self.assertGreater(b["uniform_B_reduction_factor_needed_at_largest_witness"], 1.0e72)
        self.assertIn("Simply enlarging q_design", b["interpretation"])
        self.assertIn("staged outer-H bridge", self.d["required_next_certificate"])


if __name__ == "__main__":
    unittest.main()
