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
            "W_inner": cls.d.get("P4_H_inner_level_W"),
            "W_capture": cls.d.get("P4_H_strict_decrease_W_threshold_lower"),
            "q_design": cls.d.get("P4_H_nonlinear_design_canonical_norm_radius"),
            "weakest_axis_witness": cls.d.get("weakest_axis_witness"),
            "largest_axis_witness": cls.d.get("largest_axis_witness"),
        }, sort_keys=True), flush=True)

    def test_current_P5_first_obstruction_is_source_certified(self):
        d = self.d
        self.assertEqual(P5.validate(d), [])
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertEqual(d["P1_STARTUP_CERTIFICATE"], "PASS")
        self.assertEqual(d["P4_EXACT_NONLINEAR_WORD_CERTIFICATE"], "PASS")
        self.assertEqual(d["P5_OBSTRUCTION_IDENTIFIED"], "PASS")
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
        self.assertGreaterEqual(weak["axis_witness_canonical_norm"], 0.01)
        self.assertGreater(weak["W_lower_over_strict_decrease_threshold"], 1.0e100)
        strong = d["largest_axis_witness"]
        self.assertEqual(strong["group"], "S")
        self.assertGreater(strong["W_lower_over_strict_decrease_threshold"], weak["W_lower_over_strict_decrease_threshold"])

    def test_P5_does_not_fake_a_word_count_outside_P4_domain(self):
        d = self.d
        self.assertIn("would extrapolate a local certificate outside its proof domain", d["reason_iteration_is_not_permitted"])
        self.assertIn("outer capture bridge", d["required_next_certificate"])
        self.assertIn("do not proceed", d["next_obligation"])


if __name__ == "__main__":
    unittest.main()
