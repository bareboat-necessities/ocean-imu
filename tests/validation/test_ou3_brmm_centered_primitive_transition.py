import copy
import pathlib
import sys
import unittest
from fractions import Fraction as F

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_brmm_centered_primitive_transition as P


class CenteredPrimitiveTransitionTests(unittest.TestCase):
    def test_exact_transition_and_status(self):
        d = P.build()
        self.assertEqual(P.validate(d), [])
        self.assertTrue(d["primitive_transition_operator_materialized"])
        self.assertFalse(d["same_history_acceleration_moments"]["may_be_selected_independently"])
        self.assertTrue(d["cross_word_primitive_out_is_next_word_primitive_in_required"])
        self.assertTrue(d["centered_S_origin_id_is_invariant_under_advance"])
        self.assertIsNone(d["D_S_numeric_value_used"])
        self.assertFalse(d["P4_PASS"])
        self.assertEqual(d["P3_delta"], 1e-18)

    def test_zero_moment_map_reproduces_translation_kernel_polynomials(self):
        z = (F(0), F(0), F(0))
        root = P.PrimitiveState("a", "origin", (F(1), F(0), F(0)), z, z)
        m = P.MomentWitness("tr", z, z, z)
        child = P.advance(root, m, next_primitive_id="b", h=F(2))
        self.assertEqual(child.v[0], F(1))
        self.assertEqual(child.p[0], F(2))
        self.assertEqual(child.S_L[0], F(2))

        root = P.PrimitiveState("a", "origin", z, (F(1), F(0), F(0)), z)
        child = P.advance(root, m, next_primitive_id="b", h=F(2))
        self.assertEqual(child.p[0], F(1))
        self.assertEqual(child.S_L[0], F(2))

    def test_origin_witness_cannot_change(self):
        z = (F(0), F(0), F(0))
        root = P.PrimitiveState("a", "origin", z, z, z)
        child = P.advance(root, P.MomentWitness("tr", z, z, z), next_primitive_id="b")
        self.assertEqual(child.centered_S_origin_witness_id, "origin")

    def test_mutated_transition_matrix_is_rejected(self):
        d = P.build()
        d["exact_scalar_axis_transition_matrix"][2][0] = "0"
        self.assertIn("exact primitive transition matrix changed", P.validate(d))

    def test_false_promotion_and_independent_moments_are_rejected(self):
        d = P.build()
        d["P4_PASS"] = True
        d["same_history_acceleration_moments"]["may_be_selected_independently"] = True
        f = P.validate(d)
        self.assertIn("P4_PASS not false", f)
        self.assertIn("moments became independent ports", f)


if __name__ == "__main__":
    unittest.main()
