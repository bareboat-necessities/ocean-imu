from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

from ou3_interval import matrix_point  # noqa: E402
import ou3_p4_complete_brmm_finite_map_mean_value as mod  # noqa: E402


class CompleteBrmmFiniteMapMeanValueTest(unittest.TestCase):
    def test_contract_keeps_complete_brmm_and_paper_finite_state_theorem(self):
        d = mod.build()
        self.assertEqual(mod.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_BRMM_NORMAL_LIVE_WORD")
        self.assertTrue(d["paper_finite_state_quadratic_theorem_retained"])
        self.assertTrue(d["finite_map_not_differential_metric_replacement"])
        self.assertTrue(d["same_source_correlated_generalized_Jacobian_required"])
        self.assertTrue(d["actual_applied_RS_required_inside_same_word"])
        self.assertTrue(d["Clarke_generalized_Jacobian_handles_A21_projection"])
        self.assertFalse(d["A21_projection_assumed_inactive"])
        self.assertFalse(d["independent_P_H_R_K_boxes_authorized"])
        self.assertFalse(d["independent_tuner_RS_schedule_authorized"])
        self.assertFalse(d["packetwise_remainder_sum_authorized"])
        self.assertFalse(d["finite_harmonic_or_replay_source_authorized"])
        self.assertFalse(d["P4_promoted_here"])

    def test_strict_endpoint_ldlt_accepts_a_contracting_finite_map(self):
        M0 = matrix_point([[1.0, 0.0], [0.0, 1.0]])
        MN = matrix_point([[1.0, 0.0], [0.0, 1.0]])
        J = matrix_point([[0.5, 0.0], [0.0, 0.7]])
        ok, pivots = mod.certify_strict_endpoint_contraction(M0, MN, J, 0.6)
        self.assertTrue(ok)
        self.assertEqual(len(pivots), 2)
        self.assertGreater(min(pivots), 0.0)

    def test_strict_endpoint_ldlt_rejects_noncontracting_map(self):
        M0 = matrix_point([[1.0, 0.0], [0.0, 1.0]])
        MN = matrix_point([[1.0, 0.0], [0.0, 1.0]])
        J = matrix_point([[1.0, 0.0], [0.0, 1.0]])
        ok, _ = mod.certify_strict_endpoint_contraction(M0, MN, J, 0.9)
        self.assertFalse(ok)

    def test_rectangular_H18_to_A21_style_map_is_supported(self):
        M0 = matrix_point([[1.0, 0.0], [0.0, 1.0]])
        M1 = matrix_point([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ])
        J = matrix_point([
            [0.5, 0.0],
            [0.0, 0.5],
            [0.0, 0.0],
        ])
        ok, pivots = mod.certify_strict_endpoint_contraction(M0, M1, J, 0.5)
        self.assertTrue(ok)
        self.assertGreater(min(pivots), 0.0)

    def test_prefix_gain_uses_same_full_matrix_test(self):
        M0 = matrix_point([[1.0, 0.0], [0.0, 1.0]])
        Mell = matrix_point([[1.0, 0.0], [0.0, 1.0]])
        Jell = matrix_point([[1.1, 0.0], [0.0, 1.05]])
        ok, pivots = mod.certify_prefix_gain(M0, Mell, Jell, 1.5)
        self.assertTrue(ok)
        self.assertGreater(min(pivots), 0.0)

    def test_bad_dimensions_and_invalid_gains_fail_closed(self):
        M0 = matrix_point([[1.0]])
        M1 = matrix_point([[1.0, 0.0], [0.0, 1.0]])
        Jbad = matrix_point([[1.0]])
        with self.assertRaises(ValueError):
            mod.finite_quadratic_margin(M0, M1, Jbad, 1.0)
        with self.assertRaises(ValueError):
            mod.certify_strict_endpoint_contraction(M0, M0, M0, 1.0)
        with self.assertRaises(ValueError):
            mod.certify_prefix_gain(M0, M0, M0, 0.5)


if __name__ == "__main__":
    unittest.main()
