from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

from ou3_interval import Interval, matrix_point  # noqa: E402
import ou3_p4_complete_brmm_joint_sector_master as mod  # noqa: E402


class CompleteBrmmJointSectorMasterTests(unittest.TestCase):
    def test_master_matrix_is_exact_endpoint_energy_identity(self):
        D = matrix_point([[0.75]])
        M = matrix_point([[0.5]])
        J = matrix_point([[1.0]])
        B = matrix_point([[0.1]])
        master = mod.master_quadratic_matrix(D, M, J, B)

        x = 0.4
        w = -0.3
        direct = -0.75 * x * x + 2.0 * (0.5 * x) * (0.1 * w) + (0.1 * w) ** 2
        value = mod.quadratic_value(master, [Interval.point(x), Interval.point(w)])
        self.assertTrue(value.contains(direct))
        self.assertTrue(master[0][0].contains(-0.75))
        self.assertTrue(master[0][1].contains(0.05))
        self.assertTrue(master[1][1].contains(0.01))

    def test_full_augmented_sprocedure_closes_strict_fixture(self):
        master = mod.master_quadratic_matrix(
            matrix_point([[0.75]]),
            matrix_point([[0.5]]),
            matrix_point([[1.0]]),
            matrix_point([[0.1]]),
        )
        # The graph sector is 0.01*x^2-w^2 >= 0, i.e. |w| <= 0.1|x|.
        sector = mod.quadratic_graph_sector(
            matrix_point([[1.0, 0.0]]),
            matrix_point([[0.01]]),
            matrix_point([[0.0, 1.0]]),
            matrix_point([[1.0]]),
        )
        closed, pivots = mod.certify_strict_joint_sector_domination(
            master, [sector], [0.02]
        )
        self.assertTrue(closed)
        self.assertGreater(min(pivots), 0.0)

    def test_sprocedure_rejects_invalid_multipliers_or_dimensions(self):
        master = matrix_point([[-1.0, 0.0], [0.0, -1.0]])
        sector = matrix_point([[1.0, 0.0], [0.0, -1.0]])
        with self.assertRaises(ValueError):
            mod.joint_sector_sprocedure_matrix(master, [sector], [-0.1])
        with self.assertRaises(ValueError):
            mod.joint_sector_sprocedure_matrix(master, [matrix_point([[1.0]])], [0.1])
        with self.assertRaises(ValueError):
            mod.joint_sector_sprocedure_matrix(master, [sector], [])

    def test_contract_retains_source_p3_full_state_and_actual_rs(self):
        channel = mod.CHANNEL.build()
        residual = mod.RESIDUAL.build()
        signed = {"joint_complete_word_signed_information_composition_available": True}
        p3 = {
            "P3_CONDITIONAL_BRMM_PASS": True,
            "mems_bias_preconditions": mod.BIAS.build(),
            "conditional_composition": {
                "A21_finite_bias_correlation_route_consumed": True,
                "A21_detectability_completion_closed": True,
                "A21_paper_UES_hypotheses_closed": True,
                "A21_comparison_observer_is_proof_only_not_alternate_estimator": True,
                "A21_uses_eta9_packet_shortcut": False,
                "A21_detectability_asymptotic_word_energy_gap_lower": 1.0e-18,
                "A21_bias_homogeneous_contraction_gap_lower": 0.001,
            },
            "modes": {
                "H18": {
                    "Omega_minus_delta_P_full_matrix_closed": True,
                    "relative_Riccati_injection_margin_lower": 1.0e-18,
                },
                "A21": {
                    "Omega_minus_delta_P_full_matrix_closed": True,
                    "relative_Riccati_injection_margin_lower": 1.0e-18,
                },
            },
        }
        with (
            patch.object(mod.SIGNED, "validate", return_value=[]),
            patch.object(mod.P3, "validate", return_value=[]),
        ):
            d = mod.build(
                channel_contract=channel,
                residual_contract=residual,
                signed_contract=signed,
                p3_contract=p3,
            )
        self.assertEqual(mod.validate(d), [])
        self.assertEqual(d["mems_bias_preconditions"], p3["mems_bias_preconditions"])
        self.assertEqual(d["canonical_source"], "COMPLETE_BRMM_NORMAL_LIVE_WORD")
        self.assertEqual(d["P3_delta_required"], 1.0e-18)
        self.assertEqual(d["full_state_dimensions"], {"H18": 18, "A21": 21})
        self.assertTrue(d["full_A21_finite_taub_P3_margin_required"])
        self.assertTrue(d["canonical_full_matrix_P3_consumed"])
        self.assertTrue(d["finite_taub_A21_detectability_consumed"])
        self.assertTrue(d["detectability_comparison_observer_only"])
        self.assertFalse(d["finite_taub_detectability_used_as_standalone_P4_promotion"])
        self.assertTrue(d["actual_RS_provenance_preserved"])
        self.assertTrue(d["all_21_state_cross_terms_retained"])
        self.assertTrue(d["same_history_quadratic_graph_sector_assembler_available"])
        self.assertTrue(d["arbitrary_joint_nondiagonal_sector_weights_retained"])
        self.assertFalse(d["per_event_scalar_sector_required"])
        self.assertFalse(d["source_uniform_same_history_joint_sector_closed"])
        self.assertFalse(d["source_uniform_full_augmented_LDLT_closed"])
        self.assertFalse(d["P4_promoted_here"])


if __name__ == "__main__":
    unittest.main()
