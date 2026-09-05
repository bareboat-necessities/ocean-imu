from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_h18_process_forgetting as mod  # noqa: E402


class Sea3H18ProcessForgettingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_complete_source_and_no_replacement_route(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(self.d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(self.d["theorem_prefix_starts_from_shipping_Live_seed"])
        self.assertTrue(self.d["two_piece_tau_values_are_same_SEA3_path_coordinates"])
        self.assertFalse(self.d["independent_tau_sigma_RS_source_created"])
        self.assertFalse(self.d["trajectory_replay_used"])
        self.assertFalse(self.d["source_history_graph_used"])

    def test_full_H18_interval_LDLT_closes_at_useful_gate(self):
        self.assertEqual(self.d["useful_gate"], 1.0e-18)
        self.assertTrue(self.d["full_H18_Omega_minus_delta_P_LDLT_closed"])
        self.assertTrue(self.d["H18_PROCESS_FORGETTING_PASS"])
        pivot = self.d["parameter_cover"]["worst_full_H18_LDLT_pivot_lower"]
        self.assertTrue(math.isfinite(pivot))
        self.assertGreater(pivot, 0.0)
        self.assertGreater(self.d["parameter_cover"]["accepted_cells"], 0)

    def test_measurement_conditioning_retains_applied_spectral_mse_rs(self):
        m = self.d["measurement_conditioning"]
        self.assertTrue(self.d["actual_applied_SpectralMSE_R_S_consumed"])
        self.assertGreater(m["actual_applied_R_S_base_lower"], 0.0)
        self.assertEqual(m["actual_applied_R_S_axis_factors"], [0.72, 0.72, 1.0])
        self.assertGreater(m["posterior_attenuation_lower"], 0.0)
        self.assertLessEqual(m["posterior_attenuation_lower"], 1.0)
        self.assertTrue(m["all_accelerometer_packets_admitted"])
        self.assertTrue(m["S_packet_admitted_at_every_IMU_slot"])

    def test_no_retired_scalar_promotion_shortcuts(self):
        self.assertFalse(self.d["blockwise_minimum_ratio_used"])
        self.assertFalse(self.d["scalar_information_beta_used_for_contraction"])
        self.assertFalse(self.d["determinant_trace_eigenvalue_conversion_used"])
        self.assertFalse(self.d["D_W_L_W_product_used"])
        self.assertFalse(self.d["one_step_process_Q_used_as_contraction_ratio"])
        self.assertFalse(self.d["P3_promoted"])

    def test_reset_is_only_quotiented_as_exact_congruence(self):
        self.assertTrue(self.d["left_error_resets_quotiented_as_nonsingular_congruences"])


if __name__ == "__main__":
    unittest.main()
