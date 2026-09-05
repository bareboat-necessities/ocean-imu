from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_sea3_rs_tau_lag_envelope as mod  # noqa: E402


class Sea3RsTauLagEnvelopeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_validated_low_dimensional_invariant(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertTrue(self.d["source_generated_not_trajectory_fit"])
        self.assertFalse(self.d["old_P2_800_state_graph_consumed"])
        self.assertFalse(self.d["source_history_graph_consumed"])
        self.assertFalse(self.d["predecessor_path_enumeration_consumed"])
        self.assertTrue(self.d["candidate_snapshot_commit_preserves_invariant"])

    def test_fractional_mse_power_removed_from_pass_decision(self):
        self.assertTrue(
            self.d["SpectralMSE_fractional_power_removed_by_14th_power_identity"]
        )
        self.assertFalse(self.d["ordinary_libm_fractional_power_used_in_pass_decision"])
        self.assertTrue(self.d["validated_exponential_arithmetic"])
        self.assertTrue(self.d["target_lower_curve"]["pass"])

    def test_applied_curve_improves_high_tau_floor(self):
        floor = self.d["R_S_hard_floor"]
        applied = self.d["applied_invariant_lower_curve"]
        self.assertTrue(applied["pass"])
        self.assertTrue(applied["initial_state_inside"])
        self.assertGreater(applied["R_at_tau_max_lower"], floor)
        self.assertGreater(self.d["target_lower_curve"]["R_at_tau_max_lower"], applied["R_at_tau_max_lower"])

    def test_does_not_promote_p3_by_itself(self):
        self.assertFalse(self.d["P3_promoted"])
        self.assertFalse(self.d["P4_promoted"])


if __name__ == "__main__":
    unittest.main()
