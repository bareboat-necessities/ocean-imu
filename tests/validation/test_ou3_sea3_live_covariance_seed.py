from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_live_covariance_seed as mod  # noqa: E402


class Sea3LiveCovarianceSeedTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_seed_is_shipping_generated_and_valid(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertTrue(self.d["live_entry_seed_is_source_generated_not_arbitrary_PSD"])
        self.assertFalse(self.d["bootstrap_mekf_covariance_propagated_before_live"])
        self.assertEqual(self.d["source_parity_failures"], [])

    def test_translation_seed_is_not_an_arbitrary_covariance_box(self):
        seed = self.d["translation_seed"]
        self.assertEqual(seed["P_v"], 1.0)
        self.assertEqual(seed["P_p"], 400.0)
        self.assertEqual(seed["P_S"], 2500.0)
        self.assertTrue(seed["all_translation_cross_covariances_zero_before_first_prediction"])

    def test_aw_seed_uses_same_committed_tuner_state(self):
        aw = self.d["aw_live_seed"]
        self.assertTrue(aw["reset_to_committed_stationary_covariance_before_first_prediction"])
        self.assertTrue(aw["isotropic_default"])
        self.assertEqual(aw["S_factor"], 1.0)
        lo, hi = aw["committed_vertical_std_interval_mps2"]
        self.assertGreater(lo, 0.0)
        self.assertGreaterEqual(hi, lo)

    def test_held_bias_source_semantics_are_retained(self):
        held = self.d["held_ba"]
        self.assertTrue(held["excluded_from_H18"])
        self.assertTrue(held["identity_homogeneous_dynamics"])
        self.assertTrue(held["no_process_injection_while_held"])
        self.assertTrue(held["cross_covariances_zero"])
        self.assertTrue(held["measurement_rows_frozen"])
        self.assertGreater(held["seed_variance"], 0.0)

    def test_h_to_a_is_hybrid_not_hidden_in_same_mode_word(self):
        release = self.d["H_to_A_release"]
        self.assertTrue(release["hybrid_transition_not_inside_same_mode_word"])
        self.assertGreater(release["bias_diagonal_floor_variance"], 0.0)
        self.assertFalse(self.d["P3_promoted"])


if __name__ == "__main__":
    unittest.main()
