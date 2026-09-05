from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_word_process_memory as mod  # noqa: E402


class Sea3WordProcessMemoryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_complete_sea3_source_is_retained(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(self.d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(self.d["component_of_complete_SEA3_full_word"])
        self.assertFalse(self.d["source_family_replaced"])
        self.assertFalse(self.d["trajectory_replay_used"])
        self.assertFalse(self.d["independent_tau_sigma_RS_box_used_as_source"])

    def test_shipping_commit_cadence_supplies_nontrivial_process_memory(self):
        self.assertEqual(self.d["shipping_staged_commit_parity_failures"], [])
        self.assertGreaterEqual(self.d["guaranteed_constant_active_prediction_steps"], 4)
        self.assertGreater(self.d["guaranteed_constant_active_process_horizon_s"], 0.0)
        self.assertGreater(self.d["integrated_OU_validated_leaf_count"], 0)
        self.assertGreater(self.d["integrated_OU_scaled_lambda_min_lower"], 0.0)

    def test_full_H_and_A_premeasurement_memory_is_strict(self):
        for mode in ("H", "A"):
            q = self.d["modes"][mode]["pre_measurement_process_memory_lambda_min_lower"]
            self.assertTrue(math.isfinite(q))
            self.assertGreater(q, 0.0)

    def test_process_memory_cannot_replace_actual_rs_or_promote(self):
        self.assertFalse(self.d["one_step_full_state_Q_used_as_contraction_ratio"])
        self.assertFalse(self.d["actual_R_S_information_replaced_by_process_strictness"])
        self.assertTrue(self.d["actual_applied_SpectralMSE_R_S_must_be_composed_separately"])
        self.assertFalse(self.d["P3_promoted"])


if __name__ == "__main__":
    unittest.main()
