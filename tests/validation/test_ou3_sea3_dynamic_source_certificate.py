from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_sea3_dynamic_source_certificate as mod  # noqa: E402


class Sea3DynamicSourceCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_dynamic_source_is_canonical_not_800_state_graph(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(self.d["P2_DYNAMIC_SOURCE_CERTIFICATE"], "PASS")
        self.assertFalse(self.d["old_P2_800_state_graph_consumed"])
        self.assertFalse(self.d["old_P2_history_word_enumeration_consumed"])

    def test_shipping_smoothing_gives_strict_motion_bounds(self):
        r = self.d["validated_rate_and_jump_bounds"]
        self.assertGreater(r["tau_sigma_alpha_per_sample_upper"], 0.0)
        self.assertLess(r["tau_sigma_alpha_per_sample_upper"], 1.0)
        self.assertGreater(r["R_S_alpha_per_sample_upper"], 0.0)
        self.assertLess(r["R_S_alpha_per_sample_upper"], 1.0)
        self.assertGreater(r["active_commit_gap_samples_upper"], 0)
        self.assertTrue(r["proof_relies_on_implemented_smoothing_not_unproved_sea_parameter_derivatives"])

    def test_normal_live_has_no_rejected_accelerometer_branch(self):
        live = self.d["normal_live_contract"]
        self.assertTrue(live["accelerometer_update_required_each_valid_sample"])
        self.assertFalse(live["accelerometer_rejection_in_scope"])


if __name__ == "__main__":
    unittest.main()
