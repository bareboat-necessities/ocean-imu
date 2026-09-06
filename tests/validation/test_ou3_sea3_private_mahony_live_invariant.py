from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_sea3_private_mahony_live_invariant as mod  # noqa: E402


class Sea3PrivateMahonyLiveInvariantTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_continuous_all_live_invariant_closes(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertTrue(self.d["continuous_all_live_PI_invariant_closed"])
        self.assertTrue(self.d["initial_set_inside_invariant"])
        self.assertTrue(self.d["invariant_strictly_inside_60deg_chart"])
        self.assertGreater(
            self.d["boundary_validation"]["strict_inward_margin_lower"], 0.0
        )

    def test_same_sea3_forcing_not_an_alternate_source(self):
        self.assertTrue(self.d["same_SEA3_specific_force_direction_required"])
        self.assertTrue(self.d["same_SEA3_gyro_bias_forcing_required"])
        self.assertFalse(self.d["source_generator"])
        self.assertFalse(self.d["trajectory_replay_used"])
        self.assertFalse(self.d["arbitrary_bounded_input_source_used"])

    def test_discrete_shipping_composition_remains_fail_closed(self):
        self.assertFalse(self.d["shipping_binary32_discrete_invariant_closed"])
        self.assertFalse(self.d["complete_SEA3_family_materialized_here"])
        self.assertFalse(self.d["P3_promoted"])


if __name__ == "__main__":
    unittest.main()
