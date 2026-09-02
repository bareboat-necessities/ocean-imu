import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p4_sample_clock_source_refinement as R


class TestOu3P4SampleClockSourceRefinement(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = R.build()

    def test_certificate_validates(self):
        self.assertEqual(R.validate(self.d), [])
        self.assertEqual(self.d["P2_SAMPLE_CLOCK_REFINEMENT_CERTIFICATE"], "PASS")

    def test_shipping_clock_nominal_and_infinite_time_envelope(self):
        c = self.d["clock"]
        self.assertEqual(c["nominal_stage_spacing_valid_samples"], 21)
        self.assertEqual(c["pending_apply_delay_valid_samples"], 1)
        self.assertAlmostEqual(c["nominal_commit_to_commit_elapsed_s"], 0.10499999765306711, places=15)
        self.assertEqual(c["finite_stage_spacing_valid_samples_lower"], 13)
        self.assertEqual(c["finite_stage_spacing_valid_samples_upper"], 26)
        self.assertTrue(c["boundary_enumeration_exceeds_max_stage_spacing"])
        self.assertTrue(c["floating_clock_stagnation_verified"])

    def test_refinement_removes_all_to_all_jump(self):
        self.assertEqual(self.d["partition"]["states"], 800)
        self.assertLess(self.d["transition_edges"], self.d["base_transition_edges"])
        self.assertFalse(self.d["source_graph_all_to_all"])
        self.assertTrue(self.d["arbitrary_late_commit_jump_removed"])

    def test_scope_remains_source_only_and_fail_closed(self):
        self.assertTrue(self.d["source_only"])
        self.assertFalse(self.d["trajectory_replay_used"])
        self.assertFalse(self.d["filter_changed"])
        self.assertTrue(self.d["EMA_composed_sample_by_sample"])
        self.assertTrue(self.d["sample_varying_target_and_horizon_boxes_admitted"])
        self.assertTrue(self.d["RS_target_full_deployed_clamp_overapprox"])
        self.assertFalse(self.d["RS_target_powf_tightening_used"])
        self.assertFalse(self.d["P4_USABLE_CERTIFICATE_PROMOTED"])


if __name__ == "__main__":
    unittest.main()
