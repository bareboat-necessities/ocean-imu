from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_startup_capture_gate as mod  # noqa: E402


class StartupCaptureGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_gate_is_structurally_valid_and_fail_closed(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertTrue(self.d["shipping_handoff_fiber_validated"])
        self.assertTrue(self.d["physical_first_sample_tilt_lemma_validated"])
        self.assertTrue(self.d["mahony_one_sample_binary32_map_validated"])
        self.assertFalse(self.d["P4_PASS"])
        self.assertFalse(self.d["P5_PASS"])

    def test_declared_chart_is_not_capture(self):
        self.assertTrue(self.d["existing_mahony_timeout_argument_has_circular_chart_dependency"])
        self.assertFalse(self.d["legacy_declared_mahony_chart_may_establish_capture"])
        self.assertFalse(self.d["mahony_arithmetic_boundedness_is_capture"])

    def test_quality_gate_leaky_memory_is_retained(self):
        q = self.d["quality_gate_memory"]
        self.assertEqual(q["good_counter_cap_s"], 10.0)
        self.assertEqual(q["required_hold_s"], 2.0)
        self.assertEqual(q["bad_sample_credit_decay_rate"], 2.0)
        self.assertEqual(q["max_bad_tail_after_last_good_sample_s"], 4.0)
        self.assertFalse(q["current_align_sin_threshold_implied_at_quality_handoff"])

    def test_timeout_is_live_handoff_not_p4_capture(self):
        t = self.d["shipping_timeout_semantics"]
        self.assertEqual(t["live_handoff_upper_bound_s_under_declared_aligned_branch"], 150.0)
        self.assertFalse(t["requires_tuner_ready"])
        self.assertFalse(t["requires_magnetic_north_when_with_mag"])
        self.assertTrue(t["requires_aligned_gravity_branch"])
        self.assertFalse(t["goLive_equals_P4_capture"])
        self.assertTrue(t["postlive_first_north_path_validated"])

    def test_raw_magnetic_obstruction_and_capture_admission_are_separate(self):
        m = self.d["magnetic_observability"]
        self.assertTrue(m["structural_yaw_obstruction_proved"])
        self.assertFalse(m["raw_COMPLETE_BRMM_source_forces_finite_north_acquisition"])
        self.assertFalse(m["normal_live_PE_alone_forces_first_north"])
        self.assertTrue(m["explicit_capture_admission_added"])
        self.assertFalse(m["capture_admission_changes_BRMM_motion_caps"])
        self.assertFalse(m["capture_admission_changes_BIAS_family"])
        self.assertTrue(m["finite_north_event_under_explicit_capture_admission_closed"])
        self.assertTrue(m["ungauged_gravity_quotient_capture_is_still_meaningful"])
        self.assertTrue(m["full_attitude_capture_requires_north_event"])
        self.assertEqual(m["shipping_min_accepted_mag_samples"], 128)
        self.assertEqual(m["shipping_min_accepted_mag_window_s"], 15.0)
        self.assertEqual(m["counterexample_minimax_full_attitude_error_deg_without_progress_admission"], 45.0)
        self.assertTrue(self.d["FINITE_NORTH_EVENT_UNDER_EXPLICIT_CAPTURE_ADMISSION_CLOSED"])
        self.assertFalse(self.d["UNCONDITIONAL_FINITE_NORTH_FROM_RAW_COMPLETE_BRMM_CLOSED"])

    def test_translation_entry_is_correlated_handoff_fiber(self):
        e = self.d["entry_translation_fiber"]
        self.assertEqual(e["e_p_at_handoff"], 0.0)
        self.assertEqual(e["e_S_at_handoff"], 0.0)
        self.assertEqual(e["e_v_at_handoff"], "v_true(t_h)")
        self.assertFalse(e["independent_300_m_s_S_entry_allowed"])


if __name__ == "__main__":
    unittest.main()
