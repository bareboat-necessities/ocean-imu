from __future__ import annotations

from pathlib import Path
import sys
import unittest

TOOLS = Path(__file__).resolve().parents[2] / "tools" / "stability"
sys.path.insert(0, str(TOOLS))

import ou3_p4_uniform_closure as C


class UniformClosureTests(unittest.TestCase):
    def test_audit_closes_admissions_but_not_unproved_uniform_map(self):
        d = C.build()
        self.assertEqual(C.validate(d), [])
        self.assertTrue(d["qualified_hard_entry_error_set_as_declared_theorem_assumption"])
        self.assertTrue(d["hard_entry_is_shipping_reachable_fiber"])
        self.assertFalse(d["legacy_independent_S_entry_ball_consumed"])
        self.assertFalse(d["shipping_covariance_used_as_entry_membership_test"])
        self.assertTrue(d["BIAS1_admission"]["conditional_BIAS1_SOURCE_ADMISSION_PASS"])
        self.assertFalse(d["BIAS1_admission"]["assembled_sensor_BIAS0_deployment_qualification_pass"])
        self.assertTrue(d["P3_execution_admission"]["conditional_execution_premises_admitted"])
        self.assertFalse(d["P3_execution_admission"]["physical_execution_admission_proved"])
        self.assertTrue(d["projection_global_nonlinear_sector_consumed"])
        self.assertFalse(d["source_uniform_Kalman_reset_coefficient_family_enclosed"])
        self.assertFalse(d["consecutive_compatible_storage_inequality_closed"])
        self.assertFalse(d["finite_precision_enclosure_closed"])
        self.assertFalse(d["P4_MOTION_PASS"])
        self.assertFalse(d["P4_PASS"])
        self.assertFalse(d["P5_MAY_START"])
        self.assertGreaterEqual(len(d["remaining_blockers"]), 3)

    def test_hard_entry_is_correlated_handoff_fiber_not_legacy_translation_box(self):
        d = C.build()
        entry = d["hard_entry_error_set"]
        self.assertAlmostEqual(entry["attitude_cayley_norm_upper"], 2.0 * __import__("math").tan(__import__("math").pi / 8.0))
        self.assertEqual(entry["velocity_error_norm_upper_mps"], 5.0)
        self.assertEqual(entry["position_error_norm_upper_m"], 0.0)
        self.assertEqual(entry["integral_displacement_error_norm_upper_m_s"], 0.0)
        self.assertGreater(entry["working_prefix_position_norm_upper_m"], 0.0)
        self.assertGreaterEqual(entry["working_prefix_integral_displacement_norm_upper_m_s"], 300.0)
        self.assertEqual(entry["accelerometer_bias_error_norm_upper_mps2"], 0.4)
        self.assertTrue(entry["entry_is_correlated_fiber_not_cartesian_box"])
        self.assertFalse(entry["legacy_independent_S_entry_ball_promotable"])
        self.assertFalse(entry["covariance_ellipsoid_used_for_membership"])


if __name__ == "__main__":
    unittest.main()
