from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_private_mahony_state_step as mod  # noqa: E402


class Sea3PrivateMahonyStateStepTest(unittest.TestCase):
    def test_contract_is_source_neutral_and_shipping_bound(self):
        d = mod.build()
        self.assertEqual(mod.validate(d), [])
        self.assertTrue(d["shipping_source_parity_pass"])
        self.assertTrue(d["requires_same_SEA3_gyro_and_specific_force"])
        self.assertTrue(d["requires_same_SEA3_private_observer_state"])
        self.assertTrue(d["actual_fast_inverse_sqrt_used"])
        self.assertFalse(d["source_generator"])
        self.assertFalse(d["independent_vertical_acceleration_source"])
        self.assertFalse(d["independent_quaternion_or_integral_box_promotable"])
        self.assertFalse(d["P3_promoted"])

    def test_declared_sea3_chart_closes_outer_timeout_branch(self):
        d = mod.build()
        self.assertTrue(d["live_entry_integral_state_starts_from_reset_zero"])
        self.assertTrue(d["low_level_TunerReady_can_wait_for_external_bootstrap"])
        self.assertTrue(d["deployed_outer_wrapper_has_timeout_logic"])
        self.assertEqual(d["deployed_proxy_startup_timeout_s"], 150.0)
        self.assertEqual(d["deployed_mag_acquire_deadline_s"], 60.0)
        self.assertEqual(d["deployed_timeout_s"], 150.0)
        self.assertTrue(d["timeout_path_requires_gravity_aligned_branch"])
        self.assertEqual(d["declared_startup_mahony_chart_deg"], 60.0)
        self.assertGreater(d["gravity_aligned_branch_margin_rad"], 0.0)
        self.assertTrue(d["declared_startup_chart_implies_gravity_aligned_branch"])
        self.assertTrue(d["declared_domain_live_entry_upper_bound_closed"])
        self.assertEqual(d["declared_domain_live_entry_upper_bound_s"], 150.0)
        self.assertFalse(d["unconditional_live_entry_upper_bound_closed"])
        self.assertFalse(d["live_entry_private_observer_invariant_closed"])
        self.assertIn("<=150 s declared-domain startup horizon", d["next_obligation"])

    def test_point_step_is_finite_but_cannot_promote(self):
        d = mod.build()
        self.assertTrue(d["smoke"]["finite"])
        self.assertTrue(d["point_smoke_only_not_P3"])
        self.assertFalse(d["complete_SEA3_family_materialized_here"])

    def test_component_box_that_loses_nonzero_norm_is_rejected(self):
        st = mod.State(
            mod.I(1.0), mod.I(0.0), mod.I(0.0), mod.I(0.0),
            mod.I(0.0), mod.I(0.0), mod.I(0.0), mod.I(0.0),
        )
        wide = mod.Interval(-10.0, 10.0)
        with self.assertRaisesRegex(ValueError, "SEA3 nonzero-norm coupling"):
            mod.advance_initialized_live(
                st,
                dt=mod.I(0.005),
                gyro=mod.Vec3(mod.I(0.0), mod.I(0.0), mod.I(0.0)),
                acc_specific_force=mod.Vec3(wide, wide, wide),
                gravity_ms2=mod.I(9.80665),
                two_kp=mod.I(0.2),
                two_ki=mod.I(0.02),
            )


if __name__ == "__main__":
    unittest.main()
