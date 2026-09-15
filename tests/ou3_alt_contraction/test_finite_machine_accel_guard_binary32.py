"""Binary32 shipping AccelVibrationGuard source regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_machine_accel_guard_binary32 as X
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B


class Tests(unittest.TestCase):
    def test_exact_zero_detector_rms_is_a_valid_machine_root(self):
        self.assertEqual(X._sqrt(0,0,'RMS'),0)
        with self.assertRaisesRegex(ValueError,'zero radicand requires zero result'):
            X._sqrt(0,B.rn32(F(1,1000)),'RMS')

    def test_update_API_adapter_requires_exact_binary32_RNE_image(self):
        exact=(F(1,10),F(-2,10),F(3,10)); stored=tuple(B.rn32(x) for x in exact)
        self.assertEqual(X.api_vec(exact,stored,'acc'),stored)
        bad=(stored[0],stored[1],B.rn32(F(31,100)))
        with self.assertRaisesRegex(ValueError,'API rounding'):
            X.api_vec(exact,bad,'acc')

    def test_first_enabled_sample_seeds_all_shipping_guard_states_and_is_transparent(self):
        s=X.State(); c=X.Config(); g=(B.rn32(F(1,100)),B.rn32(F(2,100)),B.rn32(F(3,100)))
        a=(B.rn32(F(1,5)),B.rn32(F(-1,10)),B.rn32(F(49,5)))
        out=X.step(s,c,raw_gyro=g,raw_acc=a,dt=B.rn32(F(1,200)))
        self.assertTrue(out.seeded); self.assertEqual(out.conditioned_acc,a)
        self.assertEqual(out.state.stages,(a,a,a,a)); self.assertEqual(out.state.detect_stages[0],a)
        self.assertEqual(out.state.detect_stages[1],(X.ZERO,X.ZERO,X.ZERO))
        self.assertEqual(out.removed_rms,X.ZERO); self.assertEqual(out.excess_rms,X.ZERO)

    def test_disabled_guard_is_bit_transparent_without_initializing_filter_memory(self):
        s=X.State(); c=X.Config(cutoff_hz=B.rn32(0)); a=(B.rn32(1),B.rn32(2),B.rn32(3)); g=(B.rn32(0),)*3
        out=X.step(s,c,raw_gyro=g,raw_acc=a,dt=B.rn32(F(1,200)))
        self.assertFalse(out.enabled); self.assertFalse(out.state.initialized); self.assertEqual(out.conditioned_acc,a)

    def test_readiness_keeps_runtime_and_native_correspondence_open(self):
        r=X.readiness()
        for k in ('shipping_guard_source_shape_matches','binary32_update_API_boundary_materialized',
                  'first_sample_guard_seed_and_exact_transparency_materialized','conditioning_and_detector_cascades_materialized',
                  'removed_RMS_EMA_and_sqrt_materialized','engagement_slew_rails_and_mixed_output_materialized'):
            self.assertTrue(r[k])
        for k in ('target_libm_and_Eigen_expression_correspondence_closed','mutable_guard_runtime_config_ancestry_closed',
                  'source_uniform_guard_supply_bound_closed','storage_search_allowed','ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[k])

if __name__=='__main__': unittest.main()
