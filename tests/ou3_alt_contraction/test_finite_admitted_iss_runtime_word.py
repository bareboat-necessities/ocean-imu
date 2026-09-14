"""Master admitted runtime + bounded ISS-history regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_imu_disturbance as DIST
from tools.stability.ou3_alt_contraction import finite_admitted_iss_runtime_word as X
import test_finite_admitted_runtime_word as BASE
import test_finite_source_bound_live_word as LOWER


class Tests(unittest.TestCase):
    def root(self,bound=1000):
        base=BASE.Tests().root()
        d=DIST.BoundedHistory(base.admitted.live_word.sensor_root,'temperature-history',F(bound))
        return X.bind_startup(base,d)

    def test_master_IMU_charges_same_kth_exact_forcing_to_W(self):
        s=self.root(); witness,r,b,raw,runtime=BASE.Tests().next_operands(s.runtime)
        out=X.imu_step(s,restricted=r,bias_restricted=b,witness=witness,
                       raw=raw,packet_id='iss-runtime-1',**runtime)
        self.assertEqual(out.restriction.ordinal,1)
        self.assertEqual(out.restriction.history,s.disturbance)
        self.assertEqual(out.restriction.forcing,out.forcing)
        self.assertEqual(out.state.disturbance,s.disturbance)
        self.assertEqual(out.state.runtime.origin,s.runtime.origin)

    def test_MAG_and_HOLD_do_not_replace_bounded_history(self):
        s=self.root(); hist=s.disturbance
        s,m=X.mag_step(s,**LOWER.mag_kwargs(s.runtime.admitted.live_word))
        self.assertIs(s.disturbance,hist)
        before=s.runtime.admitted.live_word.source.next_ordinal
        s=X.set_hold(s,hold=False)
        self.assertIs(s.disturbance,hist)
        self.assertEqual(s.runtime.admitted.live_word.source.next_ordinal,before)

    def test_master_state_rejects_detached_sensor_history(self):
        base=BASE.Tests().root(); sr=base.admitted.live_word.sensor_root
        other=type(sr)(sr.source_root,'other-gyro','other-accel')
        d=DIST.BoundedHistory(other,'temperature-history',1000)
        with self.assertRaisesRegex(ValueError,'detached from carried sensor residual histories'):
            X.bind_startup(base,d)

    def test_readiness_closes_only_bounded_input_structural_blocker(self):
        r=X.readiness()
        self.assertTrue(r['bounded_input_history_qualified'])
        self.assertTrue(r['bounded_IMU_ISS_history_is_part_of_master_runtime_state'])
        self.assertTrue(r['same_kth_executed_IMU_supply_charged_to_symbolic_W'])
        self.assertTrue(r['Racc_covariance_not_used_as_pathwise_noise_bound'])
        self.assertFalse(r['finite_horizon_probability_used_to_prune_disturbances'])
        self.assertFalse(r['deployment_exp_expm1_trig_Eigen_LDLT_closed'])
        self.assertFalse(r['shipping_signed_mag_counter_lifetime_closed'])
        self.assertFalse(r['complete_600_step_shipping_word_composed_from_restrictions'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
