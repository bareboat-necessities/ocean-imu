"""Strong admitted-source + bounded-ISS IMU edge regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_imu_disturbance as DIST
from tools.stability.ou3_alt_contraction import finite_admitted_disturbance_imu_word as X
import test_finite_admitted_source_imu_word as BASE
import test_finite_source_bound_prediction_word as PBASE


def root(bound=1000):
    live=BASE.root()
    hist=DIST.BoundedHistory(live.live_word.sensor_root,'temperature-history',F(bound))
    return X.State(live,hist)


class Tests(unittest.TestCase):
    def test_same_event_forcing_is_charged_to_persistent_history(self):
        s=root(); witness,segment,raw,r,b,dynamic=BASE.operands(s.live)
        out=X.imu_step(s,restricted=r,bias_restricted=b,witness=witness,raw=raw,
                       packet_id='imu-1',**PBASE.root_args(),**dynamic)
        self.assertEqual(out.state.disturbance,s.disturbance)
        self.assertEqual(out.restriction.ordinal,1)
        self.assertEqual(out.restriction.forcing,out.forcing)
        self.assertEqual(out.state.live.admitted_history,s.live.admitted_history)
        self.assertEqual(out.state.live.bias_history,s.live.bias_history)

    def test_detached_sensor_history_is_rejected_at_product_construction(self):
        live=BASE.root()
        other=replace(live.live_word.sensor_root,accel_residual_history_id='other-accel-history')
        hist=DIST.BoundedHistory(other,'temperature-history',1000)
        with self.assertRaisesRegex(ValueError,'detached from carried sensor residual histories'):
            X.State(live,hist)

    def test_small_symbolic_bound_rejects_actual_executed_supply(self):
        s=root(0); witness,segment,raw,r,b,dynamic=BASE.operands(s.live)
        # Perturb the same-packet accelerometer residual, not gyro.  This keeps
        # the attitude prediction on the original branch while still making the
        # executed ISS forcing nonzero, so this regression isolates the carried
        # disturbance-bound guard rather than requiring unrelated trig witnesses.
        residual=(F(1),0,0)
        raw2=replace(raw,accel_residual_internal=residual,
                     raw_accel_body=(raw.raw_accel_body[0]+1,raw.raw_accel_body[1],raw.raw_accel_body[2]))
        with self.assertRaisesRegex(ValueError,'exceeds carried theorem bound'):
            X.imu_step(s,restricted=r,bias_restricted=b,witness=witness,raw=raw2,
                       packet_id='imu-1',**PBASE.root_args(),**dynamic)
        self.assertEqual(len(s.live.live_word.source.steps),0)

    def test_readiness_does_not_promote_storage(self):
        q=X.readiness()
        self.assertTrue(q['bounded_IMU_ISS_history_carried_as_theorem_state'])
        self.assertTrue(q['same_kth_executed_forcing_charged_to_history'])
        self.assertTrue(q['Racc_covariance_not_used_as_pathwise_bound'])
        self.assertFalse(q['complete_600_step_shipping_word_composed'])
        self.assertFalse(q['storage_search_allowed'])
        self.assertFalse(q['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
