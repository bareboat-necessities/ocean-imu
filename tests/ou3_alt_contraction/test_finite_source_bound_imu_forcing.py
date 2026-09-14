"""Same-packet IMU forcing regressions; not contraction/storage evidence."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_source_bound_imu_forcing as X
import test_finite_source_bound_prediction_word as BASE


class Tests(unittest.TestCase):
    def test_forcing_is_derived_from_executed_same_packet_event(self):
        s=BASE.BASE.root_state(); witness,segment,raw,dynamic=BASE.operands(s)
        out=X.imu_step(s,witness=witness,segment=segment,raw=raw,packet_id='imu-force-1',
                       **BASE.root_args(),**dynamic)
        prefix=out.word.event.event.live
        f=out.forcing
        self.assertEqual(f.gyro_residual_internal,raw.gyro_residual_internal)
        self.assertEqual(f.guarded_accel_residual_internal,
                         prefix.guarded.effective_accel_residual_internal)
        conditioning=BASE.X._accel_conditioning(F(35))
        expected_thermal=tuple(conditioning.k_a_hat_internal[i]*conditioning.temperature_delta
                               for i in range(3))
        self.assertEqual(f.thermal_model_internal,expected_thermal)
        expected_nu=tuple(f.guarded_accel_residual_internal[i]-expected_thermal[i]
                          for i in range(3))
        self.assertEqual(f.nu_acc_internal,expected_nu)
        self.assertEqual(out.state.source.steps[-1].witness.ordinal,1)

    def test_nonreference_temperature_is_retained_as_supply_not_coefficient_change(self):
        s=BASE.BASE.root_state(); witness,segment,raw,dynamic=BASE.operands(s)
        args=BASE.root_args(); args['temperature_c']=F(37)
        out=X.imu_step(s,witness=witness,segment=segment,raw=raw,packet_id='imu-force-hot',
                       **args,**dynamic)
        f=out.forcing
        self.assertEqual(f.temperature_delta_c,F(2))
        self.assertEqual(f.thermal_model_internal,(F(1,250),)*3)
        self.assertEqual(f.nu_acc_internal,
                         tuple(f.guarded_accel_residual_internal[i]-F(1,250)
                               for i in range(3)))

    def test_forcing_record_rejects_detached_accelerometer_nu(self):
        with self.assertRaisesRegex(ValueError,'detached'):
            X.ImuForcing((0,0,0),(1,2,3),(F(1,10),)*3,(1,2,3),0)

    def test_supply_vector_contains_actual_gyro_nu_and_thermal_coordinates(self):
        f=X.ImuForcing((1,2,3),(4,5,6),(F(1,2),)*3,
                       (F(7,2),F(9,2),F(11,2)),F(1))
        self.assertEqual(f.supply_vector,
                         (F(1),F(2),F(3),F(7,2),F(9,2),F(11,2),
                          F(1,2),F(1,2),F(1,2)))

    def test_readiness_retains_supply_without_inventing_input_bound(self):
        r=X.readiness()
        self.assertTrue(r['same_packet_gyro_residual_retained_as_ISS_supply'])
        self.assertTrue(r['post_guard_accel_residual_derived_not_free'])
        self.assertTrue(r['temperature_model_term_retained_separately'])
        self.assertTrue(r['accelerometer_nu_equals_guarded_residual_minus_thermal_term'])
        self.assertTrue(r['forcing_supply_derived_from_executed_source_owned_IMU_event'])
        self.assertTrue(r['lower_model_roots_source_bound'])
        self.assertFalse(r['sensor_or_temperature_amplitude_bound_invented'])
        self.assertFalse(r['bounded_input_history_qualified'])
        self.assertFalse(r['deployment_roundoff_supply_attached'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
