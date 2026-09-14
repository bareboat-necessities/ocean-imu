"""Same-event full machine prediction-supply carrier regressions."""
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_machine_prediction_supply_interleaved_prefix as X
import test_finite_admitted_machine_scheduler_interleaved_prefix as BASE
import test_finite_admitted_machine_prediction_interleaved_prefix as PBASE
import test_finite_source_bound_live_word as LBASE


class Tests(unittest.TestCase):
    def test_begin_roots_supply_counter_at_actual_admitted_ordinal(self):
        s=X.begin(BASE.state())
        self.assertEqual(s.supply_steps,0)
        self.assertEqual(X._imu_steps(s.base),s.entry_imu_steps)

    def test_counter_cannot_detach_from_lower_word(self):
        base=BASE.state(); s=X.begin(base)
        with self.assertRaisesRegex(ValueError,'count detached'):
            X.State(base,s.entry_imu_steps,1)

    def test_exact_root_witness_extractor_fails_closed(self):
        with self.assertRaisesRegex(TypeError,'witnesses missing'):
            X._exact_root_kwargs({'ou_alpha':1})

    def test_MAG_and_HOLD_preserve_prediction_supply_count(self):
        s=X.begin(BASE.state()); n=s.supply_steps
        word=PBASE._word(s.base.base.base)
        s,_=X.mag_step(s,**LBASE.mag_kwargs(word))
        self.assertEqual(s.supply_steps,n)
        s,_=X.set_hold(s,hold=False)
        self.assertEqual(s.supply_steps,n)

    def test_complete_cannot_skip_full_prediction_supply(self):
        with self.assertRaises((ValueError,TypeError)):
            X.complete(X.begin(BASE.state()))

    def test_readiness_stops_before_postprediction_measurement_and_storage(self):
        r=X.readiness()
        self.assertTrue(r['machine_root_effect_injected_into_joint24_prediction_relation'])
        self.assertTrue(r['rebuilt_exact_prediction_required_equal_executed_shipping_prediction'])
        self.assertTrue(r['complete_word_requires_full_prediction_supply_on_all_600_IMU_edges'])
        for k in ('post_prediction_floor_and_S_service_propagate_machine_prediction_supply',
                  'machine_RS_measurement_effect_attached','accelerometer_measurement_propagates_machine_prediction_supply',
                  'source_uniform_machine_prediction_supply_bound_closed','all_target_libm_and_Eigen_correspondence_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed',
                  'ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
