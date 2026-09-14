"""Admitted whole-machine prediction-root interleaver regressions."""
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_machine_prediction_interleaved_prefix as X
import test_finite_admitted_machine_tunestate_interleaved_prefix as BASE
import test_finite_source_bound_live_word as LBASE


class Tests(unittest.TestCase):
    def test_begin_counts_from_actual_admitted_IMU_ordinal(self):
        base=BASE.state(); s=X.begin(base)
        self.assertEqual(s.prediction_steps,0)
        self.assertEqual(X._imu_steps(s.base),s.entry_imu_steps)

    def test_counter_cannot_detach_from_lower_admitted_word(self):
        base=BASE.state(); s=X.begin(base)
        with self.assertRaisesRegex(ValueError,'count detached'):
            X.State(base,s.entry_imu_steps,1)

    def test_MAG_and_HOLD_preserve_prediction_obligation_counter(self):
        s=X.begin(BASE.state()); n=s.prediction_steps
        word=X._word(s.base)
        s,_=X.mag_step(s,**LBASE.mag_kwargs(word))
        self.assertEqual(s.prediction_steps,n)
        s,_=X.set_hold(s,hold=False)
        self.assertEqual(s.prediction_steps,n)

    def test_complete_word_cannot_skip_machine_prediction_roots(self):
        s=X.begin(BASE.state())
        with self.assertRaises((ValueError,TypeError)):
            X.complete(s)

    def test_machine_root_witness_dicts_cannot_override_event_ancestry(self):
        s=X.begin(BASE.state())
        # This must fail before any lower shipping event is attempted.
        with self.assertRaisesRegex(TypeError,'cannot override source/event ancestry'):
            X.imu_step(s,separate_machine_root_kwargs={'raw':'detached'},
                       fma_machine_root_kwargs={},restricted=object(),witness=object(),raw=object())

    def test_readiness_promotes_only_prediction_root_attachment(self):
        r=X.readiness()
        self.assertTrue(r['machine_prediction_root_relation_attached_to_admitted_event'])
        self.assertTrue(r['post_boundary_exact_ActiveParameters_taken_from_executed_prediction'])
        self.assertTrue(r['complete_word_requires_machine_prediction_roots_on_all_600_IMU_edges'])
        for k in ('machine_pseudo_period_scheduler_effect_attached','machine_RS_measurement_effect_attached',
                  'OU_Qaxis_target_libm_and_Eigen_correspondence_closed',
                  'source_uniform_machine_coefficient_supply_bound_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed',
                  'ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
