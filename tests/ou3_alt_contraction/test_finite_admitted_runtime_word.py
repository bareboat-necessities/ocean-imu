"""Admitted-history + strong-runtime composition regressions."""
from dataclasses import replace
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as A
from tools.stability.ou3_alt_contraction import finite_admitted_bias_history as B
from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as ADMITTED
from tools.stability.ou3_alt_contraction import finite_admitted_runtime_word as X
import test_finite_source_bound_prediction_word as RUNTIME
import test_finite_source_bound_live_word as LOWER


class Tests(unittest.TestCase):
    def root(self):
        lower=LOWER.root_state()
        history=A.AdmittedHistory(lower.source.root.history_id)
        bias=B.AdmittedBiasHistory(lower.bias_history_id,lower.source.root.bias_family,LOWER.PHI)
        return ADMITTED.State(lower,history,bias)

    def next_operands(self,s):
        witness,segment,raw,dynamic=RUNTIME.operands(s.live_word)
        restricted=A.RestrictedSegment(s.admitted_history,witness.ordinal,segment)
        bias_restricted=B.RestrictedBiasStep(s.bias_history,witness.ordinal,segment)
        runtime=RUNTIME.root_args(); runtime.update(dynamic)
        return witness,restricted,bias_restricted,raw,runtime

    def test_IMU_requires_admitted_histories_and_executes_strong_runtime(self):
        s=self.root(); witness,r,b,raw,runtime=self.next_operands(s)
        out=X.imu_step(s,restricted=r,bias_restricted=b,witness=witness,
                       raw=raw,packet_id='admitted-runtime-1',**runtime)
        self.assertEqual(out.state.admitted_history,s.admitted_history)
        self.assertEqual(out.state.bias_history,s.bias_history)
        self.assertEqual(out.state.live_word.source.steps[-1].witness.ordinal,1)
        self.assertEqual(out.state.live_word.live.live.live.mekf.reference,r.segment.after)
        self.assertEqual(out.forcing.gyro_residual_internal,raw.gyro_residual_internal)
        self.assertIsNotNone(out.runtime_word)

    def test_detached_BRMM_history_cannot_reach_strong_runtime(self):
        s=self.root(); witness,r,b,raw,runtime=self.next_operands(s)
        bad=A.RestrictedSegment(A.AdmittedHistory('other-history'),r.ordinal,r.segment)
        with self.assertRaisesRegex(ValueError,'detached from carried admitted history'):
            X.imu_step(s,restricted=bad,bias_restricted=b,witness=witness,
                       raw=raw,packet_id='bad',**runtime)
        self.assertEqual(len(s.live_word.source.steps),0)

    def test_detached_BIAS_history_cannot_reach_strong_runtime(self):
        s=self.root(); witness,r,b,raw,runtime=self.next_operands(s)
        other=B.AdmittedBiasHistory('other-bias',s.bias_history.family,s.bias_history.phi_true)
        bad=B.RestrictedBiasStep(other,b.ordinal,b.segment)
        with self.assertRaisesRegex(ValueError,'detached from carried admitted BIAS history'):
            X.imu_step(s,restricted=r,bias_restricted=bad,witness=witness,
                       raw=raw,packet_id='bad',**runtime)
        self.assertEqual(len(s.live_word.source.steps),0)

    def test_BRMM_BIAS_source_ordinal_mismatch_fails_before_execution(self):
        s=self.root(); witness,r,b,raw,runtime=self.next_operands(s)
        bad=replace(witness,ordinal=2)
        with self.assertRaisesRegex(ValueError,'ordinals differ'):
            X.imu_step(s,restricted=r,bias_restricted=b,witness=bad,
                       raw=raw,packet_id='bad',**runtime)
        self.assertEqual(len(s.live_word.source.steps),0)

    def test_magnetic_event_preserves_admitted_histories_and_source_ordinal(self):
        s=self.root(); before=s.live_word.source.next_ordinal
        out=X.mag_step(s,**LOWER.mag_kwargs(s.live_word))
        self.assertEqual(out.state.admitted_history,s.admitted_history)
        self.assertEqual(out.state.bias_history,s.bias_history)
        self.assertEqual(out.state.live_word.source.next_ordinal,before)

    def test_readiness_records_join_without_promotion(self):
        r=X.readiness()
        for key in ('admitted_COMPLETE_BRMM_history_and_strong_runtime_joined',
                    'admitted_BIAS_history_and_strong_runtime_joined',
                    'each_IMU_event_requires_same_admitted_BRMM_BIAS_and_source_ordinal',
                    'same_admitted_segment_drives_source_bound_prediction_and_measurement',
                    'prediction_model_roots_cannot_bypass_admitted_history_edge',
                    'same_event_IMU_ISS_supply_attached_to_admitted_step',
                    'same_event_magnetic_ISS_supply_attached_to_admitted_product',
                    'BIAS_generating_history_attached'):
            self.assertTrue(r[key])
        for key in ('finite_tokens_used_as_source_membership_oracle',
                    'sensor_or_temperature_amplitude_bound_invented',
                    'startup_sample_zero_equal_to_admitted_history_restriction_proved',
                    'bounded_input_history_qualified',
                    'one_radian_attitude_guard_closed_for_every_admitted_prefix',
                    'deployment_exp_expm1_trig_Eigen_LDLT_closed',
                    'complete_600_step_shipping_word_composed_from_restrictions',
                    'storage_search_allowed','ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[key])


if __name__=='__main__': unittest.main()
