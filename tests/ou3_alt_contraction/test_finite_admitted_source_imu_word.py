"""Strong admitted-source IMU edge regressions; not stability evidence."""
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as ABRMM
from tools.stability.ou3_alt_contraction import finite_admitted_bias_history as ABIAS
from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as ADLIVE
from tools.stability.ou3_alt_contraction import finite_admitted_source_imu_word as X
import test_finite_source_bound_live_word as BASE
import test_finite_source_bound_prediction_word as PBASE


def root():
    lower=BASE.root_state()
    physical=ABRMM.AdmittedHistory(lower.source.root.history_id)
    bias=ABIAS.AdmittedBiasHistory(lower.bias_history_id,lower.source.root.bias_family,BASE.PHI)
    return ADLIVE.State(lower,physical,bias)


def operands(state):
    witness,segment,raw,dynamic=PBASE.operands(state.live_word)
    restricted=ABRMM.RestrictedSegment(state.admitted_history,witness.ordinal,segment)
    bias_restricted=ABIAS.RestrictedBiasStep(state.bias_history,witness.ordinal,segment)
    return witness,segment,raw,restricted,bias_restricted,dynamic


class Tests(unittest.TestCase):
    def test_one_edge_composes_admission_prediction_and_executed_forcing(self):
        s=root(); witness,segment,raw,r,b,dynamic=operands(s)
        out=X.imu_step(s,restricted=r,bias_restricted=b,witness=witness,raw=raw,
                       packet_id='imu-1',**PBASE.root_args(),**dynamic)
        self.assertEqual(out.state.admitted_history,s.admitted_history)
        self.assertEqual(out.state.bias_history,s.bias_history)
        self.assertEqual(out.state.live_word.source.steps[-1].witness.ordinal,1)
        self.assertEqual(out.state.live_word.live.live.live.mekf.reference,segment.after)
        self.assertEqual(out.forcing.gyro_residual_internal,raw.gyro_residual_internal)
        self.assertEqual(out.forcing.temperature_delta_c,0)

    def test_detached_bias_restriction_fails_before_shipping_event(self):
        s=root(); witness,_,raw,r,b,dynamic=operands(s)
        other=ABIAS.AdmittedBiasHistory('other-bias',b.history.family,b.history.phi_true)
        detached=ABIAS.RestrictedBiasStep(other,witness.ordinal,b.segment)
        with self.assertRaisesRegex(ValueError,'detached from carried admitted BIAS history'):
            X.imu_step(s,restricted=r,bias_restricted=detached,witness=witness,raw=raw,
                       packet_id='bad',**PBASE.root_args(),**dynamic)
        self.assertEqual(len(s.live_word.source.steps),0)

    def test_prediction_roots_still_cannot_be_injected_through_stronger_entry(self):
        s=root(); witness,_,raw,r,b,dynamic=operands(s)
        dynamic=dict(dynamic); dynamic['ou']='detached'
        with self.assertRaisesRegex(TypeError,'cannot be overridden'):
            X.imu_step(s,restricted=r,bias_restricted=b,witness=witness,raw=raw,
                       packet_id='bad',**PBASE.root_args(),**dynamic)

    def test_readiness_keeps_complete_word_and_storage_closed(self):
        q=X.readiness()
        self.assertTrue(q['admitted_BRMM_and_BIAS_histories_carried_jointly'])
        self.assertTrue(q['same_kth_BRMM_BIAS_restriction_consumed_before_shipping_execution'])
        self.assertTrue(q['same_restriction_drives_source_bound_prediction_coefficients'])
        self.assertTrue(q['same_executed_packet_forcing_retained'])
        self.assertTrue(q['physical_moments_sensor_forcing_and_model_roots_share_one_event'])
        self.assertFalse(q['sensor_residual_admissibility_attached'])
        self.assertFalse(q['startup_sample_zero_admitted_history_equality_closed'])
        self.assertFalse(q['complete_600_step_shipping_word_composed'])
        self.assertFalse(q['storage_search_allowed'])
        self.assertFalse(q['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
