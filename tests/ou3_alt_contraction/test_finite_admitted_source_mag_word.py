"""Admitted-history magnetic composition regressions; not stability evidence."""
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_source_imu_word as IMU
from tools.stability.ou3_alt_contraction import finite_admitted_source_mag_word as X
import test_finite_admitted_source_imu_word as IBASE
import test_finite_source_bound_live_word as BASE
import test_finite_source_bound_prediction_word as PBASE


class Tests(unittest.TestCase):
    def first_imu(self):
        s=IBASE.root(); witness,segment,raw,r,b,dynamic=IBASE.operands(s)
        return IMU.imu_step(s,restricted=r,bias_restricted=b,witness=witness,raw=raw,
                            packet_id='imu-1',**PBASE.root_args(),**dynamic)

    def test_post_first_IMU_mag_preserves_both_admitted_histories(self):
        first=self.first_imu(); s=first.state
        before_source=s.live_word.source
        out=X.mag_step(s,**BASE.mag_kwargs(s.live_word))
        self.assertEqual(out.state.admitted_history,s.admitted_history)
        self.assertEqual(out.state.bias_history,s.bias_history)
        self.assertIs(out.state.live_word.source,before_source)
        self.assertEqual(len(out.state.live_word.source.steps),1)

    def test_sample_zero_fails_closed_until_startup_history_equality_bridge(self):
        s=IBASE.root()
        with self.assertRaisesRegex(ValueError,'sample-zero admitted-history magnetic equality'):
            X.mag_step(s,**BASE.mag_kwargs(s.live_word))
        self.assertEqual(len(s.live_word.source.steps),0)

    def test_readiness_preserves_missing_startup_and_sensor_obligations(self):
        r=X.readiness()
        self.assertTrue(r['post_first_IMU_magnetic_edge_uses_carried_admitted_source_endpoint'])
        self.assertTrue(r['admitted_BRMM_and_BIAS_histories_preserved_across_magnetic_event'])
        self.assertTrue(r['same_event_correlated_magnetic_forcing_retained'])
        self.assertFalse(r['fresh_sample_zero_magnetic_admitted_history_equality_closed'])
        self.assertFalse(r['sensor_residual_admissibility_attached'])
        self.assertFalse(r['complete_interleaved_600_step_word_composed'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
