"""Interleaved admitted-source + bounded-ISS regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_imu_disturbance as DIST
from tools.stability.ou3_alt_contraction import finite_admitted_iss_interleaved_prefix as X
import test_finite_admitted_interleaved_prefix as BASE
import test_finite_admitted_source_imu_word as IBASE
import test_finite_source_bound_live_word as LBASE
import test_finite_source_bound_prediction_word as PBASE


def begin(bound=1000):
    p=BASE.begin()
    d=DIST.BoundedHistory(p.live.live_word.sensor_root,'temperature-history',F(bound))
    return X.begin(p,d)


class Tests(unittest.TestCase):
    def test_MAG_IMU_MAG_HOLD_preserves_one_bounded_history(self):
        t=begin(); hist=t.disturbance
        t,m0=X.mag_step(t,**LBASE.mag_kwargs(t.prefix.live.live_word))
        self.assertIs(t.disturbance,hist); self.assertEqual(t.imu_steps,0)

        witness,segment,raw,r,b,dynamic=IBASE.operands(t.prefix.live)
        t,i1=X.imu_step(t,restricted=r,bias_restricted=b,witness=witness,raw=raw,
                        packet_id='imu-1',**PBASE.root_args(),**dynamic)
        self.assertIs(t.disturbance,hist); self.assertEqual(t.imu_steps,1)
        self.assertEqual(i1.restriction.ordinal,1)
        self.assertEqual(i1.restriction.history,hist)

        t,m1=X.mag_step(t,**LBASE.mag_kwargs(t.prefix.live.live_word))
        self.assertIs(t.disturbance,hist); self.assertEqual(t.imu_steps,1)
        t,h=X.set_hold(t,hold=False)
        self.assertIs(t.disturbance,hist); self.assertEqual(t.imu_steps,1)
        self.assertEqual(tuple(e.kind for e in t.prefix.events),('mag','imu','mag','hold'))

    def test_wrong_sensor_root_rejected_before_interleaving(self):
        p=BASE.begin()
        other=DIST.BoundedHistory(
            type(p.live.live_word.sensor_root)(p.live.live_word.source.root,'other-gyro','other-accel'),
            'temperature-history',1000)
        with self.assertRaisesRegex(ValueError,'detached from carried sensor residual histories'):
            X.begin(p,other)

    def test_complete_word_cannot_be_claimed_from_short_prefix(self):
        t=begin()
        with self.assertRaisesRegex(ValueError,'exactly 600 physical source transitions'):
            X.complete(t)
        witness,segment,raw,r,b,dynamic=IBASE.operands(t.prefix.live)
        t,_=X.imu_step(t,restricted=r,bias_restricted=b,witness=witness,raw=raw,
                       packet_id='imu-1',**PBASE.root_args(),**dynamic)
        with self.assertRaisesRegex(ValueError,'exactly 600 physical source transitions'):
            X.complete(t)

    def test_readiness_keeps_complete_word_and_storage_closed(self):
        r=X.readiness()
        self.assertTrue(r['one_bounded_IMU_ISS_history_persists_through_entire_prefix'])
        self.assertTrue(r['only_IMU_consumes_kth_bounded_forcing_restriction'])
        self.assertTrue(r['MAG_and_HOLD_preserve_bounded_IMU_history_without_consumption'])
        self.assertTrue(r['bounded_input_history_structurally_attached'])
        self.assertTrue(r['complete_600_transition_strengthening_requires_actual_full_horizon'])
        self.assertTrue(r['complete_600_transition_strengthening_checks_exact_1_to_600_ordinals'])
        self.assertTrue(r['complete_600_transition_strengthening_checks_IMU_ledger_consecutivity'])
        self.assertTrue(r['complete_600_transition_strengthening_preserves_same_BRMM_BIAS_ISS_product'])
        self.assertFalse(r['all_event_arithmetic_witnesses_source_uniformly_qualified'])
        self.assertFalse(r['magnetic_counter_lifetime_closed'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
