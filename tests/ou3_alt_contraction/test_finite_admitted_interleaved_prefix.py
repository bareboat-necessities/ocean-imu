"""Strong admitted interleaving prefix regressions; not stability evidence."""
from dataclasses import replace
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as ABRMM
from tools.stability.ou3_alt_contraction import finite_admitted_interleaved_prefix as X
import test_finite_admitted_source_imu_word as IBASE
import test_finite_source_bound_live_word as BASE
import test_finite_source_bound_prediction_word as PBASE


def begin():
    live=IBASE.root()
    ref=live.live_word.live.live.live.mekf.reference
    origin=ABRMM.RestrictedOrigin(live.admitted_history,ref)
    return X.begin(live,origin)


class Tests(unittest.TestCase):
    def test_sample0_MAG_then_IMU_then_MAG_then_HOLD_is_one_successor_chain(self):
        t=begin(); origin=t.origin
        t,m0=X.mag_step(t,**BASE.mag_kwargs(t.live.live_word))
        self.assertEqual(t.imu_steps,0); self.assertEqual(t.origin,origin)

        witness,segment,raw,r,b,dynamic=IBASE.operands(t.live)
        t,i1=X.imu_step(t,restricted=r,bias_restricted=b,witness=witness,raw=raw,
                        packet_id='imu-1',**PBASE.root_args(),**dynamic)
        self.assertEqual(t.imu_steps,1); self.assertEqual(t.origin,origin)
        self.assertEqual(t.live.live_word.source.steps[0].segment.before,origin.endpoint)
        before=t.live.live_word.source
        t,m1=X.mag_step(t,**BASE.mag_kwargs(t.live.live_word))
        self.assertIs(t.live.live_word.source,before); self.assertEqual(t.origin,origin)
        t,h=X.set_hold(t,hold=False)
        self.assertEqual(t.imu_steps,1); self.assertEqual(t.origin,origin)
        self.assertEqual(tuple(e.kind for e in t.events),('mag','imu','mag','hold'))
        self.assertEqual(tuple((e.source_steps_before,e.source_steps_after) for e in t.events),
                         ((0,0),(0,1),(1,1),(1,1)))

    def test_caller_cannot_replace_persistent_origin_on_mag_edge(self):
        t=begin()
        with self.assertRaisesRegex(TypeError,'owns the admitted origin'):
            X.mag_step(t,origin=t.origin,**BASE.mag_kwargs(t.live.live_word))

    def test_first_IMU_must_start_at_persistent_origin(self):
        t=begin(); witness,segment,raw,r,b,dynamic=IBASE.operands(t.live)
        # Mutating the physical predecessor alone now violates the lower-level
        # moment/bias same-history recurrence.  That is the earliest valid
        # rejection point and is stronger than waiting for the interleave origin guard.
        with self.assertRaisesRegex(ValueError,'physical moment/bias recurrence is not one history'):
            replace(segment,before=replace(segment.before,position=(1,0,0)))
        self.assertEqual(t.imu_steps,0)
        self.assertEqual(t.origin.endpoint,segment.before)

    def test_event_record_rejects_non_IMU_source_advance(self):
        with self.assertRaisesRegex(ValueError,'ordinal accounting'):
            X.EventRecord('mag',2,3)
        with self.assertRaisesRegex(ValueError,'ordinal accounting'):
            X.EventRecord('imu',2,2)

    def test_readiness_distinguishes_composable_horizon_from_qualified_word(self):
        r=X.readiness()
        for key in ('finite_successor_induction_over_IMU_MAG_HOLD_closed',
                    'canonical_admitted_tL_origin_persistent_in_interleaved_product',
                    'startup_sample_zero_identity_consumed_before_interleaving',
                    'first_IMU_forced_to_start_at_same_admitted_tL_origin',
                    'only_IMU_advances_admitted_source_ordinal',
                    'sample_zero_MAG_uses_persistent_admitted_origin',
                    'both_admitted_histories_and_origin_persist_through_all_event_types',
                    'canonical_600_transition_source_horizon_representable'):
            self.assertTrue(r[key])
        self.assertFalse(r['all_event_arithmetic_witnesses_source_uniformly_qualified'])
        self.assertFalse(r['sensor_disturbance_admissibility_attached'])
        self.assertFalse(r['startup_reachability_to_admitted_fresh_state_closed'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
