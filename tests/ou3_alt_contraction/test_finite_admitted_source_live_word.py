"""Admitted-history restriction composition regressions; not stability evidence."""
from dataclasses import replace
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as A
from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as X
import test_finite_source_bound_live_word as BASE


class Tests(unittest.TestCase):
    def root(self):
        lower=BASE.root_state()
        history=A.AdmittedHistory(lower.source.root.history_id)
        return X.State(lower,history)

    def test_next_IMU_consumes_kth_restriction_of_same_admitted_history(self):
        s=self.root(); witness,segment,raw,dynamic=BASE.next_imu_operands(s.live_word)
        restricted=A.RestrictedSegment(s.admitted_history,witness.ordinal,segment)
        out=X.imu_step(s,restricted=restricted,witness=witness,raw=raw,
                       packet_id='imu-1',**dynamic)
        self.assertEqual(out.state.admitted_history,s.admitted_history)
        self.assertEqual(out.state.live_word.source.steps[-1].witness.ordinal,1)
        self.assertEqual(out.state.live_word.live.live.live.mekf.reference,segment.after)

    def test_different_admitted_history_cannot_supply_next_segment(self):
        s=self.root(); witness,segment,raw,dynamic=BASE.next_imu_operands(s.live_word)
        other=A.AdmittedHistory('different-primary-history')
        restricted=A.RestrictedSegment(other,witness.ordinal,segment)
        with self.assertRaisesRegex(ValueError,'detached from carried admitted history'):
            X.imu_step(s,restricted=restricted,witness=witness,raw=raw,
                       packet_id='bad',**dynamic)

    def test_restriction_ordinal_must_match_source_witness(self):
        s=self.root(); witness,segment,_,_=BASE.next_imu_operands(s.live_word)
        restricted=A.RestrictedSegment(s.admitted_history,1,segment)
        bad=replace(witness,ordinal=2)
        with self.assertRaisesRegex(ValueError,'ordinal detached'):
            A.qualify_step(s.live_word.source.root,restricted,bad)

    def test_history_is_theorem_quantifier_not_runtime_admission_flag(self):
        h=A.AdmittedHistory('H')
        self.assertEqual(h.canonical_source,A.CANONICAL_SOURCE)
        with self.assertRaises(TypeError):
            A.AdmittedHistory('H',admitted=True)
        r=A.readiness()
        self.assertTrue(r['universal_theorem_quantifier_over_admitted_primary_history_explicit'])
        self.assertTrue(r['finite_word_derived_by_restriction_not_runtime_membership_inference'])
        self.assertTrue(r['arbitrary_runtime_tokens_do_not_prove_COMPLETE_BRMM_membership'])
        self.assertFalse(r['complete_600_step_shipping_word_composed_from_restrictions'])

    def test_readiness_does_not_promote_source_or_storage(self):
        r=X.readiness()
        self.assertTrue(r['universal_admitted_COMPLETE_BRMM_history_carried_in_Live_product'])
        self.assertTrue(r['every_IMU_event_requires_same_history_restriction_ordinal'])
        self.assertFalse(r['startup_sample_zero_equal_to_admitted_history_restriction_proved'])
        self.assertFalse(r['BIAS_generating_history_attached'])
        self.assertFalse(r['sensor_disturbance_admissibility_attached'])
        self.assertFalse(r['complete_600_step_shipping_word_composed_from_restrictions'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
