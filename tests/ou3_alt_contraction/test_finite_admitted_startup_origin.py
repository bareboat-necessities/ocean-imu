"""Admitted sample-zero startup identity regressions; not reachability proof."""
from dataclasses import replace
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as A
from tools.stability.ou3_alt_contraction import finite_admitted_startup_origin as X
import test_finite_admitted_source_live_word as BASE


class Tests(unittest.TestCase):
    def test_fresh_Live_reference_is_exact_admitted_sample_zero(self):
        s=BASE.Tests().root()
        ref=s.live_word.live.live.live.mekf.reference
        origin=X.RestrictedOrigin(s.admitted_history,ref)
        bound=X.bind(s,origin)
        self.assertEqual(bound.origin.reference,ref)
        self.assertEqual(bound.admitted.live_word.source.next_ordinal,1)

    def test_different_admitted_history_cannot_claim_same_sample_zero(self):
        s=BASE.Tests().root(); ref=s.live_word.live.live.live.mekf.reference
        other=A.AdmittedHistory('other-history')
        origin=X.RestrictedOrigin(other,ref)
        with self.assertRaisesRegex(ValueError,'detached from carried admitted history'):
            X.bind(s,origin)

    def test_modified_physical_reference_is_not_startup_truth(self):
        s=BASE.Tests().root(); ref=s.live_word.live.live.live.mekf.reference
        altered=replace(ref,position=(1,0,0))
        origin=X.RestrictedOrigin(s.admitted_history,altered)
        with self.assertRaisesRegex(ValueError,'not admitted-history sample zero'):
            X.bind(s,origin)

    def test_origin_requires_live_time_and_zero_centered_S(self):
        s=BASE.Tests().root(); ref=s.live_word.live.live.live.mekf.reference
        with self.assertRaisesRegex(ValueError,'one-time Live origin'):
            X.RestrictedOrigin(s.admitted_history,replace(ref,time=ref.time+1))
        with self.assertRaisesRegex(ValueError,'centered S'):
            X.RestrictedOrigin(s.admitted_history,replace(ref,centered_S=(1,0,0)))

    def test_readiness_closes_identity_not_reachability_or_storage(self):
        r=X.readiness()
        self.assertTrue(r['admitted_history_sample_zero_restriction_datum_explicit'])
        self.assertTrue(r['fresh_joint24_physical_reference_equals_admitted_sample_zero'])
        self.assertTrue(r['startup_sample_zero_equal_to_admitted_history_restriction_proved'])
        self.assertFalse(r['finite_endpoint_checks_used_as_membership_oracle'])
        self.assertFalse(r['startup_reachability_for_every_admitted_history_proved'])
        self.assertFalse(r['deployment_binary32_correspondence_closed'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_STARTUP_PASS'])


if __name__=='__main__': unittest.main()
