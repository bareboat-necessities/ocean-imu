"""Binary32 SeaStateAutoTuner frequency-storage regressions; not WPE qualification."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_frequency_binary32 as X
from tools.stability.ou3_alt_contraction import finite_tuner_tau_binary32 as TAU
import test_finite_tuner_candidate as TC


class Tests(unittest.TestCase):
    def test_shipping_store_clamps_binary32_input_and_getter_is_identity(self):
        lo=B.rn32(F(1,20)); hi=B.rn32(5)
        mid=B.rn32(F(3,10))
        s=X.store(mid,lo,hi)
        self.assertEqual(s.stored_hz,mid)
        self.assertEqual(X.get_frequency_hz(s),mid)
        self.assertEqual(X.store(B.rn32(F(1,100)),lo,hi).stored_hz,lo)
        self.assertEqual(X.store(B.rn32(6),lo,hi).stored_hz,hi)

    def test_store_rejects_exact_real_value_not_already_binary32(self):
        with self.assertRaisesRegex(ValueError,'already be binary32'):
            X.store(F(1,10),B.rn32(F(1,20)),B.rn32(5))

    def test_tau_strong_entry_consumes_stored_float_without_requantizing(self):
        c=TC.cfg(); stored=X.store(B.rn32(F(1,2)),B.rn32(F(1,10)),B.rn32(2))
        # Conservative binary32 exp witness inside the exact-real enclosure.
        e=B.rn32(F(99,100))
        out=TAU.step_from_stored_frequency(B.rn32(F(1,2)),stored,c,dt=F(1,100),exp_decay=e)
        self.assertEqual(out.frequency,stored.stored_hz)

    def test_detached_or_nonbinary_frequency_cannot_reach_strong_tau_edge(self):
        c=TC.cfg(); e=B.rn32(F(99,100))
        with self.assertRaises(TypeError):
            TAU.step_from_stored_frequency(B.rn32(F(1,2)),F(1,2),c,dt=F(1,100),exp_decay=e)

    def test_readiness_closes_store_to_tau_edge_but_not_upstream_wpe(self):
        r=X.readiness(); t=TAU.readiness()
        self.assertTrue(r['shipping_frequency_store_source_shape_matches'])
        self.assertTrue(r['frequency_clamp_and_store_exact_binary32'])
        self.assertTrue(r['getFrequencyHz_is_identity_on_stored_binary32'])
        self.assertFalse(r['upstream_WPE_binary32_frequency_production_closed'])
        self.assertTrue(t['source_tuner_frequency_binary32_store_to_tau_edge_closed'])
        self.assertTrue(t['source_frontend_frequency_binary32_storage_correspondence_closed'])
        self.assertFalse(t['upstream_WPE_binary32_frequency_production_closed'])
        self.assertFalse(t['tuner_exp_libm_binary32_correspondence_closed'])
        self.assertFalse(t['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
