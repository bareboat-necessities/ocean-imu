"""Exact-real SpectralMSE irrational-root enclosure regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as X


def cfg(*,shipping_cadence=False):
    ratio=F(1)
    pmin=F(1,100); pmax=F(12)
    if shipping_cadence:
        # Literal compiled float operands for 0.015f/1.1f, 0.005f and 0.15f.
        ratio=B.div(B.rn32(F(3,200)),B.rn32(F(11,10)))
        pmin=B.rn32(F(1,200)); pmax=B.rn32(F(3,20))
    return C.CandidateConfig(
        TARGET.FLOOR,TARGET.CEIL,F(1),F(1),TARGET.TAU_MIN,TARGET.TAU_MAX,F(4),
        ratio,pmin,pmax,F(3,20),F(100),F(1),F(1,2),F(1),
        B.rn32(F(9,5)),B.rn32(F(2,5)),F(3,2),F(0),F(1,10),True)


def sample(freq):
    return C.WaveBandSample(F(freq),True,F(1),F(0),False,F(0),F(1),F(1))


class Tests(unittest.TestCase):
    def test_legacy_exact_rational_cell_embeds_in_interval_theorem(self):
        c=cfg(); t=C.targets(sample(F(1,2)),c)
        self.assertEqual(t.tau_target,F(1)); self.assertEqual(t.sigma_target,F(1))
        box=X.enclose(c,t)
        self.assertTrue(X.contains_legacy_exact_witness(c,t,C.SpectralWitness(1,1),box))
        self.assertLess(box.sqrt_TS_hi-box.sqrt_TS_lo,F(1,1<<90))
        self.assertLess(box.u_pow_6_7_hi-box.u_pow_6_7_lo,F(1,1<<90))

    def test_real_shipping_prior_cell_no_longer_needs_rational_root_equality(self):
        c=cfg(shipping_cadence=True); t=C.targets(sample(F(1,5)),c)
        self.assertEqual(t.frequency,F(1,5)); self.assertEqual(t.tau_target,F(5,2))
        box=X.enclose(c,t)
        expected_TS=C.clamp(c.pseudo_tau_ratio*t.tau_target,c.pseudo_min,c.pseudo_max)
        self.assertEqual(box.TS,expected_TS)
        self.assertLessEqual(box.sqrt_TS_lo**2,box.TS)
        self.assertGreaterEqual(box.sqrt_TS_hi**2,box.TS)
        self.assertLessEqual(box.u_pow_6_7_lo**7,box.u**6)
        self.assertGreaterEqual(box.u_pow_6_7_hi**7,box.u**6)
        self.assertLess(box.sqrt_TS_hi-box.sqrt_TS_lo,F(1,1<<90))
        self.assertLess(box.u_pow_6_7_hi-box.u_pow_6_7_lo,F(1,1<<90))
        self.assertLessEqual(box.target_RS_lo,box.target_RS_hi)

    def test_enclosure_is_monotone_and_fail_closed(self):
        c=cfg(shipping_cadence=True); t=C.targets(sample(F(1,5)),c)
        box=X.enclose(c,t,bits=64)
        self.assertLessEqual(box.raw_RS_lo,box.raw_RS_hi)
        self.assertLessEqual(box.target_RS_lo,box.target_RS_hi)
        with self.assertRaisesRegex(ValueError,'at least 16'):
            X.enclose(c,t,bits=8)

    def test_readiness_advances_real_roots_not_binary32_libm_or_master_word(self):
        r=X.readiness()
        for k in ('sqrt_TS_real_root_enclosed_by_exact_rationals',
                  'u_pow_6_7_real_root_enclosed_by_exact_rationals',
                  'SpectralMSE_real_RS_target_interval_propagated_monotonically',
                  'irrational_roots_do_not_require_fake_rational_equalities',
                  'legacy_exact_rational_spectral_cells_embed_in_interval_theorem'):
            self.assertTrue(r[k])
        self.assertFalse(r['binary32_sqrt_target_libm_correspondence_closed'])
        self.assertFalse(r['binary32_pow_target_libm_correspondence_closed'])
        self.assertFalse(r['finite_tuner_candidate_interval_recurrence_composed'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
