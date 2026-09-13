"""Machine band-noise-floor binary32 regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_band_binary32_contraction as C
from tools.stability.ou3_alt_contraction import finite_band_machine_ledger as L
from tools.stability.ou3_alt_contraction import finite_band_noise_floor_binary32 as X
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as RNE
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT


def sqrtw(x):
    lo,hi=ROOT.sqrt_enclosure(x); return B.rn32((lo+hi)/2)


class Tests(unittest.TestCase):
    def test_unready_band_returns_bench_sigma_exactly(self):
        s=L.initial(); bench=B.rn32(F(3,100))
        out=X.evaluate(s,bench_sigma=bench)
        self.assertEqual(out.noise_sigma,bench); self.assertIsNone(out.sqrt_gain)

    def test_ready_band_uses_same_stored_p11_and_binary32_multiply(self):
        p11=B.rn32(F(1,4)); machine=C.State(B.rn32(0),B.rn32(0),B.rn32(F(1,5)),B.rn32(0),p11,True)
        s=L.State(machine,7); bench=B.rn32(F(3,100)); sg=sqrtw(p11)
        out=X.evaluate(s,bench_sigma=bench,sqrt_gain=sg)
        self.assertEqual(out.sqrt_gain,sg); self.assertEqual(out.noise_sigma,B.mul(bench,sg))
        lo,hi=ROOT.sqrt_enclosure(p11); self.assertTrue(RNE._interval_hits_rne_cell(lo,hi,sg))

    def test_detached_ready_sqrt_fails_closed(self):
        p11=B.rn32(F(1,4)); machine=C.State(B.rn32(0),B.rn32(0),B.rn32(0),B.rn32(0),p11,True)
        with self.assertRaisesRegex(ValueError,'sqrt witness detached'):
            X.evaluate(L.State(machine,1),bench_sigma=B.rn32(F(1,100)),sqrt_gain=B.rn32(F(3,4)))

    def test_readiness_closes_consumer_not_platform_membership(self):
        r=X.readiness()
        for k in ('shipping_band_noise_floor_source_shape_matches','unready_machine_band_returns_bench_sigma_by_identity',
                  'ready_machine_band_consumes_same_stored_p11_gain','sqrt_gain_related_to_exact_root_by_binary32_RNE_cell',
                  'bench_times_sqrt_gain_binary32_multiply_materialized','common_TuneState_boundary_can_consume_this_noise_floor'):
            self.assertTrue(r[k])
        for k in ('target_sqrt_libm_correspondence_closed','machine_band_execution_membership_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
