"""Binary32 WPE log-period global compiler-track regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as X


class Tests(unittest.TestCase):
    def test_nonproducing_sample_preserves_both_log_tracks(self):
        s=X.initial(); out=X.hold_sample(s)
        self.assertEqual(out.state.separate,s.separate); self.assertEqual(out.state.fma,s.fma)
        self.assertEqual(out.state.samples,1); self.assertFalse(out.produced_period)

    def test_first_valid_period_initializes_each_global_track_from_its_own_log_result(self):
        a=X.InitWitness(B.rn32(F(7,10))); b=X.InitWitness(B.rn32(F(7001,10000)))
        out=X.first_valid(X.initial(),separate=a,fma=b)
        self.assertEqual(out.state.separate.log_period,a.log_raw)
        self.assertEqual(out.state.fma.log_period,b.log_raw)
        self.assertNotEqual(out.state.separate.log_period,out.state.fma.log_period)
        self.assertEqual(out.state.separate.valid_updates,1); self.assertEqual(out.state.fma.valid_updates,1)

    def test_smoothed_update_keeps_mode_specific_libm_and_arithmetic_histories(self):
        s=X.first_valid(X.initial(),separate=X.InitWitness(B.rn32(F(7,10))),
                        fma=X.InitWitness(B.rn32(F(7001,10000)))).state
        sw=X.SmoothWitness(B.rn32(F(3,4)),B.rn32(2),B.rn32(F(99,100)))
        fw=X.SmoothWitness(B.rn32(F(751,1000)),B.rn32(F(201,100)),B.rn32(F(989,1000)))
        out=X.smooth_valid(s,separate=sw,fma=fw)
        self.assertFalse(out.separate.contracted); self.assertTrue(out.fma.contracted)
        self.assertEqual(out.separate.sea_period_exp,sw.sea_period_exp)
        self.assertEqual(out.fma.sea_period_exp,fw.sea_period_exp)
        self.assertEqual(out.state.separate.valid_updates,2); self.assertEqual(out.state.fma.valid_updates,2)
        self.assertNotEqual(out.state.separate.log_period,out.state.fma.log_period)

    def test_cross_track_identity_is_not_required_and_invalid_branch_fails_closed(self):
        s=X.initial()
        with self.assertRaisesRegex(ValueError,'initialized WPE track'):
            X._init(X.Track(B.rn32(F(1,2)),1),X.InitWitness(B.rn32(F(3,5))),False)
        with self.assertRaisesRegex(ValueError,'smoothing requires initialized'):
            X.smooth_valid(s,separate=X.SmoothWitness(B.rn32(1),B.rn32(2),B.rn32(F(99,100))),
                           fma=X.SmoothWitness(B.rn32(1),B.rn32(2),B.rn32(F(99,100))))

    def test_readiness_closes_topology_not_libm_or_complete_word(self):
        r=X.readiness()
        self.assertTrue(r['shipping_WPE_log_update_source_shape_matches'])
        self.assertTrue(r['global_separate_and_FMA_WPE_log_tracks_persist_without_branch_explosion'])
        self.assertTrue(r['global_compiler_track_coherence_includes_WPE_log_state'])
        self.assertFalse(r['WPE_log_std_log_target_libm_correspondence_closed'])
        self.assertFalse(r['WPE_log_exp_target_libm_correspondence_closed'])
        self.assertFalse(r['WPE_raw_period_binary32_production_closed'])
        self.assertFalse(r['storage_search_allowed']); self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
