"""Guard-persistent interval tuner frontend regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as G
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_interval_prefix as X
from tools.stability.ou3_alt_contraction import finite_tuner_candidate_interval as IC
from tools.stability.ou3_alt_contraction import finite_tuner_frontend_interval_prefix as T
import test_finite_guarded_tuner_prefix as BASE
import test_finite_complete_word_tau_qualification as QBASE


def tuner_state(stage='Live'):
    p=BASE.tuner_state()
    return T.State(p.vertical,p.wpe,p.band,p.stats,p.tracker_lpf,p.stillness,
                   IC.IntervalTuneState.point(p.tune),p.last_adapt_time,p.pending,
                   p.sample_index,p.time,stage,p.stage_time,p.warmup_sec)


def kwargs():
    k=BASE.kwargs(); k.pop('spectral',None)
    k['candidate_cfg']=QBASE.shipping_runtime().candidate_cfg
    return k


class Tests(unittest.TestCase):
    def test_first_guarded_postCold_sample_reaches_real_interval_candidate(self):
        s=X.State(G.State(),tuner_state())
        out=X.step(s,BASE.packet(),dt=BASE.DT,guard_cfg=G.Config(),**kwargs())
        self.assertTrue(out.state.guard.initialized)
        self.assertIs(out.tuner.raw_sample,out.guarded_sample)
        self.assertEqual(out.tuner.candidate.target.frequency,F(1,5))
        self.assertEqual(out.tuner.candidate.target.tau_target,F(5,2))
        self.assertLessEqual(out.state.tuner.tune.RS_lo,out.state.tuner.tune.RS_hi)

    def test_second_sample_cannot_restart_guard(self):
        s=X.State(G.State(),tuner_state())
        first=X.step(s,BASE.packet(),dt=BASE.DT,guard_cfg=G.Config(),**kwargs())
        second=X.step(first.state,BASE.packet(),dt=BASE.DT,guard_cfg=G.Config(),
                      guard_decay=G.DecayWitness(1,1,0,0),guard_rms=G.RmsWitness(0),**kwargs())
        self.assertTrue(second.state.guard.initialized)
        self.assertEqual(second.state.tuner.sample_index,2)
        self.assertIs(second.tuner.raw_sample,second.guarded_sample)

    def test_Cold_guard_path_preserves_interval_tune(self):
        s=X.State(G.State(),tuner_state('Cold')); before=s.tuner.tune
        k=kwargs()
        for n in ('candidate_cfg','sigma_wave_sqrt','ema'): k.pop(n,None)
        out=X.step(s,BASE.packet(),dt=BASE.DT,guard_cfg=G.Config(),**k)
        self.assertIs(out.state.tuner.tune,before)
        self.assertIsNone(out.tuner.candidate)
        self.assertTrue(out.state.guard.initialized)

    def test_readiness_closes_guarded_real_path_not_machine_or_MEKF(self):
        r=X.readiness()
        for k in ('raw_accel_guard_state_persists_across_interval_tuner_samples',
                  'one_guard_successor_feeds_private_vertical_band_tuner_and_WPE',
                  'guarded_sample_preserves_raw_COMPLETE_BRMM_sensor_ancestry',
                  'postCold_general_SpectralMSE_interval_candidate_composed',
                  'actual_0p2Hz_prior_postCold_candidate_representable_on_guarded_path'):
            self.assertTrue(r[k])
        self.assertFalse(r['guard_exp_sqrt_binary32_ancestry_attached'])
        self.assertFalse(r['Racc_inflation_from_same_guard_excess_attached'])
        self.assertFalse(r['finite_MEKF_event_interleaved_in_same_sample_word'])
        self.assertFalse(r['binary32_sqrt_pow_exp_correspondence_closed'])
        self.assertFalse(r['goLive_interval_RS_to_actual_MEKF_commit_closed'])
        self.assertFalse(r['source_uniform_complete_startup_reachability_closed'])
        self.assertFalse(r['storage_search_allowed'])


if __name__=='__main__': unittest.main()
