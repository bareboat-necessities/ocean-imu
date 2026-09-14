"""Interval-valued tracker-free tuner frontend regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate_interval as IC
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_frontend_interval_prefix as X
import test_finite_guarded_tuner_prefix as BASE


def state(stage='Live'):
    p=BASE.tuner_state()
    return X.State(p.vertical,p.wpe,p.band,p.stats,p.tracker_lpf,p.stillness,
                   IC.IntervalTuneState.point(p.tune),p.last_adapt_time,p.pending,
                   p.sample_index,p.time,stage,p.stage_time,p.warmup_sec)


def deployment_cfg(cache_result=1):
    # Explicit machine-cache witness; exact-real frontend consumes only the
    # source-rooted qeff interval, not this unqualified binary32 result.
    return D.shipping_defaults(qeff_pow_result=B.rn32(cache_result))


def kwargs(cache_result=1):
    k=BASE.kwargs(); k.pop('spectral',None)
    k['candidate_cfg']=deployment_cfg(cache_result)
    return k


class Tests(unittest.TestCase):
    def test_real_postCold_prior_cell_executes_with_deployment_interval_RS(self):
        s=state('Live'); out=X.step(s,BASE.packet(),dt=BASE.DT,**kwargs())
        self.assertIsNotNone(out.candidate)
        self.assertEqual(out.candidate.target.frequency,F(1,5))
        self.assertEqual(out.candidate.target.tau_target,F(5,2))
        self.assertLess(out.candidate.spectral.sqrt_TS_lo,out.candidate.spectral.sqrt_TS_hi)
        self.assertLess(out.candidate.spectral.u_pow_6_7_lo,out.candidate.spectral.u_pow_6_7_hi)
        self.assertIsNotNone(out.candidate.spectral.qeff_root_lo)
        self.assertLessEqual(out.state.tune.RS_lo,out.state.tune.RS_hi)
        self.assertEqual(out.preupdate_wpe.state,s.wpe)
        self.assertEqual(out.state.sample_index,s.sample_index+1)

    def test_machine_qeff_cache_change_cannot_move_exact_real_frontend_interval(self):
        s=state('Live')
        a=X.step(s,BASE.packet(),dt=BASE.DT,**kwargs(1))
        b=X.step(s,BASE.packet(),dt=BASE.DT,**kwargs(2))
        self.assertEqual(a.candidate.spectral.raw_RS_lo,b.candidate.spectral.raw_RS_lo)
        self.assertEqual(a.candidate.spectral.raw_RS_hi,b.candidate.spectral.raw_RS_hi)

    def test_current_sample_WPE_cannot_feed_its_own_interval_candidate(self):
        s=state('Live'); out=X.step(s,BASE.packet(),dt=BASE.DT,**kwargs())
        self.assertEqual(out.candidate.target.frequency,F(1,5))
        self.assertIs(out.preupdate_wpe.state,s.wpe)
        self.assertEqual(out.wpe.state,out.state.wpe)

    def test_Cold_branch_preserves_interval_tune_and_consumes_no_candidate_operands(self):
        p=state('Cold'); k=kwargs()
        for n in ('candidate_cfg','sigma_wave_sqrt','ema'): k.pop(n,None)
        before=p.tune
        out=X.step(p,BASE.packet(),dt=BASE.DT,**k)
        self.assertIs(out.state.tune,before); self.assertIsNone(out.candidate)
        bad=dict(k); bad['candidate_cfg']=deployment_cfg()
        with self.assertRaisesRegex(ValueError,'Cold tuner branch'):
            X.step(p,BASE.packet(),dt=BASE.DT,**bad)

    def test_pending_state_cannot_skip_commit_boundary(self):
        s=replace(state('Live'),pending=True)
        with self.assertRaisesRegex(ValueError,'pending interval tuner state'):
            X.step(s,BASE.packet(),dt=BASE.DT,**kwargs())

    def test_readiness_advances_real_frontend_not_goLive_or_machine_libm(self):
        r=X.readiness()
        for k in ('same_private_vertical_WPE_band_stillness_order_as_shipping_frontend',
                  'sample_entry_WPE_view_precedes_current_sample_WPE_update',
                  'postCold_general_SpectralMSE_interval_candidate_composed',
                  'tau_sigma_frontend_state_remains_exact',
                  'RS_frontend_state_is_rigorous_exact_real_interval',
                  'actual_0p2Hz_prior_postCold_candidate_representable'):
            self.assertTrue(r[k])
        self.assertFalse(r['point_RS_representative_selected'])
        self.assertFalse(r['goLive_interval_RS_to_actual_MEKF_commit_closed'])
        self.assertFalse(r['binary32_sqrt_pow_exp_correspondence_closed'])
        self.assertFalse(r['source_uniform_complete_startup_reachability_closed'])
        self.assertFalse(r['storage_search_allowed']); self.assertFalse(r['ALT_STARTUP_PASS'])


if __name__=='__main__': unittest.main()
