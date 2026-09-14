"""Persistent dual-compiler sigma deployment-ledger regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_common_alpha_qualification as A
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as S
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_deployment_ledger as X
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT
from tools.stability.ou3_alt_contraction import finite_tuner_tau_binary32 as T
import test_finite_complete_word_tau_qualification as QBASE


def sqrt_witness(x):
    lo,hi=ROOT.sqrt_enclosure(x); return B.rn32((lo+hi)/2)


def mode(freq,av,bn):
    c=QBASE.shipping_runtime().candidate_cfg; d=D.shipping_defaults(qeff_pow_result=B.rn32(1)); dt=B.rn32(F(1,200)); f=B.rn32(F(freq))
    _,_,_,adapt=T._floats_from_frequency(f,c,dt); x=B.div(dt,adapt); elo,ehi,_,_=EXP.enclosure(x); e=B.rn32((elo+ehi)/2)
    tau=T._step_from_frequency(B.rn32(F(11,10)),f,c,dt=dt,exp_decay=e)
    alpha=A.qualify(tau,c,d,dt=dt)
    av=B.rn32(F(av)); bn=B.rn32(F(bn)); vn=B.mul(bn,bn); vw=max(B.rn32(0),B.sub(av,vn)); vw=max(vw,S.VAR_FLOOR)
    target=S.target(d,var_ready=True,accel_variance=av,band_noise_sigma=bn,still=False,sqrt_result=sqrt_witness(vw))
    return target,alpha


class Tests(unittest.TestCase):
    def test_initial_state_is_literal_shipping_sigma_seed(self):
        s=X.initial(); self.assertEqual(s.separate,B.rn32(F(1,100))); self.assertEqual(s.fma,B.rn32(F(1,100)))
        self.assertEqual(s.updates,0); self.assertIs(X.hold(s),s)

    def test_distinct_mode_coherent_targets_and_alphas_advance_global_tracks(self):
        st,sa=mode(F(1,5),F(1,4),F(1,10)); ft,fa=mode(F(1,4),F(9,25),F(3,25))
        out=X.step(X.initial(),separate_target=st,separate_alpha=sa,fma_target=ft,fma_alpha=fa)
        self.assertEqual(out.state.separate,B.ema(X.INITIAL,st.sigma_target,sa.alpha,contracted=False))
        self.assertEqual(out.state.fma,B.ema(X.INITIAL,ft.sigma_target,fa.alpha,contracted=True))
        self.assertEqual(out.state.updates,1)
        self.assertNotEqual(st.sigma_target,ft.sigma_target)

    def test_successive_updates_use_each_track_previous_state(self):
        st,sa=mode(F(1,5),F(1,4),F(1,10)); ft,fa=mode(F(1,4),F(9,25),F(3,25))
        first=X.step(X.initial(),separate_target=st,separate_alpha=sa,fma_target=ft,fma_alpha=fa)
        second=X.step(first.state,separate_target=st,separate_alpha=sa,fma_target=ft,fma_alpha=fa)
        self.assertEqual(second.separate.previous,first.state.separate)
        self.assertEqual(second.fma.previous,first.state.fma)
        self.assertEqual(second.state.updates,2)

    def test_alpha_from_other_deployment_config_is_rejected(self):
        st,sa=mode(F(1,5),F(1,4),F(1,10)); ft,fa=mode(F(1,4),F(9,25),F(3,25))
        # Build a sigma target with a distinct qeff cache; common scalar fields
        # are unchanged, but one persistent ledger step still requires one exact
        # deployment config object across both tracks.
        other=D.shipping_defaults(qeff_pow_result=B.rn32(F(9,10)))
        bad=S.Target(other,ft.var_ready,ft.accel_variance,ft.band_noise_sigma,ft.still,ft.still_time,
                     ft.attenuation,ft.var_noise,ft.var_total,ft.var_wave_pre_attenuation,
                     ft.var_wave_attenuated,ft.var_wave,ft.sqrt_result,ft.sigma_wave,ft.scaled_sigma,ft.sigma_target)
        with self.assertRaisesRegex(ValueError,'common deployment config'):
            X.step(X.initial(),separate_target=st,separate_alpha=sa,fma_target=bad,fma_alpha=fa)

    def test_readiness_closes_persistence_not_boundary_or_upstream_libm(self):
        r=X.readiness()
        for k in ('shipping_sigma_seed_and_update_source_shape_matches','source_locked_binary32_sigma_seed_0p01',
                  'persistent_separate_and_FMA_sigma_histories_carried','sigma_each_track_consumes_qualified_common_tau_sigma_alpha',
                  'compiler_histories_may_consume_distinct_sigma_targets_and_alphas','cold_or_nonadapting_sigma_identity_branch_materialized'):
            self.assertTrue(r[k])
        for k in ('pending_common_TuneState_boundary_attached','startup_frontend_sigma_machine_history_attached',
                  'Live_600_step_sigma_machine_history_attached','sigma_sqrt_and_still_exp_libm_correspondence_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
