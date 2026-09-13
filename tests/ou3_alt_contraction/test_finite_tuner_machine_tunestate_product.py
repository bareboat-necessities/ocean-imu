"""Coherent tau/sigma/R_S machine TuneState product regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_common_alpha_qualification as A
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_frequency_binary32 as FREQ
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as X
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_binary32 as RAB
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_machine_real_join as RA
from tools.stability.ou3_alt_contraction import finite_tuner_rs_deployment_ledger as RS
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as S
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_deployment_ledger as SIG
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_machine_real_join as ST
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as TAU
import test_finite_complete_word_tau_qualification as QBASE


def sqrt_witness(x):
    lo,hi=ROOT.sqrt_enclosure(x); return B.rn32((lo+hi)/2)


def exp_witness(x):
    lo,hi,_,_=EXP.enclosure(x); return B.rn32((lo+hi)/2)


def mode_inputs(freq,av,bn):
    c=QBASE.shipping_runtime().candidate_cfg
    d=D.shipping_defaults(qeff_pow_result=B.rn32(1))
    stored=FREQ.store(B.rn32(F(freq)),B.rn32(c.min_freq),B.rn32(c.max_freq))
    # tau common-alpha exp witness from the same stored frequency/config.
    _,_,_,adapt=__import__('tools.stability.ou3_alt_contraction.finite_tuner_tau_binary32',fromlist=['x'])._floats_from_frequency(stored.stored_hz,c,TAU.DT)
    e=exp_witness(B.div(TAU.DT,adapt))
    return c,d,stored,e,B.rn32(F(av)),B.rn32(F(bn))


def build_sample(previous:X.State, *, sep=(F(1,5),F(1,4),F(1,10)), fma=(F(1,4),F(9,25),F(3,25)), pending=True):
    sc,sd,sf,se,sav,sbn=mode_inputs(*sep); fc,fd,ff,fe,fav,fbn=mode_inputs(*fma)
    if sc!=fc or sd!=fd: raise AssertionError('test fixture must share persistent configs')
    tau=TAU.step_tracks(previous.tau,separate_frequency=sf,fma_frequency=ff,cfg=sc,dt=TAU.DT,
                        separate_exp_decay=se,fma_exp_decay=fe)
    sa=A.qualify(tau.separate_step,sc,sd,dt=TAU.DT); fa=A.qualify(tau.fma_step,fc,fd,dt=TAU.DT)

    def sigma_target(av,bn):
        vn=B.mul(bn,bn); vw=max(S.VAR_FLOOR,max(B.rn32(0),B.sub(av,vn)))
        return S.target(sd,var_ready=True,accel_variance=av,band_noise_sigma=bn,still=False,
                        sqrt_result=sqrt_witness(vw))
    st=sigma_target(sav,sbn); ft=sigma_target(fav,fbn)
    sigma=SIG.step(previous.sigma,separate_target=st,separate_alpha=sa,fma_target=ft,fma_alpha=fa)

    def spectral(tau_step,sigma_target,freq):
        t=C.TargetState(F(freq),F(1),tau_step.tau_target,sigma_target.sigma_target)
        return ST.join(sd,t,pow_result=B.rn32(16),sqrt_result=B.rn32(F(1,5)))
    srs=spectral(tau.separate_step,st,sep[0]); frs=spectral(tau.fma_step,ft,fma[0])
    def rs_alpha(t):
        tau_t=t.target.tau_target; safe=min(max(tau_t,RAB.TIME_MIN),RAB.TIME_MAX)
        requested=B.mul(sd.adapt_RS_mult,safe); lo=min(max(TAU.DT,RAB.HORIZON_MIN),RAB.HORIZON_MAX)
        rssec=min(max(requested,lo),RAB.HORIZON_MAX); e=exp_witness(B.div(TAU.DT,rssec))
        return RA.join(sd,tau_target=tau_t,dt=TAU.DT,exp_decay=e)
    rs=RS.step(previous.rs,separate_target=srs,separate_alpha=rs_alpha(srs),
               fma_target=frs,fma_alpha=rs_alpha(frs))
    return X.compose_after_sample(previous,tau=tau,sigma=sigma,rs=rs,pending_after=pending)


class Tests(unittest.TestCase):
    def test_initial_product_uses_three_literal_shipping_seeds(self):
        s=X.initial()
        self.assertEqual(s.tau.separate,B.rn32(F(11,10)))
        self.assertEqual(s.sigma.separate,B.rn32(F(1,100)))
        self.assertEqual(s.rs.separate,B.rn32(F(1,2)))
        self.assertEqual(s.tau.updates,s.sigma.updates); self.assertEqual(s.sigma.updates,s.rs.updates)
        self.assertFalse(s.pending); self.assertIs(X.hold(s),s)

    def test_one_sample_advances_all_three_coherent_mode_tracks(self):
        out=build_sample(X.initial())
        self.assertEqual(out.state.tau.updates,1); self.assertEqual(out.state.sigma.updates,1); self.assertEqual(out.state.rs.updates,1)
        self.assertTrue(out.state.pending)
        self.assertIs(out.sigma.separate.alpha.step,out.tau.separate_step)
        self.assertIs(out.sigma.fma.alpha.step,out.tau.fma_step)
        self.assertEqual(out.rs.separate.target.target.sigma_target,out.sigma.separate.target.sigma_target)
        self.assertEqual(out.rs.fma.target.target.sigma_target,out.sigma.fma.target.sigma_target)

    def test_successive_sample_uses_whole_prior_machine_tunestate(self):
        first=build_sample(X.initial()).state; second=build_sample(first,pending=False)
        self.assertEqual(second.tau.separate_step.previous,first.tau.separate)
        self.assertEqual(second.sigma.separate.previous,first.sigma.separate)
        self.assertEqual(second.rs.separate.ema.previous,first.rs.separate)
        self.assertEqual(second.state.tau.updates,2); self.assertFalse(second.state.pending)

    def test_cross_mode_sigma_splice_is_rejected(self):
        out=build_sample(X.initial())
        crossed=replace(out.sigma,separate=out.sigma.fma,fma=out.sigma.separate)
        with self.assertRaisesRegex(ValueError,'compiler modes crossed|common-alpha source'):
            X.compose_after_sample(X.initial(),tau=out.tau,sigma=crossed,rs=out.rs,pending_after=True)

    def test_rs_sigma_target_splice_is_rejected(self):
        out=build_sample(X.initial())
        bad_sep=replace(out.rs.separate.target.target,sigma_target=B.rn32(F(7,10)))
        bad_join=replace(out.rs.separate.target,target=bad_sep)
        bad_track=replace(out.rs.separate,target=bad_join)
        bad_rs=replace(out.rs,separate=bad_track)
        with self.assertRaisesRegex(ValueError,'SpectralMSE sigma target'):
            X.compose_after_sample(X.initial(),tau=out.tau,sigma=out.sigma,rs=bad_rs,pending_after=True)

    def test_readiness_closes_joint_persistence_not_boundary_or_master(self):
        r=X.readiness()
        for k in ('shipping_tau_sigma_RS_persistent_machine_states_composed','all_three_scalar_ledgers_advance_with_one_common_update_index',
                  'sigma_common_alpha_is_same_mode_tau_step','SpectralMSE_RS_uses_same_mode_tau_and_sigma_targets',
                  'sigma_and_RS_share_same_deployment_config_per_mode','one_common_pending_bit_attached_after_whole_TuneState_successor'):
            self.assertTrue(r[k])
        for k in ('next_boundary_common_machine_commit_attached','startup_frontend_machine_TuneState_product_attached',
                  'Live_600_step_machine_TuneState_product_attached','all_target_libm_correspondence_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
