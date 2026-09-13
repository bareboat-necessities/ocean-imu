"""Whole tau/sigma/R_S machine-candidate transition regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_frequency_binary32 as FREQ
from tools.stability.ou3_alt_contraction import finite_tuner_tau_binary32 as TB
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as TL
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as SM
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_machine_real_join as SJ
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_binary32 as RA
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as P
from tools.stability.ou3_alt_contraction import finite_tuner_machine_candidate_step as X
import test_finite_complete_word_tau_qualification as QBASE

DT=B.rn32(F(1,200))


def cfgs():
    c=QBASE.shipping_runtime().candidate_cfg
    d=D.shipping_defaults(qeff_pow_result=B.rn32(1))
    return c,d


def exp_for(x):
    lo,hi,_,_=EXP.enclosure(x); return B.rn32((lo+hi)/2)


def tau_result(prev,c):
    f=FREQ.store(B.rn32(F(1,5)),B.rn32(c.min_freq),B.rn32(c.max_freq))
    _,_,_,adapt=TB._floats_from_frequency(f.stored_hz,c,DT)
    e=exp_for(B.div(DT,adapt))
    return TL.step_tracks(prev,separate_frequency=f,fma_frequency=f,cfg=c,dt=DT,
                          separate_exp_decay=e,fma_exp_decay=e)


def exact_sample():
    return C.WaveBandSample(F(1,5),True,F(1,4),0,False,0,1,F(1,2))


def exact_candidate(c):
    t=C.targets(exact_sample(),c)
    return C.CandidateResult(t.frequency,t.variance_wave,t.tau_target,t.sigma_target,F(1),
        C.TuneState(F(11,10),F(1,100),F(1,2)),True,DT,F(1),F(1))


def sigma_join(c,d):
    m=SM.target(d,var_ready=True,accel_variance=B.rn32(F(1,4)),
        band_noise_sigma=B.rn32(0),still=False,still_time=B.rn32(0),
        sqrt_result=B.rn32(F(1,2)))
    return SJ.join(exact_sample(),c,d,m)


def rs_exp(d,tau_target):
    safe=RA.clamp(tau_target,RA.TIME_MIN,RA.TIME_MAX)
    requested=B.mul(d.adapt_RS_mult,safe)
    lo=min(max(DT,RA.HORIZON_MIN),RA.HORIZON_MAX)
    sec=RA.clamp(requested,lo,RA.HORIZON_MAX)
    return exp_for(B.div(DT,sec))


def spectral_sqrt_for(tau,d):
    requested=B.mul(d.pseudo_tau_ratio,tau)
    ts=min(max(requested,d.pseudo_min),d.pseudo_max)
    # Pick the correctly rounded positive sqrt witness through exact enclosure.
    from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT
    lo,hi=ROOT.sqrt_enclosure(ts); return B.rn32((lo+hi)/2)


def spectral_pow_for(tau,sigma,d):
    sdiv=B.div(sigma,d.sigma_coeff)
    sab=max(sdiv,B.rn32(F(1,10**6)))
    tau2=B.mul(tau,tau); u=B.mul(B.mul(sab,tau2),tau2)
    from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT
    lo,hi=ROOT.pow_6_7_enclosure(u); return B.rn32((lo+hi)/2)


class Tests(unittest.TestCase):
    def test_one_candidate_advances_all_three_ledgers_once_and_sets_one_pending(self):
        c,d=cfgs(); prev=P.initial(); tr=tau_result(prev.tau,c); sj=sigma_join(c,d)
        st=tr.separate_step; ft=tr.fma_step
        out=X.step(prev,exact_candidate(c),c,d,tr,dt=DT,
            separate_sigma_join=sj,fma_sigma_join=sj,
            separate_spectral_pow=spectral_pow_for(st.tau_target,sj.machine.sigma_target,d),
            separate_spectral_sqrt=spectral_sqrt_for(st.tau_target,d),
            fma_spectral_pow=spectral_pow_for(ft.tau_target,sj.machine.sigma_target,d),
            fma_spectral_sqrt=spectral_sqrt_for(ft.tau_target,d),
            separate_rs_exp_decay=rs_exp(d,st.tau_target),
            fma_rs_exp_decay=rs_exp(d,ft.tau_target))
        self.assertEqual(out.product.state.tau.updates,1)
        self.assertEqual(out.product.state.sigma.updates,1)
        self.assertEqual(out.product.state.rs.updates,1)
        self.assertTrue(out.product.state.pending)
        self.assertIs(out.product.tau,tr)
        self.assertEqual(out.separate.common_alpha.step,tr.separate_step)
        self.assertEqual(out.fma.common_alpha.step,tr.fma_step)

    def test_sigma_join_from_different_exact_candidate_is_rejected_before_RS(self):
        c,d=cfgs(); prev=P.initial(); tr=tau_result(prev.tau,c); sj=sigma_join(c,d)
        bad=replace(exact_candidate(c),variance_wave=F(1,3))
        with self.assertRaisesRegex(ValueError,'detached from common exact candidate'):
            X.step(prev,bad,c,d,tr,dt=DT,separate_sigma_join=sj,fma_sigma_join=sj,
                separate_spectral_pow=B.rn32(1),separate_spectral_sqrt=B.rn32(1),
                fma_spectral_pow=B.rn32(1),fma_spectral_sqrt=B.rn32(1),
                separate_rs_exp_decay=B.rn32(1),fma_rs_exp_decay=B.rn32(1))

    def test_tau_result_from_other_whole_predecessor_is_rejected(self):
        c,d=cfgs(); prev=P.initial(); other=replace(prev.tau,separate=B.rn32(1))
        tr=tau_result(other,c); sj=sigma_join(c,d)
        with self.assertRaisesRegex(ValueError,'detached from whole machine TuneState predecessor'):
            X.step(prev,exact_candidate(c),c,d,tr,dt=DT,separate_sigma_join=sj,fma_sigma_join=sj,
                separate_spectral_pow=B.rn32(1),separate_spectral_sqrt=B.rn32(1),
                fma_spectral_pow=B.rn32(1),fma_spectral_sqrt=B.rn32(1),
                separate_rs_exp_decay=B.rn32(1),fma_rs_exp_decay=B.rn32(1))

    def test_readiness_stays_before_startup_and_storage_promotion(self):
        r=X.readiness()
        self.assertTrue(r['tau_sigma_RS_scalar_ledgers_advance_once_before_pending_bit_is_attached'])
        self.assertTrue(r['machine_input_and_libm_supplies_retained_not_collapsed'])
        for k in ('upstream_frontend_binary32_correspondence_closed','all_target_libm_correspondence_closed',
                  'source_uniform_machine_supply_bounds_closed','startup_frontend_machine_TuneState_product_attached',
                  'Live_600_step_machine_TuneState_product_attached','source_uniform_complete_600_step_word_qualified',
                  'storage_search_allowed','ALT_STARTUP_PASS','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
