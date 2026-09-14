"""Admitted Live whole-machine TuneState interleaver regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_machine_tunestate_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_admitted_wpe_tau_interleaved_prefix as LOWER
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as PRODUCT
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_deployment_ledger as SIG
from tools.stability.ou3_alt_contraction import finite_tuner_rs_deployment_ledger as RS
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
import test_finite_admitted_wpe_tau_interleaved_prefix as BASE
import test_finite_source_bound_live_word as LBASE
import test_finite_admitted_source_imu_word as IBASE
import test_finite_source_bound_prediction_word as PBASE
import test_finite_tuner_machine_candidate_step as MBASE
from tools.stability.ou3_alt_contraction import finite_machine_frontend_sigma_source as MF
from test_finite_machine_frontend_sigma_source import source_step, ready_pair
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as CAND
from tools.stability.ou3_alt_contraction import finite_tuner_projection_bridge as PROJECTION
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as SM
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND


def dcfg(): return D.shipping_defaults(qeff_pow_result=B.rn32(1))

def state(*,usable=False,pending=False):
    b0=BASE.base_state(usable=usable)
    if pending:
        # Isolate a carried pending boundary; this is not startup reachability.
        p=b0.prefix; word=p.prefix.live.live_word; imu=word.live.live.live
        imu=replace(imu,tuner=replace(imu.tuner,pending=True))
        word=replace(word,live=replace(word.live,live=replace(word.live.live,live=imu)))
        admitted=replace(p.prefix.live,live_word=word)
        b0=replace(b0,prefix=replace(p,prefix=replace(p.prefix,live=admitted)))
    b=LOWER.begin(b0,BASE.machine_wpe_for(b0))
    n=b.base.tau.updates
    exact=X._entry_live(b)
    m=PRODUCT.State(b.base.tau,SIG.State(updates=n),RS.State(updates=n),exact.tuner.pending)
    # The inherited rational-root component fixture is not the full shipping
    # source. Compile its carried sigma scale and maximum consistently.
    c=b.base.prefix.prefix.live.live_word.runtime.candidate_cfg
    d=replace(dcfg(),sigma_coeff=B.rn32(c.sigma_coeff),max_sigma=B.rn32(c.max_sigma))
    return X.begin(b,m,d,separate_active=exact.active,fma_active=exact.active,frontends=MF.initial_pair())


def event_operands(s):
    """Conditional rational-root cell; no claim of universal libm qualification."""
    witness,_,raw,r,br,dynamic=IBASE.operands(s.base.base.prefix.prefix.live)
    dynamic=dict(dynamic,preupdate_period=F(2),preupdate_frequency=F(1,2),
                 spectral=CAND.SpectralWitness(1,1))
    if s.machine.pending: dynamic['boundary_noise_sqrt']=BAND.NoiseSqrtWitness(0)
    sg,fg=BASE.half_getters(s.base)
    c=s.base.base.prefix.prefix.live.live_word.runtime.candidate_cfg
    e=BASE.TBASE.machine_decay(F(1,2),c)
    kwargs=dict(separate_getter=sg,fma_getter=fg,shadow_frequency=F(1,2),
        separate_tau_exp_decay=e,fma_tau_exp_decay=e,
        restricted=r,bias_restricted=br,witness=witness,raw=raw,
        packet_id='imu-whole-machine',**PBASE.root_args(),**dynamic)
    # Derive sigma operands from this same event's projected frontend, not a
    # second free candidate. The explicit machine/shadow supplies remain open.
    lower=LOWER.imu_step(s.base,**kwargs); suffix=X._live_result(lower).tuner_suffix
    sample=PROJECTION.sample_from_projection(suffix.band,suffix.stillness,
                                             sigma_wave_sqrt=kwargs['sigma_wave_sqrt'])
    assert sample.accel_variance==1 and sample.band_noise_sigma==0 and sample.still
    runtime=s.base.base.prefix.prefix.live.live_word.runtime
    st=B.rn32(sample.still_time); lo,hi=SM.exp_minus_enclosure(st)
    atten=B.rn32((lo+hi)/2)
    machine={}
    for mode,fr in (('separate',lower.separate_frequency),('fma',lower.fma_frequency)):
        source=source_step(getattr(s.frontends,mode),band_cfg=runtime.band_cfg,stats_cfg=runtime.stats_cfg,
            frequency=fr.external.stored.input_hz,bench_noise_sigma=runtime.bench_noise_sigma,last=(mode=='fma'))
        vn=B.mul(source.band_noise_sigma,source.band_noise_sigma)
        pre=max(F(0),B.sub(source.accel_variance,vn))
        var=max(SM.VAR_FLOOR,B.mul(pre,atten))
        lo,hi=ROOT.sqrt_enclosure(var)
        sm=SM.target(s.deployment_cfg,var_ready=source.state.stats.var_ready,accel_variance=source.accel_variance,
            band_noise_sigma=source.band_noise_sigma,still=True,still_time=st,
            still_exp_result=atten,sqrt_result=B.rn32((lo+hi)/2))
        machine[mode+'_sigma_machine']=sm; machine[mode+'_frontend']=source
    for mode,step in (('separate',lower.tau_step.separate_step),('fma',lower.tau_step.fma_step)):
        sm=machine[mode+'_sigma_machine']
        machine[mode+'_spectral_pow']=MBASE.spectral_pow_for(step.tau_target,sm.sigma_target,s.deployment_cfg)
        machine[mode+'_spectral_sqrt']=MBASE.spectral_sqrt_for(step.tau_target,s.deployment_cfg)
        machine[mode+'_rs_exp_decay']=MBASE.rs_exp(s.deployment_cfg,step.tau_target)
    return kwargs,machine


class Tests(unittest.TestCase):
    def test_product_rejects_tau_pending_and_config_splices(self):
        s=state()
        with self.assertRaisesRegex(ValueError,'tau ledger detached'):
            X.State(s.base,replace(s.machine,tau=replace(s.machine.tau,separate=B.rn32(1))),
                    s.deployment_cfg,s.separate_active,s.fma_active,s.live_entry_machine_updates)
        with self.assertRaisesRegex(ValueError,'pending bit detached'):
            X.State(s.base,replace(s.machine,pending=not s.machine.pending),
                    s.deployment_cfg,s.separate_active,s.fma_active,s.live_entry_machine_updates)
        bad=replace(s.deployment_cfg,sigma_coeff=B.rn32(F(4,5)))
        with self.assertRaisesRegex(ValueError,'config detached'):
            X.State(s.base,s.machine,bad,s.separate_active,s.fma_active,s.live_entry_machine_updates)

    def test_applied_parameters_are_distinct_state_from_candidate_memory(self):
        s=state()
        self.assertEqual(s.separate_active_join.supply.tau,0)
        self.assertEqual(s.fma_active_join.supply.tau,0)
        changed=replace(s.machine,sigma=replace(s.machine.sigma,separate=B.rn32(F(1,5))))
        # A candidate-memory change alone does not rewrite already-applied active
        # parameters; it would only take effect after a later pending boundary.
        altered=X.State(s.base,changed,s.deployment_cfg,s.separate_active,s.fma_active,
                        s.live_entry_machine_updates)
        self.assertEqual(altered.separate_active,s.separate_active)

    def test_MAG_and_HOLD_preserve_candidate_and_applied_machine_state(self):
        s=state(); m=s.machine; a=s.separate_active; w=s.base.wpe
        s,_=X.mag_step(s,**LBASE.mag_kwargs(s.base.base.prefix.prefix.live.live_word))
        self.assertIs(s.machine,m); self.assertIs(s.separate_active,a); self.assertIs(s.base.wpe,w)
        s,_=X.set_hold(s,hold=False)
        self.assertIs(s.machine,m); self.assertIs(s.separate_active,a); self.assertIs(s.base.wpe,w)

    def test_complete_word_cannot_be_claimed_without_600_common_machine_updates(self):
        s=state()
        with self.assertRaises((ValueError,TypeError)):
            X.complete(s)

    def test_executed_Live_event_advances_candidate_but_holds_active_without_pending(self):
        s=state(usable=True); kwargs,machine=event_operands(s)
        out=X.imu_step(s,**machine,**kwargs)
        self.assertEqual(out.state.machine.tau.updates,s.machine.tau.updates+1)
        self.assertEqual(out.state.machine.sigma.updates,out.state.machine.tau.updates)
        self.assertEqual(out.state.machine.rs.updates,out.state.machine.tau.updates)
        self.assertEqual(out.state.base.wpe.samples,s.base.wpe.samples+1)
        self.assertFalse(out.machine_boundary.consumed)
        self.assertIs(out.state.separate_active,s.separate_active)
        self.assertEqual(out.separate_sigma_join.exact_target.sigma_target,1)

    def test_executed_Live_pending_event_installs_rounded_mode_outputs(self):
        s=state(usable=True,pending=True)
        # Scalar boundary fixture with different carried machine gains; not an
        # assertion these sample-zero states arise from real startup.
        s=replace(s,frontends=ready_pair())
        prefix=s.base.base.prefix; word=prefix.prefix.live.live_word
        bench=B.rn32(F(1,10))
        word=replace(word,runtime=replace(word.runtime,boundary_bench_noise_sigma=bench,bench_noise_sigma=bench))
        admitted=replace(prefix.prefix.live,live_word=word)
        prefix=replace(prefix,prefix=replace(prefix.prefix,live=admitted))
        s=replace(s,base=replace(s.base,base=replace(s.base.base,prefix=prefix)))
        kwargs,machine=event_operands(s)
        sb,fb=B.mul(bench,11),B.mul(bench,12)
        out=X.imu_step(s,**machine,**kwargs,
            separate_boundary_noise_sqrt_gain=11,fma_boundary_noise_sqrt_gain=12)
        self.assertTrue(out.machine_boundary.consumed)
        self.assertIs(out.machine_boundary.arithmetic.before,s.machine)
        self.assertEqual(out.state.separate_active.Sigma_aw[2][2],B.mul(sb,sb))
        self.assertEqual(out.state.fma_active.Sigma_aw[2][2],B.mul(fb,fb))
        self.assertNotEqual(out.state.separate_active.Sigma_aw[2][2],sb*sb)
        self.assertNotEqual(out.state.separate_active,out.state.fma_active)
        self.assertIs(out.noise_floors[0].band,s.frontends.separate.band)
        self.assertEqual(out.state.frontends.samples,s.frontends.samples+1)
        self.assertEqual(out.state.machine.pending,X._entry_live(out.state.base).tuner.pending)

    def test_readiness_attaches_full_tuner_state_but_keeps_live_coefficients_open(self):
        r=X.readiness()
        self.assertTrue(r['whole_tau_sigma_RS_machine_TuneState_carried_in_same_Live_product'])
        self.assertTrue(r['candidate_memory_and_applied_machine_parameters_carried_separately'])
        self.assertTrue(r['exact_vs_machine_applied_parameter_supplies_exposed_every_Live_state'])
        self.assertTrue(r['Live_600_step_machine_TuneState_product_attached'])
        self.assertTrue(r['MAG_and_HOLD_preserve_WPE_TuneState_and_applied_machine_parameters_by_identity'])
        self.assertFalse(r['machine_active_parameter_displacement_injected_into_Live_coefficients'])
        self.assertFalse(r['source_uniform_machine_supply_bounds_closed'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS']); self.assertFalse(r['ALT_END_TO_END_PASS'])


if __name__=='__main__': unittest.main()
