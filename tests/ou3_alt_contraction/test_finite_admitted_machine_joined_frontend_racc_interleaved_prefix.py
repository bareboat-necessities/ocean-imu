"""Same-event guard/Mahony/frontend/sigma -> Racc measurement join regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_machine_joined_frontend_racc_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
import test_finite_admitted_machine_racc_supply_interleaved_prefix as RBASE
import test_finite_admitted_machine_vertical_stillness_interleaved_prefix as VBASE
import test_finite_admitted_machine_tunestate_interleaved_prefix as TBASE
import test_finite_live_magnetic_word as MAG
import test_finite_source_bound_live_word as LBASE


def fixture():
    r=RBASE.state(); mt=X._mtune_state(r); runtime=X._runtime(r)
    # Reproduce only the arithmetic witnesses for the already-defined same MTUNE
    # event so the joined source can be checked against the nested event result.
    tkw,_=TBASE.event_operands(mt); lower0=TBASE.LOWER.imu_step(mt.base,**tkw)
    live=TBASE.X._live_result(lower0); h=tkw['restricted'].segment.h
    api=VBASE.machine_api(live,h)
    gcfg=X.CONFIG._machine_guard_cfg(runtime); guard0=X.GUARD.State()
    guarded=X.GUARD.step(guard0,gcfg,raw_gyro=api['machine_gyro_body'],raw_acc=api['machine_acc_body'],dt=api['machine_dt'])
    sep0=VBASE.source_state(mt); fma0=VBASE.source_state(mt)
    sep,sw=VBASE.source_witness(sep0,guarded,runtime,h,last=False)
    fma,fw=VBASE.source_witness(fma0,guarded,runtime,h,last=True)
    assert sep.mahony==fma.mahony
    s=X.begin(r,guard=guard0,guard_cfg=gcfg,separate_source=sep0,fma_source=fma0)
    join=dict(api,
        separate_lpf_alpha_exp=sw['lpf_alpha_exp'],separate_lpf_successor=sw['lpf_successor'],
        separate_still_energy_successor=sw['still_energy_successor'],separate_still_attenuation_exp=sw['still_attenuation_exp'],
        fma_lpf_alpha_exp=fw['lpf_alpha_exp'],fma_lpf_successor=fw['lpf_successor'],
        fma_still_energy_successor=fw['still_energy_successor'],fma_still_attenuation_exp=fw['still_attenuation_exp'])
    # Feed the real same-event Mahony output into the lower band/sigma graph;
    # the old fixture supplied an unrelated zero vertical input here.
    kw=RBASE.event_operands(r)
    for mode,src,freq,tau in (('separate',sep,lower0.separate_frequency,lower0.tau_step.separate_step),
                             ('fma',fma,lower0.fma_frequency,lower0.tau_step.fma_step)):
        front=VBASE.MFBASE.source_step(getattr(mt.frontends,mode),band_cfg=runtime.band_cfg,
            stats_cfg=runtime.stats_cfg,frequency=freq.external.stored.input_hz,
            bench_noise_sigma=runtime.bench_noise_sigma,x=src.band_input,last=mode=='fma')
        sigma=VBASE.sigma_target(mt,front,src)
        kw[mode+'_frontend']=front; kw[mode+'_sigma_machine']=sigma
        kw[mode+'_spectral_pow']=VBASE.MBASE.spectral_pow_for(tau.tau_target,sigma.sigma_target,mt.deployment_cfg)
        kw[mode+'_spectral_sqrt']=VBASE.MBASE.spectral_sqrt_for(tau.tau_target,mt.deployment_cfg)
        kw[mode+'_rs_exp_decay']=VBASE.MBASE.rs_exp(mt.deployment_cfg,tau.tau_target)
    return s,kw,join,guarded,sep,fma


class Tests(unittest.TestCase):
    def test_same_Racc_event_consumes_frontend_and_sigma_from_joined_machine_source(self):
        s,kw,join,guarded,sep,fma=fixture()
        out=X.imu_step(s,separate_racc_accel_ldlt=MAG.REJECT,fma_racc_accel_ldlt=MAG.REJECT,**join,**kw)
        mt=X.LOWER._mtune_result(out.lower.lower)
        self.assertEqual(out.guard,guarded)
        self.assertEqual(out.separate_source,sep); self.assertEqual(out.fma_source,fma)
        self.assertEqual(mt.separate_frontend.band.envelope.x,out.separate_source.band_input)
        self.assertEqual(mt.fma_frontend.band.envelope.x,out.fma_source.band_input)
        self.assertEqual(mt.separate_sigma_join.machine.still_time,out.separate_source.stillness.state.still_time)
        self.assertEqual(out.state.source_steps,1)
        self.assertEqual(out.state.base.racc_steps,s.base.racc_steps+1)

    def test_frontend_input_splice_is_rejected_against_nested_Racc_TuneState_event(self):
        s,kw,join,_,_,_=fixture()
        # The lower machine frontend is fixed by its own arithmetic witnesses;
        # changing only the joined Mahony API source must fail at the same-event
        # frontend equality rather than creating a second accepted history.
        bad=dict(join); a=list(bad['machine_acc_body']); a[2]=B.rn32(F(a[2])+F(1,100)); bad['machine_acc_body']=tuple(a)
        with self.assertRaises((ValueError,TypeError)):
            X.imu_step(s,separate_racc_accel_ldlt=MAG.REJECT,fma_racc_accel_ldlt=MAG.REJECT,**bad,**kw)

    def test_detached_guard_runtime_config_is_rejected_before_join(self):
        s,_,_,_,_,_=fixture()
        with self.assertRaisesRegex(ValueError,'guard configuration detached'):
            replace(s,guard_cfg=replace(s.guard_cfg,cutoff_hz=B.rn32(13)))

    def test_MAG_and_HOLD_preserve_joined_machine_source_history(self):
        s,_,_,_,_,_=fixture(); g=s.guard; a=s.separate_source; n=s.source_steps
        word=X._runtime(s.base)  # establish runtime path before async event
        pre=X._mtune_state(s.base).base.base.prefix.prefix.live.live_word
        s,_=X.mag_step(s,**LBASE.mag_kwargs(pre)); self.assertIs(s.guard,g); self.assertIs(s.separate_source,a); self.assertEqual(s.source_steps,n)
        s,_=X.set_hold(s,hold=False); self.assertIs(s.guard,g); self.assertIs(s.separate_source,a); self.assertEqual(s.source_steps,n)

    def test_complete_cannot_skip_join_on_600_measurement_edges(self):
        s,_,_,_,_,_=fixture()
        with self.assertRaises((ValueError,TypeError)):
            X.complete(s)

    def test_readiness_closes_parallel_product_gap_only(self):
        r=X.readiness()
        self.assertTrue(r['startup_joined_machine_history_attached'])
        self.assertTrue(r['startup_attachment_is_conditional_not_universal_reachability'])
        for k in ('same_executed_TuneState_event_supplies_frontend_and_sigma_join',
                  'no_second_physical_or_filter_event_executed_for_frontend_ancestry',
                  'common_machine_guard_and_private_Mahony_history_joined_to_Racc_word',
                  'machine_band_input_bound_to_frontend_consumed_by_Racc_word',
                  'machine_sigma_stillness_bound_to_sigma_target_consumed_by_Racc_word',
                  'machine_guard_runtime_config_bound_in_joined_word',
                  'tracker_LPF_default_cutoff_bound_in_joined_word',
                  'complete_word_requires_join_on_all_600_Racc_measurement_IMU_edges'):
            self.assertTrue(r[k])
        for k in ('tracker_LPF_mutable_setter_ancestry_closed','private_Mahony_mutable_config_setter_ancestry_closed',
                  'target_libm_Eigen_and_compiler_profile_correspondence_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
