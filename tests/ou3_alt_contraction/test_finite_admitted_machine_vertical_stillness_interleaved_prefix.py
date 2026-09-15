"""Admitted machine guard -> Mahony -> LPF/stillness ancestry regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest
from tools.stability.ou3_alt_contraction import finite_admitted_machine_vertical_stillness_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_machine_vertical_stillness_source as VS
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_binary32_mahony_startup as MAH
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
from tools.stability.ou3_alt_contraction import finite_stillness_sigma_binary32 as STILL
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as SM
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT
import test_finite_admitted_machine_tunestate_interleaved_prefix as TBASE
import test_finite_machine_frontend_sigma_source as MFBASE
import test_finite_tuner_machine_candidate_step as MBASE
import test_finite_source_bound_live_word as LBASE

def expw(x):
    lo,hi=SM.exp_minus_enclosure(F(x)); return B.rn32((lo+hi)/2)
def sqrtw(x):
    lo,hi=ROOT.sqrt_enclosure(F(x)); return B.rn32((lo+hi)/2)
def source_state(mt):
    n=mt.frontends.samples
    return VS.begin(V.State(initialized=True),0,n>0,STILL.State(samples=n),cutoff_hz=B.rn32(6),samples=n)
def machine_api(live,h):
    return dict(machine_dt=B.rn32(h),machine_gyro_body=tuple(B.rn32(F(x)) for x in live.guarded.raw_gyro_body),machine_acc_body=tuple(B.rn32(F(x)) for x in live.guarded.raw_accel_body))
def source_witness(source,guarded,runtime,h,*,last=False):
    hq=B.rn32(h); vcfg=X._machine_vertical_cfg(runtime)
    mah=MAH.step(source.vertical,vcfg,dt=hq,gyro=guarded.raw_gyro,acc=guarded.conditioned_acc)
    x=B.rn32(mah.vertical.vertical_accel)
    if source.lpf.initialized:
        mag=B.mul(B.mul(B.mul(VS.TWO,VS.PI_F),source.lpf.cutoff_hz),hq); alpha=expw(mag)
        om=B.sub(VS.ONE,alpha); ax=B.mul(om,x); ap=B.mul(alpha,source.lpf.value)
        vals=tuple(sorted(set((B.add(ax,ap),B.fma(om,x,ap),B.fma(alpha,source.lpf.value,ax)))))
        lpf_successor=vals[-1] if last else vals[0]
    else: alpha=None; lpf_successor=None
    lp=x if not source.lpf.initialized else lpf_successor
    cfg=runtime.still_cfg; g=B.rn32(cfg.gravity); a=B.rn32(cfg.energy_alpha)
    an=B.div(lp,g); inst=B.mul(an,an); decay=B.sub(B.rn32(1),a)
    energies=STILL._sum_products(decay,source.stillness.energy,a,inst); en=energies[-1] if last else energies[0]
    is_still=en < B.rn32(cfg.energy_thresh); st=min(B.add(source.stillness.still_time,hq),B.rn32(60)) if is_still else B.rn32(0)
    atten=expw(st) if is_still else None
    lps=VS.lpf_step(source.lpf,x=x,dt=hq,alpha_exp=alpha,successor=lpf_successor)
    sts=STILL.step(source.stillness,cfg,vertical_lp=lps.state.value,dt=hq,energy_successor=en,attenuation_exp=atten)
    out=VS.Result(source,VS.State(mah.vertical.state,lps.state,sts.state,source.samples+1),mah,lps,sts,x)
    return out,dict(lpf_alpha_exp=alpha,lpf_successor=lpf_successor,still_energy_successor=en,still_attenuation_exp=atten)
def sigma_target(mt,front,src):
    vn=B.mul(front.band_noise_sigma,front.band_noise_sigma); vt=max(B.rn32(0),front.accel_variance) if front.state.stats.var_ready else vn
    pre=max(B.rn32(0),B.sub(vt,vn)); atten=src.stillness.attenuation if src.stillness.state.is_still else B.rn32(1)
    vw=max(B.mul(pre,atten) if src.stillness.state.is_still else pre,SM.VAR_FLOOR)
    return SM.target(mt.deployment_cfg,var_ready=front.state.stats.var_ready,accel_variance=front.accel_variance,band_noise_sigma=front.band_noise_sigma,still=src.stillness.state.is_still,still_time=src.stillness.state.still_time,still_exp_result=src.stillness.exp_result,sqrt_result=sqrtw(vw))
def executed_fixture():
    mt=TBASE.state(usable=True); kwargs,_=TBASE.event_operands(mt); lower0=TBASE.LOWER.imu_step(mt.base,**kwargs); live=TBASE.X._live_result(lower0)
    runtime=mt.base.base.prefix.prefix.live.live_word.runtime; h=kwargs['restricted'].segment.h; api=machine_api(live,h)
    guard0=X.GUARD.State(); gcfg=X.GUARD.Config(); guarded=X.GUARD.step(guard0,gcfg,raw_gyro=api['machine_gyro_body'],raw_acc=api['machine_acc_body'],dt=api['machine_dt'])
    sep0=source_state(mt); fma0=source_state(mt); sep_src,sw=source_witness(sep0,guarded,runtime,h,last=False); fma_src,fw=source_witness(fma0,guarded,runtime,h,last=True); assert sep_src.mahony==fma_src.mahony
    sep_front=MFBASE.source_step(mt.frontends.separate,band_cfg=runtime.band_cfg,stats_cfg=runtime.stats_cfg,frequency=lower0.separate_frequency.external.stored.input_hz,bench_noise_sigma=runtime.bench_noise_sigma,x=sep_src.band_input,last=False)
    fma_front=MFBASE.source_step(mt.frontends.fma,band_cfg=runtime.band_cfg,stats_cfg=runtime.stats_cfg,frequency=lower0.fma_frequency.external.stored.input_hz,bench_noise_sigma=runtime.bench_noise_sigma,x=fma_src.band_input,last=True)
    sep_sigma=sigma_target(mt,sep_front,sep_src); fma_sigma=sigma_target(mt,fma_front,fma_src); machine=dict(separate_frontend=sep_front,fma_frontend=fma_front,separate_sigma_machine=sep_sigma,fma_sigma_machine=fma_sigma)
    for mode,step,sm in (('separate',lower0.tau_step.separate_step,sep_sigma),('fma',lower0.tau_step.fma_step,fma_sigma)):
        machine[mode+'_spectral_pow']=MBASE.spectral_pow_for(step.tau_target,sm.sigma_target,mt.deployment_cfg); machine[mode+'_spectral_sqrt']=MBASE.spectral_sqrt_for(step.tau_target,mt.deployment_cfg); machine[mode+'_rs_exp_decay']=MBASE.rs_exp(mt.deployment_cfg,step.tau_target)
    state=X.begin(mt,guard=guard0,guard_cfg=gcfg,separate_source=sep0,fma_source=fma0)
    witness=dict(api,separate_lpf_alpha_exp=sw['lpf_alpha_exp'],separate_lpf_successor=sw['lpf_successor'],separate_still_energy_successor=sw['still_energy_successor'],separate_still_attenuation_exp=sw['still_attenuation_exp'],fma_lpf_alpha_exp=fw['lpf_alpha_exp'],fma_lpf_successor=fw['lpf_successor'],fma_still_energy_successor=fw['still_energy_successor'],fma_still_attenuation_exp=fw['still_attenuation_exp'])
    return state,kwargs,machine,witness,guarded,sep_src,fma_src

class Tests(unittest.TestCase):
    def test_executed_event_binds_float_API_guard_band_and_sigma_on_one_history(self):
        s,kwargs,machine,witness,guarded,sep_src,fma_src=executed_fixture(); out=X.imu_step(s,**witness,**machine,**kwargs)
        self.assertEqual(out.guard,guarded); self.assertEqual(out.separate_source,sep_src); self.assertEqual(out.fma_source,fma_src); self.assertEqual(out.state.source_steps,1)
        self.assertEqual(out.lower.separate_frontend.band.envelope.x,out.separate_source.band_input); self.assertEqual(out.lower.separate_sigma_join.machine.still_time,out.separate_source.stillness.state.still_time)
    def test_API_float_splice_is_rejected_before_machine_guard(self):
        s,kwargs,machine,witness,_,_,_=executed_fixture(); witness=dict(witness); bad=list(witness['machine_acc_body']); bad[0]=B.rn32(F(bad[0])+F(1,100)); witness['machine_acc_body']=tuple(bad)
        with self.assertRaisesRegex(ValueError,'API rounding'): X.imu_step(s,**witness,**machine,**kwargs)
    def test_state_rejects_private_Mahony_divergence(self):
        mt=TBASE.state(usable=True); a=source_state(mt); b=replace(source_state(mt),vertical=replace(a.vertical,elapsed=F(1,200)))
        with self.assertRaisesRegex(ValueError,'private-Mahony state'): X.begin(mt,guard=X.GUARD.State(),guard_cfg=X.GUARD.Config(),separate_source=a,fma_source=b)
    def test_MAG_and_HOLD_preserve_common_guard_and_sources(self):
        mt=TBASE.state(usable=True); a=source_state(mt); g=X.GUARD.State(); s=X.begin(mt,guard=g,guard_cfg=X.GUARD.Config(),separate_source=a,fma_source=a); word=mt.base.base.prefix.prefix.live.live_word
        s,_=X.mag_step(s,**LBASE.mag_kwargs(word)); self.assertIs(s.guard,g); s,_=X.set_hold(s,hold=False); self.assertIs(s.guard,g)
    def test_complete_cannot_skip_600_machine_source_edges(self):
        mt=TBASE.state(usable=True); a=source_state(mt)
        with self.assertRaises((ValueError,TypeError)): X.complete(X.begin(mt,guard=X.GUARD.State(),guard_cfg=X.GUARD.Config(),separate_source=a,fma_source=a))
    def test_readiness_closes_composition_but_not_config_startup_or_native(self):
        r=X.readiness()
        for k in ('admitted_machine_TuneState_word_consumed',
                  'one_common_binary32_guard_feeds_private_Mahony_and_band_history',
                  'shipping_float_API_dt_gyro_accel_projection_is_mandatory',
                  'same_machine_private_Mahony_vertical_feeds_band_and_tracker_LPF',
                  'machine_band_input_bound_to_same_frontend_consumed_by_sigma_target',
                  'machine_tracker_LPF_and_stillness_bound_to_same_sigma_target',
                  'complete_word_requires_machine_source_on_all_600_IMU_edges'):
            self.assertTrue(r[k])
        for k in ('guard_runtime_config_ancestry_closed',
                  'private_Mahony_runtime_config_and_startup_ancestry_closed',
                  'tracker_LPF_runtime_config_ancestry_closed',
                  'guard_exp_sqrt_compiler_correspondence_closed',
                  'private_Mahony_compiler_profile_qualified',
                  'machine_guard_displacement_already_injected_into_Racc',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed',
                  'ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[k])
if __name__=='__main__': unittest.main()
