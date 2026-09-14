"""Admitted machine TuneState + vertical/LPF/stillness ancestry regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_machine_vertical_stillness_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_machine_vertical_stillness_source as VS
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_binary32_mahony as MAH
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


def source_state(mt,*,fma=False):
    n=mt.frontends.samples
    # This is a carried Live-entry fixture, not a startup-reachability claim.
    vertical=V.State(initialized=True)
    still=STILL.State(samples=n)
    return VS.begin(vertical,0,n>0,still,cutoff_hz=B.rn32(6),samples=n)


def source_witness(source,guarded,runtime,h,*,last=False):
    hq=B.rn32(h)
    mah=MAH.step_initialized(source.vertical,runtime.vertical_cfg,dt=hq,
        gyro=guarded.raw_gyro_body,acc=guarded.conditioned_accel_body)
    x=B.rn32(mah.vertical.vertical_accel)
    if source.lpf.initialized:
        mag=B.mul(B.mul(B.mul(VS.TWO,VS.PI_F),source.lpf.cutoff_hz),hq)
        alpha=expw(mag); om=B.sub(VS.ONE,alpha)
        ax=B.mul(om,x); ap=B.mul(alpha,source.lpf.value)
        vals=tuple(sorted(set((B.add(ax,ap),B.fma(om,x,ap),B.fma(alpha,source.lpf.value,ax)))))
        lpf_successor=vals[-1] if last else vals[0]
    else:
        alpha=None; lpf_successor=None
        # first sample is the exact input seed
        vals=(x,)
    lp=x if not source.lpf.initialized else lpf_successor
    cfg=runtime.still_cfg; g=B.rn32(cfg.gravity); a=B.rn32(cfg.energy_alpha)
    an=B.div(lp,g); inst=B.mul(an,an); decay=B.sub(B.rn32(1),a)
    energies=STILL._sum_products(decay,source.stillness.energy,a,inst)
    en=energies[-1] if last else energies[0]
    still=en < B.rn32(cfg.energy_thresh)
    st=min(B.add(source.stillness.still_time,hq),B.rn32(60)) if still else B.rn32(0)
    atten=expw(st) if still else None
    out=VS.step(source,guarded,runtime.vertical_cfg,runtime.still_cfg,dt=hq,
        lpf_alpha_exp=alpha,lpf_successor=lpf_successor,
        still_energy_successor=en,still_attenuation_exp=atten)
    return out,dict(lpf_alpha_exp=alpha,lpf_successor=lpf_successor,
                    still_energy_successor=en,still_attenuation_exp=atten)


def sigma_target(mt,front,src):
    vn=B.mul(front.band_noise_sigma,front.band_noise_sigma)
    vt=max(B.rn32(0),front.accel_variance) if front.state.stats.var_ready else vn
    pre=max(B.rn32(0),B.sub(vt,vn))
    atten=src.stillness.attenuation if src.stillness.state.is_still else B.rn32(1)
    vw=max(B.mul(pre,atten) if src.stillness.state.is_still else pre,SM.VAR_FLOOR)
    return SM.target(mt.deployment_cfg,var_ready=front.state.stats.var_ready,
        accel_variance=front.accel_variance,band_noise_sigma=front.band_noise_sigma,
        still=src.stillness.state.is_still,still_time=src.stillness.state.still_time,
        still_exp_result=src.stillness.exp_result,sqrt_result=sqrtw(vw))


def executed_fixture():
    mt=TBASE.state(usable=True)
    kwargs,_=TBASE.event_operands(mt)
    # Reexecute only the existing lower WPE/tau relation to expose the exact
    # same-event guarded sample and machine frequency objects used by MTUNE.
    lower0=TBASE.LOWER.imu_step(mt.base,**kwargs)
    live=TBASE.X._live_result(lower0)
    runtime=mt.base.base.prefix.prefix.live.live_word.runtime
    sep0=source_state(mt); fma0=source_state(mt,fma=True)
    sep_src,sw=source_witness(sep0,live.guarded,runtime,kwargs['restricted'].segment.h,last=False)
    fma_src,fw=source_witness(fma0,live.guarded,runtime,kwargs['restricted'].segment.h,last=True)
    # Mahony must be common even if later contraction choices diverge.
    assert sep_src.mahony==fma_src.mahony
    sep_front=MFBASE.source_step(mt.frontends.separate,band_cfg=runtime.band_cfg,stats_cfg=runtime.stats_cfg,
        frequency=lower0.separate_frequency.external.stored.input_hz,
        bench_noise_sigma=runtime.bench_noise_sigma,x=sep_src.band_input,last=False)
    fma_front=MFBASE.source_step(mt.frontends.fma,band_cfg=runtime.band_cfg,stats_cfg=runtime.stats_cfg,
        frequency=lower0.fma_frequency.external.stored.input_hz,
        bench_noise_sigma=runtime.bench_noise_sigma,x=fma_src.band_input,last=True)
    sep_sigma=sigma_target(mt,sep_front,sep_src); fma_sigma=sigma_target(mt,fma_front,fma_src)
    machine=dict(separate_frontend=sep_front,fma_frontend=fma_front,
                 separate_sigma_machine=sep_sigma,fma_sigma_machine=fma_sigma)
    for mode,step,sm in (('separate',lower0.tau_step.separate_step,sep_sigma),
                         ('fma',lower0.tau_step.fma_step,fma_sigma)):
        machine[mode+'_spectral_pow']=MBASE.spectral_pow_for(step.tau_target,sm.sigma_target,mt.deployment_cfg)
        machine[mode+'_spectral_sqrt']=MBASE.spectral_sqrt_for(step.tau_target,mt.deployment_cfg)
        machine[mode+'_rs_exp_decay']=MBASE.rs_exp(mt.deployment_cfg,step.tau_target)
    state=X.begin(mt,separate_source=sep0,fma_source=fma0)
    witness=dict(separate_lpf_alpha_exp=sw['lpf_alpha_exp'],separate_lpf_successor=sw['lpf_successor'],
                 separate_still_energy_successor=sw['still_energy_successor'],separate_still_attenuation_exp=sw['still_attenuation_exp'],
                 fma_lpf_alpha_exp=fw['lpf_alpha_exp'],fma_lpf_successor=fw['lpf_successor'],
                 fma_still_energy_successor=fw['still_energy_successor'],fma_still_attenuation_exp=fw['still_attenuation_exp'])
    return state,kwargs,machine,witness,sep_src,fma_src


class Tests(unittest.TestCase):
    def test_executed_event_binds_band_and_sigma_to_same_machine_vertical_stillness_history(self):
        s,kwargs,machine,witness,sep_src,fma_src=executed_fixture()
        out=X.imu_step(s,**witness,**machine,**kwargs)
        self.assertEqual(out.state.source_steps,1)
        self.assertEqual(out.separate_source,sep_src); self.assertEqual(out.fma_source,fma_src)
        self.assertEqual(out.lower.separate_frontend.band.envelope.x,out.separate_source.band_input)
        self.assertEqual(out.lower.separate_sigma_join.machine.still,out.separate_source.stillness.state.is_still)
        self.assertEqual(out.lower.separate_sigma_join.machine.still_time,out.separate_source.stillness.state.still_time)
        self.assertEqual(out.lower.separate_sigma_join.machine.attenuation,out.separate_source.stillness.attenuation)

    def test_frontend_band_input_splice_is_rejected(self):
        s,kwargs,machine,witness,_,_=executed_fixture()
        runtime=X._runtime(s.base)
        fr=machine['separate_frontend']
        bad=MFBASE.source_step(s.base.frontends.separate,band_cfg=runtime.band_cfg,stats_cfg=runtime.stats_cfg,
            frequency=fr.stats_coefficients.frequency,bench_noise_sigma=runtime.bench_noise_sigma,
            x=B.rn32(F(1,2)),last=False)
        machine=dict(machine); machine['separate_frontend']=bad
        machine['separate_sigma_machine']=sigma_target(s.base,bad,VS.step(s.separate_source,
            TBASE.X._live_result(TBASE.LOWER.imu_step(s.base.base,**kwargs)).guarded,
            runtime.vertical_cfg,runtime.still_cfg,dt=B.rn32(kwargs['restricted'].segment.h),
            lpf_alpha_exp=witness['separate_lpf_alpha_exp'],lpf_successor=witness['separate_lpf_successor'],
            still_energy_successor=witness['separate_still_energy_successor'],still_attenuation_exp=witness['separate_still_attenuation_exp']))
        with self.assertRaisesRegex(ValueError,'band input detached'):
            X.imu_step(s,**witness,**machine,**kwargs)

    def test_state_rejects_private_Mahony_divergence(self):
        mt=TBASE.state(usable=True); a=source_state(mt); b=replace(source_state(mt),vertical=replace(a.vertical,elapsed=F(1,200)))
        with self.assertRaisesRegex(ValueError,'private-Mahony state'):
            X.begin(mt,separate_source=a,fma_source=b)

    def test_MAG_and_HOLD_preserve_machine_source_histories(self):
        mt=TBASE.state(usable=True); a=source_state(mt); s=X.begin(mt,separate_source=a,fma_source=a)
        word=mt.base.base.prefix.prefix.live.live_word
        s,_=X.mag_step(s,**LBASE.mag_kwargs(word)); self.assertIs(s.separate_source,a)
        s,_=X.set_hold(s,hold=False); self.assertIs(s.separate_source,a)

    def test_complete_cannot_skip_600_machine_source_edges(self):
        mt=TBASE.state(usable=True); a=source_state(mt)
        with self.assertRaises((ValueError,TypeError)): X.complete(X.begin(mt,separate_source=a,fma_source=a))

    def test_readiness_closes_ancestry_only(self):
        r=X.readiness()
        for k in ('persistent_separate_and_FMA_machine_vertical_stillness_histories_attached',
                  'private_Mahony_source_is_common_across_compiler_histories',
                  'same_guarded_shipping_sample_drives_machine_private_Mahony',
                  'machine_band_input_bound_to_same_private_Mahony_successor',
                  'machine_sigma_stillness_bound_to_same_tracker_LPF_successor',
                  'complete_word_requires_machine_source_on_all_600_IMU_edges'):
            self.assertTrue(r[k])
        for k in ('tracker_LPF_cutoff_runtime_ancestry_closed','machine_guard_binary32_history_attached',
                  'startup_machine_vertical_stillness_history_attached','target_libm_and_compiler_profile_correspondence_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed',
                  'ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
