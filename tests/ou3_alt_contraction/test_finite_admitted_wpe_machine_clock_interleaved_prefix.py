import unittest
from fractions import Fraction as F
from tools.stability.ou3_alt_contraction import finite_admitted_wpe_machine_clock_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_admitted_machine_clock_qualified_interleaved_prefix as CLOCK
from tools.stability.ou3_alt_contraction import finite_admitted_machine_joined_frontend_racc_interleaved_prefix as JOIN
from tools.stability.ou3_alt_contraction import finite_wpe_machine_binary32 as WM
from tools.stability.ou3_alt_contraction import finite_wpe_moment_binary32 as M
import test_finite_admitted_machine_joined_frontend_racc_interleaved_prefix as BASE
import test_finite_live_magnetic_word as MAG


def expw(arg):
    lo,hi=M.EXP.exp_minus_enclosure(arg); return M.B.rn32((lo+hi)/2)


def witness(wpe,vertical,dt):
    cfg=wpe.cfg; mom=wpe.separate; h=M.B.rn32(dt)
    d=expw(M.B.mul(cfg.lambda_,h)); gain=M.B.div(M.B.sub(M.ONE,d),cfg.lambda_)
    hp1=M.B.mul(d,M.B.sub(M.B.add(mom.hp1,vertical),mom.accel_prev))
    hp2=M.B.mul(d,M.B.sub(M.B.add(mom.hp2,hp1),mom.hp1_prev))
    vel=M._sum_products(d,mom.velocity,gain,hp2)[0]; elev=M._sum_products(d,mom.elevation,gain,vel)[0]
    period=WM.HorizonPeriodWitness(wpe.logs.separate.log_period,M.B.rn32(1))
    horizon=min(max(M.B.mul(cfg.moment_horizon_periods,period.period_exp),cfg.min_horizon_sec),cfg.max_horizon_sec)
    md=expw(M.B.div(h,horizon)); a=M.B.sub(M.ONE,md)
    ms=M.MomentSuccessors(M._ema(mom.weight,M.ONE,a)[0],M._ema(mom.velocity_mean,vel,a)[0],
        M._second_moment(mom.velocity_sq,vel,a)[0],M._ema(mom.elevation_mean,elev,a)[0],
        M._second_moment(mom.elevation_sq,elev,a)[0])
    vm=M.B.div(ms.velocity_mean,ms.weight); em=M.B.div(ms.elevation_mean,ms.weight)
    vs=M.B.div(ms.velocity_sq,ms.weight); es=M.B.div(ms.elevation_sq,ms.weight)
    vv=max(M.ZERO,M.B.sub(vs,M.B.mul(vm,vm))); ev=max(M.ZERO,M.B.sub(es,M.B.mul(em,em)))
    return WM.ModeWitnesses(dict(decay_exp=d,velocity_successor=vel,elevation_successor=elev,
        moment_decay_exp=md,moment_successors=ms,velocity_var_successor=vv,elevation_var_successor=ev),horizon=period)


def general_witness(wpe,vertical,dt,mode):
    """A component arithmetic witness, not a startup-reachability fixture."""
    cfg=wpe.cfg; mom=getattr(wpe,mode); track=getattr(wpe.logs,mode); h=M.B.rn32(dt)
    d=expw(M.B.mul(cfg.lambda_,h)); gain=M.B.div(M.B.sub(M.ONE,d),cfg.lambda_)
    hp1=M.B.mul(d,M.B.sub(M.B.add(mom.hp1,vertical),mom.accel_prev))
    hp2=M.B.mul(d,M.B.sub(M.B.add(mom.hp2,hp1),mom.hp1_prev))
    vel=M._sum_products(d,mom.velocity,gain,hp2)[0]
    elev=M._sum_products(d,mom.elevation,gain,vel)[0]
    kw=dict(decay_exp=d,velocity_successor=vel,elevation_successor=elev)
    if M.B.add(mom.elapsed,h)<M.B.div(3,cfg.lambda_): return WM.ModeWitnesses(kw)
    period=None if track.log_period is None else WM.HorizonPeriodWitness(track.log_period,M.B.rn32(1))
    sea=M.SIX if period is None else period.period_exp
    horizon=min(max(M.B.mul(cfg.moment_horizon_periods,sea),cfg.min_horizon_sec),cfg.max_horizon_sec)
    md=expw(M.B.div(h,horizon)); a=M.B.sub(M.ONE,md)
    ms=M.MomentSuccessors(M._ema(mom.weight,M.ONE,a)[0],M._ema(mom.velocity_mean,vel,a)[0],
        M._second_moment(mom.velocity_sq,vel,a)[0],M._ema(mom.elevation_mean,elev,a)[0],
        M._second_moment(mom.elevation_sq,elev,a)[0])
    kw.update(moment_decay_exp=md,moment_successors=ms)
    if ms.weight<=M.WEIGHT_GATE: return WM.ModeWitnesses(kw,horizon=period)
    vm=M.B.div(ms.velocity_mean,ms.weight); em=M.B.div(ms.elevation_mean,ms.weight)
    vv=max(M.ZERO,M.B.sub(M.B.div(ms.velocity_sq,ms.weight),M.B.mul(vm,vm)))
    ev=max(M.ZERO,M.B.sub(M.B.div(ms.elevation_sq,ms.weight),M.B.mul(em,em)))
    kw.update(velocity_var_successor=vv,elevation_var_successor=ev)
    if vv<=M.VAR_GATE or ev<=M.VAR_GATE: return WM.ModeWitnesses(kw,horizon=period)
    omega=M.B.sub(M.B.div(vv,ev),M.B.mul(cfg.lambda_,cfg.lambda_))
    if omega<=M.OMEGA_GATE: return WM.ModeWitnesses(kw,horizon=period)
    root=M.SQRT.sqrt32(omega); kw['sqrt_omega']=root
    raw=M.B.div(M.TWO_PI,root); lr=M.B.rn32(2)
    if period is None: log=WM.LOG.InitWitness(lr)
    else:
        lh=min(max(M.B.mul(WM.LOG.LOG_SMOOTH_PERIODS,sea),WM.LOG.HORIZON_MIN),WM.LOG.HORIZON_MAX)
        log=WM.LOG.SmoothWitness(lr,sea,expw(M.B.div(h,lh)))
    if getattr(wpe,mode+'_usable'): usable=None
    else:
        nextlog=WM.LOG.advance_modes(wpe.logs,separate_produced=mode=='separate',fma_produced=mode=='fma',
            separate_witness=log if mode=='separate' else None,fma_witness=log if mode=='fma' else None)
        usable=WM.USABLE.PeriodWitness(getattr(nextlog.state,mode).log_period,M.B.rn32(1))
    return WM.ModeWitnesses(kw,horizon=period,raw_log=WM.RawLogBinding(raw,lr),log=log,usable=usable)


class Tests(unittest.TestCase):
    def test_one_live_edge_advances_same_full_wpe_history(self):
        base,kw,join,_,_,_=BASE.fixture(); clock=CLOCK.begin(base)
        mt=JOIN._mtune_state(base.base)
        cfg=WM.MOM.Config.from_shadow(JOIN._runtime(base.base).wpe_cfg)
        mom=M.State(elapsed=M.B.rn32(3),samples=0)
        wpe=WM.State(cfg,mom,mom,mt.base.wpe,separate_usable=True,fma_usable=True)
        s=X.begin(clock,wpe)
        probe=JOIN.imu_step(base,separate_racc_accel_ldlt=MAG.REJECT,fma_racc_accel_ldlt=MAG.REJECT,**join,**kw)
        sw=witness(wpe,probe.separate_source.band_input,join['machine_dt'])
        fw=WM.ModeWitnesses(dict(sw.moment),horizon=WM.HorizonPeriodWitness(wpe.logs.fma.log_period,M.B.rn32(1)))
        out=X.imu_step(s,separate_wpe=sw,fma_wpe=fw,separate_racc_accel_ldlt=MAG.REJECT,fma_racc_accel_ldlt=MAG.REJECT,**join,**kw)
        self.assertEqual(out.state.wpe_steps,1)
        self.assertEqual(out.state.wpe.logs,X._logs(out.state.base))
        self.assertEqual(out.separate_wpe.vertical_accel,out.lower.lower.separate_source.band_input)

    def test_distinct_production_and_takeover_flow_through_full_Live_tuner_Racc_join(self):
        from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as WF
        from tools.stability.ou3_alt_contraction import finite_admitted_tau_interleaved_prefix as TJ
        from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as TL
        import test_finite_admitted_tau_interleaved_prefix as TB
        base,kw,join,_,sep,fma=BASE.fixture(); mt=JOIN._mtune_state(base.base)
        runtime=JOIN._runtime(base.base); cfg=WM.MOM.Config.from_shadow(runtime.wpe_cfg)
        # Component predecessor: exact is usable, only the separate machine is.
        sm=M.State(elapsed=M.B.rn32(4),weight=M.ONE,velocity_sq=M.B.rn32(9),elevation_sq=M.ONE)
        fm=M.State(elapsed=M.B.rn32(4))
        wpe=WM.State(cfg,sm,fm,mt.base.wpe,separate_usable=True,fma_usable=False)
        state=X.begin(CLOCK.begin(base),wpe)
        sf,ff=WF.machine_frequencies(TJ._entry_wpe(mt.base.base),wpe,logs=mt.base.wpe,
            separate_getter=kw['separate_getter'],fma_getter=kw['fma_getter'],shadow_frequency=F(1,2),
            stats_cfg=runtime.stats_cfg,exact_min_hz=runtime.candidate_cfg.min_freq,exact_max_hz=runtime.candidate_cfg.max_freq)
        old_fma_front=kw['fma_frontend']
        es=TB.machine_decay(sf.stored.stored_hz,runtime.candidate_cfg)
        ef=TB.machine_decay(ff.stored.stored_hz,runtime.candidate_cfg)
        tau=TL.step_tracks(mt.base.base.tau,separate_frequency=sf.stored,fma_frequency=ff.stored,
            cfg=runtime.candidate_cfg,dt=TL.DT,separate_exp_decay=es,fma_exp_decay=ef)
        kw.update(separate_tau_exp_decay=es,fma_tau_exp_decay=ef)
        for mode,src,freq,t in (('separate',sep,sf,tau.separate_step),('fma',fma,ff,tau.fma_step)):
            front=BASE.VBASE.MFBASE.source_step(getattr(mt.frontends,mode),band_cfg=runtime.band_cfg,
                stats_cfg=runtime.stats_cfg,frequency=freq.external.stored.input_hz,
                bench_noise_sigma=runtime.bench_noise_sigma,x=src.band_input,last=mode=='fma')
            sigma=BASE.VBASE.sigma_target(mt,front,src)
            kw[mode+'_frontend']=front; kw[mode+'_sigma_machine']=sigma
            kw[mode+'_spectral_pow']=BASE.VBASE.MBASE.spectral_pow_for(t.tau_target,sigma.sigma_target,mt.deployment_cfg)
            kw[mode+'_spectral_sqrt']=BASE.VBASE.MBASE.spectral_sqrt_for(t.tau_target,mt.deployment_cfg)
            kw[mode+'_rs_exp_decay']=BASE.VBASE.MBASE.rs_exp(mt.deployment_cfg,t.tau_target)
        sw=general_witness(wpe,sep.band_input,join['machine_dt'],'separate')
        fw=general_witness(wpe,fma.band_input,join['machine_dt'],'fma')
        args=dict(separate_wpe=sw,fma_wpe=fw,separate_racc_accel_ldlt=MAG.REJECT,
                  fma_racc_accel_ldlt=MAG.REJECT,**join,**kw)
        out=X.imu_step(state,**args)
        mtout=JOIN.LOWER._mtune_result(out.lower.lower.lower.lower)
        self.assertEqual((out.separate_wpe.branch,out.fma_wpe.branch),('valid-period','degenerate-moments'))
        self.assertEqual((mtout.lower.separate_frequency.external.branch,mtout.lower.fma_frequency.external.branch),('wpe','prior'))
        self.assertEqual(mtout.lower.fma_frequency.external.stored.input_hz,WF.PRIOR)
        self.assertNotEqual(mtout.lower.fma_supply.frequency,0)
        self.assertIs(out.state.wpe.logs.fma,wpe.logs.fma)
        self.assertNotEqual(out.state.wpe.logs.separate,wpe.logs.separate)
        self.assertEqual(out.state.wpe_steps,1)
        # The next event retains every machine branch state and the same source.
        self.assertEqual(out.state.wpe.logs,X._logs(out.state.base))
        with self.assertRaisesRegex(ValueError,'low-corner exp witness detached'):
            X.imu_step(state,**dict(args,fma_frontend=old_fma_front))
        with self.assertRaisesRegex(TypeError,'owned by the full WPE'):
            X.imu_step(state,**dict(args,machine_wpe_entry=wpe))

    def test_complete_cannot_skip_wpe_edges(self):
        base,_,_,_,_,_=BASE.fixture(); clock=CLOCK.begin(base); mt=JOIN._mtune_state(base.base)
        cfg=WM.MOM.Config.from_shadow(JOIN._runtime(base.base).wpe_cfg); mom=M.State(elapsed=M.B.rn32(3),samples=0)
        s=X.begin(clock,WM.State(cfg,mom,mom,mt.base.wpe,separate_usable=True,fma_usable=True))
        with self.assertRaises((ValueError,TypeError)): X.complete(s)

    def test_readiness_attaches_all_600_edges_and_MEMS_WPE_supplies_not_complete_word(self):
        r=X.readiness()
        self.assertTrue(r['Live_600_step_WPE_machine_history_attached'])
        self.assertTrue(r['complete_word_requires_full_WPE_machine_step_on_all_600_IMU_edges'])
        self.assertFalse(r['machine_vs_exact_WPE_period_branch_agreement_required'])
        self.assertTrue(r['independent_machine_WPE_production_and_frequency_branches_composed'])
        self.assertTrue(r['source_uniform_WPE_machine_supply_bounds_closed'])
        self.assertTrue(r['WPE_supply_bound_requires_no_startup_deadline'])
        self.assertTrue(r['target_exp_log_error_bounds_fit_WPE_supply_profile'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])

if __name__=='__main__': unittest.main()
