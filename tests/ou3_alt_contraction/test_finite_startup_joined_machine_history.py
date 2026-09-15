"""Reset and conditional goLive relation tests, not source-uniform capture."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_startup_joined_machine_history as X
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_binary32_mahony_startup as SEED
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as G
from tools.stability.ou3_alt_contraction import finite_racc_runtime as RACC
from tools.stability.ou3_alt_contraction import finite_periodic_aw_sync as AWSYNC
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
from tools.stability.ou3_alt_contraction import deployment_scope as SCOPE
import test_finite_admitted_machine_tunestate_interleaved_prefix as TBASE
import test_finite_guarded_tuner_wpe_tau_deployment as BASE
import test_finite_guarded_tuner_prefix as RAW
import test_finite_admitted_machine_vertical_stillness_interleaved_prefix as VSBASE
import test_finite_machine_frontend_sigma_source as MFBASE
import test_finite_startup_live_machine_tunestate_bridge as GOBASE


def initial():
    mt=TBASE.state(usable=False)
    return X.initial(mt.base.base.prefix.prefix.live.live_word.runtime,mt.deployment_cfg)


def guard_witnesses(s,acc,h):
    """Deterministically choose allowed local arithmetic outcomes for a test.

    These ideal-RNE witnesses are not asserted to qualify a target toolchain.
    """
    c=s.guard_cfg; before=s.guard; G32=X.GUARD
    if not before.initialized: return {}
    def decay(f): return VSBASE.expw(B.mul(B.mul(B.mul(G32.TWO,G32.PI_F),f),h))
    alpha=decay(c.cutoff_hz); gamma=decay(c.detect_hz); be=decay(c.removed_rms_hz)
    def vec(a,x,old): return tuple(G32._lin_values(a,y,z)[0] for y,z in zip(x,old))
    lp=[]; low=acc
    for old in before.stages[:c.poles]:
        low=vec(alpha,low,old); lp.append(low)
    det=[]; hp=acc
    for old in before.detect_stages:
        d=vec(gamma,hp,old); det.append(d); hp=G32._sub(hp,d)
    beta=B.sub(B.rn32(1),be)
    ms=tuple(B.add(old,B.mul(beta,B.sub(B.mul(x,x),old))) for x,old in zip(hp,before.removed_ms))
    rms=VSBASE.sqrtw(max(F(0),G32._sum3_values(ms)[0]))
    target=F(1) if c.engage_lo<=0 else (B.div(B.sub(rms,c.engage_lo),B.sub(c.engage_hi,c.engage_lo))
        if c.engage_hi>c.engage_lo else F(rms>=c.engage_lo))
    target=min(F(1),max(F(0),target)); se=VSBASE.expw(B.div(h,c.slew_tau))
    w=B.add(before.weight,B.mul(B.sub(B.rn32(1),se),B.sub(target,before.weight)))
    wr=F(0) if w<c.weight_epsilon else (F(1) if w>B.sub(B.rn32(1),c.weight_epsilon) else w)
    out=None if wr in (0,1) else tuple(B.add(x,B.mul(wr,B.sub(y,x))) for x,y in zip(acc,low))
    return dict(lp_alpha_exp=alpha,lp_successors=lp,detect_gamma_exp=gamma,detect_successors=det,
        removed_beta_exp=be,removed_ms_successor=ms,removed_rms_sqrt=rms,
        slew_exp=se,weight_successor=w,output_successor=out)


def cold_operands(s):
    raw=RAW.packet(); h=RAW.DT; r=s.runtime
    raw=replace(raw,physical=replace(raw.physical,time=s.base.lower.frontend.tuner.time,
        centered_S=(0,0,0),velocity=(0,0,0),position=(0,0,0)))
    api=dict(machine_dt=B.rn32(h),machine_gyro_body=tuple(B.rn32(x) for x in raw.raw_gyro_body),
             machine_acc_body=tuple(B.rn32(x) for x in raw.raw_accel_body))
    gw=guard_witnesses(s,api['machine_acc_body'],api['machine_dt'])
    guard=X.GUARD.step(s.guard,s.guard_cfg,dt=api['machine_dt'],raw_gyro=api['machine_gyro_body'],raw_acc=api['machine_acc_body'],**gw)
    sep,sw=VSBASE.source_witness(s.separate_source,guard,r,h)
    fma,fw=VSBASE.source_witness(s.fma_source,guard,r,h,last=True)
    kw={k:v for k,v in BASE.cold_kwargs().items() if k not in X.STARTUP_CONFIG_KEYS}
    kw['seed']=None if s.separate_source.vertical.initialized else V.SeedWitness((1,0,0,0))
    if s.guard.initialized:
        kw.update(guard_decay=G.DecayWitness(1,1,0,0),guard_rms=G.RmsWitness(0))
    # The reset exact band uses rational witnesses; machine band uses its own
    # rigorous RNE-cell relation. These are not a universal libm certificate.
    a=dict(band_cfg=r.band_cfg,stats_cfg=r.stats_cfg,frequency=B.rn32(F(1,5)),dt=h,
           bench_noise_sigma=r.bench_noise_sigma,x=sep.band_input)
    kw['separate_frontend']=MFBASE.source_step(s.base.frontends.separate,**a)
    kw['fma_frontend']=MFBASE.source_step(s.base.frontends.fma,**a,last=True)
    return raw,dict(dt=h,guard_witnesses=gw,separate_source_witnesses=sw,fma_source_witnesses=fw,**api,**kw)


def edge_fixture(pending):
    """Isolate the handoff implication, not fabricate a capture trajectory."""
    base,entry,fresh,active,scheduler=GOBASE.startup(pending)
    root=initial(); runtime=replace(root.runtime,commit_cfg=GOBASE.BASE.cfg(),bench_noise_sigma=F(0),boundary_bench_noise_sigma=F(0))
    state=X.State(base,runtime,root.deployment_cfg,root.guard,root.separate_source,root.fma_source,root.racc)
    return state,entry,fresh,scheduler


def admitted_fixture(pending=False):
    """Compose the conditional handoff with real lower source-admission objects."""
    from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as ADMIT
    from tools.stability.ou3_alt_contraction import finite_admitted_bias_history as BIAS
    import test_finite_admitted_source_live_word as SOURCE
    import test_finite_live_magnetic_word as MAG
    import test_finite_source_bound_live_word as LIVE
    s,entry,fresh,scheduler=edge_fixture(pending)
    go=X.go_live(s,entry,fresh,scope=SCOPE.certified_scope(),
        noise_sqrt=BAND.NoiseSqrtWitness(0),scheduler=scheduler,aw_sync=AWSYNC.State())
    _,magnetic=SOURCE.Tests().dual_startup()
    ref=go.lower.live.state.mekf.reference
    origin=ADMIT.RestrictedOrigin(ADMIT.AdmittedHistory(ref.history_id),ref)
    bias=BIAS.AdmittedBiasHistory(ref.bias_root,ref.bias_family,LIVE.PHI)
    kw=dict(gyro_residual_history_id='gyro-start',accel_residual_history_id='acc-start',
        temperature_model_history_id='temp-start',supply_norm_upper=1,
        proxy_q_norm=MAG.X.TILT.SqrtWitness(1,1),proxy_yaw_half=MAG.zero_yaw())
    return go,magnetic,origin,bias,kw


class Tests(unittest.TestCase):
    def test_literal_reset_has_no_external_machine_seed_port(self):
        s=initial(); self.assertEqual(s.guard,X.GUARD.State())
        self.assertEqual(s.separate_source.vertical,V.State())
        self.assertEqual(s.base.frontends.samples,0)
        self.assertEqual(s.base.lower.frontend.tuner.warmup_sec,10)
        self.assertEqual(s.base.machine.sigma.separate,B.rn32(F(1,100)))
        with self.assertRaises(TypeError): X.initial(s.runtime,s.deployment_cfg,vertical=V.State(initialized=True))

    def test_two_Cold_samples_preserve_one_guard_Mahony_LPF_band_history(self):
        s=initial(); raw,kw=cold_operands(s); a=X.step(s,raw,**kw)
        self.assertEqual(a.separate_source.mahony.startup.branch,'ordinary-FromTwoVectors')
        raw,kw=cold_operands(a.state); b=X.step(a.state,raw,**kw)
        self.assertIsNone(b.separate_source.mahony.startup)
        self.assertEqual(b.state.base.frontends.samples,2)
        self.assertEqual(b.state.guard.samples,2); self.assertEqual(b.state.separate_source.samples,2)
        self.assertEqual(b.state.base.machine,s.base.machine)  # Cold no candidate
        self.assertEqual(b.state.racc,s.racc)  # bootstrap does not drive MEKF
        self.assertIs(b.separate_source.before,a.state.separate_source)
        self.assertEqual(b.separate_source.mahony,b.fma_source.mahony)
        self.assertEqual(b.lower.separate_frontend.band.envelope.x,b.separate_source.band_input)

    def test_band_input_splice_is_rejected_at_same_event_join(self):
        s=initial(); raw,kw=cold_operands(s)
        a=dict(band_cfg=s.runtime.band_cfg,stats_cfg=s.runtime.stats_cfg,frequency=B.rn32(F(1,5)),
            dt=RAW.DT,bench_noise_sigma=s.runtime.bench_noise_sigma,x=B.rn32(F(1)))
        kw['separate_frontend']=MFBASE.source_step(s.base.frontends.separate,**a)
        with self.assertRaisesRegex(ValueError,'SAME Mahony'):
            X.step(s,raw,**kw)

    def test_boundary_preserves_upstream_machine_memory(self):
        s=initial(); raw,kw=cold_operands(s); s=X.step(s,raw,**kw).state
        out,event=X.boundary(s)
        self.assertIs(out.guard,s.guard); self.assertIs(out.separate_source,s.separate_source)
        self.assertIs(out.racc,s.racc); self.assertFalse(event.machine.consumed)

    def test_goLive_preserves_pending_and_upstream_memory_without_EMA(self):
        for pending in (False,True):
            s,entry,fresh,scheduler=edge_fixture(pending)
            go=X.go_live(s,entry,fresh,scope=SCOPE.certified_scope(),
                noise_sqrt=BAND.NoiseSqrtWitness(0),scheduler=scheduler,aw_sync=AWSYNC.State())
            self.assertIs(go.startup,s); self.assertIs(go.lower.machine,s.base.machine)
            self.assertIs(go.lower.frontends,s.base.frontends)
            self.assertEqual(go.lower.machine.pending,pending)
            self.assertEqual(go.lower.live.state.racc,s.racc)
            self.assertIsNotNone(go.lower.separate_commit)  # unconditional commit

    def test_goLive_cannot_supply_new_Racc_state(self):
        s,entry,fresh,scheduler=edge_fixture(False)
        with self.assertRaisesRegex(TypeError,'cannot replace carried'):
            X.go_live(s,entry,fresh,racc=RACC.State(),scope=SCOPE.certified_scope(),scheduler=scheduler,aw_sync=AWSYNC.State())
        with self.assertRaisesRegex(ValueError,'drive_mekf=false'):
            replace(s,racc=RACC.State(True,(1,1,1)))

    def test_runtime_override_and_source_count_splices_are_rejected(self):
        s=initial(); raw,kw=cold_operands(s)
        with self.assertRaisesRegex(TypeError,'override carried'):
            X.step(s,raw,guard_cfg=s.runtime.guard_cfg,**kw)
        with self.assertRaisesRegex(ValueError,'guard count'):
            replace(s,guard=replace(s.guard,samples=1))

    def test_GoLive_object_rejects_a_different_startup_history(self):
        s,entry,fresh,scheduler=edge_fixture(False)
        go=X.go_live(s,entry,fresh,scope=SCOPE.certified_scope(),noise_sqrt=BAND.NoiseSqrtWitness(0),scheduler=scheduler,aw_sync=AWSYNC.State())
        other=replace(s,base=replace(s.base,lower=replace(s.base.lower,
            frontend=replace(s.base.lower.frontend,tuner=replace(s.base.lower.frontend.tuner,stage_time=F(9))))))
        with self.assertRaisesRegex(ValueError,'detached from joined'):
            replace(go,startup=other)

    def test_admitted_Live_factory_carries_exact_same_conditional_startup(self):
        from tools.stability.ou3_alt_contraction import finite_admitted_machine_joined_frontend_racc_interleaved_prefix as JOIN
        for pending in (False,True):
            go,magnetic,origin,bias,kw=admitted_fixture(pending)
            live=X.admitted_live(go,magnetic,origin,bias,**kw)
            mt=JOIN._mtune_state(live.base)
            self.assertIs(live.guard,go.startup.guard)
            self.assertIs(live.separate_source,go.startup.separate_source)
            self.assertIs(live.fma_source,go.startup.fma_source)
            self.assertIs(mt.machine,go.lower.machine)
            self.assertIs(mt.frontends,go.lower.frontends)
            self.assertEqual(mt.base.wpe,go.lower.wpe)
            self.assertEqual(mt.machine.pending,pending)
            self.assertEqual(live.base.separate_racc,go.startup.racc)
            self.assertEqual(live.base.fma_racc,go.startup.racc)
            self.assertEqual(live.source_steps,0)
            self.assertEqual(live.base.racc_steps,0)
            # Entry attachment is not permission to skip the literal 600 edges.
            with self.assertRaises((ValueError,TypeError)):
                JOIN.complete(live)

    def test_admitted_factory_rejects_new_runtime_or_scheduler_predecessor(self):
        go,magnetic,origin,bias,kw=admitted_fixture()
        with self.assertRaisesRegex(TypeError,'owned by joined startup'):
            X.admitted_live(go,magnetic,origin,bias,runtime=go.startup.runtime,**kw)
        with self.assertRaisesRegex(ValueError,'scheduler detached'):
            replace(go,scheduler_before=replace(go.scheduler_before,elapsed=go.scheduler_before.period/2))

    def test_admitted_hold_does_not_restart_startup_machine_memory(self):
        from tools.stability.ou3_alt_contraction import finite_admitted_machine_joined_frontend_racc_interleaved_prefix as JOIN
        go,magnetic,origin,bias,kw=admitted_fixture()
        live=X.admitted_live(go,magnetic,origin,bias,**kw)
        after,_=JOIN.set_hold(live,hold=False)
        self.assertIs(after.guard,live.guard)
        self.assertIs(after.separate_source,live.separate_source)
        self.assertIs(after.fma_source,live.fma_source)
        self.assertEqual(after.source_steps,0)
        self.assertEqual(JOIN._mtune_state(after.base).machine,go.lower.machine)

    def test_startup_raw_packet_history_and_clock_splices_are_rejected(self):
        s=initial(); raw,kw=cold_operands(s)
        s=X.step(s,raw,**kw).state
        raw,kw=cold_operands(s)
        with self.assertRaisesRegex(ValueError,'sample-entry clock'):
            X.step(s,replace(raw,physical=replace(raw.physical,time=F(99))),**kw)
        with self.assertRaisesRegex(ValueError,'history identity changed'):
            X.step(s,replace(raw,physical=replace(raw.physical,history_id='other-history')),**kw)

    def test_readiness_and_seed_gap_do_not_unlock_storage(self):
        r=X.readiness()
        for k in ('literal_reset_roots_guard_observer_LPF_stillness_band_and_TuneState',
                  'Cold_and_postCold_machine_sources_joined_to_one_executed_frontend_event',
                  'conditional_goLive_preserves_joined_machine_source_history'): self.assertTrue(r[k])
        for k in ('startup_source_membership_and_clock_reachability_closed','nearly_antiparallel_seed_SVD_closed',
                  'ungauged_timeout_machine_handoff_closed','source_uniform_complete_600_step_word_qualified',
                  'storage_search_allowed','ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'): self.assertFalse(r[k])
        computed=SEED.seed((0,0,B.rn32(10)))
        self.assertIsNotNone(computed.solver)
        self.assertEqual(computed.branch,'near-antiparallel-JacobiSVD')


if __name__=='__main__': unittest.main()
