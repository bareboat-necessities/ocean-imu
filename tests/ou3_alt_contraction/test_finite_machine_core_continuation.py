"""Consecutive local supplies must not impersonate a deployed machine word."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_machine_core_continuation as X
from tools.stability.ou3_alt_contraction import finite_machine_prediction_displacement as DISP
from tools.stability.ou3_alt_contraction import finite_machine_active_prediction_roots as MROOT
from tools.stability.ou3_alt_contraction import finite_source_bound_prediction_word as EROOT
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_admitted_wpe_machine_clock_interleaved_prefix as LIVE
from tools.stability.ou3_alt_contraction import finite_mag_bias_gate as GATE
import test_finite_source_bound_live_word as WORDTEST
import test_finite_admitted_machine_prediction_supply_interleaved_prefix as BASE
import test_finite_source_bound_prediction_word as ROOTBASE
import test_finite_admitted_wpe_machine_clock_interleaved_prefix as WPEBASE
import test_finite_admitted_machine_joined_frontend_racc_interleaved_prefix as JOINBASE
import test_finite_tilt_reset_runtime as RESETBASE
import test_finite_startup_joined_machine_history as GUARDTEST
import test_finite_ungauged_live_word as UNGAUGED


def joined_operands(state):
    """Extend the existing conditional fixture from its actual prior history."""
    J=JOINBASE; V=J.VBASE; T=J.TBASE
    r=state.base; mt=J.X._mtune_state(r); runtime=J.X._runtime(r)
    tkw,_=T.event_operands(mt); lower=T.LOWER.imu_step(mt.base,**tkw)
    live=T.X._live_result(lower); h=tkw['restricted'].segment.h
    api=V.machine_api(live,h)
    gw=GUARDTEST.guard_witnesses(state,api['machine_acc_body'],api['machine_dt'])
    guarded=J.X.GUARD.step(state.guard,state.guard_cfg,
        raw_gyro=api['machine_gyro_body'],raw_acc=api['machine_acc_body'],dt=api['machine_dt'],**gw)
    sep,sw=V.source_witness(state.separate_source,guarded,runtime,h,last=False)
    fma,fw=V.source_witness(state.fma_source,guarded,runtime,h,last=True)
    join=dict(api,**{'guard_'+name:value for name,value in gw.items()})
    for mode,w in (('separate',sw),('fma',fw)):
        for name in ('lpf_alpha_exp','lpf_successor','still_energy_successor','still_attenuation_exp'):
            join[mode+'_'+name]=w[name]
    kw=J.RBASE.event_operands(r)
    for mode,src,freq,tau in (('separate',sep,lower.separate_frequency,lower.tau_step.separate_step),
                            ('fma',fma,lower.fma_frequency,lower.tau_step.fma_step)):
        front=V.MFBASE.source_step(getattr(mt.frontends,mode),band_cfg=runtime.band_cfg,
            stats_cfg=runtime.stats_cfg,frequency=freq.external.stored.input_hz,
            bench_noise_sigma=runtime.bench_noise_sigma,x=src.band_input,last=mode=='fma')
        sigma=V.sigma_target(mt,front,src)
        kw[mode+'_frontend']=front; kw[mode+'_sigma_machine']=sigma
        kw[mode+'_spectral_pow']=V.MBASE.spectral_pow_for(tau.tau_target,sigma.sigma_target,mt.deployment_cfg)
        kw[mode+'_spectral_sqrt']=V.MBASE.spectral_sqrt_for(tau.tau_target,mt.deployment_cfg)
        kw[mode+'_rs_exp_decay']=V.MBASE.rs_exp(mt.deployment_cfg,tau.tau_target)
    return kw,join,sep,fma


class Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.before, cls.first = BASE.executed_supply(pending=True)

    def initial(self):
        return X.begin(self.first.separate.predecessor, self.first.fma.predecessor)

    def first_observation(self):
        return X.observe_imu(self.initial(),
            separate_predecessor=self.first.separate.predecessor,
            fma_predecessor=self.first.fma.predecessor,
            separate_after_accel=self.first.separate.machine,
            fma_after_accel=self.first.fma.machine)

    def test_two_local_predictions_do_not_chain_the_machine_covariance(self):
        # This is a conditional algebra fixture, never source-admission evidence.
        # All subsequent measurement branches in the lower fixture reject, so
        # its stored CORE is exactly the shadow prediction.
        word = BASE.X._preword(self.first.state.base)
        core = word.live.live.live.mekf
        self.assertEqual(core, self.first.separate.exact)
        self.assertNotEqual(core.covariance, self.first.separate.machine.covariance)
        residual = max(abs(x-y) for a,b in zip(core.covariance,
                       self.first.separate.machine.covariance) for x,y in zip(a,b))
        self.assertGreater(residual, F(14, 1000))
        self.assertLess(residual, F(15, 1000))

        witness, segment, raw, _ = ROOTBASE.operands(word)
        physical = SOURCE.QualifiedPhysicalSegment(word.source.root, witness, segment)
        args = ROOTBASE.root_args(); args.pop('temperature_c')
        exact = EROOT.build(word, physical, raw, **args)
        relations = []
        for mode in ('separate', 'fma'):
            old = getattr(self.first, mode)
            machine_args = dict(args, ou_alpha=old.machine_roots.ou.alpha,
                                ou_em1=old.machine_roots.ou.em1)
            machine = MROOT.build(word, physical, raw, old.machine_roots.active,
                mode=mode, exact_active=word.live.live.live.active, **machine_args)
            relations.append(DISP.compare(core, segment, raw, exact, machine,
                                           Qbase=word.runtime.Qbase))
        carried = X.observe_imu(self.first_observation(),
            separate_predecessor=relations[0].predecessor,
            fma_predecessor=relations[1].predecessor,
            separate_after_accel=relations[0].machine,
            fma_after_accel=relations[1].machine)
        self.assertFalse(carried.observations[-1].separate_predecessor_connected)
        self.assertFalse(carried.observations[-1].fma_predecessor_connected)
        self.assertIn('separate_machine_CORE_successor_not_consumed', carried.missing)
        self.assertEqual(carried.separate, relations[0].machine)

    def test_MAG_HOLD_do_not_silently_use_the_shadow_successor(self):
        initial = self.first_observation()
        state = X.observe_async(initial, kind='mag')
        state = X.observe_async(state, kind='hold')
        self.assertIs(state.separate, initial.separate)
        self.assertEqual(state.imu_steps, 1)
        self.assertIn('mag_machine_CORE_and_control_successor_not_attached', state.missing)
        self.assertIn('hold_machine_CORE_and_control_successor_not_attached', state.missing)

    def test_600_counters_cannot_replace_event_correspondence(self):
        events = tuple(X.Observation('imu', n, True, True) for n in range(1, 601))
        forged = replace(self.initial(), observations=events)
        with self.assertRaisesRegex(ValueError, 'successor algorithms not attached'):
            X.require_complete(forged)

    def test_covariance_atlas_bias_and_source_are_part_of_predecessor_identity(self):
        state = self.first_observation()
        self.assertEqual(len(state.separate.z), 24)
        self.assertEqual(len(state.separate.covariance), 21)
        self.assertEqual(len(state.separate.covariance[0]), 21)
        self.assertEqual(state.separate.reference.beta, state.separate.z[21:])
        self.assertIn('machine_post_accelerometer_event_and_control_suffix_not_attached',
                      state.missing)

    def test_strongest_word_requires_this_continuation(self):
        r = LIVE.readiness()
        self.assertTrue(r['complete_word_requires_machine_CORE_and_control_continuation'])
        self.assertTrue(r['machine_CORE_local_successors_persist_and_next_predecessors_checked'])
        self.assertTrue(r['machine_CORE_and_control_successor_algorithms_attached'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])

    def test_distinct_predecessor_relation_keeps_prior_covariance_discrepancy(self):
        word=BASE.X._preword(self.first.state.base)
        witness,segment,raw,_=ROOTBASE.operands(word)
        physical=SOURCE.QualifiedPhysicalSegment(word.source.root,witness,segment)
        args=ROOTBASE.root_args(); args.pop('temperature_c')
        exact=EROOT.build(word,physical,raw,**args)
        previous=self.first.separate.machine
        old=self.first.separate.machine_roots
        machine=MROOT.build(word,physical,raw,old.active,mode='separate',
            machine_predecessor=previous,exact_active=word.live.live.live.active,
            **dict(args,ou_alpha=old.ou.alpha,ou_em1=old.ou.em1))
        out=DISP.compare(word.live.live.live.mekf,segment,raw,exact,machine,
            machine_predecessor=previous,Qbase=word.runtime.Qbase)
        direct=DISP.PRED.prediction_from_raw(previous,segment,raw,
            angular=machine.angular,ou=machine.ou,bias=machine.bias,qaxis=machine.qaxis,
            Qbase=word.runtime.Qbase)
        self.assertEqual(out.machine,direct)
        self.assertIs(out.machine_predecessor,previous)
        self.assertNotEqual(out.machine_predecessor.covariance,out.predecessor.covariance)
        # Covariance hygiene may take a different solver branch in each mode.
        retry={'attitude_first_ldlt_success':False,'attitude_second_ldlt_success':False}
        retried=DISP.compare(word.live.live.live.mekf,segment,raw,exact,machine,
            machine_predecessor=previous,machine_attitude_solver=retry,Qbase=word.runtime.Qbase)
        retry_direct=DISP.PRED.prediction_from_raw(previous,segment,raw,
            angular=machine.angular,ou=machine.ou,bias=machine.bias,qaxis=machine.qaxis,
            Qbase=word.runtime.Qbase,**retry)
        self.assertEqual(retried.exact,out.exact)
        self.assertEqual(retried.machine,retry_direct)
        self.assertNotEqual(retried.machine.covariance,out.machine.covariance)
        with self.assertRaisesRegex(TypeError,'only same-mode attitude LDLT'):
            DISP.compare(word.live.live.live.mekf,segment,raw,exact,machine,
                machine_predecessor=previous,machine_attitude_solver={'Qbase':word.runtime.Qbase},
                Qbase=word.runtime.Qbase)
        with self.assertRaisesRegex(ValueError,'roots detached from carried predecessor'):
            DISP.compare(word.live.live.live.mekf,segment,raw,exact,machine,
                         Qbase=word.runtime.Qbase)

    def _actual_runtime(self,*,unlocked=False):
        runtime=BASE.X._preword(self.first.state.base).live
        core=self.first.separate.machine
        magnetic=runtime.magnetic
        if unlocked:
            core=GATE._active(core,magnetic.memory.cfg.gate.initial_ba_std)
            cov=[list(row) for row in core.covariance]
            cov[6][18]=cov[18][6]=F(1,100)
            core=replace(core,covariance=tuple(tuple(row) for row in cov))
            magnetic=replace(magnetic,control=replace(magnetic.control,locked=False,hold=False))
        runtime=replace(runtime,live=replace(runtime.live,
            live=replace(runtime.live.live,mekf=core)),magnetic=magnetic)
        return X.begin_from_interleaved(runtime)

    def test_carried_gyro_bias_error_drives_its_own_angular_roots_and_attitude(self):
        word=BASE.X._preword(self.first.state.base)
        witness,segment,raw,_=ROOTBASE.operands(word)
        physical=SOURCE.QualifiedPhysicalSegment(word.source.root,witness,segment)
        args=ROOTBASE.root_args(); args.pop('temperature_c')
        exact=EROOT.build(word,physical,raw,**args)
        # A nonzero bias error inside the existing small-rate branch needs no
        # caller-selected trigonometric witness, and still changes attitude.
        z=list(self.first.separate.machine.z); z[3]+=F(1,10**8)
        previous=replace(self.first.separate.machine,z=tuple(z))
        old=self.first.separate.machine_roots
        machine=MROOT.build(word,physical,raw,old.active,mode='separate',
            machine_predecessor=previous,exact_active=word.live.live.live.active,
            **dict(args,ou_alpha=old.ou.alpha,ou_em1=old.ou.em1))
        self.assertNotEqual(machine.angular.w,exact.angular.w)
        self.assertEqual(machine.angular.w,raw.required_bias_corrected_relation(previous.z[3:6]))
        out=DISP.compare(word.live.live.live.mekf,segment,raw,exact,machine,
            machine_predecessor=previous,Qbase=word.runtime.Qbase)
        self.assertNotEqual(out.machine.q_hat,out.exact.q_hat)
        self.assertEqual(out.machine.reference,out.exact.reference)

    def test_machine_prediction_retains_nonzero_source_API_rounding_defect(self):
        word=BASE.X._preword(self.first.state.base)
        witness,segment,raw,_=ROOTBASE.operands(word)
        rate=F(1,10**8)+F(1,10**18)
        raw=replace(raw,raw_gyro_body=(rate,0,0),gyro_residual_internal=(rate,0,0))
        packet=X.INPUT.check_packet(raw)
        self.assertNotEqual(packet.machine_gyro,raw.raw_gyro_body)
        physical=SOURCE.QualifiedPhysicalSegment(word.source.root,witness,segment)
        args=ROOTBASE.root_args(); args.pop('temperature_c')
        exact=EROOT.build(word,physical,raw,**args)
        core=word.live.live.live.mekf; old=self.first.separate.machine_roots
        roots=MROOT.build(word,physical,raw,old.active,mode='separate',machine_predecessor=core,
            machine_packet=packet,exact_active=word.live.live.live.active,
            **dict(args,ou_alpha=old.ou.alpha,ou_em1=old.ou.em1))
        self.assertEqual(roots.angular.w,packet.machine_gyro)
        self.assertNotEqual(roots.angular.w,exact.angular.w)
        out=DISP.compare(core,segment,raw,exact,roots,machine_predecessor=core,Qbase=word.runtime.Qbase)
        direct=DISP.PRED.prediction(core,segment,gyro_body=packet.machine_gyro,
            angular=roots.angular,ou=roots.ou,bias=roots.bias,qaxis=roots.qaxis,Qbase=word.runtime.Qbase)
        self.assertEqual(out.machine,direct)
        self.assertNotEqual(out.machine.q_hat,out.exact.q_hat)
        other=replace(raw,raw_gyro_body=(rate/2,0,0),gyro_residual_internal=(rate/2,0,0))
        with self.assertRaisesRegex(ValueError,'another source packet'):
            MROOT.build(word,physical,raw,old.active,mode='separate',machine_predecessor=core,
                machine_packet=X.INPUT.check_packet(other),exact_active=word.live.live.live.active,
                **dict(args,ou_alpha=old.ou.alpha,ou_em1=old.ou.em1))

    def test_ungauged_gravity_reads_same_rounded_API_accel_and_gyro(self):
        runtime=UNGAUGED.root(south=False); core=runtime.live.live.mekf
        raw,_,_=UNGAUGED.I.imu_operands(runtime)
        rate=F(1,10**8)+F(1,10**18)
        raw=replace(raw,raw_gyro_body=(rate,0,0),gyro_residual_internal=(rate,0,0))
        packet=X.INPUT.check_packet(raw)
        sample=X.SENSOR.guarded_sample(raw,X.SENSOR.GUARD.State(),
                                      X.SENSOR.GUARD.Config(cutoff_hz=0),dt=SOURCE.DT)
        G=X.LIVE.MAG.GRAVITY; mg=packet.machine_gyro[0]
        gravity={'lpf_exp':None,'gyro_norm':G.SqrtWitness(mg*mg,abs(mg))}
        out=X._imu_suffix(runtime,runtime,core,runtime.live.live.tuner.vertical,
                          sample,SOURCE.DT,reset={},gravity=gravity,machine_packet=packet)
        self.assertEqual(out.magnetic.startup.word.gate.lpf,packet.machine_accel)
        self.assertNotEqual(packet.machine_accel,raw.raw_accel_body)
        with self.assertRaisesRegex(ValueError,'gyro norm witness detached'):
            X._imu_suffix(runtime,runtime,core,runtime.live.live.tuner.vertical,
                sample,SOURCE.DT,reset={},gravity={'lpf_exp':None,'gyro_norm':G.SqrtWitness(rate*rate,rate)},
                machine_packet=packet)

    def test_startup_entry_seats_aw_covariance_from_each_actual_machine_commit(self):
        # Conditional handoff fixture; no claim that this component predecessor
        # is reached by every startup source, or has target-qualified q/P.
        startup,entry,fresh,scheduler=GUARDTEST.edge_fixture(False)
        machine=startup.base.machine
        sigma=replace(machine.sigma,separate=X.INPUT.B.rn32(F(3,10)),
                      fma=X.INPUT.B.rn32(F(2,5)))
        startup=replace(startup,base=replace(startup.base,machine=replace(machine,sigma=sigma)))
        go=GUARDTEST.X.go_live(startup,entry,fresh,scope=GUARDTEST.SCOPE.certified_scope(),
            noise_sqrt=GUARDTEST.BAND.NoiseSqrtWitness(0),scheduler=scheduler,
            aw_sync=GUARDTEST.AWSYNC.State())
        _,magnetic,origin,bias,kw=GUARDTEST.admitted_fixture(False)
        wpe=LIVE.WPE.initial(startup.runtime.wpe_cfg)
        full=LIVE.STARTWPE.State(startup,wpe)
        full_go=LIVE.STARTWPE.GoLive(full,go,wpe)
        out=LIVE.begin_from_startup(full_go,magnetic,origin,bias,**kw)
        carry=out.machine_core; exact=LIVE._preword(out.base).live.live.live.mekf
        self.assertIs(carry.startup_go_live,go)
        self.assertEqual(carry.imu_steps,0)
        self.assertNotEqual(carry.separate.covariance,carry.fma.covariance)
        for mode in ('separate','fma'):
            core=getattr(carry,mode); active=getattr(go.lower,mode+'_active')
            self.assertEqual(tuple(tuple(row[15:18]) for row in core.covariance[15:18]),active.Sigma_aw)
            self.assertEqual(core.z,exact.z); self.assertEqual(core.q_hat,exact.q_hat)
            self.assertTrue(all(core.covariance[a][i]==core.covariance[i][a]==0
                for a in range(15,18) for i in range(21) if i<15 or i>=18))
            self.assertIs(getattr(carry,mode+'_runtime').live.live.tuner.vertical,
                          getattr(startup,mode+'_source').vertical)
        held,_=LIVE.set_hold(out,hold=True)
        self.assertIs(held.machine_core.startup_go_live,go)
        self.assertFalse(X.readiness()['machine_startup_quaternion_and_remaining_covariance_producer_attached'])

    def test_concrete_HOLD_release_preserves_machine_mean_and_full_covariance(self):
        state=self._actual_runtime(unlocked=True)
        core=state.separate
        held,events=X.advance_hold(state,hold=True)
        self.assertEqual(held.separate.mode,'H')
        self.assertEqual(held.separate.covariance[6][18],0)
        self.assertEqual(held.separate.z,core.z)
        self.assertEqual(held.separate.q_hat,core.q_hat)
        self.assertEqual(held.separate.attitude_chart,core.attitude_chart)
        released,events=X.advance_hold(held,hold=False)
        self.assertEqual(released.separate.mode,'A')
        self.assertEqual(released.separate.z,core.z)
        self.assertEqual(released.separate_runtime.live.live.mekf,released.separate)
        self.assertEqual(released.missing,())
        self.assertEqual(released.imu_steps,0)

    def test_firing_watchdog_uses_carried_attitude_and_its_own_timer(self):
        state=self._actual_runtime(); runtime=state.separate_runtime
        q=(F(3,5),F(4,5),F(0),F(0)); z=list(state.separate.z)
        attitude=X.CORE.ATLAS.encode(X.LIVE.P.quat_mul(state.separate.reference.q_world_to_body,
                                                    X.LIVE.P.quat_conj(q)))
        z[:3]=attitude.coordinates
        tilted=replace(state.separate,q_hat=q,z=tuple(z),attitude_chart=attitude.chart)
        runtime=replace(runtime,live=replace(runtime.live,
            watchdog=X.WATCH.State(X.WATCH.HOLD_SEC-SOURCE.DT)))
        sample=RESETBASE.guarded(tilted); T=RESETBASE.T; R=X.RESET
        witnesses=dict(old_q_norm=T.SqrtWitness(1,1),old_yaw_half=RESETBASE.atan_zero(),
            acc_tilt=R.AccTiltWitness(T.SqrtWitness(RESETBASE.GRAV**2,RESETBASE.GRAV)),
            pitch_cos=T.SqrtWitness(1,1),pitch_half=R.SignedHalfWitness(0,1,1,0),
            roll_half=RESETBASE.atan_zero())
        out=X._imu_suffix(runtime,runtime,tilted,runtime.live.live.tuner.vertical,
                          sample,SOURCE.DT,reset=witnesses,gravity=None)
        self.assertEqual(out.live.live.mekf.q_hat,(1,0,0,0))
        self.assertEqual(out.live.live.mekf.z[3:],tilted.z[3:])
        self.assertEqual(out.live.watchdog,X.WATCH.State(0,X.WATCH.COOLDOWN_SEC))
        self.assertTrue(all(out.live.live.mekf.covariance[i][j]==0
                            for i in range(3) for j in range(3,21)))
        # The same after-accel attitude with another carried timer does not fire.
        fresh=replace(runtime,live=replace(runtime.live,watchdog=X.WATCH.State()))
        noreset=X._imu_suffix(fresh,fresh,tilted,fresh.live.live.tuner.vertical,
                             sample,SOURCE.DT,reset={},gravity=None)
        self.assertEqual(noreset.live.live.mekf,tilted)
        self.assertEqual(noreset.live.watchdog.over_limit,SOURCE.DT)

    def test_concrete_MAG_consumes_persisted_machine_covariance_and_control(self):
        state=self._actual_runtime()
        word=BASE.X._preword(self.first.state.base)
        kwargs=WORDTEST.mag_kwargs(word)
        out,events=X.advance_mag(state,**kwargs)
        expected=X.LIVE.mag_step(state.separate_runtime,**kwargs)
        self.assertEqual(out.separate,expected.state.live.live.mekf)
        self.assertEqual(out.separate_runtime.magnetic,expected.state.magnetic)
        self.assertEqual(out.separate.reference,state.separate.reference)
        self.assertEqual(out.imu_steps,0)
        self.assertEqual(out.missing,())
        with self.assertRaisesRegex(TypeError,'override common source operands'):
            X.advance_mag(state,separate={'packet_id':'spliced'},**kwargs)

    def test_strongest_two_IMU_word_consumes_actual_prior_compiler_CORE(self):
        base,kw,join,_,sep,fma=JOINBASE.fixture()
        cfg=WPEBASE.WM.MOM.Config.from_shadow(JOINBASE.X._runtime(base.base).wpe_cfg)
        mom=WPEBASE.M.State(elapsed=WPEBASE.M.B.rn32(3),samples=0)
        logs=JOINBASE.X._mtune_state(base.base).base.wpe
        wpe=WPEBASE.WM.State(cfg,mom,mom,logs,separate_usable=True,fma_usable=True)
        state=LIVE.begin(WPEBASE.CLOCK.begin(base),wpe)
        for n in (1,2):
            if n==2:
                kw,join,sep,fma=joined_operands(state.base.base)
            before=state.machine_core
            sw=WPEBASE.general_witness(state.wpe,sep.band_input,join['machine_dt'],'separate')
            fw=WPEBASE.general_witness(state.wpe,fma.band_input,join['machine_dt'],'fma')
            out=LIVE.imu_step(state,separate_wpe=sw,fma_wpe=fw,
                separate_racc_accel_ldlt=WPEBASE.MAG.REJECT,
                fma_racc_accel_ldlt=WPEBASE.MAG.REJECT,**join,**kw)
            observation=out.state.machine_core.observations[-1]
            for mode in ('separate','fma'):
                measurement,racc,runtime=getattr(observation,mode+'_event')
                self.assertEqual(measurement.prediction.machine_predecessor,getattr(before,mode))
                self.assertEqual(measurement.prediction.machine_roots.machine_packet.raw,kw['raw'])
                self.assertEqual(measurement.prediction.machine_roots.machine_packet.machine_gyro,
                                 join['machine_gyro_body'])
                self.assertEqual(runtime.live.live.mekf,getattr(out.state.machine_core,mode))
            self.assertEqual(out.state.machine_core.imu_steps,n)
            self.assertEqual(out.state.machine_core.missing,())
            state=out.state


if __name__ == '__main__':
    unittest.main()
