"""Successive finite product-state identities, never source/stability evidence."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_live_interleave as X
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
from tools.stability.ou3_alt_contraction import finite_magnetic_wrapper_clock as WCLOCK
import test_finite_live_magnetic_word as MAG
import test_finite_startup_first_live_step as FIRST


def root(**kwargs):
    bridge, magnetic = MAG.startup(**kwargs)
    return X.from_startup(bridge, magnetic, proxy_q_norm=MAG.X.TILT.SqrtWitness(1,1),
                          proxy_yaw_half=MAG.zero_yaw())


def imu_operands(state):
    """Exact stationary physical segment and conditional finite runtime witnesses.

    This is a component fixture, NOT a deployment profile or arbitrary-word
    source qualification. The covariance and all persistent memory are carried,
    never refitted to this fixture. No storage/feasibility experiment occurs.
    """
    s=state.live.live; ref=s.mekf.reference; h=F(1,200); z=(0,0,0)
    segment=PHYS.PhysicalSegment(ref,replace(ref,time=ref.time+h),z,z,z,1,z)
    raw=FIRST.SENSOR.RawImuSample(ref,z,z,z,z,(0,0,-FIRST.G),(0,0,FIRST.G))
    kw=dict(
      dt=h,commit_cfg=FIRST.SB.cfg(),boundary_bench_noise_sigma=0,
      guard_cfg=FIRST.GUARD.Config(),guard_decay=FIRST.GUARD.DecayWitness(1,1,0,0),
      guard_rms=FIRST.GUARD.RmsWitness(0),
      vertical_cfg=FIRST.V.Config(0,0,FIRST.G,20),
      accel_invnorm=FIRST.V.InvSqrtWitness(FIRST.G**2,1/FIRST.G),
      quat_invnorm=FIRST.V.InvSqrtWitness(1,1),
      band_cfg=FIRST.B.BandConfig(F(1,2),4,F(1,100),6,F(3,100),F(6,5),F(1,5)),
      racc_cfg=FIRST.RACC.Config(),nominal_racc_std=(1,1,1),
      angular=FIRST.A.AngularRuntime(raw.bias_corrected_internal_gyro(s.mekf.z[3:6]),h),
      Qbase=FIRST.M.zeros(6,6),ou=FIRST.O.OUDecay(h,s.active.tau,F(199,200)),
      bias=FIRST.O.BiasDecay(s.mekf.mode=='A',s.active.tau,
            F(199,200) if s.mekf.mode=='A' else 1,FIRST.M.zeros(3,3)),
      qaxis=FIRST.PR.QAxisBranch(False,s.active.Sigma_aw,
            (FIRST.PASS,)*3,(FIRST.PASS,)*3,F(1,10**7)),
      accel_conditioning=FIRST.SENSOR.AccelConditioning(0,z),accel_ldlt=MAG.REJECT,
      wpe_cfg=FIRST.W.WPEConfig(1,4,F(1,2),1,180),wpe_decay=FIRST.W.ExpWitness(1),
      stats_cfg=FIRST.B.StatsConfig(4,F(3,10),60,F(1,20),5),
      band_decay=FIRST.B.BandDecayWitness(1,1),variance_decay=FIRST.B.VarianceDecayWitness(1),
      bench_noise_sigma=0,noise_sqrt=FIRST.B.NoiseSqrtWitness(0),
      tracker_lpf_decay=FIRST.FRONT.LPFDecayWitness(1),still_cfg=FIRST.SF.Config(),
      still_attenuation=FIRST.SF.AttenuationWitness(1),candidate_cfg=FIRST.LP.candidate_cfg(),
      sigma_wave_sqrt=1,spectral=FIRST.C.SpectralWitness(1,1),ema=FIRST.C.EmaWitness(1,1))
    if s.scheduler.elapsed+h >= s.scheduler.period:
        kw['S_ldlt']=MAG.REJECT
    if s.tuner.pending:
        kw['boundary_noise_sqrt']=FIRST.B.NoiseSqrtWitness(0)
    return raw,segment,kw


def imu(state):
    raw,segment,kw=imu_operands(state)
    return X.imu_step(state,raw,segment,**kw)


def mag(state,**overrides):
    kw=MAG.live_kwargs(state.magnetic)
    memory=state.magnetic.memory; physical=state.live.live.mekf.reference.time
    if memory.cfg.continuous_enabled:
        ts=WCLOCK.at_physical_time(physical)
        dt=WCLOCK.shipping_elapsed(ts,memory.last_hi_time,fallback_dt=memory.cfg.sample_dt)
        kw['hi_decay']=MAG.HI.Decay(dt,memory.cfg.continuous.memory,1)
    kw.update(overrides)
    return X.mag_step(state,**kw)


class Tests(unittest.TestCase):
    def test_successive_IMU_mag_IMU_consumes_actual_magnetic_covariance(self):
        s=root(); first=imu(s); before=first.state
        event=mag(before,ldlt=MAG.PASS)
        self.assertTrue(event.event.measurement.measurement_accepted)
        self.assertNotEqual(event.state.live.live.mekf.covariance,before.live.live.mekf.covariance)
        self.assertIs(event.state.live.live.tuner,before.live.live.tuner)
        self.assertIs(event.state.live.live.scheduler,before.live.live.scheduler)
        second=imu(event.state)
        raw,segment,kw=imu_operands(event.state)
        direct=X.LIVE.step_from_shipping_operands(event.state.live,raw,segment,**kw)
        self.assertEqual(second.state.live,direct.state)
        without_mag=imu(before)
        self.assertNotEqual(second.state.live.live.mekf.covariance,
                            without_mag.state.live.live.mekf.covariance)
        self.assertEqual(second.state.live.live.tuner.sample_index,s.live.live.tuner.sample_index+2)
        self.assertEqual(second.state.live.live.mekf.reference.live_origin,s.clock.live_time)
        self.assertIs(second.state.magnetic,event.state.magnetic)
        self.assertEqual(second.state.clock.calls,1)
        self.assertEqual(len(second.state.live.live.mekf.z),24)
        self.assertEqual(len(second.state.live.live.mekf.covariance),21)

    def test_rejected_repeated_timestamp_calls_keep_memory_and_count(self):
        s=imu(root()).state
        a=mag(s); b=mag(a.state)
        self.assertFalse(a.event.measurement.measurement_accepted)
        self.assertFalse(b.event.measurement.measurement_accepted)
        self.assertEqual(b.state.clock.calls,2)
        self.assertEqual(b.state.magnetic.control.updates,2)
        self.assertEqual(b.state.clock.last_time,a.state.clock.last_time)
        self.assertEqual(b.state.magnetic.memory.continuous.weight,3)
        self.assertIs(b.state.live.live.mekf,s.live.live.mekf)

    def test_refinement_rewrite_and_hold_release_reach_following_IMU(self):
        s=imu(root(refinement_start=0)).state
        out=mag(s,**MAG.refine_kwargs())
        self.assertTrue(out.state.magnetic.refinement_done)
        self.assertFalse(out.state.magnetic.control.hold)
        self.assertTrue(out.state.magnetic.control.locked)
        self.assertEqual(out.state.live.live.mekf.mode,'H')
        self.assertEqual(out.state.magnetic.active.generation,1)
        nxt=imu(out.state)
        self.assertIs(nxt.state.magnetic,out.state.magnetic)
        self.assertEqual(nxt.state.live.live.mekf.mode,'H')

    def test_external_hold_edge_preserves_frontend_and_physical_origin(self):
        s=root()
        out=X.set_hold(s,hold=False)
        self.assertFalse(out.state.magnetic.control.hold)
        self.assertEqual(out.state.live.live.mekf.mode,'H')
        self.assertIs(out.state.live.live.tuner,s.live.live.tuner)
        self.assertIs(out.state.live.live.mekf.reference,s.live.live.mekf.reference)
        self.assertIs(out.state.clock,s.clock)
        self.assertEqual(X.set_hold(out.state,hold=True).state.live.live.mekf.mode,'H')

    def test_clock_count_mode_and_history_cannot_be_restarted(self):
        s=mag(imu(root()).state).state
        with self.assertRaisesRegex(ValueError,'counted magnetic calls'):
            replace(s,clock=X.SCHEDULE.Clock(s.clock.live_time))
        bad=replace(s.magnetic.memory,history_id='detached')
        with self.assertRaisesRegex(ValueError,'physical history'):
            replace(s,magnetic=replace(s.magnetic,memory=bad))
        with self.assertRaisesRegex(ValueError,'H18/A21'):
            replace(s,magnetic=replace(s.magnetic,control=X.GATE.State(1,s.clock.last_time,False,False)))

    def test_missed_schedule_deadline_fails_without_inventing_a_call(self):
        s=replace(root(),schedule=X.SCHEDULE.Schedule(F(1,1000),F(1,25)))
        with self.assertRaisesRegex(ValueError,'deadline missed'):
            imu(s)
        self.assertEqual(s.magnetic.control.updates,0)

    def test_second_physical_endpoint_or_proxy_cannot_be_supplied_to_mag(self):
        s=root()
        for key,value in (('core',s.live.live.mekf),('proxy',s.live.live.tuner.vertical)):
            with self.assertRaises(TypeError):
                X.mag_step(s,**{key:value})

    def test_startup_proxy_and_pending_gauge_must_match_fresh_handoff(self):
        bridge,s=MAG.startup()
        changed=replace(bridge.frontend_before,tuner=replace(bridge.frontend_before.tuner,
            vertical=replace(bridge.frontend_before.tuner.vertical,elapsed=1)))
        with self.assertRaisesRegex(ValueError,'same startup private observer'):
            X.from_startup(replace(bridge,frontend_before=changed),s,
                proxy_q_norm=MAG.X.TILT.SqrtWitness(1,1),proxy_yaw_half=MAG.zero_yaw())

    def test_all_three_bias_labels_survive_interleaving_without_admission_claim(self):
        for family in ('BIAS0','BIAS1','BIAS2'):
            s=root(); core=s.live.live.mekf
            core=replace(core,reference=replace(core.reference,bias_family=family))
            s=replace(s,live=replace(s.live,live=replace(s.live.live,mekf=core)))
            s=mag(imu(s).state).state
            self.assertEqual(s.live.live.mekf.reference.bias_family,family)
            self.assertEqual(s.live.live.mekf.reference.bias_root,core.reference.bias_root)
            self.assertEqual(s.live.live.mekf.z[21:24],core.reference.beta)

    def test_zero_heel_scope_is_enforced_at_each_IMU_boundary(self):
        s=root(); raw,segment,kw=imu_operands(s)
        bad=replace(raw,deheel_body_to_internal=((-1,0,0),(0,-1,0),(0,0,1)))
        with self.assertRaisesRegex(ValueError,'zero-wind-heel'):
            X.imu_step(s,bad,segment,**kw)

    def test_free_tilt_and_reset_outputs_are_forbidden_in_product_word(self):
        s=root(); raw,segment,kw=imu_operands(s)
        with self.assertRaisesRegex(TypeError,'free tilt/reset'):
            X.imu_step(s,raw,segment,tilt_deg=0,**kw)

    def test_future_calibration_clock_is_not_a_valid_current_product_state(self):
        s=root()
        memory=replace(s.magnetic.memory,last_hi_time=F(1))
        with self.assertRaisesRegex(ValueError,'future physical endpoint'):
            replace(s,magnetic=replace(s.magnetic,memory=memory))

    def test_readiness_cannot_promote_finite_prefixes_to_universal_word(self):
        r=X.readiness()
        self.assertTrue(r['interleaved_IMU_uses_same_operand_tilt_reset_entry'])
        self.assertTrue(r['free_watchdog_angle_and_reset_quaternion_forbidden'])
        self.assertTrue(r['live_magnetic_outer_inner_dual_clock_composed'])
        self.assertTrue(r['live_magnetic_wrapper_clock_prefix_arithmetic_closed'])
        self.assertFalse(r['startup_magnetic_dual_clock_history_required_at_handoff'])
        for key in ('infinite_schedule_qualified_by_finite_prefix',
                    'source_uniform_complete_600_step_word_qualified','storage_search_allowed',
                    'ALT_STARTUP_PASS','ALT_LIVE_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(r[key])


if __name__=='__main__': unittest.main()
