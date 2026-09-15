"""Conditional ungauged runtime composition; no all-source capture claim."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

import test_finite_live_interleave as I
import test_finite_live_magnetic_word as T
from tools.stability.ou3_alt_contraction import finite_live_interleave as X
from tools.stability.ou3_alt_contraction import finite_live_magnetic_word as MAG
from tools.stability.ou3_alt_contraction import finite_startup_handoff_seed as SEED
from tools.stability.ou3_alt_contraction import finite_startup_handoff_init as INIT
from tools.stability.ou3_alt_contraction import finite_startup_live_entry as ENTRY
from tools.stability.ou3_alt_contraction import finite_fresh_joint24_entry as FRESH
from tools.stability.ou3_alt_contraction import deployment_scope as SCOPE


def root(*,samples=2,refine=False,field=(30,0,0),south=True):
    bridge=I.FIRST.startup(); core=bridge.state.mekf; proxy=bridge.frontend_before.tuner.vertical
    seed=SEED.seed(proxy,None)
    hand=INIT.initialize_from_seed_zero_heel(seed,(0,)*21,MAG.M.eye(21),
        scope=SCOPE.certified_scope(),q_norm=MAG.TILT.SqrtWitness(1,1))
    entry=ENTRY.enter_live(hand,bridge.active,scope=SCOPE.certified_scope())
    ref=replace(core.reference,q_world_to_body=(0,0,0,1) if south else (1,0,0,0))
    core=FRESH.build(entry,ref,scope=SCOPE.certified_scope())
    bridge=replace(bridge,state=replace(bridge.state,mekf=core))
    cfg=MAG.Config(gate=MAG.GATE.Config(mag_delay=0),gravity=MAG.GRAVITY.Config(mag_delay=0),
        tuner=MAG.TUNER.Config(min_samples=samples,min_window=0),
        refinement_enabled=refine,refinement_start=0,refinement_window=0)
    word=MAG.START.State(MAG.GRAVITY.State(gravity_good=10,aligned_branch=True),T.PREFIX.State(proxy))
    magnetic=MAG.begin_startup(word,cfg,MAG.SOURCE.Model(field,(0,0,0),'fixed-field'),
        core.reference.history_id,'same-reference-producer')
    return X.from_startup(bridge,magnetic,proxy_q_norm=MAG.TILT.SqrtWitness(1,1))


def packet(state,*,ready=False,norm=30,**kwargs):
    ts=I.WCLOCK.at_physical_time(state.live.live.mekf.reference.time)
    dt=I.WCLOCK.shipping_elapsed(ts,state.magnetic.memory.last_hi_time,
                                 fallback_dt=state.magnetic.memory.cfg.sample_dt)
    initial={'q_norm':MAG.TILT.SqrtWitness(1,1),'yaw_half':T.zero_yaw(),
             'mag_norm':MAG.TUNER.SqrtWitness(norm*norm,norm)}
    if ready:
        initial.update(mean_norm=MAG.TUNER.SqrtWitness(norm*norm,norm),
            horizontal_sqrt=T.GAUGE.HorizontalSqrt(norm*norm,norm),
            gauge_half=MAG.TILT.YawHalfWitness(-norm,0,MAG.TILT.SqrtWitness(norm*norm,norm),0,1))
        kwargs.setdefault('ldlt',T.REJECT)
    return X.mag_step(state,residual_body=(0,0,0),packet_id='same-packet',
        proxy_q_norm=MAG.TILT.SqrtWitness(1,1),proxy_yaw_half=T.zero_yaw(),
        hi_decay=MAG.HI.Decay(dt,state.magnetic.memory.cfg.continuous.memory,1),
        initial=initial,**kwargs)


class Tests(unittest.TestCase):
    def test_ungauged_IMU_waiting_north_and_next_IMU_preserve_one_history(self):
        s=root(); self.assertEqual(s.live.live.mekf.attitude_chart,3)
        first=packet(s)
        self.assertIsInstance(first.state.magnetic,MAG.UngaugedLiveState)
        self.assertIs(first.event.filter,s.live.live.mekf)
        self.assertEqual(first.state.clock.calls,0)
        raw,seg,kw=I.imu_operands(first.state)
        with self.assertRaisesRegex(ValueError,'gravity-gate arithmetic'):
            X.imu_step(first.state,raw,seg,**kw)
        middle=X.imu_step(first.state,raw,seg,gravity_witnesses={
            'lpf_exp':None,'gyro_norm':MAG.GRAVITY.SqrtWitness(0,0)},**kw)
        north=packet(middle.state,ready=True)
        self.assertIsInstance(north.state.magnetic,MAG.LiveState)
        self.assertEqual(north.state.live.live.mekf.attitude_chart,0)
        self.assertEqual(north.state.live.live.mekf.z[:3],(0,0,0))
        self.assertEqual(north.state.live.live.mekf.reference.live_origin,0)
        self.assertEqual(north.state.clock.live_time,seg.after.time)
        self.assertEqual(north.state.clock.calls,1)
        self.assertEqual(north.state.magnetic.control.updates,1)
        self.assertEqual(north.state.magnetic.memory.continuous.weight,2)
        self.assertEqual(north.state.live.live.mekf.covariance,middle.state.live.live.mekf.covariance)
        self.assertEqual(north.state.live.live.mekf.z[3:],middle.state.live.live.mekf.z[3:])
        after=I.imu(north.state)
        self.assertEqual(after.state.live.live.mekf.reference.live_origin,0)
        self.assertIs(after.state.magnetic,north.state.magnetic)

    def test_Live_acquisition_uses_MEKF_tilt_continuous_uses_private_tilt(self):
        s=root(field=(15,20,0),south=False)
        core=s.live.live.mekf; q=(F(4,5),F(3,5),0,0)
        attitude=MAG.CORE.ATLAS.encode(MAG.P.quat_conj(q))
        core=replace(core,q_hat=q,z=attitude.coordinates+core.z[3:],attitude_chart=attitude.chart)
        s=replace(s,live=replace(s.live,live=replace(s.live.live,mekf=core)))
        out=packet(s,norm=25)
        self.assertEqual(out.state.magnetic.startup.word.mag.tuner.last_world_sample,(15,F(28,5),F(-96,5)))
        self.assertEqual(out.state.magnetic.memory.continuous.level_sum,(15,20,0))
        self.assertIs(out.state.live.live.tuner,s.live.live.tuner)

    def test_same_packet_initial_north_and_refinement_accumulate_continuous_once(self):
        s=root(samples=1,refine=True)
        south=MAG.TILT.YawHalfWitness(-1,0,MAG.TILT.SqrtWitness(1,1),0,1)
        gauge=MAG.TILT.YawHalfWitness(-30,0,MAG.TILT.SqrtWitness(900,30),0,1)
        out=packet(s,ready=True,mag_norm=MAG.TUNER.SqrtWitness(900,30),
            mean_norm=MAG.TUNER.SqrtWitness(900,30),horizontal_sqrt=T.GAUGE.HorizontalSqrt(900,30),
            gauge_half=gauge,mekf_q_norm=MAG.TILT.SqrtWitness(1,1),mekf_yaw_half=south)
        self.assertTrue(out.state.magnetic.refinement_done)
        self.assertEqual(out.state.magnetic.memory.continuous.weight,1)
        self.assertEqual(out.state.magnetic.control.updates,1)
        self.assertEqual(out.state.live.live.mekf.attitude_chart,0)

    def test_waiting_call_cannot_consume_an_inner_measurement(self):
        with self.assertRaisesRegex(ValueError,'measurement suffix'):
            packet(root(),ldlt=T.REJECT)

    def test_later_north_clock_cannot_reanchor_physical_S(self):
        s=root(samples=1); out=packet(s,ready=True)
        with self.assertRaisesRegex(ValueError,'detached'):
            replace(out.state,clock=X.SCHEDULE.Clock(F(1),F(1),1))


if __name__=='__main__': unittest.main()
