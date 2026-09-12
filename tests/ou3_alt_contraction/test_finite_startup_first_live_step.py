"""First literal Live sample from exact startup ancestry."""
from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_startup_first_live_step as X
from tools.stability.ou3_alt_contraction import finite_startup_live_runtime_bridge as BRIDGE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as GUARD
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
from tools.stability.ou3_alt_contraction import finite_racc_runtime as RACC
from tools.stability.ou3_alt_contraction import finite_periodic_aw_sync as AWSYNC
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as B
from tools.stability.ou3_alt_contraction import finite_frontend_runtime as FRONT
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as SF
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W
from tools.stability.ou3_alt_contraction import finite_measurement_runtime as MR
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_runtime as PR
from tools.stability.ou3_alt_contraction import finite_qaxis_runtime as Q
from tools.stability.ou3_alt_contraction import finite_attitude_runtime as A
from tools.stability.ou3_alt_contraction import finite_ou_runtime_primitives as O
from tools.stability.ou3_alt_contraction import deployment_scope as SCOPE
import test_finite_startup_live_runtime_bridge as SB
import test_finite_live_imu_prefix as LP
import test_finite_core as FC

G=F(196133,20000); PASS=Q.PSDWitness(True)


def startup():
    front,entry,fresh,active,scheduler=SB.objects(False)
    level=(F(0),F(0),-G); zero=(F(0),F(0),F(0))
    guard=GUARD.State(stages=(level,level,level,level),detect_stages=(zero,zero),
                      removed_ms=zero,weight=1,initialized=True)
    front=replace(front,guard=guard)
    return BRIDGE.bridge(entry,fresh,front,scope=SCOPE.certified_scope(),
        commit_cfg=SB.cfg(),bench_noise_sigma=0,noise_sqrt=B.NoiseSqrtWitness(0),
        scheduler=scheduler,racc=RACC.State(False,(1,1,1)),aw_sync=AWSYNC.State())


def physical_step(s):
    segment,_=FC.physical_successor(s.state.mekf)
    ref=s.state.mekf.reference
    # Fresh reference is zero-rate/zero-bias/zero-beta/zero-acceleration.  One
    # exact held level accelerometer sample therefore satisfies the raw source
    # identity and is reused after prediction by the existing held-sample lemma.
    raw=SENSOR.RawImuSample(ref,(0,0,0),(0,0,0),(0,0,0),
                            (0,0,0),(0,0,-G),(0,0,G))
    h=segment.h
    angular=A.AngularRuntime((0,0,0),h)
    ou=O.OUDecay(h,s.active.tau,F(199,200))
    bias=O.BiasDecay(False,s.active.tau,1,M.zeros(3,3))
    qaxis=PR.QAxisBranch(False,s.active.Sigma_aw,(PASS,PASS,PASS),(PASS,PASS,PASS),F(1,10**7))
    return raw,segment,angular,ou,bias,qaxis


def run(s=None,**overrides):
    s=s or startup(); raw,segment,angular,ou,bias,qaxis=physical_step(s)
    kw=dict(
      dt=segment.h,commit_cfg=SB.cfg(),boundary_bench_noise_sigma=0,
      guard_cfg=GUARD.Config(),guard_decay=GUARD.DecayWitness(1,1,0,0),guard_rms=GUARD.RmsWitness(0),
      vertical_cfg=V.Config(0,0,G,20),accel_invnorm=V.InvSqrtWitness(G*G,1/G),quat_invnorm=V.InvSqrtWitness(1,1),seed=None,
      band_cfg=B.BandConfig(F(1,2),4,F(1,100),6,F(3,100),F(6,5),F(1,5)),
      racc_cfg=RACC.Config(),nominal_racc_std=(1,1,1),
      angular=angular,Qbase=M.zeros(6,6),ou=ou,bias=bias,qaxis=qaxis,
      accel_conditioning=SENSOR.AccelConditioning(0,(0,0,0)),
      accel_ldlt=MR.SafeLDLT(False,False,1,F(1,10**7)),tilt_reset_due=False,
      wpe_cfg=W.WPEConfig(1,4,F(1,2),1,180),wpe_decay=W.ExpWitness(1),
      stats_cfg=B.StatsConfig(4,F(3,10),60,F(1,20),5),
      band_decay=B.BandDecayWitness(1,1),variance_decay=B.VarianceDecayWitness(1),
      bench_noise_sigma=0,noise_sqrt=B.NoiseSqrtWitness(0),
      tracker_lpf_decay=FRONT.LPFDecayWitness(1),still_cfg=SF.Config(),
      still_attenuation=SF.AttenuationWitness(1),candidate_cfg=LP.candidate_cfg(),
      sigma_wave_sqrt=1,spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(1,1))
    kw.update(overrides)
    return X.step(s,raw,segment,**kw)


class Tests(unittest.TestCase):
    def test_first_Live_prefix_really_starts_from_startup_fresh_H18(self):
        s=startup(); out=run(s)
        self.assertIs(out.startup,s)
        self.assertEqual(s.state.mekf.mode,'H')
        self.assertEqual(out.first_live.state.mekf.mode,'H')
        self.assertEqual(out.first_live.prediction.state.reference,s.state.mekf.reference)
        self.assertEqual(out.first_live.state.tuner.sample_index,s.state.tuner.sample_index+1)
        self.assertEqual(out.first_live.state.tuner.time,out.first_live.state.mekf.reference.time)
        self.assertFalse(out.first_live.accelerometer.accepted)

    def test_frontend_memory_is_not_restarted_on_first_Live_sample(self):
        s=startup(); before=s.state.tuner
        out=run(s).first_live
        self.assertEqual(out.vertical.state.elapsed,before.vertical.elapsed+F(1,200))
        self.assertEqual(out.tuner_suffix.state.sample_index,before.sample_index+1)
        self.assertEqual(out.tuner_suffix.stage_before,'Live')

    def test_detached_first_packet_or_segment_fails_closed(self):
        s=startup(); raw,segment,angular,ou,bias,qaxis=physical_step(s)
        other=replace(raw,physical=segment.after)
        with self.assertRaisesRegex(ValueError,'detached'):
            X.step(s,other,segment)
        bad=replace(segment,before=segment.after)
        with self.assertRaisesRegex((ValueError,TypeError),'detached|physical'):
            X.step(s,raw,bad)

    def test_readiness_keeps_async_reset_precision_and_indefinite_word_open(self):
        r=X.readiness()
        self.assertTrue(r['synthetic_Live_root_removed_at_first_sample'])
        self.assertTrue(r['first_prediction_S_accel_tuner_WPE_prefix_composed_from_startup'])
        for k in ('async_mag_at_first_prefix_composed','firing_tilt_reset_at_first_prefix_composed',
                  'deployment_finite_precision_closed','source_uniform_indefinite_continuation_closed',
                  'complete_word_finite_identity','ALT_STARTUP_PASS','ALT_LIVE_PASS'):
            self.assertFalse(r[k])

if __name__=='__main__': unittest.main()
