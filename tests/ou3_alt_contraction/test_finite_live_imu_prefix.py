"""Composed default-policy Live IMU prefix regressions; not stability admission."""
from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_live_imu_prefix as X
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as GUARD
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
from tools.stability.ou3_alt_contraction import finite_tuner_frontend_prefix as T
from tools.stability.ou3_alt_contraction import finite_tuner_commit as TC
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import finite_post_prediction as POST
from tools.stability.ou3_alt_contraction import finite_racc_runtime as RACC
from tools.stability.ou3_alt_contraction import finite_periodic_aw_sync as AWSYNC
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as B
from tools.stability.ou3_alt_contraction import finite_frontend_runtime as FRONT
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as SF
from tools.stability.ou3_alt_contraction import finite_tuner_stillness_projection as SP
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W
from tools.stability.ou3_alt_contraction import finite_measurement_runtime as MR
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_runtime as PR
from tools.stability.ou3_alt_contraction import finite_qaxis_runtime as Q
from tools.stability.ou3_alt_contraction import finite_attitude_runtime as A
from tools.stability.ou3_alt_contraction import finite_ou_runtime_primitives as O
import test_finite_core as FC
import test_finite_runtime_parameters as RP

G=F(196133,20000); PASS=Q.PSDWitness(True)


def candidate_cfg():
    return C.CandidateConfig(F(1,10),2,F(2,5),1,F(1,10),2,2,1,F(1,100),2,
                             F(1,10),2,1,F(1,2),1,1,F(2,5),F(3,2),0,F(1,10),True)


def commit_cfg():
    return TC.CommitConfig(F(3,2),F(1,100),2,F(3,200),True,True,F(2,5),5,F(4,5),F(6,5),F(7,5))


def fixture():
    base=FC.root('A')
    z=list(base.z)
    for i in range(3): z[3+i]=base.reference.gyro_bias[i]
    mekf=replace(base,z=tuple(z))
    segment,_=FC.physical_successor(mekf)
    h=segment.h

    # Same raw packet roots prediction at the physical predecessor.  Choosing
    # e_bg=b_g makes the estimator's corrected angular rate exactly zero while
    # retaining the raw sensor identity through omega_sample=-b_g.
    omega=tuple(-x for x in mekf.reference.gyro_bias)
    inertial=tuple(mekf.reference.acceleration[i]-(0,0,G)[i] for i in range(3))
    fbody=SENSOR.q_rotate(mekf.reference.q_world_to_body,inertial)
    acc=tuple(fbody[i]+mekf.reference.beta[i] for i in range(3))
    raw=SENSOR.RawImuSample(mekf.reference,omega,(0,0,0),(0,0,0),(0,0,0),acc,(0,0,G))

    # A valid persistent guard state can make acc_in a simple level packet; the
    # difference from raw remains in the exact held-sample guard_delta term.
    level=(F(0),F(0),-G); zero=(F(0),F(0),F(0))
    guard=GUARD.State(stages=(level,level,level,level),detect_stages=(zero,zero),
                      removed_ms=zero,weight=1,initialized=True)

    band=B.BandState(band=0,p00=0,p01=0,p11=0,ready=True)
    stats=B.StatsState(frequency=F(1,2),mean_value=0,mean_weight=1,sq_value=1,sq_weight=1)
    tuner=T.State(V.State(initialized=True),W.WPEState(),band,stats,FRONT.LPFState(),
                  SP.State(),TC.TuneState(1,1,1),last_adapt_time=0,pending=False,
                  sample_index=7,time=mekf.reference.time,stage='Live')
    active=RP.committed(True)
    state=X.State(mekf,guard,tuner,RACC.State(False,(1,1,1)),active,
                  POST.Scheduler(active.pseudo_period,0),AWSYNC.State(False,mekf.reference.time,None))

    angular=A.AngularRuntime((0,0,0),h)
    ou=O.OUDecay(h,active.tau,F(199,200))
    bias=O.BiasDecay(True,active.tau,F(199,200),M.scaled(M.eye(3),F(1,100000)))
    qaxis=PR.QAxisBranch(False,active.Sigma_aw,(PASS,PASS,PASS),(PASS,PASS,PASS),F(1,10**7))
    return state,raw,segment,angular,ou,bias,qaxis


def run(state=None,**overrides):
    if state is None: state,raw,segment,angular,ou,bias,qaxis=fixture()
    else:
        _,raw,segment,angular,ou,bias,qaxis=fixture()
        # Preserve caller state but root raw/segment in its MEKF predecessor only
        # for tests that do not alter the physical reference.
    kw=dict(
      dt=segment.h,commit_cfg=commit_cfg(),boundary_bench_noise_sigma=0,
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
      still_attenuation=SF.AttenuationWitness(1),candidate_cfg=candidate_cfg(),
      sigma_wave_sqrt=1,spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(1,1))
    kw.update(overrides)
    return X.step(state,raw,segment,**kw)


class Tests(unittest.TestCase):
    def test_one_live_sample_composes_same_raw_guarded_prediction_and_tuner_history(self):
        state,_,segment,_,_,_,_=fixture(); out=run(state)
        self.assertEqual(out.state.mekf.reference,segment.after)
        self.assertEqual(out.state.tuner.time,segment.after.time)
        self.assertEqual(out.state.tuner.sample_index,state.tuner.sample_index+1)
        self.assertIs(out.tuner_suffix.raw_sample,out.guarded)
        self.assertEqual(out.vertical.vertical_accel,0)
        self.assertFalse(out.S_service.S_service_due if hasattr(out.S_service,'S_service_due') else out.S_service.prefix.S_service_due)
        self.assertIsNone(out.S_service.measurement)
        self.assertFalse(out.accelerometer.accepted)
        self.assertEqual(out.state.scheduler.period,state.active.pseudo_period)
        self.assertFalse(out.state.aw_sync.pending)

    def test_guard_delta_and_held_physical_evolution_are_both_nontrivial(self):
        state,_,segment,_,_,_,_=fixture(); out=run(state)
        from tools.stability.ou3_alt_contraction import finite_held_accel_runtime as H
        held=H.observation(out.guarded,segment,SENSOR.AccelConditioning(0,(0,0,0)))
        self.assertNotEqual(held.guard_delta,(0,0,0))
        self.assertNotEqual(held.physical_hold_forcing,(0,0,0))

    def test_current_sample_does_not_commit_candidate_or_retarget_period(self):
        state,_,_,_,_,_,_=fixture(); out=run(state)
        self.assertIsNone(out.boundary.commit)
        self.assertFalse(out.tuner_suffix.state.pending)
        self.assertEqual(out.state.active,state.active)
        self.assertEqual(out.state.scheduler.period,state.scheduler.period)

    def test_firing_tilt_watchdog_fails_closed_until_reset_branch_is_composed(self):
        state,_,_,_,_,_,_=fixture()
        with self.assertRaises(NotImplementedError): run(state,tilt_reset_due=True)

    def test_readiness_does_not_promote_live_or_end_to_end_theorem(self):
        r=X.readiness()
        self.assertTrue(r['one_guarded_accel_feeds_Mahony_and_postprediction_accel_update'])
        self.assertTrue(r['held_accel_physical_evolution_forcing_retained'])
        self.assertTrue(r['post_Mahony_MEKF_interleave_before_tuner_suffix'])
        for k in ('live_tilt_reset_branch_attached','async_magnetometer_branch_attached',
                  'source_uniform_COMPLETE_BRMM_bounds_attached','deployment_finite_precision_closed',
                  'complete_word_finite_identity','ALT_LIVE_PASS'):
            self.assertFalse(r[k])

if __name__=='__main__': unittest.main()
