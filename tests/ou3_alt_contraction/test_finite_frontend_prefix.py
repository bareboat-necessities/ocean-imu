"""Temporal raw-IMU -> frontend prefix regressions; not source qualification."""
from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_frontend_prefix as X
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as RAW
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as B
from tools.stability.ou3_alt_contraction import finite_frontend_runtime as FRONT
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as S
import test_finite_core as FC

G=F(196133,20000)
DT=F(1,200)


def configs():
    return (
        V.Config(0,0,G,20),
        W.WPEConfig(1,4,F(1,2),1,180),
        B.BandConfig(F(1,2),4,F(1,100),6,F(3,100),F(6,5),F(1,5)),
        B.StatsConfig(4,F(3,10),60,F(1,20),5),
        S.Config(),
    )


def runtime_state():
    return X.State(V.State(initialized=True),W.WPEState(),B.BandState(),B.StatsState(),FRONT.LPFState(),S.State())


def packet():
    state=FC.root('A')
    phys=replace(state.reference,acceleration=(0,0,0),beta=(0,0,0),gyro_bias=(0,0,0))
    z=list(state.z); z[21:24]=[F(0),F(0),F(0)]
    state=replace(state,z=tuple(z),reference=phys)
    sample=RAW.RawImuSample(phys,(0,0,0),(0,0,0),(0,0,0),(0,0,0),(0,0,-G),(0,0,G))
    return state,sample


def do_step(state,sample):
    vc,wc,bc,sc,stc=configs()
    return X.step(
        state,sample,dt=DT,
        vertical_cfg=vc,accel_invnorm=V.InvSqrtWitness(G*G,1/G),
        quat_invnorm=V.InvSqrtWitness(1,1),seed=None,
        wpe_cfg=wc,wpe_decay=W.ExpWitness(1),
        band_cfg=bc,stats_cfg=sc,band_decay=B.BandDecayWitness(1,1),
        variance_decay=B.VarianceDecayWitness(1),bench_noise_sigma=F(3,25),
        noise_sqrt=B.NoiseSqrtWitness(0),
        tracker_lpf_decay=FRONT.LPFDecayWitness(1),
        tracker=FRONT.TrackerOutputWitness(F(3,10)),still_cfg=stc,
        still_attenuation=S.AttenuationWitness(1))


class Tests(unittest.TestCase):
    def test_one_raw_packet_drives_all_frontend_consumers(self):
        mekf,sample=packet(); out=do_step(runtime_state(),sample)
        self.assertIs(out.raw_sample,sample)
        self.assertEqual(out.vertical.vertical_accel,0)
        self.assertEqual(out.wpe.state.accel_prev,0)
        self.assertEqual(out.band.filtered_accel,0)
        self.assertEqual(out.tracker_lpf.output,0)
        self.assertEqual(out.band.external_tuner_frequency,F(1,5))
        self.assertFalse(out.wpe.state.usable_period)
        self.assertEqual(X.assert_prediction_packet(out,mekf,sample.raw_gyro_body),sample.internal_gyro)
        self.assertEqual(X.assert_measurement_packet(out,sample.raw_accel_body),sample.internal_accel)

    def test_two_samples_compose_without_frontend_state_restart(self):
        _,sample=packet(); first=do_step(runtime_state(),sample); second=do_step(first.state,sample)
        self.assertEqual(first.state.sample_index,1); self.assertEqual(second.state.sample_index,2)
        self.assertEqual(first.state.time,DT); self.assertEqual(second.state.time,2*DT)
        self.assertEqual(second.state.wpe.elapsed,2*DT)
        self.assertEqual(second.state.vertical.elapsed,2*DT)
        self.assertTrue(second.state.tracker_lpf.initialized)
        self.assertGreater(second.state.stillness.still_time,first.state.stillness.still_time)
        self.assertEqual(second.band.external_tuner_frequency,F(1,5))

    def test_detached_MEKF_packet_is_rejected_after_frontend_consumed_sample(self):
        mekf,sample=packet(); out=do_step(runtime_state(),sample)
        bad=list(sample.raw_gyro_body); bad[0]+=1
        with self.assertRaises(ValueError): X.assert_prediction_packet(out,mekf,bad)
        bad_a=list(sample.raw_accel_body); bad_a[0]+=1
        with self.assertRaises(ValueError): X.assert_measurement_packet(out,bad_a)

    def test_readiness_remains_fail_closed(self):
        r=X.readiness()
        self.assertTrue(r['same_raw_packet_private_vertical_and_MEKF_bindable'])
        self.assertTrue(r['successive_frontend_samples_compose_without_state_restart'])
        self.assertFalse(r['tracker_algorithm_attached'])
        self.assertFalse(r['sensor_residual_source_bounds_attached'])
        self.assertFalse(r['complete_word_finite_identity'])
        self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
