"""Startup magnetic gravity-gate to tuner composition regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_startup_word as X
from tools.stability.ou3_alt_contraction import finite_mag_gravity_gate as G
from tools.stability.ou3_alt_contraction import finite_mag_startup_prefix as P
from tools.stability.ou3_alt_contraction import finite_mag_tuner_default as T
from tools.stability.ou3_alt_contraction import finite_mag_tilt_frame as TF
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V


def groot(v): return G.SqrtWitness(v*v,v)
def troot(v): return TF.SqrtWitness(v*v,v)
def mroot(v): return T.SqrtWitness(v*v,v)
def proxy(): return V.State((1,0,0,0),(0,0,0),True,F(2),0)


class Tests(unittest.TestCase):
    def test_IMU_event_advances_gate_only_and_uses_same_proxy_quaternion(self):
        s=X.State(G.State(),P.State(proxy()))
        out=X.imu_gate_step(s,G.Config(world_warmup=5),acc_body=(0,0,-10),gyro_body=(0,0,0),dt=F(1,200),
                            lpf_exp=None,gyro_norm=groot(0))
        self.assertEqual(out.state.mag,s.mag)
        self.assertEqual(out.state.gate.lpf,(0,0,-10))
        self.assertTrue(out.state.gate.aligned_branch)

    def test_nonadmitted_updateMag_cannot_advance_tuner_or_last_mag_clock(self):
        s=X.State(G.State(),P.State(proxy()))
        packet=P.Packet(6,(3,4,0),'pre-delay')
        out=X.update_mag_call(s,G.Config(mag_delay=7),T.Config(min_samples=10,min_window=0),packet,
                              begun=True,have_last_imu=True,sample_dt=F(1,200))
        self.assertFalse(out.admission.admitted); self.assertIsNone(out.magnetic)
        self.assertEqual(out.state.mag,s.mag); self.assertIsNone(out.state.mag.last_mag_time)
        with self.assertRaisesRegex(ValueError,'nonadmitted'):
            X.update_mag_call(s,G.Config(mag_delay=7),T.Config(),packet,begun=True,have_last_imu=True,
                              sample_dt=F(1,200),boat_q_norm=troot(1))

    def test_gravity_trusted_admission_advances_tuner_and_clock(self):
        gate=G.State((0,0,-10),True,5,2,True,0)
        s=X.State(gate,P.State(proxy()))
        packet=P.Packet(7,(3,4,0),'admitted')
        out=X.update_mag_call(s,G.Config(mag_delay=0,hold_sec=2),T.Config(min_samples=10,min_window=0),packet,
                              begun=True,have_last_imu=True,sample_dt=F(1,200),
                              boat_q_norm=troot(1),yaw_half=None,mag_norm=mroot(5))
        self.assertTrue(out.admission.admitted); self.assertIsNotNone(out.magnetic)
        self.assertEqual(out.state.mag.last_mag_time,7)
        self.assertEqual(out.state.mag.tuner.accumulator.accepted_count,1)
        self.assertEqual(out.state.gate,out.admission.state)

    def test_fallback_admission_uses_same_latched_t0(self):
        gate=G.State(eligible_t0=5)
        s=X.State(gate,P.State(proxy()))
        packet=P.Packet(35,(3,4,0),'fallback')
        out=X.update_mag_call(s,G.Config(mag_delay=0,fallback_sec=30),T.Config(min_samples=10,min_window=0),packet,
                              begun=True,have_last_imu=True,sample_dt=F(1,200),
                              boat_q_norm=troot(1),yaw_half=None,mag_norm=mroot(5))
        self.assertTrue(out.admission.fallback_ok); self.assertTrue(out.admission.admitted)
        self.assertEqual(out.state.gate.eligible_t0,5)

    def test_no_last_imu_blocks_tuner_even_when_gravity_trusted(self):
        gate=G.State(gravity_good=3,aligned_branch=True,eligible_t0=0)
        s=X.State(gate,P.State(proxy()))
        out=X.update_mag_call(s,G.Config(mag_delay=0),T.Config(),P.Packet(10,(3,4,0),'noimu'),
                              begun=True,have_last_imu=False,sample_dt=F(1,200))
        self.assertFalse(out.admission.admitted); self.assertEqual(out.admission.reason,'no_last_imu')
        self.assertEqual(out.state.mag,s.mag)

    def test_readiness_keeps_source_schedule_and_roundoff_open(self):
        r=X.readiness()
        self.assertTrue(r['admitted_updateMag_composes_literal_wrapper_gate_to_tuner'])
        self.assertFalse(r['startup_raw_mag_physical_source_relation_attached'])
        self.assertFalse(r['startup_mag_call_schedule_source_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_STARTUP_PASS'])

if __name__=='__main__': unittest.main()
