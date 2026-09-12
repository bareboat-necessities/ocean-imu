"""Post-eligibility startup magnetic prefix regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_startup_prefix as X
from tools.stability.ou3_alt_contraction import finite_mag_tuner_default as T
from tools.stability.ou3_alt_contraction import finite_mag_tilt_frame as TF
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V


def troot(v): return TF.SqrtWitness(v*v,v)
def mroot(v): return T.SqrtWitness(v*v,v)

def proxy_yaw():
    return V.State((F(3,5),0,0,F(4,5)),(0,0,0),True,F(2),0)

def yaw_witness():
    return TF.YawHalfWitness(F(-7,25),F(24,25),troot(1),F(3,5),F(4,5))


class Tests(unittest.TestCase):
    def test_first_eligible_mag_uses_configured_dt_and_same_Mahony_quaternion(self):
        s=X.State(proxy_yaw())
        p=X.Packet(F(7),'raw-1' and (3,4,0),'raw-1')
        out=X.eligible_update(s,T.Config(min_samples=10,min_window=0),p,sample_dt=F(1,100),
                              boat_q_norm=troot(1),yaw_half=yaw_witness(),mag_norm=mroot(5))
        self.assertEqual(out.dt_mag,F(1,100)); self.assertEqual(out.state.last_mag_time,7)
        self.assertEqual(out.state.proxy,s.proxy)
        # Pure proxy yaw is removed before accumulation.
        self.assertEqual(out.state.tuner.last_world_sample,(3,4,0))

    def test_strictly_advancing_time_uses_elapsed_mag_interval(self):
        s=X.State(proxy_yaw(),last_mag_time=F(7))
        p=X.Packet(F(29,4),(3,4,0),'raw-2')
        out=X.eligible_update(s,T.Config(min_samples=10,min_window=0),p,sample_dt=F(1,200),
                              boat_q_norm=troot(1),yaw_half=yaw_witness(),mag_norm=mroot(5))
        self.assertEqual(out.dt_mag,F(1,4)); self.assertEqual(out.state.last_mag_time,F(29,4))

    def test_equal_or_reversed_time_uses_fallback_but_stores_current_time(self):
        for now in (F(7),F(27,4)):
            with self.subTest(now=now):
                s=X.State(proxy_yaw(),last_mag_time=F(7))
                p=X.Packet(now,(3,4,0),'raw-fallback')
                out=X.eligible_update(s,T.Config(min_samples=10,min_window=0),p,sample_dt=F(1,200),
                                      boat_q_norm=troot(1),yaw_half=yaw_witness(),mag_norm=mroot(5))
                self.assertEqual(out.dt_mag,F(1,200)); self.assertEqual(out.state.last_mag_time,now)

    def test_clock_update_survives_tuner_rejection(self):
        # Seed an accepted running norm of 5, then reject norm 7 as >35% outlier.
        first=X.eligible_update(X.State(proxy_yaw()),T.Config(min_samples=10,min_window=0),
              X.Packet(F(7),(3,4,0),'raw-a'),sample_dt=F(1,200),
              boat_q_norm=troot(1),yaw_half=yaw_witness(),mag_norm=mroot(5))
        second=X.eligible_update(first.state,T.Config(min_samples=10,min_window=0),
              X.Packet(F(8),(7,0,0),'raw-b'),sample_dt=F(1,200),
              boat_q_norm=troot(1),yaw_half=yaw_witness(),mag_norm=mroot(7))
        self.assertEqual(second.tuner_step.rejection,'running_norm_outlier')
        self.assertEqual(second.state.last_mag_time,8)
        self.assertEqual(second.state.tuner.accumulator,first.state.tuner.accumulator)
        self.assertEqual(second.state.proxy,first.state.proxy)

    def test_packet_id_and_body_are_persistent_event_ancestry(self):
        p=X.Packet(3,(1,2,3),'mag-42')
        self.assertEqual(p.packet_id,'mag-42'); self.assertEqual(p.raw_body,(1,2,3))
        with self.assertRaises(ValueError): X.Packet(3,(1,2,3),'')

    def test_readiness_keeps_outer_gate_and_physical_source_open(self):
        r=X.readiness()
        self.assertTrue(r['startup_mag_reads_persistent_Mahony_quaternion'])
        self.assertTrue(r['clock_update_survives_tuner_rejection'])
        self.assertFalse(r['gravity_alignment_eligibility_gate_attached'])
        self.assertFalse(r['startup_raw_mag_physical_source_relation_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_STARTUP_PASS'])

if __name__=='__main__': unittest.main()
