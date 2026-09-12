"""Startup magnetic yaw-stripped tilt-frame regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_tilt_frame as X
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P


def root(v): return X.SqrtWitness(v*v,v)


class Tests(unittest.TestCase):
    def test_pure_yaw_is_removed_exactly(self):
        # Half-yaw (3/5,4/5) => full heading cos=-7/25, sin=24/25.
        q=(F(3,5),0,0,F(4,5))
        w=X.YawHalfWitness(F(-7,25),F(24,25),root(1),F(3,5),F(4,5))
        out=X.yaw_removed(q,q_norm=root(1),yaw_half=w)
        self.assertEqual(out.heading_c,F(-7,25)); self.assertEqual(out.heading_s,F(24,25))
        self.assertEqual(out.q_tilt,(1,0,0,0))

    def test_yaw_removed_from_roll_tilt_preserves_tilt(self):
        q_yaw=(F(3,5),0,0,F(4,5))
        q_roll=(F(5,13),F(12,13),0,0)
        q=tuple(P.quat_mul(q_yaw,q_roll))
        w=X.YawHalfWitness(F(-7,25),F(24,25),root(1),F(3,5),F(4,5))
        out=X.yaw_removed(q,q_norm=root(1),yaw_half=w)
        self.assertEqual(out.q_tilt,q_roll)

    def test_witness_must_belong_to_same_normalized_quaternion(self):
        q=(F(3,5),0,0,F(4,5))
        bad=X.YawHalfWitness(1,0,root(1),1,0)
        with self.assertRaisesRegex(ValueError,'detached'):
            X.yaw_removed(q,q_norm=root(1),yaw_half=bad)
        with self.assertRaisesRegex(ValueError,'norm witness detached'):
            X.yaw_removed(q,q_norm=X.SqrtWitness(4,2),yaw_half=None)

    def test_tiny_quaternion_matches_shipping_identity_fallback_structure(self):
        out=X.yaw_removed((0,0,0,0),q_norm=root(0),yaw_half=None)
        self.assertEqual(out.q_tilt,(1,0,0,0)); self.assertEqual(out.q_normalized,(1,0,0,0))

    def test_readiness_keeps_native_transcendentals_open(self):
        r=X.readiness()
        self.assertTrue(r['startup_tilt_quaternion_is_not_free_operand'])
        self.assertTrue(r['heading_c_R00_and_s_R10_from_same_normalized_quaternion'])
        self.assertFalse(r['atan2_AngleAxis_binary32_attached'])
        self.assertFalse(r['startup_mag_sensor_packet_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_STARTUP_PASS'])

if __name__=='__main__': unittest.main()
