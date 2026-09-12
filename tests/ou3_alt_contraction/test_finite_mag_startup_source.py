"""Startup raw magnetic physical-source regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_startup_source as X


class Tests(unittest.TestCase):
    def test_raw_packet_is_true_rotated_field_plus_body_hard_iron_plus_residual(self):
        p=X.PhysicalEndpoint(7,(1,0,0,0),'physical-startup-1')
        m=X.Model((20,0,40),(1,2,3),'model-1')
        s=X.make_sample(p,m,(F(1,10),F(-1,5),F(3,10)),'mag-1')
        self.assertEqual(s.predicted_body_field,(20,0,40))
        self.assertEqual(s.raw_body,(F(211,10),F(9,5),F(433,10)))
        packet=s.packet(); self.assertEqual(packet.wrapper_time,7); self.assertEqual(packet.raw_body,s.raw_body)

    def test_nontrivial_true_attitude_rotates_same_world_field(self):
        # 180 deg about z: x,y reverse while z stays fixed.
        p=X.PhysicalEndpoint(8,(0,0,0,1),'physical-startup-1')
        m=X.Model((20,5,40),(0,0,0),'model-1')
        s=X.make_sample(p,m,(0,0,0),'mag-2')
        self.assertEqual(s.raw_body,(-20,-5,40))

    def test_detached_raw_packet_is_rejected(self):
        p=X.PhysicalEndpoint(7,(1,0,0,0),'physical-startup-1')
        m=X.Model((20,0,40),(1,0,0),'model-1')
        with self.assertRaisesRegex(ValueError,'detached'):
            X.Sample(p,m,(0,0,0),(20,0,40),'bad')

    def test_field_hard_iron_and_physical_roots_cannot_restart_mid_history(self):
        p1=X.PhysicalEndpoint(7,(1,0,0,0),'hist')
        p2=X.PhysicalEndpoint(8,(1,0,0,0),'hist')
        model=X.Model((20,0,40),(1,0,0),'root')
        a=X.make_sample(p1,model,(0,0,0),'a'); b=X.make_sample(p2,model,(0,0,0),'b')
        self.assertTrue(X.assert_same_model_root((a,b)))
        changed=X.make_sample(p2,X.Model((21,0,40),(1,0,0),'root'),(0,0,0),'changed')
        with self.assertRaisesRegex(ValueError,'model restarted'):
            X.assert_same_model_root((a,changed))
        other=X.make_sample(X.PhysicalEndpoint(8,(1,0,0,0),'other'),model,(0,0,0),'other')
        with self.assertRaisesRegex(ValueError,'history restarted'):
            X.assert_same_model_root((a,other))

    def test_readiness_keeps_numerical_source_bounds_open(self):
        r=X.readiness()
        self.assertTrue(r['startup_stream_is_uncorrected_before_MagAutoTuner'])
        self.assertTrue(r['world_field_and_body_hard_iron_have_persistent_model_root'])
        self.assertFalse(r['startup_mag_noise_bound_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_STARTUP_PASS'])

if __name__=='__main__': unittest.main()
