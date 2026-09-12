"""Startup magnetic bridge to main finite physical endpoint regressions."""
from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_startup_physical_bridge as X
from tools.stability.ou3_alt_contraction import finite_mag_startup_source as M
from tools.stability.ou3_alt_contraction import finite_physical_prediction as P


def physical():
    return P.PhysicalKinematics(
        time=F(7),live_origin=F(0),q_world_to_body=(1,0,0,0),
        velocity=(0,0,0),position=(0,0,0),centered_S=(0,0,0),
        acceleration=(0,0,0),gyro_bias=(0,0,0),beta=(0,0,0))


class Tests(unittest.TestCase):
    def test_endpoint_time_and_true_attitude_are_derived_from_same_finite_physical_object(self):
        p=physical(); e=X.endpoint(p,'hist')
        self.assertEqual(e.time,p.time); self.assertEqual(e.q_world_to_body,p.q_world_to_body)
        self.assertEqual(e.history_id,'hist')

    def test_source_sample_is_constructed_from_main_endpoint_not_free_duplicate(self):
        p=physical(); model=M.Model((3,4,0),(0,0,0),'model')
        s=X.sample(p,'hist',model,(0,0,0),'mag-1')
        self.assertTrue(X.assert_same_endpoint(p,s,'hist'))
        self.assertEqual(s.raw_body,(3,4,0))

    def test_time_attitude_or_history_detachment_fails(self):
        p=physical(); model=M.Model((3,4,0),(0,0,0),'model')
        s=X.sample(p,'hist',model,(0,0,0),'mag-1')
        p_time=replace(p,time=p.time+F(1,200))
        with self.assertRaisesRegex(ValueError,'detached'):
            X.assert_same_endpoint(p_time,s,'hist')
        p_att=replace(p,q_world_to_body=(0,0,0,1))
        with self.assertRaisesRegex(ValueError,'detached'):
            X.assert_same_endpoint(p_att,s,'hist')
        with self.assertRaisesRegex(ValueError,'detached'):
            X.assert_same_endpoint(p,s,'other')

    def test_object_bridge_is_not_BRMM_admission(self):
        r=X.readiness()
        self.assertTrue(r['duplicate_free_physical_endpoint_bridge'])
        self.assertFalse(r['COMPLETE_BRMM_admission_proved_by_object_identity'])
        self.assertFalse(r['magnetic_source_envelope_declared'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_STARTUP_PASS'])

if __name__=='__main__': unittest.main()
