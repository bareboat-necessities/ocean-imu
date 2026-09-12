"""Persistent active magnetic reference regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_reference_runtime as X
from tools.stability.ou3_alt_contraction import finite_mag_runtime as MAG
import test_finite_core as FC


class Tests(unittest.TestCase):
    def test_sample_inherits_active_model_exactly(self):
        s=FC.root('H'); model=MAG.Model((22,0,43),(F(1,5),F(1,4),F(1,3)))
        active=X.State(model,3,'one-mag-root')
        n=(F(1,100),0,F(-1,100)); p=MAG.SENSOR.q_rotate(s.reference.q_world_to_body,model.world_reference)
        packet=X.sample(active,s.reference,tuple(p[i]+n[i] for i in range(3)),n)
        self.assertEqual(packet.model,active.model); self.assertIs(X.require_same(active,packet),packet)

    def test_world_reference_or_sigma_change_requires_new_active_generation(self):
        s=FC.root('H'); model=MAG.Model((22,0,43),(1,1,1)); active=X.State(model,0,'mag-root')
        for detached in (MAG.Model((23,0,43),(1,1,1)), MAG.Model((22,0,43),(2,1,1))):
            n=(0,0,0); p=MAG.SENSOR.q_rotate(s.reference.q_world_to_body,detached.world_reference)
            packet=MAG.Sample(s.reference,tuple(p),n,detached)
            with self.assertRaisesRegex(ValueError,'persistent active magnetic model'):
                X.require_same(active,packet)

    def test_generation_and_root_are_persistent_identity_metadata(self):
        model=MAG.Model((22,0,43),(1,1,1))
        with self.assertRaises(ValueError): X.State(model,-1,'root')
        with self.assertRaises(ValueError): X.State(model,0,'')
        a=X.State(model,0,'root'); b=X.State(model,1,'root')
        self.assertNotEqual(a,b)

    def test_readiness_keeps_reference_generation_sources_open(self):
        r=X.readiness()
        self.assertTrue(r['updateMag_cannot_choose_free_world_reference'])
        self.assertTrue(r['updateMag_cannot_choose_free_Rmag'])
        self.assertFalse(r['startup_mag_reference_generation_attached'])
        self.assertFalse(r['second_stage_mag_refinement_attached'])
        self.assertFalse(r['continuous_hard_iron_reference_update_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
