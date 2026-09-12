"""Startup MagAutoTuner-ready transition regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_startup_ready as X
from tools.stability.ou3_alt_contraction import finite_mag_tuner_default as T
from tools.stability.ou3_alt_contraction import finite_mag_accumulator as A
from tools.stability.ou3_alt_contraction import finite_mag_tilt_frame as TF


def troot(v): return TF.SqrtWitness(v*v,v)

def ready_state():
    # mean=(3,4,12), gauge-fixed reference=(5,0,12)
    acc=A.State((3,4,12),13,1,F(1,200),1)
    return T.State(acc,0,True,(3,4,12),(5,0,12),1,(3,4,12))


class Tests(unittest.TestCase):
    def test_ready_transition_writes_generation_zero_reference_and_pending_inverse_yaw(self):
        tuner=ready_state(); cfg=X.Config((F(1,5),F(1,4),F(1,3)),'startup-producer')
        # Direction (3,4), radial 5; half-angle cos/sin=(2/sqrt5,1/sqrt5) is irrational,
        # so use an exact 3-4-5 direction whose half angle is 2/sqrt5? Not rational.
        # Instead choose mean direction with full-angle cos=3/5,sin=4/5 and exact
        # rational half-angle from the 3-4-5 half-angle triple: cos_half=2/sqrt5
        # is not rational either. Build an alternative ready state below with
        # full-angle cos=-7/25,sin=24/25, half=(3/5,4/5).
        acc=A.State((F(-7),24,0),25,1,F(1,200),1)
        tuner=T.State(acc,0,True,(F(-7),24,0),(25,0,0),1,(F(-7),24,0))
        yaw=TF.YawHalfWitness(F(-7),24,troot(25),F(3,5),F(4,5))
        out=X.transition(tuner,cfg,yaw_half=yaw)
        self.assertEqual(out.active_reference.generation,0)
        self.assertEqual(out.active_reference.root,'startup-producer')
        self.assertEqual(out.active_reference.model.world_reference,(25,0,0))
        self.assertEqual(out.active_reference.model.sigma_internal,cfg.sigma_internal)
        self.assertEqual(out.pending_yaw.q_abs,(F(3,5),0,0,F(-4,5)))
        self.assertEqual(out.pending_yaw.mean_xy,(F(-7),24))

    def test_transition_requires_ready_tuner_and_same_mean_yaw_witness(self):
        cfg=X.Config((1,1,1),'root')
        with self.assertRaisesRegex(ValueError,'ready'):
            X.transition(T.State(),cfg,yaw_half=TF.YawHalfWitness(1,0,troot(1),1,0))
        acc=A.State((F(-7),24,0),25,1,F(1,200),1)
        tuner=T.State(acc,0,True,(F(-7),24,0),(25,0,0),1,(F(-7),24,0))
        bad=TF.YawHalfWitness(1,0,troot(1),1,0)
        with self.assertRaisesRegex(ValueError,'detached'):
            X.transition(tuner,cfg,yaw_half=bad)

    def test_reference_write_does_not_apply_yaw_to_any_filter_state(self):
        tuner=ready_state()
        # Use zero-yaw mean for an exact trivial pending lock.
        acc=A.State((5,0,12),13,1,F(1,200),1)
        tuner=T.State(acc,0,True,(5,0,12),(5,0,12),1,(5,0,12))
        out=X.transition(tuner,X.Config((1,1,1),'root'),
                         yaw_half=TF.YawHalfWitness(5,0,troot(5),1,0))
        self.assertIs(out.tuner_state,tuner)
        self.assertEqual(out.pending_yaw.q_abs,(1,0,0,0))
        self.assertEqual(out.active_reference.model.world_reference,tuner.world_reference)

    def test_readiness_keeps_transcendental_and_source_envelope_open(self):
        r=X.readiness()
        self.assertTrue(r['pending_absolute_yaw_tied_to_same_accepted_mean'])
        self.assertTrue(r['pending_yaw_not_applied_to_MEKF_before_handoff'])
        self.assertFalse(r['atan2_wrapPi_AngleAxis_binary32_attached'])
        self.assertFalse(r['magnetic_source_envelope_declared'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_STARTUP_PASS'])

if __name__=='__main__': unittest.main()
