"""Held accelerometer post-prediction ancestry regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_held_accel_runtime as H
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as S
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as G
from tools.stability.ou3_alt_contraction import finite_racc_runtime as RACC
from tools.stability.ou3_alt_contraction import finite_measurement_runtime as MR
from tools.stability.ou3_alt_contraction.finite_tuner_commit import TuneState
import test_finite_core as FC

GRAV=F(196133,20000)


def packet_and_segment():
    state=FC.root('A'); segment,kw=FC.physical_successor(state); b=segment.before
    gyro=tuple(b.gyro_bias)
    inertial=[b.acceleration[i]-(GRAV if i==2 else 0) for i in range(3)]
    f=S.q_rotate(b.q_world_to_body,inertial)
    acc=tuple(f[i]+b.beta[i] for i in range(3))
    raw=S.RawImuSample(b,(0,0,0),(0,0,0),(0,0,0),gyro,acc,(0,0,GRAV))
    guarded=S.guarded_sample(raw,G.State(),G.Config(cutoff_hz=0),dt=segment.h)
    return state,segment,kw,guarded


class Tests(unittest.TestCase):
    def test_nontrivial_segment_retains_physical_hold_forcing(self):
        _,segment,_,guarded=packet_and_segment()
        out=H.observation(guarded,segment,S.AccelConditioning(0,(0,0,0)))
        self.assertEqual(out.guard_delta,(0,0,0))
        self.assertNotEqual(out.physical_hold_forcing,(0,0,0))
        self.assertEqual(out.nu_acc,out.physical_hold_forcing)
        after=segment.after
        inertial=[after.acceleration[i]-(GRAV if i==2 else 0) for i in range(3)]
        fa=S.q_rotate(after.q_world_to_body,inertial)
        self.assertEqual(out.observed,tuple(fa[i]+after.beta[i]+out.nu_acc[i] for i in range(3)))

    def test_guard_and_temperature_terms_are_retained_separately(self):
        _,segment,_,guarded=packet_and_segment()
        # Use an exact artificial conditioned descendant from the same raw packet.
        gr=G.Result(guarded.guard.state,
                    tuple(x+F(1,100) for x in guarded.raw.raw_accel_body),
                    guarded.guard.low_passed,guarded.guard.detector_high_passed,
                    guarded.guard.removed_rms,guarded.guard.excess_rms,
                    guarded.guard.target,guarded.guard.raw_weight)
        changed=S.GuardedImuSample(guarded.raw,gr)
        cond=S.AccelConditioning(2,(F(1,100),0,0))
        out=H.observation(changed,segment,cond)
        self.assertNotEqual(out.guard_delta,(0,0,0))
        self.assertEqual(out.nu_acc[0],out.physical_hold_forcing[0]+out.guard_delta[0]-F(1,50))

    def test_shipping_measurement_requires_exact_segment_endpoints(self):
        state,segment,kw,guarded=packet_and_segment()
        predicted=FC.C.prediction(state,segment,**kw)
        rr=RACC.step(RACC.State(False,(1,1,1)),RACC.Config(),guarded.guard,
                     nominal_std=(1,1,1),tune=TuneState(1,1,1),
                     preupdate_frequency=F(1,5),live=True)
        reject=MR.SafeLDLT(False,False,1,F(1,10**7))
        out=MR.accelerometer_from_held_guarded_racc(predicted,segment,guarded,
                S.AccelConditioning(0,(0,0,0)),rr,ldlt=reject)
        self.assertFalse(out.accepted); self.assertIs(out.state,predicted)
        from dataclasses import replace
        with self.assertRaises(ValueError):
            MR.accelerometer_from_held_guarded_racc(predicted,replace(segment,before=segment.after),guarded,
                S.AccelConditioning(0,(0,0,0)),rr,ldlt=reject)

    def test_readiness_keeps_source_bound_open(self):
        r=H.readiness()
        self.assertTrue(r['pre_to_post_physical_specific_force_hold_term_retained'])
        self.assertTrue(r['pre_to_post_true_bias_change_retained'])
        self.assertFalse(r['held_accel_COMPLETE_BRMM_bound_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
