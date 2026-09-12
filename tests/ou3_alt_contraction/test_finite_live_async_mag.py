"""Async updateMag composition regressions."""
from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_live_async_mag as X
from tools.stability.ou3_alt_contraction import finite_mag_runtime as MAG
from tools.stability.ou3_alt_contraction import finite_mag_bias_gate as GATE
from tools.stability.ou3_alt_contraction import finite_measurement_runtime as MR
import test_finite_core as FC


def mag_sample(state,*,zero=False):
    if zero:
        model=MAG.Model((0,0,0),(1,1,1)); return MAG.Sample(state.reference,(0,0,0),(0,0,0),model)
    model=MAG.Model((22,0,43),(1,1,1)); n=(F(1,100),0,F(-1,100))
    p=MAG.SENSOR.q_rotate(state.reference.q_world_to_body,model.world_reference)
    return MAG.Sample(state.reference,tuple(p[i]+n[i] for i in range(3)),n,model)


class Tests(unittest.TestCase):
    def test_delay_gate_consumes_no_packet_or_ldlt(self):
        s=FC.root('H'); st=X.State(s,GATE.State()); cfg=GATE.Config(mag_delay=7)
        out=X.update_mag_call(st,cfg,time=s.reference.time,live=True)
        self.assertFalse(out.wrapper_attempted); self.assertEqual(out.state,st)
        with self.assertRaisesRegex(ValueError,'consumes no magnetic proof operands'):
            X.update_mag_call(st,cfg,time=s.reference.time,live=True,sample=mag_sample(s))

    def test_sanity_rejection_still_increments_updateMag_count(self):
        s=FC.root('H'); st=X.State(s,GATE.State()); cfg=GATE.Config(mag_delay=0)
        out=X.update_mag_call(st,cfg,time=s.reference.time,live=True,sample=mag_sample(s,zero=True))
        self.assertTrue(out.wrapper_attempted); self.assertFalse(out.measurement_accepted)
        self.assertFalse(out.magnetic.attempted_measurement)
        self.assertEqual(out.state.control.updates,1); self.assertEqual(out.state.filter,s)

    def test_ldlt_rejection_still_increments_count(self):
        s=FC.root('H'); st=X.State(s,GATE.State()); cfg=GATE.Config(mag_delay=0)
        out=X.update_mag_call(st,cfg,time=s.reference.time,live=True,sample=mag_sample(s),
                              ldlt=MR.SafeLDLT(False,False,1,F(1,10**7)))
        self.assertTrue(out.magnetic.attempted_measurement); self.assertFalse(out.measurement_accepted)
        self.assertEqual(out.state.control.updates,1); self.assertEqual(out.state.filter,s)

    def test_unlock_and_H18_A21_release_do_not_require_mag_acceptance(self):
        s=FC.root('H'); ref=replace(s.reference,time=F(2)); s=replace(s,reference=ref)
        control=GATE.State(249,0,True,False); cfg=GATE.Config(mag_delay=0,unlock_count=250)
        out=X.update_mag_call(X.State(s,control),cfg,time=2,live=True,sample=mag_sample(s),
                              ldlt=MR.SafeLDLT(False,False,1,F(1,10**7)))
        self.assertFalse(out.measurement_accepted)
        self.assertTrue(out.gate.unlocked_now); self.assertFalse(out.state.control.locked)
        self.assertEqual(out.state.filter.mode,'A')

    def test_accepted_event_keeps_same_physical_endpoint(self):
        s=FC.root('A'); st=X.State(s,GATE.State()); cfg=GATE.Config(mag_delay=0)
        out=X.update_mag_call(st,cfg,time=s.reference.time,live=True,sample=mag_sample(s),
                              ldlt=MR.SafeLDLT(True,None,1,F(1,10**7)))
        self.assertTrue(out.measurement_accepted)
        self.assertEqual(out.state.filter.reference,s.reference)
        self.assertNotEqual(out.state.filter,s)

    def test_clock_and_packet_ancestry_fail_closed(self):
        s=FC.root('H'); cfg=GATE.Config(mag_delay=0); st=X.State(s,GATE.State())
        with self.assertRaisesRegex(ValueError,'clock detached'):
            X.update_mag_call(st,cfg,time=s.reference.time+F(1,200),live=True,sample=mag_sample(s),
                              ldlt=MR.SafeLDLT(False,False,1,F(1,10**7)))
        other=FC.root('H','BIAS2')
        with self.assertRaisesRegex(ValueError,'detached'):
            X.update_mag_call(st,cfg,time=s.reference.time,live=True,sample=mag_sample(other),
                              ldlt=MR.SafeLDLT(False,False,1,F(1,10**7)))

    def test_readiness_remains_fail_closed(self):
        r=X.readiness()
        self.assertTrue(r['H18_to_A21_release_independent_of_mag_acceptance'])
        self.assertFalse(r['async_mag_call_schedule_source_attached'])
        self.assertFalse(r['mag_noise_source_bound_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
