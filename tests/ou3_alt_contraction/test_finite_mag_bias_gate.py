"""Asynchronous magnetometer control / BA release regressions."""
from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_bias_gate as X
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
import test_finite_core as FC


def with_ba_cross(state):
    cov=[list(r) for r in state.covariance]
    cov[18][0]=cov[0][18]=F(1,7)
    cov[19][7]=cov[7][19]=F(1,9)
    cov[20][15]=cov[15][20]=F(1,11)
    return replace(state,covariance=tuple(tuple(r) for r in cov))


class Tests(unittest.TestCase):
    def test_delay_or_disabled_gate_consumes_no_measurement_successor(self):
        s=FC.root('H'); c=X.State(); cfg=X.Config(with_mag=True,mag_delay=7)
        out=X.update_mag_call(c,s,cfg,time=F(699,100),live=True)
        self.assertFalse(out.attempted); self.assertEqual(out.state,c); self.assertEqual(out.filter_state,s)
        with self.assertRaisesRegex(ValueError,'gated-off'):
            X.update_mag_call(c,s,cfg,time=1,live=True,measurement_state=s)
        out=X.update_mag_call(c,s,X.Config(with_mag=False),time=100,live=True)
        self.assertFalse(out.attempted)

    def test_attempt_count_is_independent_of_measurement_acceptance(self):
        s=FC.root('H'); c=X.State(); cfg=X.Config(mag_delay=7)
        # A rejected MEKF magnetic measurement returns the unchanged state.  The
        # wrapper still increments count and stamps first attempt after the call.
        out=X.update_mag_call(c,s,cfg,time=7,live=True,measurement_state=s)
        self.assertTrue(out.attempted); self.assertEqual(out.state.updates,1)
        self.assertEqual(out.state.first_time,7); self.assertEqual(out.filter_state,s)
        self.assertTrue(out.state.locked)

    def test_unlock_requires_count_and_strictly_more_than_one_second(self):
        s=FC.root('H'); cfg=X.Config(mag_delay=7,unlock_count=250)
        c=X.State(updates=249,first_time=7,locked=True,hold=False)
        # Equality at one second is not enough.
        equal=X.update_mag_call(c,s,cfg,time=8,live=True,measurement_state=s)
        self.assertFalse(equal.unlocked_now); self.assertTrue(equal.state.locked)
        # One later attempted call beyond the guard unlocks and performs H->A.
        later_control=replace(c,updates=249)
        later=X.update_mag_call(later_control,s,cfg,time=F(801,100),live=True,measurement_state=s)
        self.assertTrue(later.unlocked_now); self.assertFalse(later.state.locked)
        self.assertEqual(later.filter_state.mode,'A')
        floor=F(4,1000)**2
        for i in range(18,21): self.assertGreaterEqual(later.filter_state.covariance[i][i],floor)

    def test_external_hold_clears_lock_without_enabling_and_zeroes_BA_crosses(self):
        s=with_ba_cross(FC.root('A')); cfg=X.Config()
        held=X.set_hold(X.State(10,7,False,False),s,cfg,hold=True,live=True)
        self.assertEqual(held.filter_state.mode,'H'); self.assertTrue(held.state.hold)
        for i in range(18,21):
            for j in range(18): self.assertEqual(held.filter_state.covariance[i][j],0)
        released=X.set_hold(held.state,held.filter_state,cfg,hold=False,live=True)
        self.assertEqual(released.filter_state.mode,'A'); self.assertFalse(released.state.hold)

        # Unlock while an external hold exists clears only the wrapper lock.
        h=FC.root('H'); c=X.State(249,7,True,True)
        out=X.update_mag_call(c,h,cfg,time=F(801,100),live=True,measurement_state=h)
        self.assertTrue(out.unlocked_now); self.assertFalse(out.state.locked)
        self.assertEqual(out.filter_state.mode,'H'); self.assertFalse(out.enabled_bias_now)

    def test_readiness_keeps_mag_source_and_precision_open(self):
        r=X.readiness()
        self.assertTrue(r['attempt_count_independent_of_measurement_acceptance'])
        self.assertTrue(r['H_to_A_BA_variance_floor_materialized'])
        self.assertFalse(r['magnetic_measurement_same_source_packet_attached'])
        self.assertFalse(r['Rmag_runtime_ancestry_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
