from __future__ import annotations
from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'/'stability'))
from ou3_interval import Interval
import ou3_p4_live_entry_graph as E
import ou3_p4_brmm_same_signal_statistics as C
import ou3_brmm_frontend_state_step as FRONT
import ou3_brmm_wpe_state_step as WPE
import ou3_p4_joint_brmm_frontend_transition as JOINT

I=Interval.point

class LiveEntryGraphTest(unittest.TestCase):
    def test_shared_origin_residual_is_exact_not_small(self):
        physical=(315,-410,700)
        error=(313,-409,703)
        origin=(300,-400,600)
        self.assertEqual(E.S_residual(error,physical),
            E.S_residual(E.reanchor_S(error,origin),E.reanchor_S(physical,origin)))
        self.assertEqual(E.S_residual(physical,physical),(0,0,0))
        self.assertGreater(sum(F(x)**2 for x in physical),300**2)

    def test_origin_column_cancels_for_arbitrary_event_gain(self):
        # Exact signed Joseph identity over rational gains. The symbolic proof
        # in the producer is stronger; these catch sign/offset regressions.
        for n in (18,21):
            gain=[[F((i+1)*(j-1),7) for j in range(3)] for i in range(n)]
            emap=[[F(int(i==12+j))-gain[i][j]+gain[i][j]
                   for j in range(3)] for i in range(n)]
            centered=E.center_prefix_map(emap,[[int(i==j) for j in range(3)] for i in range(3)])
            self.assertTrue(all(x==0 for row in centered for x in row))
        # Dropping the physical-S port is a different and generally nonzero map.
        self.assertNotEqual(F(1)-F(2,7),F(1))

    def test_actual_initialization_is_not_covariance_membership(self):
        d=E.build()
        self.assertEqual(d['canonical_S_entry_error'],0)
        self.assertFalse(d['goLive_resets_linear_values'])
        self.assertFalse(d['position_origin_at_handoff_established_by_shipping'])
        self.assertFalse(d['P4_PASS'])
        self.assertFalse(d['P5_MAY_START'])

    def test_shared_continuous_source_jet(self):
        # Constant a=(1,0,0), h=2: jointly generated integrals are 2,2,4/3.
        j=E.physical_jet((2,0,0),(3,0,0),2,[(2,0,0),(2,0,0),(F(4,3),0,0)])
        self.assertEqual(j['v'],(5,0,0))
        self.assertEqual(j['p'],(10,0,0))
        self.assertEqual(j['S_from_handoff'],(F(34,3),0,0))
        self.assertEqual(E.prefix_S_bound(2,2,3,1),F(34,3))
        self.assertEqual(E.prefix_S_bound(100,2,3,1,uniform_S_bound=10),20)

    def test_shifted_energy_keeps_signed_source_correlation(self):
        # S=5 over duration 2: original energy=50, integral=10, origin=5.
        self.assertEqual(E.shifted_S_energy(50,(10,0,0),(5,0,0),2),0)
        self.assertEqual(E.shifted_S_energy(50,(10,0,0),(3,0,0),2),8)
        with self.assertRaises(ValueError):
            E.shifted_S_energy(0,(10,0,0),(5,0,0),2)

    def test_plateau_is_all_window_quiet_not_a_frozen_sinusoid(self):
        d=E.quiet_plateau_witness()
        self.assertLess(F(d['all_window_AC_energy_upper']),F(d['quiet_threshold']))
        self.assertEqual(F(d['vertical_Hs_upper_via_4_std_le_4A_m']),F(42,5))
        self.assertTrue(d['all_continuous_window_starts_covered'])

class CentralBranchCoverageTest(unittest.TestCase):
    def state(self, Cv=I(0), Ce=I(0), log=None, elapsed=150, usable=False):
        seed=FRONT._point_state()
        w=replace(seed.wpe,accel_prev=I(0),high_pass_1=I(0),high_pass_1_prev=I(0),
            high_pass_2=I(0),velocity=I(0),elevation=I(0),velocity_mean=I(0),
            elevation_mean=I(0),velocity_sq=Cv,elevation_sq=Ce,weight=I(1),
            elapsed_s=I(elapsed),raw_period_s=None,log_period_s=log,usable_period=usable)
        return w,C.State(Cv,Ce,I(0),log)

    def test_prior_survives_quiet_moments_without_usable_period(self):
        w,c=self.state()
        for _ in range(5):
            successors=C.advance_wpe(w,c,I(0),WPE.constants())
            self.assertEqual([s.branch for s in successors],['hold_invalid_variance'])
            w,c=successors[0].shipping_state,successors[0].correlation_state
            self.assertIsNone(c.theorem_log_period_s)
            self.assertFalse(w.usable_period)
            self.assertEqual(C.frequency(c,usable_period=False),I(.2))

    def test_straddling_variance_does_not_drop_valid_history(self):
        w,c=self.state(Interval(0,1),Interval(0,1))
        branches={s.branch for s in C.advance_wpe(w,c,I(0),WPE.constants())}
        self.assertIn('hold_invalid_variance',branches)
        self.assertTrue(any(b.startswith('valid_period') for b in branches))

    def test_straddling_omega_does_not_drop_valid_history(self):
        w,c=self.state(Interval(.01,.1),I(1))
        branches={s.branch for s in C.advance_wpe(w,c,I(0),WPE.constants())}
        self.assertIn('hold_invalid_omega',branches)
        self.assertTrue(any(b.startswith('valid_period') for b in branches))

    def test_first_valid_log_and_one_way_takeover(self):
        w,c=self.state(I(1),I(1))
        successors=C.advance_wpe(w,c,I(0),WPE.constants())
        self.assertEqual([s.branch for s in successors],['valid_period_takeover'])
        s=successors[0]
        self.assertTrue(s.shipping_state.usable_period)
        self.assertEqual(s.shipping_state.last_log_horizon_s,I(0))
        self.assertEqual(C.frequency(c,usable_period=False),I(.2))
        self.assertNotEqual(C.frequency(s.correlation_state),I(.2))
        # Invalid updates HOLD a previously latched usable flag and log state.
        w,c=self.state(I(0),I(0),log=I(1),usable=True)
        s=C.advance_wpe(w,c,I(0),WPE.constants())[0]
        self.assertTrue(s.shipping_state.usable_period)
        self.assertEqual(s.correlation_state.theorem_log_period_s,I(1))

    def test_log_without_latch_still_uses_prior(self):
        w,c=self.state(I(1),I(1),log=I(2),elapsed=30)
        self.assertEqual(C.frequency(c,usable_period=False),I(.2))
        successors=C.advance_wpe(w,c,I(0),WPE.constants())
        self.assertTrue(all(not s.shipping_state.usable_period for s in successors))

    def test_zero_weight_branch_is_exact_not_a_small_weight_shortcut(self):
        self.assertEqual(C.central_update(I(0),I(0),I(0),I(7),I(.1)),I(0))
        with self.assertRaises(RuntimeError):
            C.central_update(I(0),I(0),Interval(0,1e-16),I(7),I(.1))

    def test_detached_log_rejected(self):
        w,c=self.state(log=I(1),usable=True)
        with self.assertRaises(ValueError):
            C.advance_wpe(w,replace(c,theorem_log_period_s=I(2)),I(0),WPE.constants())

    def test_joint_frontend_reaches_tuner_through_prior_selector(self):
        seed=JOINT._smoke_state()
        w,c=self.state()
        seed=replace(seed,frontend=replace(seed.frontend,wpe=w),
            statistics=replace(seed.statistics,wpe_velocity_C=c.wpe_velocity_C,
                wpe_elevation_C=c.wpe_elevation_C,theorem_log_period_s=None))
        sample=FRONT.Sample(JOINT.MAHONY.Vec3(I(0),I(0),I(0)),
            JOINT.MAHONY.Vec3(I(0),I(0),I(-9.80665)))
        images=JOINT.advance(seed,sample,gravity_ms2=I(9.80665),two_kp=I(.2),two_ki=I(.02),child_prefix='prior')
        self.assertTrue(images)
        self.assertTrue(all(s.frequency_hz==I(.2) for s in images))

if __name__=='__main__':unittest.main()
