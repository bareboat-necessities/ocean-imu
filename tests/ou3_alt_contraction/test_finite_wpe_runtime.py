"""Finite WPE recurrence regressions; no source/stability promotion."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W


def cfg(log_periods=F(1,2)):
    return W.WPEConfig(1,4,log_periods,1,180)


class Tests(unittest.TestCase):
    def test_two_hp_stages_and_leaky_integrators_share_same_decay(self):
        out=W.update(W.WPEState(last_moment_horizon=1),cfg(),dt=F(1,10),vertical_accel=2,
                     decay=W.ExpWitness(F(1,2)))
        s=out.state
        self.assertEqual(s.hp1,1)
        self.assertEqual(s.hp2,F(1,2))
        self.assertEqual(s.velocity,F(1,4))
        self.assertEqual(s.elevation,F(1,8))
        self.assertFalse(out.produced_period)

    def test_pre_moment_start_consumes_no_moment_witnesses(self):
        with self.assertRaises(ValueError):
            W.update(W.WPEState(last_moment_horizon=1),cfg(),dt=F(1,10),vertical_accel=0,
                     decay=W.ExpWitness(1),moment_decay=W.ExpWitness(1))

    def moment_state(self, log_period=None):
        return W.WPEState(weight=1,elapsed=4,velocity_mean=0,velocity_sq=5,
                          elevation_mean=0,elevation_sq=1,log_period=log_period,
                          last_moment_horizon=4)

    def test_valid_moment_ratio_builds_first_canonical_log_period(self):
        s=self.moment_state()
        out=W.update(s,cfg(),dt=F(1,10),vertical_accel=0,decay=W.ExpWitness(1),
                     moment_decay=W.ExpWitness(1),
                     period_witness=W.PeriodWitness(2,3,7))
        self.assertTrue(out.produced_period)
        self.assertEqual(out.state.raw_period,3)
        self.assertEqual(out.state.log_period,7)
        self.assertEqual(out.state.last_log_horizon,0)
        self.assertFalse(out.state.usable_period)

    def test_sqrt_period_witness_must_belong_to_same_moments(self):
        with self.assertRaises(ValueError):
            W.update(self.moment_state(),cfg(),dt=F(1,10),vertical_accel=0,
                     decay=W.ExpWitness(1),moment_decay=W.ExpWitness(1),
                     period_witness=W.PeriodWitness(3,3,7))

    def test_existing_canonical_state_drives_same_log_horizon_and_ema(self):
        s=self.moment_state(log_period=2)
        out=W.update(s,cfg(),dt=F(1,10),vertical_accel=0,decay=W.ExpWitness(1),
                     moment_decay=W.ExpWitness(1),period_witness=W.PeriodWitness(2,3,4),
                     log_witness=W.LogUpdateWitness(5,F(1,2),3),
                     current_period=5,current_frequency=F(1,5))
        self.assertEqual(out.state.log_period,3)
        self.assertEqual(out.state.last_log_horizon,F(5,2))
        self.assertTrue(out.state.usable_period)

    def test_period_and_frequency_are_one_canonical_reciprocal_pair(self):
        with self.assertRaises(ValueError):
            W.update(self.moment_state(log_period=2),cfg(),dt=F(1,10),vertical_accel=0,
                     decay=W.ExpWitness(1),moment_decay=W.ExpWitness(1),
                     period_witness=W.PeriodWitness(2,3,4),
                     log_witness=W.LogUpdateWitness(5,F(1,2),3),
                     current_period=5,current_frequency=F(1,4))

    def test_degenerate_moments_do_not_invent_period(self):
        s=W.WPEState(weight=1,elapsed=4,velocity_sq=0,elevation_sq=1,last_moment_horizon=4)
        out=W.update(s,cfg(),dt=F(1,10),vertical_accel=0,decay=W.ExpWitness(1),
                     moment_decay=W.ExpWitness(1))
        self.assertFalse(out.produced_period)
        self.assertIsNone(out.state.log_period)

    def test_readiness_stays_fail_closed_at_transcendentals_and_frontend(self):
        r=W.readiness()
        self.assertTrue(r['two_high_pass_stages_materialized'])
        self.assertTrue(r['moment_ratio_period_branch_materialized'])
        self.assertFalse(r['exp_log_sqrt_pi_binary32_ancestry_attached'])
        self.assertFalse(r['vertical_accel_frontend_same_history_attached'])
        self.assertFalse(r['complete_word_finite_identity'])
        self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
