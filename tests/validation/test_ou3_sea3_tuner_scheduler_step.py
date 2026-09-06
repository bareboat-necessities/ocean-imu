from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_sea3_tuner_scheduler_step as mod  # noqa: E402
from ou3_interval import Interval  # noqa: E402


class Sea3TunerSchedulerStepTest(unittest.TestCase):
    def test_contract_closes_only_tuner_scheduler_not_source(self):
        d = mod.build()
        self.assertEqual(mod.validate(d), [])
        self.assertTrue(d["tuner_scheduler_step_closed"])
        self.assertTrue(d["requires_same_SEA3_Mahony_and_WPE_state"])
        self.assertFalse(d["source_generator"])
        self.assertFalse(d["Mahony_step_closed_here"])
        self.assertFalse(d["WavePeriodEstimator_step_closed_here"])
        self.assertFalse(d["complete_SEA3_family_materialized_here"])
        self.assertFalse(d["P3_promoted"])

    def test_pending_commit_affects_only_next_sample(self):
        c = mod.constants()
        st = mod.TunerState(
            mod.BandState(mod.I(0), mod.I(0.1), mod.I(0.1), mod.I(0), mod.I(0.1), True),
            mod.MomentState(mod.I(0), mod.I(0.5), mod.I(0.1), mod.I(0.5), mod.I(0.2)),
            mod.CandidateState(mod.I(1.1), mod.I(0.5), mod.I(2.0)),
            mod.ActiveSchedule(mod.I(1.0), mod.I(0.4), mod.I(1.5), mod.pseudo_period(mod.I(1.0), c)),
            mod.SchedulerState(mod.I(0.099), False),
        )
        out = mod.advance_after_measurement(
            st, a_vertical=mod.I(0.2), f_wave_previous_wpe=mod.I(0.2), c=c
        )
        self.assertTrue(out)
        for s in out:
            self.assertEqual(s.active, st.active)
        pending = [s for s in out if s.scheduler.pending_commit]
        self.assertTrue(pending)
        committed = mod.commit_if_pending(pending[0], c)
        self.assertFalse(committed.scheduler.pending_commit)
        self.assertEqual(committed.active.tau, pending[0].candidate.tau)
        self.assertEqual(committed.active.rs_base, pending[0].candidate.rs)

    def test_actual_rs_uses_spectral_mse_base_without_extra_cadence_scale(self):
        c = mod.constants()
        active = mod.ActiveSchedule(mod.I(2.0), mod.I(0.8), mod.I(10.0), mod.pseudo_period(mod.I(2.0), c))
        xyz = mod.active_rs_std_xyz(active, c)
        self.assertEqual(xyz[0], mod.I(c.rs_x_factor) * mod.I(10.0))
        self.assertEqual(xyz[1], mod.I(c.rs_y_factor) * mod.I(10.0))
        self.assertEqual(xyz[2], mod.I(10.0))

    def test_rational_power_is_verified_not_trusted_libm(self):
        x = Interval.outward_bounds(0.01, 100.0)
        y = mod.rational_power_positive(x, 6, 7)
        self.assertGreater(y.lo, 0.0)
        x6 = mod.ipow(x, 6)
        self.assertLessEqual(mod.ipow(Interval.point(y.lo), 7).hi, x6.lo)
        self.assertGreaterEqual(mod.ipow(Interval.point(y.hi), 7).lo, x6.hi)


if __name__ == "__main__":
    unittest.main()
