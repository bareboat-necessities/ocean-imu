from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

import ou3_p3_pseudo_scheduler_progress_certificate as P


class PseudoSchedulerProgressCertificateTests(unittest.TestCase):
    def test_unexpired_retarget_preserves_elapsed_exactly(self):
        elapsed = P._f32(0.12)
        period = P._f32(0.14)
        self.assertEqual(P._retarget(elapsed, period), elapsed)

    def test_overdue_retarget_arms_next_sample(self):
        dt = P._f32(0.005)
        period = P._f32(0.13)
        armed = P._retarget(P._f32(0.14), period)
        self.assertEqual(armed, P._nextafterf_down_positive(period))
        fire, _ = P._due(dt, period, armed)
        self.assertTrue(fire)

    def test_fixed_max_period_fires_on_deployed_bound(self):
        sched = P.BASE.source_schedule()
        tau_hi = float(sched["tau_applied_invariant_s"][1])
        period = P._period_from_tau(tau_hi, sched)
        sample, _ = P._first_fire_from_zero(P._f32(sched["dt_s"]), period)
        self.assertEqual(sample, P.DEPLOYED_MAX_GAP_SAMPLES)
        self.assertEqual(sample, 30)

    def test_former_635_sample_starvation_cycle_is_broken(self):
        sched = P.BASE.source_schedule()
        fires, worst, _ = P._replay_former_starvation(P._f32(sched["dt_s"]))
        self.assertGreater(len(fires), 0)
        self.assertLessEqual(worst, P.DEPLOYED_MAX_GAP_SAMPLES)

    def test_full_certificate_validates(self):
        d = P.build()
        self.assertEqual(P.validate(d), [])
        self.assertEqual(d["certified_uniform_max_gap_samples"], P.DEPLOYED_MAX_GAP_SAMPLES)
        self.assertTrue(d["scheduler_recurrence_certificate"])
        self.assertFalse(d["P3_PROMOTED"])


if __name__ == "__main__":
    unittest.main()
