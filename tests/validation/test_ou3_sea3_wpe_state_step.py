import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from ou3_interval import Interval
import ou3_sea3_wpe_state_step as WPE
import ou3_validated_log as VLOG


class Sea3WpeStateStepTest(unittest.TestCase):
    def test_contract_is_source_neutral_and_fail_closed(self):
        d = WPE.build()
        self.assertEqual(d["qualification"], WPE.QUALIFICATION)
        self.assertTrue(d["shipping_source_parity_pass"])
        self.assertTrue(d["requires_same_SEA3_vertical_acceleration"])
        self.assertTrue(d["validity_boundaries_are_branched_not_selected"])
        self.assertTrue(d["validated_log_used"])
        self.assertTrue(d["mathematical_WPE_state_step_closed"])
        self.assertFalse(d["source_generator"])
        self.assertFalse(d["trajectory_replay_used"])
        self.assertFalse(d["target_binary32_libm_roundoff_closed"])
        self.assertFalse(d["private_Mahony_step_closed_here"])
        self.assertFalse(d["complete_SEA3_family_materialized_here"])
        self.assertFalse(d["P3_promoted"])
        self.assertEqual(WPE.validate(d), [])

    def _state(self):
        I = WPE.I
        return WPE.WPEState(
            accel_prev=I(0.0), high_pass_1=I(0.0), high_pass_1_prev=I(0.0),
            high_pass_2=I(0.0), velocity=I(0.1), elevation=I(0.1),
            velocity_mean=I(0.0), velocity_sq=I(0.4),
            elevation_mean=I(0.0), elevation_sq=I(0.1), weight=I(0.8),
            elapsed_s=I(60.0), raw_period_s=I(3.2),
            log_period_s=VLOG.log_interval(I(3.2)), usable_period=True,
            last_moment_horizon_s=I(20.0), last_log_horizon_s=I(0.16),
        )

    def test_normal_live_step_keeps_usable_latch(self):
        successors = WPE.advance(self._state(), a_vertical=WPE.I(0.2))
        self.assertGreaterEqual(len(successors), 1)
        self.assertTrue(all(s.usable_period for s in successors))
        self.assertTrue(all(s.weight.lo > 1e-3 for s in successors))

    def test_requires_already_usable_live_state(self):
        s = self._state()
        bad = WPE.replace(s, usable_period=False)
        with self.assertRaises(ValueError):
            WPE.advance(bad, a_vertical=WPE.I(0.2))

    def test_invalid_ratio_branch_holds_log_period(self):
        s = self._state()
        # Very broad moment boxes deliberately make invalid-ratio source
        # members possible. The executor must retain a hold branch rather than
        # assuming the period update succeeds.
        broad = WPE.replace(
            s,
            velocity_mean=Interval(-1.0, 1.0),
            velocity_sq=Interval(0.0, 2.0),
            elevation_mean=Interval(-1.0, 1.0),
            elevation_sq=Interval(0.0, 2.0),
        )
        successors = WPE.advance(broad, a_vertical=WPE.I(0.0))
        self.assertTrue(any(x.log_period_s == broad.log_period_s for x in successors))


if __name__ == "__main__":
    unittest.main()
