"""Physical-regime contracts; quiet packet evidence never promotes a theorem."""
from fractions import Fraction as F
import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.stability.ou3_theorem.regimes import (
    QuietEvidenceLimits, QuietEvidenceMonitor, bridge_prefixes, bridge_retained,
    certificate, moving_window_requirement, squared_bridge,
    stationary_gyro_average_radius,
)


class RegimeTests(unittest.TestCase):
    def window(self, start, span=0.1, end=100):
        return moving_window_requirement(episode_start=10, episode_end=end,
                                         window_start=start, window_s=20,
                                         theta_e=0.1, span=span)

    def test_departure_has_complete_window_grace(self):
        # [0,20] includes only the first ten seconds of motion; no pathological
        # requirement that its short moving suffix already have full span.
        self.assertFalse(self.window(0)["excitation_required"])
        self.assertFalse(self.window(9.99999)["excitation_required"])
        self.assertTrue(self.window(10)["excitation_required"])
        self.assertTrue(self.window(10)["complete_excited_window"])
        self.assertFalse(self.window(10)["regime_certificate"])

    def test_return_to_rest_and_short_episode(self):
        self.assertTrue(self.window(80)["excitation_required"])
        self.assertFalse(self.window(80.001)["excitation_required"])
        self.assertFalse(self.window(10, end=25)["excitation_required"])

    def test_sustained_constant_attitude_is_still_rejected(self):
        for start in (10, 30, 4000):
            row = self.window(start, span=0, end=math.inf)
            self.assertTrue(row["excitation_required"])
            self.assertFalse(row["complete_excited_window"])

    def test_span_uses_whole_window_not_endpoint_difference(self):
        self.assertTrue(self.window(10, span=.1)["complete_excited_window"])
        self.assertFalse(self.window(10, span=math.nextafter(.1, 0))["complete_excited_window"])
        for value in (None, math.nan, math.inf, -1, 4):
            self.assertFalse(self.window(10, span=value)["complete_excited_window"])

    def test_invalid_parameters_fail_even_on_crossing_window(self):
        for value in (0, -1, math.nan, math.inf, -math.inf):
            for key in ("theta_e", "window_s"):
                args = dict(episode_start=10, episode_end=100, window_start=0,
                            window_s=20, theta_e=.1)
                args[key] = value
                with self.assertRaises(ValueError):
                    moving_window_requirement(**args)

    def test_bridge_checks_every_prefix_and_keeps_supplies(self):
        ops = [(2, 1), (F(1, 4), 2), (3, 4)]
        rows = bridge_prefixes(2, ops)
        self.assertEqual([x["radius"] for x in rows], [5, F(13, 4), F(55, 4)])
        self.assertTrue(bridge_retained(2, ops, [5, 4, 14]))
        self.assertFalse(bridge_retained(2, ops, [4, 4, 14]))
        with self.assertRaises(ValueError):
            bridge_retained(2, [], [])
        gain, supply = rows[-1]["gain"], rows[-1]["supply"]
        sq = squared_bridge(gain, supply, F(1, 3))
        self.assertGreaterEqual(sq["storage_gain"] * 4 + sq["storage_supply"], F(55, 4)**2)

    def test_finite_bridges_do_not_prove_repeated_retention(self):
        rows = bridge_prefixes(1, [(3, 0), (F(1, 2), 0)] * 20)
        self.assertEqual(rows[-1]["radius"], F(3, 2)**20)
        self.assertGreater(rows[-1]["radius"], 3000)

    def test_gyro_average_charges_false_entry_and_supplied_fast_action(self):
        radius = stationary_gyro_average_radius(F(1, 50), F(1, 100000),
                                               [1, 0], [F(1, 2)] * 2, slow_amplitude=F(1,50))
        self.assertEqual(radius, F(4001, 200000))
        hidden = stationary_gyro_average_radius(F(1, 50), 0, [1, 0],
                                               [F(1, 2)] * 2, F(1, 100), slow_amplitude=F(1,50))
        self.assertEqual(hidden, F(3, 100))
        self.assertEqual(stationary_gyro_average_radius(F(1, 50), 0,
                         range(1000), [F(1, 1000)] * 1000, slow_amplitude=F(1,50)), F(1, 50))

    def test_exact_bias_service_and_scope_certificate(self):
        report = certificate()
        self.assertTrue(all(F(x) > 0 for x in report["witness_bias_margin"].values()))
        self.assertGreater(F(report["witness_actual_service_lower"]), 1)
        self.assertFalse(report["finite_exit_delay_from_existing_assumptions"])
        self.assertFalse(report["stationary_observable_quotient_boundedness"])
        self.assertFalse(report["theorem_closed"])


class QuietEvidenceTests(unittest.TestCase):
    def monitor(self):
        return QuietEvidenceMonitor(QuietEvidenceLimits(
            9.80665, .22516660498395405, .001, .3, .02, .00001, .02,
            dwell_s=2, max_gap_s=.006))

    def rest(self, monitor, start, count):
        result = None
        for k in range(count):
            # Amplitude-only necessary screen: passing it never certifies
            # the new fast temporal contract or physical STILL.
            result = monitor.step(start + k * .005, (.019, .001, 0),
                                  (.08, .03, -9.70665 + .05 * math.sin(k)))
        return result

    def test_single_packet_dwell_and_long_noisy_rest(self):
        m = self.monitor()
        self.assertEqual(self.rest(m, 0, 1), "TRANSITION")
        self.assertEqual(self.rest(m, .005, 399), "TRANSITION")
        self.assertEqual(self.rest(m, 2, 1), "STILL_COMPATIBLE")
        self.assertEqual(self.rest(m, 2.005, 20000), "STILL_COMPATIBLE")
        self.assertFalse(m.certified_still)
        self.assertFalse(m.certified_fast_history)

    def test_detectable_motion_exits_then_requires_fresh_dwell(self):
        m = self.monitor()
        self.rest(m, 0, 401)
        self.assertEqual(m.step(2.005, (.2, 0, 0), (0, 0, -9.80665)), "TRANSITION")
        self.assertEqual(self.rest(m, 2.010, 400), "TRANSITION")
        self.assertEqual(self.rest(m, 4.010, 2), "STILL_COMPATIBLE")

    def test_accel_variation_rejects_quiet_looking_norm(self):
        m = self.monitor()
        self.rest(m, 0, 401)
        self.assertEqual(m.step(2.005, (0, 0, 0), (3, 0, -math.sqrt(9.80665**2-9))),
                         "TRANSITION")

    def test_invalid_packets_gaps_and_clock_reversal_clear_evidence(self):
        for event in ((2.005, (math.nan, 0, 0), (0, 0, -9.80665)),
                      (2.005, (0, 0, 0), (0, 0, math.inf)),
                      (2.1, (0, 0, 0), (0, 0, -9.80665)),
                      (1, (0, 0, 0), (0, 0, -9.80665))):
            m = self.monitor()
            self.rest(m, 0, 401)
            self.assertEqual(m.step(*event), "TRANSITION")
            self.assertFalse(m.compatible)

    def test_false_entry_hidden_rest_wave_rest_has_no_detectable_exit(self):
        m = self.monitor()
        alpha, nu, start = .001, .025, 3
        period = 2 * math.pi / nu
        moving_seen = False
        for k in range(56000):
            t = k * .005
            phase = nu * (t - start)
            moving = start < t < start + period
            phi = alpha * math.sin(phase)**3 if moving else 0
            rate = 3 * alpha * nu * math.sin(phase)**2 * math.cos(phase) if moving else 0
            ba = (0, 9.80665 * math.sin(phi), 9.80665 * (math.cos(phi)-1))
            self.assertLessEqual(math.hypot(*ba), .22516660498395405)
            self.assertLessEqual(abs(rate), .02)
            # Exact identities, not rounded cancellation, supply common inputs.
            status = m.step(t, (0, 0, 0), (0, 0, -9.80665))
            if t >= 2:
                self.assertEqual(status, "STILL_COMPATIBLE")
            moving_seen |= abs(phi) > .0009
            self.assertFalse(m.certified_still)
        self.assertFalse(m.certified_fast_history)
        self.assertTrue(moving_seen)


if __name__ == "__main__":
    unittest.main()
