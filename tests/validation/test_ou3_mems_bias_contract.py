#!/usr/bin/env python3
import copy
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/stability"))

from ou3_interval import Interval, matrix_identity
import ou3_mems_bias_contract as BIAS
import ou3_p4_complete_brmm_differential_events as EVENTS


class MemsBiasContractTests(unittest.TestCase):
    def test_unqualified_device_and_unproved_separation_cannot_promote(self):
        contract = BIAS.build()
        self.assertEqual(BIAS.validate(contract), [])
        for section, field, forged in (
            ("BIAS0", "assembled_sensor_qualification_closed", True),
            ("BIAS0", "qualified_tau_interval_s", [5000.0, 5000.0]),
            ("BIAS0", "qualified_true_residual_root_bound_mps2", 0.35),
            ("BIAS1", "bias_error_is_free_GM_history", True),
            ("BIAS1", "Kalman_corrections_charged_to_exogenous_ISS", True),
            ("BIAS2", "uniform_separation_constant_lower", 0.001),
            ("stochastic_corollary", "PSD_used_to_prune_homogeneous_BRMM", True),
        ):
            with self.subTest(field=field):
                bad = copy.deepcopy(contract)
                bad[section][field] = forged
                self.assertTrue(BIAS.validate(bad))

    def test_root_history_encloses_exact_decay_at_every_prefix(self):
        root = tuple(Interval.point(x) for x in (0.12, -0.08, 0.03))
        tau = Interval.point(5000.0)
        self.assertEqual(BIAS.homogeneous_bias_at(root, tau, Interval.point(0.0)), root)
        for t in (0.005, 0.01, 1.0, 3.0):
            output = BIAS.homogeneous_bias_at(root, tau, Interval.point(t))
            for initial, cell in zip(root, output):
                exact = initial.lo * math.exp(-t / 5000.0)
                self.assertLessEqual(cell.lo, exact)
                self.assertGreaterEqual(cell.hi, exact)

    def test_invalid_time_coordinate_cannot_generate_a_bias_history(self):
        root = (Interval.point(0.1),) * 3
        for tau, elapsed in ((0.0, 1.0), (-1.0, 1.0), (5000.0, -0.1)):
            with self.assertRaises(ValueError):
                BIAS.homogeneous_bias_at(root, Interval.point(tau), Interval.point(elapsed))

    def test_band_fraction_has_two_sided_normalization_and_no_hard_cap(self):
        # Integral from -1/tau to +1/tau contains exactly half the variance.
        self.assertAlmostEqual(BIAS.ou_band_fraction(5.0, 0.0, 0.2), 0.5)
        self.assertEqual(BIAS.ou_band_fraction(5.0, 0.2, 0.2), 0.0)
        self.assertLess(BIAS.ou_band_fraction(5000.0, 0.2, 2.0), 0.001)
        with self.assertRaises(ValueError):
            BIAS.ou_band_fraction(5000.0, 2.0, 0.2)

    def test_noiseless_measurement_changes_bias_error_without_changing_true_bias(self):
        state = [Interval.point(0.0) for _ in range(21)]
        state[18] = Interval.point(0.1)
        truth = [Interval.point(0.0) for _ in range(3)]
        event = EVENTS.source_joseph_event(
            "A", state, matrix_identity(21), matrix_identity(3), "accelerometer",
            f_hat=[Interval.point(0.0), Interval.point(0.0), Interval.point(9.80665)],
            R_hat=matrix_identity(3), bias_true=truth, bias_projection_limit=0.4,
        )
        self.assertEqual(event["bias_projection_branch"], "inactive")
        self.assertLess(event["state_out"][18].hi, state[18].lo)
        # This is a measurement event: there is no elapsed-time GM decay.
        # Substituting a free GM trajectory for corrected e_b would miss it.


if __name__ == "__main__":
    unittest.main()
