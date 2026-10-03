"""Exact source proofs and fail-closed status of the moving compatibility audit."""
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.stability.ou3_theorem.moving_quiet_compatibility import (
    certificate, centered_slow_charge, sensor_polynomial_identity,
)
from tools.stability.ou3_theorem.moving_compatibility_diagnostic import (
    HEADER, PROBE, instrument,
)

class MovingQuietCompatibilityTests(unittest.TestCase):
    def test_exact_reproduction(self):
        path=ROOT/'reports/results/ou3_stability/moving-quiet-compatibility-certificate.json'
        self.assertEqual(json.loads(path.read_text()), certificate())

    def test_all_time_physical_membership_not_magnetic_promotion(self):
        c=certificate()
        self.assertTrue(c['marine_and_slow_fast_verified'])
        self.assertTrue(all(c['physical_membership_checks'].values()))
        self.assertFalse(c['actually_applied_magnetic_service_all_time_verified'])
        self.assertFalse(c['universal_MOVING_entry_refuted_on_fully_admitted_shipping_history'])
        self.assertFalse(c['shipping_instability_counterexample'])
        self.assertFalse(c['theorem_closed'])

    def test_formal_identity_not_phase_samples(self):
        result=sensor_polynomial_identity(F(196133,20000),F(400,40001),F(39999,40001))
        self.assertTrue(result['formal_polynomial_difference_zero'])
        self.assertEqual(result['two_epoch_difference_exact'],'0')

    def test_covariance_bound_uses_full_metric_marginal_comparison(self):
        c=certificate()['same_estimator_pair_metric']
        self.assertEqual(F(c['BA_marginal_ceiling']),F(1,1600))
        self.assertGreater(F(c['max_of_pair_V_lower']),15)
        self.assertFalse(c['estimator_symmetry_needed_for_pair_lower'])
        self.assertFalse(c['BA_physical_forcing_dropped'])

    def test_centered_history_charge_and_invalid_windows(self):
        self.assertEqual(centered_slow_charge(F('.05'),F('.001'),60,14,F('.006')),
                         F('.05')/14+F('.001')*14/4+F('.001')*F('.006'))
        with self.assertRaises(ValueError):
            centered_slow_charge(F('.05'),F('.001'),60,61,F('.006'))
        with self.assertRaises(ValueError):
            centered_slow_charge(F('.05'),F('.001'),60,0,F('.006'))

    def test_finite_probe_source_binding(self):
        c=json.loads((ROOT/'reports/results/ou3_stability/moving-compatibility-carried.json').read_text())
        self.assertEqual(c['probe_sha256'],hashlib.sha256(PROBE.read_bytes()).hexdigest())
        self.assertEqual(c['shipping_header_sha256'],hashlib.sha256(HEADER.read_bytes()).hexdigest())
        self.assertEqual(c['instrumented_header_sha256'],hashlib.sha256(instrument(HEADER.read_text()).encode()).hexdigest())
        self.assertFalse(c['all_time_magnetic_service_verified'])
        self.assertFalse(c['all_placed_service_windows_verified'])
        self.assertFalse(c['full_shipping_counterexample_admitted'])
        self.assertFalse(c['same_as_world_frame_fixture'])
        self.assertEqual(c['noise_profile'],{'sigma_a':.2,'sigma_g':.00135,'sigma_m':.8})
        self.assertGreater(c['native']['service_min'],1)
        self.assertGreater(c['native']['active_step'],0)
        self.assertGreater(c['native']['tail_V_lower_min'],15)

    def test_both_replays_export_real_tail_covariance_measurements(self):
        for seconds, name in ((240, 'moving-compatibility-carried.json'),
                              (1200, 'moving-compatibility-carried-1200s.json')):
            with self.subTest(duration=seconds):
                c = json.loads((ROOT/'reports/results/ou3_stability'/name).read_text())
                self.assertEqual(c['probe_sha256'], hashlib.sha256(PROBE.read_bytes()).hexdigest())
                self.assertEqual(c['shipping_header_sha256'], hashlib.sha256(HEADER.read_bytes()).hexdigest())
                self.assertEqual(c['instrumented_header_sha256'],
                                 hashlib.sha256(instrument(HEADER.read_text()).encode()).hexdigest())
                self.assertTrue(c['compiler'])
                n = c['native']
                self.assertEqual(n['duration'], seconds)
                self.assertEqual(n['tail_cov_samples'], (seconds - 200) * 200)
                for key in ('tail_pth_min', 'tail_pth_max', 'tail_pbg_norm_max',
                            'tail_pth_bg_norm_max', 'service_min'):
                    self.assertTrue(math.isfinite(n[key]))
                self.assertGreater(n['tail_pth_min'], 0)
                self.assertGreaterEqual(n['tail_pth_max'], n['tail_pth_min'])
                self.assertGreater(n['tail_pbg_norm_max'], 0)
                self.assertGreaterEqual(n['tail_pth_bg_norm_max'], 0)
                self.assertEqual(n['sliding_service_windows'], 4001)
                self.assertGreater(n['sliding_service_min'], 1)
                self.assertEqual(n['sliding_service_mags_min'], 25)
                self.assertEqual(n['parity_off_fro_max'], 0)
                self.assertLess(n['tail_fhat_max'], 10)
                self.assertEqual(n['cycle_tau_abs_diff'], 0)
                self.assertEqual(n['cycle_sigma_abs_diff'], 0)
                self.assertEqual(n['cycle_RS_abs_diff'], 0)
                self.assertEqual(n['cycle_period_abs_diff'], 0)
                self.assertGreater(n['cycle_scheduler_elapsed_abs_diff'], 0)
                self.assertFalse(c['all_time_magnetic_service_verified'])
                self.assertFalse(c['all_placed_service_windows_verified'])
                self.assertFalse(c['full_shipping_counterexample_admitted'])
                self.assertFalse(c['theorem_closed'])

if __name__=='__main__':
    unittest.main()
