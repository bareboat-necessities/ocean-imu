from fractions import Fraction as F
from pathlib import Path
from unittest.mock import patch
import os
import shutil
import subprocess
import tempfile
import unittest

from tools.stability.ou3_alt_contraction import finite_mag_counter_obstruction as X
from tools.stability.ou3_alt_contraction import finite_mag_call_schedule as CALLS

ROOT = Path(__file__).resolve().parents[2]


class Tests(unittest.TestCase):
    def test_burst_is_compressed_at_existing_fresh_source_endpoint(self):
        b = X.Burst()
        self.assertEqual(b.size, 2147483648)
        self.assertEqual(b.first_time, 0)
        self.assertEqual(b.time(1), b.time(b.size))
        self.assertLess(b.first_time, X.CLOCK.DT_REAL)
        self.assertEqual(b.projected_count_after_defined_prefix(b.size - 1), 2147483647)
        with self.assertRaisesRegex(ValueError, 'undefined increment'):
            b.projected_count_after_defined_prefix(b.size)

    def test_burst_and_unbounded_tail_obey_literal_schedule_checker(self):
        b = X.Burst()
        clock = CALLS.Clock(0)
        for ordinal in (1, 2, b.size, b.size + 1, b.size + 2):
            clock = CALLS.record_call(clock, b.schedule, time=b.time(ordinal))
        self.assertEqual(b.time(b.size + 10**20), b.first_time + 10**20 * F(1, 25))
        self.assertEqual(b.calls_through(b.first_time - F(1, 10**9)), 0)
        self.assertEqual(b.calls_through(b.first_time), b.size)
        self.assertEqual(b.calls_through(b.first_time + F(3, 25)), b.size + 3)

    def test_arbitrary_safe_entry_count_and_positive_horizon(self):
        for c in (0, 249, 1000000, 2147483646, 2147483647):
            for h in (F(3), F(1, 10**12)):
                b = X.Burst(c, h)
                self.assertTrue(0 <= b.first_time < h)
                self.assertEqual(b.projected_count_after_defined_prefix(0), c)
                self.assertEqual(b.projected_count_after_defined_prefix(b.size - 1), 2147483647)
        for invalid in (-1, True, 2147483648):
            with self.assertRaises(ValueError): X.Burst(invalid)

    def test_actual_caller_count_does_not_replace_generic_schedule(self):
        r = X.build()
        self.assertEqual(r['actual_sketch_conditional_30602_invocation_count_upper'], 30602)
        self.assertLess(X.caller_count_bound(30602), X.GATE.SIGNED_COUNTER_MAX)
        self.assertFalse(r['actual_caller_substituted_for_generic_schedule'])
        self.assertFalse(r['actual_sketch_35ms_gate_is_unconditional_min_gap'])
        self.assertFalse(r['actual_sketch_measured_dt_is_canonical_exact_5ms'])
        self.assertFalse(r['total_defined_execution_for_current_schedule_possible'])
        self.assertFalse(r['shipping_filter_instability_claimed'])
        self.assertFalse(r['storage_search_allowed'])

    def test_semantic_source_change_requires_reaudit(self):
        changed = dict(X.source_fingerprints(), inner_updateMag='different')
        with patch.object(X, 'source_fingerprints', return_value=changed):
            with self.assertRaisesRegex(RuntimeError, 're-audit'):
                X.build()


class NativeTests(unittest.TestCase):
    def test_unchanged_outer_call_reaches_last_defined_and_undefined_increment(self):
        cxx = shutil.which(os.environ.get('CXX', 'g++'))
        eigen = Path(os.environ.get('EIGEN_INCLUDE_DIR', '/usr/include/eigen3'))
        if not cxx or not (eigen / 'Eigen/Dense').is_file():
            if os.environ.get('OU3_ALT_REQUIRE_NATIVE') == '1':
                self.fail('required native counter correspondence needs g++ and Eigen')
            self.skipTest('native counter correspondence needs g++ and Eigen')
        with tempfile.TemporaryDirectory() as td:
            exe = Path(td) / 'mag-counter'
            cmd = [cxx, '-std=c++20', '-O1', '-ffp-contract=off', '-fno-fast-math',
                   '-DEIGEN_DONT_VECTORIZE', '-fsanitize=signed-integer-overflow',
                   '-fno-sanitize-recover=signed-integer-overflow', f'-I{eigen}',
                   f'-I{ROOT / "src"}',
                   str(ROOT / 'tests/ou3_alt_contraction/mag_counter_boundary.cpp'),
                   '-o', str(exe)]
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=180)
            safe = subprocess.run([str(exe)], capture_output=True, text=True, timeout=45)
            self.assertEqual(safe.returncode, 0, safe.stderr)
            self.assertIn('last_defined_count=2147483647', safe.stdout)
            bad = subprocess.run([str(exe), 'overflow'], capture_output=True, text=True, timeout=45)
            self.assertNotEqual(bad.returncode, 0)
            self.assertIn('SeaStateFusionFilter_OU_III.h:', bad.stderr)
            self.assertIn('signed integer overflow', bad.stderr)
            self.assertIn('2147483647 + 1', bad.stderr)


if __name__ == '__main__':
    unittest.main()
