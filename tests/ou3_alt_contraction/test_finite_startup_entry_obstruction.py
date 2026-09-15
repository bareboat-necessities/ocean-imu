"""Physical startup counterfamily and native correspondence, never promotion."""
from fractions import Fraction as F
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
import unittest

from tools.stability.ou3_alt_contraction import finite_startup_entry_obstruction as X
from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_binary32_mahony as B
from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as CLOCK

ROOT = Path(__file__).resolve().parents[2]


class Tests(unittest.TestCase):
    def test_exact_physical_pole_is_not_a_zero_or_invalid_quaternion(self):
        self.assertEqual(sum(v*v for v in X.SOUTH), 1)
        packet = X.quiet_packet(X.quiet_reference())
        self.assertEqual(packet.raw_accel_body, (0, 0, -X.GRAVITY))
        self.assertEqual(packet.raw_gyro_body, (0, 0, 0))
        with self.assertRaisesRegex(ValueError, 'Cayley chart pole'):
            CORE.cayley(X.SOUTH)

    def test_same_physical_IMU_and_unbounded_exact_nearby_chart_family(self):
        # The algebraic identity is proved for all integers in the note; these
        # checks catch sign, factor-two and quaternion-convention regressions.
        for n in (1, 2, 10, 1000, 10**9):
            q = X.near_south(n)
            self.assertEqual(sum(v*v for v in q), 1)
            self.assertEqual(CORE.cayley(q), (0, 0, F(n)-F(1, n)))
            self.assertEqual(X.quiet_packet(X.quiet_reference(q)).raw_accel_body,
                             (0, 0, -X.GRAVITY))

    def test_obstruction_uses_the_same_zero_history_in_all_bias_families(self):
        for family in ('BIAS0', 'BIAS1', 'BIAS2'):
            ref = X.quiet_reference(bias_family=family)
            self.assertEqual(ref.beta, (0, 0, 0))
            self.assertEqual(ref.centered_S, (0, 0, 0))
            X.quiet_packet(ref)
        report = X.build()
        self.assertEqual(X.validate(report), [])
        self.assertFalse(report['shipping_filter_instability_claimed'])
        self.assertFalse(report['universal_single_Cayley_fresh_entry_possible'])
        self.assertFalse(report['storage_search_allowed'])


class NativeTests(unittest.TestCase):
    def test_public_wrapper_reaches_ungauged_identity_at_exact_timeout_crossing(self):
        cxx = shutil.which(os.environ.get('CXX', 'g++'))
        eigen = Path(os.environ.get('EIGEN_INCLUDE_DIR', '/usr/include/eigen3'))
        if not cxx or not (eigen / 'Eigen/Dense').is_file():
            if os.environ.get('OU3_ALT_REQUIRE_NATIVE') == '1':
                self.fail('required native startup correspondence needs g++ and Eigen')
            self.skipTest('native startup correspondence needs g++ and Eigen')
        with tempfile.TemporaryDirectory() as td:
            exe = Path(td) / 'startup-entry'
            command = [cxx, '-std=c++20', '-O1', '-ffp-contract=off', '-fno-fast-math',
                       '-DEIGEN_DONT_VECTORIZE', f'-I{eigen}', f'-I{ROOT / "src"}',
                       str(ROOT / 'tests/ou3_alt_contraction/startup_entry_obstruction.cpp'),
                       str(ROOT / 'src/util/W3dSimCommon.cpp'), '-o', str(exe)]
            subprocess.run(command, check=True, capture_output=True, text=True, timeout=180)
            result = subprocess.run([str(exe)], check=True, capture_output=True,
                                    text=True, timeout=30)
        step, time, before, after, north, *q = map(int, result.stdout.split())
        self.assertEqual(step, CLOCK.STARTUP_TIMEOUT_STEPS)
        self.assertEqual(B.value(time), CLOCK.TIMEOUT_CROSSING)
        self.assertIn(before, (0, 1, 2))
        self.assertEqual(after, 3)
        self.assertEqual(north, 0)
        self.assertEqual(tuple(B.value(v) for v in q), X.IDENTITY)


if __name__ == '__main__':
    unittest.main()
