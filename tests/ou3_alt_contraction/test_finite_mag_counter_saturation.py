from pathlib import Path
from unittest.mock import patch
import os
import shutil
import subprocess
import tempfile
import unittest

from tools.stability.ou3_alt_contraction import finite_mag_counter_saturation as X

ROOT = Path(__file__).resolve().parents[2]


class Tests(unittest.TestCase):
    def test_tuner_counts_saturate_while_statistics_and_rejection_continue(self):
        from dataclasses import replace
        from tools.stability.ou3_alt_contraction import finite_mag_tuner_default as T
        c=T.Config(min_window=1000000,min_samples=X.SIGNED_MAX)
        s=T.State(accumulator=T.ACC.State(accepted_count=X.SIGNED_MAX-1),
                  rejected_count=X.SIGNED_MAX-1)
        for _ in range(3):
            out=T.step(s,c,q_tilt_bw=(1,0,0,0),mag_body=(30,0,0),dt=1,
                q_norm=T.SqrtWitness(1,1),mag_norm=T.SqrtWitness(900,30))
            s=out.state
            self.assertEqual(s.accumulator.accepted_count,X.SIGNED_MAX)
            out=T.step(s,c,q_tilt_bw=(1,0,0,0),mag_body=(0,0,0),dt=1,
                q_norm=T.SqrtWitness(1,1),mag_norm=T.SqrtWitness(0,0))
            s=out.state
            self.assertEqual(s.rejected_count,X.SIGNED_MAX)
        self.assertEqual(s.accumulator.weight_sum,3)
        self.assertEqual(s.accumulator.world_sum,(90,0,0))
        self.assertEqual(s.accumulator.accepted_window,3)
        for invalid in (-1,True,X.SIGNED_MAX+1):
            with self.assertRaises(ValueError): replace(s,rejected_count=invalid)
            with self.assertRaises(ValueError): replace(s.accumulator,accepted_count=invalid)

    def test_compressed_arbitrary_length_prefix_stays_in_int32(self):
        for c in (0, 249, X.SIGNED_MAX-1, X.SIGNED_MAX):
            for n in (0, 1, 250, 2**31, 10**100):
                result=X.after_attempts(c,n)
                self.assertTrue(c <= result <= X.SIGNED_MAX)
                self.assertEqual(X.after_attempts(result),X.after_attempts(c,n+1))
        self.assertEqual(X.after_attempts(X.SIGNED_MAX),X.SIGNED_MAX)

    def test_saturation_preserves_every_configurable_threshold_case(self):
        for c in (0, 249, 250, X.SIGNED_MAX-1, X.SIGNED_MAX):
            for n in (0, 1, 250, X.SIGNED_MAX, 10**100):
                for u in (0, 1, 250, 500, X.SIGNED_MAX):
                    self.assertTrue(X.threshold_equivalence(c,n,u))
        # Saturating at today's unlock threshold would break later increases.
        self.assertTrue(X.after_attempts(250,250) >= 500)

    def test_source_and_counter_proof_do_not_promote_whole_word(self):
        r=X.build()
        self.assertTrue(r['counter_lifetime_closed'])
        self.assertTrue(r['no_positive_minimum_call_gap_required'])
        self.assertFalse(r['complete_deployment_arithmetic_closed'])
        self.assertFalse(r['storage_search_allowed'])
        with patch.object(X,'source_fingerprint',return_value='different'):
            with self.assertRaisesRegex(RuntimeError,'re-audit'): X.build()

    def test_invalid_machine_counts_or_thresholds_rejected(self):
        for invalid in (-1,True,X.SIGNED_MAX+1):
            with self.assertRaises(ValueError): X.after_attempts(invalid)
            with self.assertRaises(ValueError): X.threshold_equivalence(0,0,invalid)

    def test_source_audit_also_binds_reset_and_tuner_paths(self):
        with tempfile.TemporaryDirectory() as td:
            source=Path(td)/'wrapper.h'
            source.write_text(X.SOURCE.read_text().replace(
                'int  mag_updates_applied_ = 0;', 'int  mag_updates_applied_ = -1;'))
            with patch.object(X,'SOURCE',source):
                with self.assertRaisesRegex(RuntimeError,'reset/configuration'): X.audit_source()
            tuner=Path(td)/'tuner.h'
            tuner.write_text(X.TUNER_SOURCE.read_text().replace(
                'incrementCount_(rejected_count_);', '++rejected_count_;',1))
            with patch.object(X,'TUNER_SOURCE',tuner):
                with self.assertRaisesRegex(RuntimeError,'acquisition/refinement'): X.audit_source()


class NativeTests(unittest.TestCase):
    def test_outer_measurement_and_release_continue_under_overflow_sanitizer(self):
        cxx=shutil.which(os.environ.get('CXX','g++'))
        eigen=Path(os.environ.get('EIGEN_INCLUDE_DIR','/usr/include/eigen3'))
        if not cxx or not (eigen/'Eigen/Dense').is_file():
            if os.environ.get('OU3_ALT_REQUIRE_NATIVE')=='1':
                self.fail('required native saturation check needs g++ and Eigen')
            self.skipTest('native saturation check needs g++ and Eigen')
        with tempfile.TemporaryDirectory() as td:
            exe=Path(td)/'mag-counter'
            cmd=[cxx,'-std=c++20','-O1','-ffp-contract=off','-fno-fast-math',
                 '-DEIGEN_DONT_VECTORIZE','-fsanitize=signed-integer-overflow',
                 '-fno-sanitize-recover=signed-integer-overflow',f'-I{eigen}',
                 f'-I{ROOT / "src"}',str(ROOT/'tests/ou3_alt_contraction/mag_counter_boundary.cpp'),
                 '-o',str(exe)]
            subprocess.run(cmd,check=True,capture_output=True,text=True,timeout=180)
            result=subprocess.run([str(exe)],capture_output=True,text=True,timeout=45)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertIn('saturated_count=2147483647',result.stdout)
            self.assertIn('measurement_continues=1 release_continues=1',result.stdout)
            self.assertNotIn('runtime error',result.stderr)


if __name__=='__main__': unittest.main()
