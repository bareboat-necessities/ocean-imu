"""Failed native experiments must leave auditable receipts and a failing exit."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / 'tools/tangent_shared_modes_study.py'


class TangentStudyFailureTests(unittest.TestCase):
    def score_rows(self, work):
        with (work/'result/scores.csv').open() as stream:
            return list(csv.DictReader(stream))

    def run_study(self, work, binary, resume=False):
        command = [sys.executable, str(RUNNER), '--profile', 'development',
                   '--cases', 'rocking', '--seeds', '17', '--workers', '1',
                   '--binary', str(binary), '--raw', str(work/'raw'),
                   '--output', str(work/'result')]
        if resume:
            command.append('--resume')
        result = subprocess.run(command, capture_output=True, text=True, timeout=30)
        manifest = json.loads((work/'result/manifest.json').read_text())
        receipt = json.loads((work/'raw/rocking_17_native_receipt.json').read_text())
        return result, manifest, receipt

    def fake_binary(self, work, body, executable=True):
        binary = work/'fake-native'
        binary.write_text('#!/usr/bin/env python3\n' + body)
        binary.chmod(0o755 if executable else 0o644)
        return binary

    def test_early_exit_preserves_diagnostics_and_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            binary = self.fake_binary(work, "import sys\nprint('captured diagnostic', file=sys.stderr)\nsys.exit(7)\n")
            result, manifest, receipt = self.run_study(work, binary)
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertEqual(receipt['exit_code'], 7)
            self.assertIn('captured diagnostic', receipt['stderr'])
            self.assertTrue(receipt['failed'])
            self.assertFalse(receipt['outputs_complete'])
            for key in ('scores_sha256', 'trace_sha256', 'model_sha256'):
                self.assertIsNone(receipt[key])
            self.assertEqual(manifest['failed_experiments'], 1)
            self.assertEqual(manifest['incomplete_experiments'], 1)
            self.assertEqual(self.score_rows(work), [])
            resumed, again, retained = self.run_study(work, binary, resume=True)
            self.assertEqual(resumed.returncode, 2, resumed.stderr)
            self.assertEqual(again['commands'][0]['stderr'], receipt['stderr'])
            self.assertEqual(retained['exit_code'], 7)

    def test_partial_scores_survive_missing_trace_and_model(self):
        body = """import pathlib,sys
p=pathlib.Path(sys.argv[sys.argv.index('--output')+1])
pathlib.Path(str(p)+'_scores.csv').write_text(
    'scenario,seed,filter,start,end,sensor_hash,failed\\n'
    'rocking,17,QSM,150,180,paired,1\\n')
sys.exit(3)
"""
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            result, manifest, receipt = self.run_study(work, self.fake_binary(work, body))
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertEqual(receipt['exit_code'], 3)
            self.assertIsNotNone(receipt['scores_sha256'])
            self.assertIsNone(receipt['model_sha256'])
            self.assertEqual(manifest['failed_experiments'], 1)
            rows = self.score_rows(work)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]['failed'], '1')

    def test_launch_error_is_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            binary = self.fake_binary(work, 'pass\n', executable=False)
            result, manifest, receipt = self.run_study(work, binary)
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertIsNone(receipt['exit_code'])
            self.assertEqual(receipt['launch_error'], 'PermissionError')
            self.assertEqual(manifest['failed_experiments'], 1)

    def test_resume_rejects_changed_source_or_commit_without_relabelling(self):
        for field, stale in (('implementation_commit', '0'*40), ('source_fingerprint', '0'*64)):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as tmp:
                work = Path(tmp)
                binary = self.fake_binary(work, 'raise SystemExit(7)\n')
                _, manifest, receipt = self.run_study(work, binary)
                self.assertEqual(receipt['implementation_commit'], manifest['implementation_commit'])
                self.assertEqual(receipt['source_fingerprint'], manifest['source_fingerprint'])
                receipt[field] = stale
                (work/'raw/rocking_17_native_receipt.json').write_text(json.dumps(receipt))
                result, retained, _ = self.run_study(work, binary, resume=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('stale resume receipt', result.stderr)
                self.assertEqual(retained, manifest)

    def test_resume_rejects_legacy_unbound_receipt(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            binary = self.fake_binary(work, 'raise SystemExit(7)\n')
            _, manifest, receipt = self.run_study(work, binary)
            del receipt['source_fingerprint']
            (work/'raw/rocking_17_native_receipt.json').write_text(json.dumps(receipt))
            result, retained, _ = self.run_study(work, binary, resume=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('stale resume receipt', result.stderr)
            self.assertEqual(retained, manifest)

    def test_unpaired_outputs_cannot_report_success(self):
        body = """import pathlib,sys
p=pathlib.Path(sys.argv[sys.argv.index('--output')+1])
pathlib.Path(str(p)+'_scores.csv').write_text(
    'scenario,seed,filter,start,end,sensor_hash,failed\\n'
    'rocking,17,QSM,150,180,first,0\\n'
    'rocking,17,TFG,150,180,second,0\\n')
pathlib.Path(str(p)+'_trace.csv').write_text('trace')
pathlib.Path(str(p)+'_model.json').write_text('{}')
"""
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            result, manifest, receipt = self.run_study(work, self.fake_binary(work, body))
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertIn('unpaired sensors', receipt['errors'])
            self.assertEqual(manifest['failed_experiments'], 1)

    def test_plain_module_make_still_builds_and_tests(self):
        result = subprocess.run(['make', '-n', '-C', str(ROOT/'tests/kalman_carma_shared_modes')],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('bash ./run_tests.sh', result.stdout)
        explicit = subprocess.run(['make', '-n', '-C', str(ROOT/'tests/kalman_carma_shared_modes'), 'all'],
                                  capture_output=True, text=True, timeout=30)
        self.assertEqual(explicit.returncode, 0, explicit.stderr)
        self.assertEqual(result.stdout, explicit.stdout)


if __name__ == '__main__':
    unittest.main()
