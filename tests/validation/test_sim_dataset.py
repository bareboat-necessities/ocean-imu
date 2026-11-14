"""Regression checks for accidental surface-data reuse during RAO migration."""
import csv
import importlib.util
from unittest.mock import patch
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('sim_dataset', ROOT/'tools/sim_dataset.py')
D = importlib.util.module_from_spec(spec)
spec.loader.exec_module(D)


class SimDatasetTests(unittest.TestCase):
    def test_exact_clock_survives_the_1000_second_boundary(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'clock.csv'
            times = [i / 200 for i in range(199998, 200004)]
            path.write_text('time,value\n' + ''.join(f'{t:.17g},0\n' for t in times))
            D.verify_csv_clock(path, first_index=199998, samples=6)
            path.write_text('time,value\n' + ''.join(f'{t:.6g},0\n' for t in times))
            with self.assertRaisesRegex(ValueError, 'invalid timestamp'):
                D.verify_csv_clock(path, first_index=199998, samples=6)

    def test_fenton_skip_and_missing_sample_are_checked(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'clock.csv'
            path.write_text('time,value\n0.01,0\n0.014999999999999999,0\n0.02,0\n')
            D.verify_csv_clock(path, first_index=2, samples=3)
            with self.assertRaisesRegex(ValueError, 'expected 4 samples'):
                D.verify_csv_clock(path, first_index=2, samples=4)

    def test_wrong_existing_archive_is_rejected_without_reuse(self):
        with tempfile.TemporaryDirectory() as directory:
            archive = Path(directory)/D.ARCHIVE
            archive.write_bytes(b'old surface data with the new archive name')
            with self.assertRaisesRegex(ValueError, 'SHA-256 mismatch'):
                D.fetch(archive)

    def test_all_download_workflows_use_verified_vessel_archive(self):
        for path in (ROOT/'.github/workflows').glob('*.yml'):
            source = path.read_text()
            self.assertNotIn('sim-data-files.zip', source, path.name)
            self.assertNotIn('releases/download/v1.1.3', source, path.name)
            self.assertNotIn('--pattern sim-data-files-vessel-rao-28ft-5kn.zip', source, path.name)
            if 'sim-data-files-vessel-rao-28ft-5kn.zip' in source and 'ou-validation' not in path.name:
                self.assertIn('tools/sim_dataset.py', source, path.name)

    def test_make_checks_all_inputs_instead_of_one_filename(self):
        source = (ROOT/'Makefile').read_text()
        self.assertIn('tools/sim_dataset.py', source)
        self.assertIn('--dest $(TEST_DIRS)', source)
        self.assertNotIn('SIM_DATA_CHECK_FILE', source)

    def test_changed_dataset_invalidates_append_only_results(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'dataset_runs.csv'
            path.write_text('old stationary rows')
            with patch.object(D, 'input_provenance', return_value={'archive_sha256': D.SHA256}):
                D.bind_dataset_results(path, Path(directory))
                self.assertFalse(path.exists())
                path.write_text('new moving rows')
                D.bind_dataset_results(path, Path(directory))
                self.assertEqual(path.read_text(), 'new moving rows')

    def test_resume_replays_partial_or_duplicate_eight_record_cells(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'dataset_runs.csv'
            with path.open('w') as stream:
                writer = csv.DictWriter(stream, fieldnames=('family', 'seed', 'input'))
                writer.writeheader()
                for seed, names in [('1', list(D.REFERENCE_CSV_SHA256)),
                                    ('2', list(D.REFERENCE_CSV_SHA256)[:3]),
                                    ('3', [next(iter(D.REFERENCE_CSV_SHA256))] * 8)]:
                    for name in names:
                        writer.writerow({'family': 'A41', 'seed': seed, 'input': name})
            self.assertEqual(D.completed_dataset_runs(path), {('A41', '1')})
            with path.open() as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), 8)

    def test_diagnostics_reject_stale_startup_operating_points(self):
        import sys
        import hashlib
        sys.path.insert(0, str(ROOT / 'tools'))
        from ou2_dual_mse_coefficients import verify_operating_points
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'startup_runs.csv'
            path.write_text('measured operating points')
            marker = path.with_name('manifest.json')
            expected = {'archive_sha256': D.SHA256}
            data = {'simulation_provenance': {'archive_sha256': 'stationary'},
                    'files': {path.name: hashlib.sha256(path.read_bytes()).hexdigest()}}
            marker.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, 'rerun tools/startup_ablation.py'):
                verify_operating_points(path, expected)
            data['simulation_provenance'] = expected
            marker.write_text(json.dumps(data))
            verify_operating_points(path, expected)
            path.write_text('different rows with a stale receipt')
            with self.assertRaises(ValueError):
                verify_operating_points(path, expected)

    def test_committed_results_consumed_the_pinned_motion_bytes(self):
        # Older bundles name the release they were run from; their inputs must
        # still be byte-identical to the motion CSVs of the pinned release.
        checked = 0
        for path in (ROOT/'reports').rglob('*.json'):
            stack = [json.loads(path.read_text())]
            while stack:
                node = stack.pop()
                if isinstance(node, list):
                    stack.extend(node)
                elif isinstance(node, dict):
                    inputs = node.get('input_sha256')
                    if (node.get('archive') == D.ARCHIVE and node.get('release') == D.RELEASE
                            and isinstance(inputs, dict)):
                        for name, sha in inputs.items():
                            self.assertEqual(D.REFERENCE_CSV_SHA256.get(name), sha, f'{path}: {name}')
                        checked += 1
                    stack.extend(node.values())
        self.assertGreater(checked, 0)


if __name__ == '__main__':
    unittest.main()
