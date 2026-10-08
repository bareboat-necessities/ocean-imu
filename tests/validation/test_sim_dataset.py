"""Regression checks for accidental surface-data reuse during RAO migration."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('sim_dataset', ROOT/'tools/sim_dataset.py')
D = importlib.util.module_from_spec(spec)
spec.loader.exec_module(D)


class SimDatasetTests(unittest.TestCase):
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
            self.assertNotIn('--pattern sim-data-files-vessel-rao-28ft.zip', source, path.name)
            if 'sim-data-files-vessel-rao-28ft.zip' in source and 'ou-validation' not in path.name:
                self.assertIn('tools/sim_dataset.py', source, path.name)

    def test_make_checks_all_inputs_instead_of_one_filename(self):
        source = (ROOT/'Makefile').read_text()
        self.assertIn('tools/sim_dataset.py', source)
        self.assertIn('--dest $(TEST_DIRS)', source)
        self.assertNotIn('SIM_DATA_CHECK_FILE', source)

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
                    if node.get('archive') == D.ARCHIVE and isinstance(inputs, dict):
                        for name, sha in inputs.items():
                            self.assertEqual(D.REFERENCE_CSV_SHA256.get(name), sha, f'{path}: {name}')
                        checked += 1
                    stack.extend(node.values())
        self.assertGreater(checked, 0)


if __name__ == '__main__':
    unittest.main()
