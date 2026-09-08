"""Regression checks for accidental surface-data reuse during RAO migration."""
import importlib.util
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


if __name__ == '__main__':
    unittest.main()
