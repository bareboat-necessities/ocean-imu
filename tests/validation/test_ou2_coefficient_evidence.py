"""A regenerated coefficient report must replace its dataset and byte receipts."""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools'))
import ou2_dual_mse_coefficients as study  # noqa: E402


class CoefficientEvidenceTests(unittest.TestCase):
    def test_report_and_manifest_are_regenerated_together(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            points = out / 'startup.csv'
            points.write_text('current startup rows\n')
            (out / 'manifest.json').write_text('{"simulation_provenance":{"release":"obsolete"}}')
            provenance = {'release': 'current', 'input_sha256': {'wave.csv': 'verified'}}

            def calculate(args):
                print('fresh coefficient diagnostic')
                return 0

            with patch.object(study, 'calculate', side_effect=calculate), \
                 patch.object(study, 'input_provenance', return_value=provenance), \
                 contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(study.main(['--exact', '--output-dir', str(out),
                                             '--operating-points', str(points)]), 0)
            manifest = json.loads((out / 'manifest.json').read_text())
            self.assertEqual(manifest['simulation_provenance'], provenance)
            self.assertEqual(manifest['operating_points_sha256'], hashlib.sha256(points.read_bytes()).hexdigest())
            self.assertEqual(manifest['files']['coefficient_report.txt'],
                             hashlib.sha256((out / 'coefficient_report.txt').read_bytes()).hexdigest())
            self.assertTrue(manifest['exact_objective'])


if __name__ == '__main__':
    unittest.main()
