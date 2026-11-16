"""Direction identity, trajectory consistency, and matrix coverage contracts."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
import sim_dataset as data
import simulation_studies as studies
from underway_reference import track


class BidirectionalSimTests(unittest.TestCase):
    def test_cases_share_names_but_never_bytes_or_receipts(self):
        self.assertEqual(set(data.CASES), {'following', 'head'})
        self.assertEqual(data.RELEASE, 'v1.4.0')
        north, south = (data.REFERENCE_CSV_SHA256_BY_CASE[c] for c in ('following', 'head'))
        self.assertEqual(set(north), set(south))
        self.assertEqual(len(north), 8)
        self.assertTrue(all(north[n] != south[n] for n in north))
        name = next(iter(north))
        with patch.object(data, 'digest', return_value=north[name]):
            receipt = data.input_provenance([Path(name)], 'following')
            self.assertEqual(receipt['case'], 'following')
            self.assertEqual(receipt['heading_period_s'], 20)
            with self.assertRaisesRegex(ValueError, 'v1.4.0/head'):
                data.input_provenance([Path(name)], 'head')
        with self.assertRaises(ValueError):
            data.selected_case('typo')

    def test_python_and_cpp_tracks_agree_in_both_directions(self):
        source = r'''
#include "util/UnderwayReference.h"
#include <iostream>
#include <iomanip>
int main() {
  std::cout << std::setprecision(17);
  for (const char* name : {"following", "head"}) {
    setenv("OCEAN_IMU_SIM_CASE", name, 1);
    for (double t : {0., 5., 15., 20., 1000., 1199.995}) {
      const auto r = UnderwayReference::at(t);
      for (const auto* v : {&r.displacement, &r.velocity, &r.acceleration})
        for (double x : *v) std::cout << x << ' ';
      std::cout << '\n';
    }
  }
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            cpp, exe = Path(tmp)/'track.cpp', Path(tmp)/'track'
            cpp.write_text(source)
            subprocess.run(['g++', '-std=c++20', '-Wall', '-Wextra', '-Werror', '-I', str(ROOT/'src'),
                            str(cpp), '-o', str(exe)], check=True, capture_output=True)
            result = subprocess.check_output([str(exe)], text=True)
        actual = np.array([[float(v) for v in line.split()] for line in result.splitlines()])
        times = [0, 5, 15, 20, 1000, 1199.995]
        expected = np.vstack([np.column_stack(track(times, case)[:3]) for case in ('following', 'head')])
        np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=1e-12)
        np.testing.assert_allclose(actual[6:], -actual[:6], atol=1e-10)
        self.assertAlmostEqual(track([5], 'head')[3][0], 183)
        self.assertAlmostEqual(track([15], 'following')[3][0], -3)

    def test_changed_case_invalidates_same_name_resumption(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = Path(tmp)/'runs.csv'
            for case in ('following', 'head'):
                result.write_text('stale rows')
                with patch.object(data, 'input_provenance', return_value={'case': case}):
                    data.bind_dataset_results(result, Path(tmp))
                self.assertFalse(result.exists())
            marker = json.loads(result.with_suffix('.inputs.json').read_text())
            self.assertEqual(marker['dataset']['case'], 'head')

    def test_complete_study_rows_are_required(self):
        rows = [dict(input=name, disp_z_rms_m='1', roll_rms_deg='1', pitch_rms_deg='1', yaw_rms_deg='1')
                for name in data.REFERENCE_CSV_SHA256]
        data.require_complete_metrics(rows, 0)
        for bad, code in [(rows[:-1], 0), (rows + rows[:1], 0), (rows, 2)]:
            with self.assertRaises(ValueError):
                data.require_complete_metrics(bad, code)
        rows[0]['disp_z_rms_m'] = 'nan'
        with self.assertRaises(ValueError):
            data.require_complete_metrics(rows, 0)

    def test_all_studies_are_in_both_case_jobs(self):
        workflow = (ROOT/'.github/workflows/simulation-studies.yml').read_text()
        self.assertIn('case: [following, head]', workflow)
        for name in studies.STUDIES:
            self.assertIn(name, workflow)
        for name in ('build', 'simulation-matrix', 'carma-lever-arm-study', 'ou3-lever-arm-study',
                     'tangent-shared-modes', 'tfg-validation', 'heave-only'):
            self.assertIn('case: [following, head]', (ROOT/f'.github/workflows/{name}.yml').read_text())
        self.assertIn('sim-regression || status=$$?', (ROOT/'Makefile').read_text())
        self.assertIn("'tangent'", (ROOT/'tools/filter_sim_gates.py').read_text())


if __name__ == '__main__':
    unittest.main()
