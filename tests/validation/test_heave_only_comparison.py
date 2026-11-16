"""Compact publication traces preserve samples and reject incomplete replays."""
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from heave_only_comparison import COLUMNS, retain_tail


class HeaveTraceTests(unittest.TestCase):
    def test_complete_trace_keeps_original_rate_and_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, target = Path(tmp)/'full.csv', Path(tmp)/'tail.csv'
            time = np.arange(240000)/200
            frame = pd.DataFrame({name: time if name == 'time' else np.sin(time)
                                  for name in COLUMNS})
            frame.to_csv(source, index=False)
            retain_tail(source, target)
            actual = pd.read_csv(target)
            self.assertGreaterEqual(len(actual), 12000)
            self.assertLessEqual(len(actual), 12002)
            np.testing.assert_allclose(actual.time.diff().dropna(), .005, atol=1e-12)
            np.testing.assert_allclose(actual.disp_est_z, np.sin(actual.time), atol=1e-12)
            self.assertEqual(actual.time.iloc[-1], time[-1])
            self.assertFalse(source.exists())

    def test_incomplete_and_nonfinite_traces_cannot_be_published(self):
        for bad in (False, True):
            with self.subTest(nonfinite=bad), tempfile.TemporaryDirectory() as tmp:
                source, target = Path(tmp)/'bad.csv', Path(tmp)/'tail.csv'
                frame = pd.DataFrame({name: np.arange(5, dtype=float) for name in COLUMNS})
                if bad:
                    frame.loc[2, 'disp_est_z'] = np.nan
                frame.to_csv(source, index=False)
                with self.assertRaises(ValueError):
                    retain_tail(source, target)
                self.assertTrue(source.exists())
                self.assertFalse(target.exists())


if __name__ == '__main__':
    unittest.main()
