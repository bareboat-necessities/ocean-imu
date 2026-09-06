#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))
    sys.path.insert(0, str(TOOLS / "stability"))

from ou3_interval import Interval
import ou3_p3_source_uniform_translation_covariance as U
import ou3_source_reachable_matrix_p3 as BASE


class SourceUniformTranslationCovarianceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = U.build()

    def test_source_uniform_contract_passes(self):
        self.assertEqual(U.validate(self.payload), [])
        self.assertTrue(self.payload["time_varying_source_parameters_covered_by_pointwise_extrema"])
        self.assertFalse(self.payload["P3_PROMOTED"])



if __name__ == "__main__":
    unittest.main()
