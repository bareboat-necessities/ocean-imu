from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p3_clock_phase_ancestor_envelope as A


class P3ClockPhaseAncestorEnvelopeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = A.build()

    def test_validates_and_keeps_full_physical_partition(self):
        self.assertEqual(A.validate(self.d), [])
        self.assertEqual(self.d["physical_source_nodes"], 800)
        self.assertEqual(self.d["clock_phase_gap_alphabet_samples"], list(range(13, 27)))
        self.assertGreater(self.d["P3_word_horizon_upper_s"], 3.0)
        self.assertGreater(self.d["reverse_stage_transitions_admitted"], 0)

    def test_every_endpoint_has_a_nonempty_source_complete_ancestor_set(self):
        self.assertEqual(len(self.d["endpoints"]), 800)
        for row in self.d["endpoints"]:
            self.assertGreaterEqual(row["ancestor_count"], 1)
            self.assertLessEqual(row["ancestor_count"], 800)
            env = row["envelope"]
            for key in (
                "tau_s", "sigma_tuner_raw_mps2", "sigma_filter_committed_mps2",
                "R_S_filter_std", "pseudo_update_period_s",
            ):
                lo, hi = env[key]
                self.assertGreater(lo, 0.0)
                self.assertLessEqual(lo, hi)

    def test_hard_corner_and_node_zero_are_reported_without_promotion(self):
        for key in ("0", "729"):
            row = self.d["diagnostic_endpoints"][key]
            self.assertGreater(row["ancestor_count"], 0)
        self.assertFalse(self.d["P3_COVARIANCE_UPPER_ESTABLISHED_HERE"])
        self.assertFalse(self.d["P3_PROMOTED"])

    def test_validation_fails_closed_on_promotion_or_clock_drift(self):
        for key, value in (
            ("P3_PROMOTED", True),
            ("P3_COVARIANCE_UPPER_ESTABLISHED_HERE", True),
            ("clock_phase_gap_alphabet_samples", [21]),
        ):
            with self.subTest(key=key):
                d = deepcopy(self.d)
                d[key] = value
                self.assertNotEqual(A.validate(d), [])


if __name__ == "__main__":
    unittest.main()
