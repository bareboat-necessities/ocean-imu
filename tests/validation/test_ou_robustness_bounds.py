import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools"))

import ou_robustness as robustness  # noqa: E402
import ou_validation as validation  # noqa: E402


class RobustnessBoundTests(unittest.TestCase):
    def test_requested_half_scale_realizes_r_s_floor(self):
        baseline = validation.TuningPoint(
            tau_s=1.2,
            sigma_a_mps2=0.8,
            RS_ms=0.65543,
        )
        point = robustness.scaled_tuning_point(
            baseline, "sigma_aw_rs", 0.5
        )
        self.assertAlmostEqual(point.RS_ms, 0.4)
        self.assertAlmostEqual(point.sigma_a_mps2, 0.4)
        robustness.validate_tuning_point(point)

    def test_tau_coupled_half_scale_also_respects_floor(self):
        baseline = validation.TuningPoint(
            tau_s=1.2,
            sigma_a_mps2=0.8,
            RS_ms=0.65543,
        )
        point = robustness.scaled_tuning_point(baseline, "tau_rs", 0.5)
        self.assertAlmostEqual(point.tau_s, 0.6)
        self.assertAlmostEqual(point.RS_ms, 0.4)
        robustness.validate_tuning_point(point)


if __name__ == "__main__":
    unittest.main()
