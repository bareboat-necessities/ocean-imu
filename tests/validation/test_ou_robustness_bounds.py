import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools"))

import ou_robustness as robustness  # noqa: E402
import ou_validation as validation  # noqa: E402

# The floor is read from the mirror rather than written out, so lowering the
# implementation clamp moves these expectations with it instead of silently
# turning a floor test into a no-op.  test_record_conventions.py is what pins
# the mirror itself to SeaStateFusionFilter_OU_III.h.
R_S_FLOOR = robustness.R_S_BOUNDS_MS[0]


def _baseline():
    return validation.TuningPoint(tau_s=1.2, sigma_a_mps2=0.8, RS_ms=0.65543)


class RobustnessBoundTests(unittest.TestCase):
    def test_sigma_coupled_half_scale_stays_above_the_r_s_floor(self):
        # r_S scales linearly with sigma_aw on this axis: 0.65543 -> 0.327715,
        # which clears the deployed floor, so the request is realized exactly.
        # It clipped under the old 0.4 floor, which is what made the low-motion
        # seas insensitive to the tuner in the first place.
        point = robustness.scaled_tuning_point(_baseline(), "sigma_aw_rs", 0.5)
        self.assertGreater(point.RS_ms, R_S_FLOOR)
        self.assertAlmostEqual(point.RS_ms, 0.327715)
        self.assertAlmostEqual(point.sigma_a_mps2, 0.4)
        robustness.validate_tuning_point(point)

    def test_tau_coupled_half_scale_realizes_r_s_floor(self):
        # r_S scales as tau^3 on this axis: 0.65543 * 0.125 = 0.081929, which is
        # under the floor, so the clamp mirror must still engage.
        point = robustness.scaled_tuning_point(_baseline(), "tau_rs", 0.5)
        self.assertAlmostEqual(point.tau_s, 0.6)
        self.assertAlmostEqual(point.RS_ms, R_S_FLOOR)
        robustness.validate_tuning_point(point)

    def test_the_floor_is_the_only_thing_clipping_that_request(self):
        # Guards the pairing above: without the clamp the tau-coupled request
        # would land well below the floor, so the previous test really is
        # exercising the clamp and not an incidental equality.
        baseline = _baseline()
        unclipped = baseline.RS_ms * 0.5**3
        self.assertLess(unclipped, R_S_FLOOR)


if __name__ == "__main__":
    unittest.main()
