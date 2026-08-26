import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_first_accel_rotation_gauge_v2 as G


class Ou3P5FirstAccelRotationGaugeV2Tests(unittest.TestCase):
    def test_historical_v2_fails_closed_on_cross_axis_terms(self):
        with self.assertRaises(RuntimeError) as ctx:
            G.build(
                source_pieces=2,
                yaw_axis_face_pieces=4,
                force_magnitude_pieces=4,
            )
        self.assertIn(
            "first-prefix linear covariance gained cross-axis terms",
            str(ctx.exception),
        )


if __name__ == "__main__":
    unittest.main()
