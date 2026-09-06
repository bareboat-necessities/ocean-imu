import math
import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))
    sys.path.insert(0, str(TOOLS / "stability"))

from ou3_interval import Interval
import ou3_validated_log as VLOG


class ValidatedLogTest(unittest.TestCase):
    def test_exact_one_contains_zero(self):
        x = VLOG.log_point(1.0)
        self.assertLessEqual(x.lo, 0.0)
        self.assertGreaterEqual(x.hi, 0.0)

    def test_reference_values_are_enclosed(self):
        # libm is used only as an independent test oracle here, never by the
        # proof implementation itself.
        for value in (0.02, 0.5, 0.9, 1.1, 2.0, 3.2, 12.0, 180.0):
            with self.subTest(value=value):
                x = VLOG.log_point(value)
                self.assertTrue(x.contains(math.log(value)))
                self.assertGreaterEqual(x.hi, x.lo)

    def test_interval_monotonicity(self):
        src = Interval.outward_bounds(2.0, 8.0)
        out = VLOG.log_interval(src)
        self.assertTrue(out.contains(math.log(2.0)))
        self.assertTrue(out.contains(math.log(8.0)))

    def test_reject_nonpositive(self):
        with self.assertRaises(ValueError):
            VLOG.log_point(0.0)
        with self.assertRaises(ValueError):
            VLOG.log_interval(Interval(-1.0, 1.0))


if __name__ == "__main__":
    unittest.main()
