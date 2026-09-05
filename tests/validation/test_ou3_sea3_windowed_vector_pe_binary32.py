from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_windowed_vector_pe as mod  # noqa: E402


class Sea3WindowedPEBinary32Test(unittest.TestCase):
    def test_shipping_float_literals_are_quantized_before_variance(self):
        d = mod.build()
        self.assertEqual(mod.validate(d), [])
        self.assertTrue(d["source_literals_converted_to_binary32_before_variance"])
        self.assertTrue(d["variance_bounds_outward_after_binary32_conversion"])
        m = d["measurement_runtime"]
        a = mod.binary32(0.2)
        b = mod.binary32(0.3)
        self.assertEqual(m["accelerometer_std_mps2"], [a, a, a])
        self.assertEqual(m["magnetometer_std_uT"], [b, b, b])
        self.assertGreaterEqual(m["accelerometer_variance_upper"], a * a)
        self.assertGreaterEqual(m["magnetometer_variance_upper"], b * b)


if __name__ == "__main__":
    unittest.main()
