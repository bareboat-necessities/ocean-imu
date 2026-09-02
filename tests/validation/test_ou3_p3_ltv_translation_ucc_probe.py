#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p3_ltv_translation_ucc_probe as P


class LtvTranslationUccProbeTests(unittest.TestCase):
    def test_relative_floor_is_strict(self):
        d = P.ltv_relative_process_floor([1.0, 1.0, 1.0, 1.0], 0.5, 1.0/3.0, 12.0, 0.05)
        self.assertTrue(math.isfinite(d["relative_process_floor_lower"]))
        self.assertGreater(d["relative_process_floor_lower"], 0.0)

    def test_common_covariance_scaling_has_expected_inverse_effect(self):
        a = P.ltv_relative_process_floor([1.0, 2.0, 3.0, 4.0], 0.5, 1.0/3.0, 12.0, 0.05)
        b = P.ltv_relative_process_floor([10.0, 20.0, 30.0, 40.0], 0.5, 1.0/3.0, 12.0, 0.05)
        ratio = a["relative_process_floor_lower"] / b["relative_process_floor_lower"]
        self.assertGreater(ratio, 9.999999999)
        self.assertLess(ratio, 10.000000001)

    def test_probe_is_fail_closed_before_measurement_attenuation(self):
        d = P.build(source_node_indices=(729,))
        self.assertTrue(d["arbitrary_time_varying_tau_inside_window_covered"])
        self.assertTrue(d["arbitrary_time_varying_sigma_inside_window_covered_by_qc_min"])
        self.assertFalse(d["frozen_parameter_Q_Nh_identity_used"])
        self.assertFalse(d["interleaved_measurement_attenuation_enclosed_here"])
        self.assertFalse(d["P3_PROMOTED"])
        self.assertEqual(P.validate(d), [])


if __name__ == "__main__":
    unittest.main()
