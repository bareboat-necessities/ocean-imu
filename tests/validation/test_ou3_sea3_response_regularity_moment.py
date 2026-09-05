import math
import sys
from pathlib import Path
import unittest

TOOLS = Path(__file__).resolve().parents[2] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_sea3_response_regularity_moment as REG


class Sea3ResponseRegularityMomentTest(unittest.TestCase):
    def test_analytical_bound_is_finite_and_nonpromoting(self):
        d = REG.build()
        self.assertEqual(REG.validate(d), [])
        self.assertTrue(d["SEA3_parameter_domain_compact"])
        self.assertTrue(d["same_continuum_RAO_family_consumed"])
        self.assertFalse(d["source_generator"])
        self.assertFalse(d["hard_finite_window_source_materialized"])
        self.assertFalse(d["may_substitute_for_hard_window_IQC"])
        self.assertFalse(d["P3_promoted"])
        for value in d["bounds"].values():
            self.assertTrue(math.isfinite(float(value)))
            self.assertGreater(float(value), 0.0)

    def test_partition_bound_uses_coupled_energy_not_three_independent_heights(self):
        d = REG.build()
        self.assertFalse(d["independent_H_T_rectangle_used"])
        self.assertFalse(d["independent_partition_H_maxima_used"])
        c = d["constants"]
        b = d["bounds"]
        expected_sum_h = math.sqrt(c["modes_max"]) * c["Hs_upper_m"]
        self.assertGreaterEqual(b["sum_H_upper_m"], expected_sum_h)
        self.assertLess(b["sum_H_upper_m"], math.nextafter(expected_sum_h, math.inf) * 1.0000001)

    def test_closed_form_jonswap_shape_bound_is_used(self):
        d = REG.build()
        expected = 7.0 * math.sqrt(5.0 * math.pi / 4.0)
        got = d["bounds"]["JONSWAP_dimensionless_I2_over_I0_upper"]
        self.assertGreaterEqual(got, expected)
        self.assertLess(got, expected * 1.0000001)


if __name__ == "__main__":
    unittest.main()
