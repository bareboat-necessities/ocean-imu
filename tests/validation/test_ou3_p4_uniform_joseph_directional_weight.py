#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p4_uniform_joseph_directional_weight as U


class UniformJosephDirectionalWeightTests(unittest.TestCase):
    def test_P3_metric_yields_strict_uniform_directional_weights(self):
        d = U.build()
        self.assertEqual(U.validate(d), [])
        self.assertTrue(d["P4_UNIFORM_JOSEPH_DIRECTIONAL_WEIGHTS_ESTABLISHED"])
        self.assertTrue(d["P3_covariance_Loewner_upper_consumed"])
        self.assertTrue(d["P3_prefix_metric_contract_consumed"])
        self.assertFalse(d["interval_K_materialized"])
        self.assertFalse(d["condition_number_conversion_used"])
        self.assertFalse(d["per_packet_full_state_scalarization_used"])
        self.assertFalse(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"])
        self.assertFalse(d["P5_FINITE_CAPTURE_ESTABLISHED"])

        for mode, n in (("H", 18), ("A", 21)):
            m = d["modes"][mode]
            self.assertEqual(m["dimension"], n)
            self.assertEqual(len(m["P3_Sigma_diagonal_Loewner_dominator"]), n)
            self.assertTrue(m["directional_rank_preserved"])
            self.assertFalse(m["complete_word_directional_sum_emitted_here"])
            self.assertFalse(m["P4_PROMOTED"])
            for key in (
                "uniform_vector_Joseph_attenuation_lower",
                "uniform_S_Joseph_attenuation_lower",
                "P3_R_inverse_vector_gyro_bias_alpha6_lower",
                "finite_angle_vector_information_retention_lower",
                "uniform_finite_angle_Joseph_vector_gyro_bias_alpha6_lower",
            ):
                x = m[key]
                self.assertTrue(math.isfinite(x), (mode, key, x))
                self.assertGreater(x, 0.0, (mode, key, x))
            self.assertLessEqual(m["uniform_vector_Joseph_attenuation_lower"], 1.00000000000001)
            self.assertLessEqual(m["uniform_S_Joseph_attenuation_lower"], 1.00000000000001)
            self.assertLessEqual(
                m["uniform_finite_angle_Joseph_vector_gyro_bias_alpha6_lower"],
                m["P3_R_inverse_vector_gyro_bias_alpha6_lower"],
            )

    def test_validation_rejects_global_scalarization(self):
        d = U.build()
        d["per_packet_full_state_scalarization_used"] = True
        f = U.validate(d)
        self.assertTrue(any("per_packet_full_state_scalarization_used" in x for x in f))


if __name__ == "__main__":
    unittest.main()
