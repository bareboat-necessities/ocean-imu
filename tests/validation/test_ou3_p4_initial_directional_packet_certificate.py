#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p4_initial_directional_packet_certificate as P


class InitialDirectionalPacketCertificateTests(unittest.TestCase):
    def test_source_certificate_is_positive_but_not_scalarized(self):
        d = P.build()
        self.assertEqual(P.validate(d), [])
        self.assertTrue(d["P4_INITIAL_DIRECTIONAL_JOSEPH_PACKET_ESTABLISHED"])
        self.assertFalse(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"])
        self.assertFalse(d["P5_FINITE_CAPTURE_ESTABLISHED"])
        self.assertTrue(d["same_information_metric_route"])
        self.assertFalse(d["standalone_eta_norm_penalty_used"])
        self.assertFalse(d["condition_number_conversion_used"])
        self.assertFalse(d["per_packet_full_state_scalarization_used"])
        self.assertGreaterEqual(d["declared_sector_angle_rad"], 0.80)

        for mode, n in (("H", 18), ("A", 21)):
            m = d["modes"][mode]
            self.assertEqual(m["dimension"], n)
            self.assertEqual(m["vector_packet_rank_exact"], 5)
            self.assertFalse(m["instantaneous_full_state_scalar_margin_valid"])
            self.assertFalse(m["complete_word_prefix_covariances_propagated_here"])
            self.assertFalse(m["complete_word_directional_accumulation_closed_here"])
            self.assertFalse(m["P4_PROMOTED"])

            a0 = m["P3_R_inverse_vector_gyro_bias_alpha6_lower"]
            k = m["finite_angle_vector_information_retention_lower"]
            c = m["initial_vector_Joseph_attenuation_lower"]
            a_sector = m["finite_angle_R_inverse_directional_alpha6_lower"]
            a_j = m["initial_finite_angle_Joseph_directional_alpha6_lower"]
            for x in (a0, k, c, a_sector, a_j):
                self.assertTrue(math.isfinite(x))
                self.assertGreater(x, 0.0)
            self.assertLessEqual(k, 1.0)
            self.assertLessEqual(c, 1.00000000000001)
            self.assertLessEqual(a_sector, a0)
            self.assertLessEqual(a_j, a_sector)

    def test_validation_rejects_premature_scalarization_or_promotion(self):
        d = P.build()
        d["per_packet_full_state_scalarization_used"] = True
        d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"] = True
        failures = P.validate(d)
        self.assertTrue(any("per_packet_full_state_scalarization_used" in x for x in failures))
        self.assertTrue(any("P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED" in x for x in failures))


if __name__ == "__main__":
    unittest.main()
