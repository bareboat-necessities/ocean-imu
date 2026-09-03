from pathlib import Path
import sys
import unittest
from unittest import mock

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))

import ou3_p4_signed_vector_directional_coefficients as D


class Ou3P4SignedVectorDirectionalCoefficientsTests(unittest.TestCase):
    def _inputs(self, acc_ret=0.8, mag_ret=0.9):
        rows = []
        for i in range(3200):
            rows.append({
                "mode": "H" if i % 2 == 0 else "A",
                "source_node": (i // 4) % 800,
                "phase_envelope": "stage_boundary_0" if (i // 2) % 2 == 0 else "positive_1_25",
                "accelerometer_innovation_retention_R_over_S_lower": acc_ret,
                "magnetometer_innovation_retention_R_over_S_lower": mag_ret,
            })
        audit = {"expected_source_phase_mode_classes": 3200, "rows": rows}
        cayley = {
            "chart_sigma_min_lower": 0.8,
            "exact_vector_information_retention_factor_lower": 0.64,
            "outer_angle_rad": 0.8,
        }
        remainder = {"acc_eta_force_rotation_quadratic_coefficient_upper": 0.1}
        mag = {
            "effective_tangent_coordinate_gain_lower": 0.75,
            "effective_vs_linear_tangent_defect_ratio_upper": 0.2,
        }
        return audit, cayley, remainder, mag

    def test_coefficients_are_against_linear_tangent_forms(self):
        audit, cayley, remainder, mag = self._inputs()
        with mock.patch.object(D.AUDIT, "validate", return_value=[]), \
             mock.patch.object(D.CAYLEY, "validate", return_value=[]), \
             mock.patch.object(D.REMAINDER, "validate", return_value=[]), \
             mock.patch.object(D.MAG, "validate", return_value=[]):
            d = D.evaluate(audit, cayley, remainder, mag)
        self.assertEqual([], D.validate(d))
        # Accelerometer: 0.8*0.8 - 0.1 = 0.54, rounded downward.
        self.assertLessEqual(d["worst_by_mode"]["H"]["accelerometer_tangent_signed_coefficient_lower"], 0.54)
        self.assertGreater(d["worst_by_mode"]["H"]["accelerometer_tangent_signed_coefficient_lower"], 0.539999999999)
        # Magnetometer: 0.9*(0.75^2) - 0.2^2 = 0.46625, outward.
        self.assertLessEqual(d["worst_by_mode"]["A"]["magnetometer_tangent_signed_coefficient_lower"], 0.46625)
        self.assertGreater(d["worst_by_mode"]["A"]["magnetometer_tangent_signed_coefficient_lower"], 0.466249999999)
        self.assertTrue(d["local_vector_tangent_signed_coefficients_positive_everywhere"])
        self.assertFalse(d["word_level_directional_accumulation_required"])

    def test_negative_local_coefficient_requires_word_accumulation(self):
        audit, cayley, remainder, mag = self._inputs(acc_ret=0.01, mag_ret=0.01)
        with mock.patch.object(D.AUDIT, "validate", return_value=[]), \
             mock.patch.object(D.CAYLEY, "validate", return_value=[]), \
             mock.patch.object(D.REMAINDER, "validate", return_value=[]), \
             mock.patch.object(D.MAG, "validate", return_value=[]):
            d = D.evaluate(audit, cayley, remainder, mag)
        self.assertEqual([], D.validate(d))
        self.assertFalse(d["local_vector_tangent_signed_coefficients_positive_everywhere"])
        self.assertTrue(d["word_level_directional_accumulation_required"])
        self.assertFalse(d["instantaneous_full_state_scalarization_attempted"])
        self.assertFalse(d["P4_USABLE_CERTIFICATE_PROMOTED"])


if __name__ == "__main__":
    unittest.main()