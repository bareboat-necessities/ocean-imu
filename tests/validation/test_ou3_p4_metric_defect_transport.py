from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_metric_defect_transport as TRANSPORT
import ou3_p4_nonlinear_word_certificate as P4


class Ou3P4MetricDefectTransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = P4.build()

    def test_certificate_keeps_the_smaller_of_the_two_valid_defect_gains(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            iso = m["transported_word_defect_B_isotropic_upper"]
            met = m["transported_word_defect_B_metric_consistent_upper"]
            chosen = m["transported_word_defect_B_upper"]
            self.assertGreater(iso, 0.0)
            self.assertGreater(met, 0.0)
            self.assertEqual(chosen, min(iso, met))
            self.assertEqual(
                m["transported_word_defect_route"],
                "METRIC_CONSISTENT_STRUCTURED_DEFECT_TRANSPORT"
                if met <= iso
                else "ISOTROPIC_EUCLIDEAN_DEFECT_ENVELOPE",
            )

    def test_structured_transport_is_the_binding_route_and_widens_the_level(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            self.assertEqual(
                m["transported_word_defect_route"],
                "METRIC_CONSISTENT_STRUCTURED_DEFECT_TRANSPORT",
            )
            # The structured route must stay far below the isotropic envelope it
            # replaces; a regression that silently reverts it is caught here.
            self.assertLess(
                m["transported_word_defect_B_metric_consistent_upper"],
                1.0e-20 * m["transported_word_defect_B_isotropic_upper"],
            )
            # The isotropic route certified 3.3e-141; the structured one must stay
            # many decades above that without ever reaching a practical level.
            self.assertGreater(m["certified_level_W"], 1.0e-100)
            self.assertLess(m["certified_level_W"], 1.0e-40)

    def test_defect_inputs_exclude_the_translation_block(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            # The nonlinear defects are quadratic in attitude, gyro bias and a_w
            # only.  The translation diagonals carry the whole covariance spread
            # and must not enter the chart scale.
            self.assertGreaterEqual(
                m["Sigma_defect_input_block_upper"], m["Sigma_attitude_block_upper"]
            )
            self.assertLess(m["Sigma_defect_input_block_upper"], m["Sigma_lambda_max_upper"])
            t = m["metric_consistent_defect_transport"]
            self.assertLess(t["attitude_chart_scale"], 1.0)
            self.assertAlmostEqual(
                t["attitude_chart_scale"],
                math.sqrt(m["Sigma_attitude_block_upper"] / m["Sigma_lambda_max_upper"]),
                delta=1.0e-12,
            )

    def test_gain_transport_uses_the_half_bound_not_the_isotropic_gain_norm(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            t = m["metric_consistent_defect_transport"]
            rmin = m["measurement_bounds"]["all_correction_R_lambda_min_lower"]
            s = m["Sigma_lambda_max_upper"]
            self.assertAlmostEqual(
                t["gain_metric_transport_upper"],
                math.sqrt(s) / (2.0 * math.sqrt(rmin)),
                delta=1.0e-6 * t["gain_metric_transport_upper"],
            )
            # sqrt(s)/(2 sqrt(Rmin)) must beat the retired sqrt(m_+)*||K|| chain.
            self.assertLess(
                t["gain_metric_transport_upper"],
                math.sqrt(m["metric_lambda_max_upper"]) * m["full_gain_norm_upper"],
            )

    def test_module_validates_and_closes_its_own_prefix_bootstrap(self):
        for mode in ("H", "A"):
            t = self.d["modes"][mode]["metric_consistent_defect_transport"]
            self.assertEqual(TRANSPORT.validate(t), [])
            self.assertLessEqual(t["prefix_metric_norm_upper"], t["design_metric_radius"])
            self.assertLess(t["prefix_attitude_cayley_norm_upper"], 1.0)
            self.assertGreater(t["certified_level_W"], 0.0)

    def test_module_fails_closed_on_a_non_positive_source_bound(self):
        base = dict(self.d["modes"]["H"]["metric_consistent_defect_transport"])
        inputs = {
            "metric_scale": self.d["modes"]["H"]["Sigma_lambda_max_upper"],
            "word_endpoint_delta_lower": 0.0,
            "correction_R_lambda_min_lower": 1.0,
            "Sigma_attitude_upper": 1.0,
            "Sigma_gyro_bias_upper": 1.0,
            "Sigma_defect_input_upper": 1.0,
            "rho_attitude_scaled_lower": 1.0,
            "Q_theta_diagonal_lower": 1.0,
            "H_attitude_norm_upper": 1.0,
            "vector_residual_quadratic_constant_upper": 1.0,
            "prediction_increment_gain_upper": 1.0,
            "state_operation_count_upper": 1.0,
        }
        with self.assertRaises(RuntimeError):
            TRANSPORT.build(inputs)
        self.assertGreater(base["certified_level_W"], 0.0)


if __name__ == "__main__":
    unittest.main()
