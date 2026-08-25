import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_cayley_eta_geometry as E


class Ou3P5CayleyEtaGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = E.build()

    def test_exact_eta_geometry_passes_on_widened_first_S_chart(self):
        d = self.d
        self.assertEqual(E.validate(d), [])
        self.assertEqual(d["P5_CAYLEY_ETA_GEOMETRY_CERTIFICATE"], "PASS")
        self.assertGreater(d["widened_cayley_norm_upper"], 1.0)
        self.assertGreater(d["widened_exact_vector_residual_factor_lower"], 0.0)
        self.assertGreater(d["widened_eta_to_rotational_residual_information_ratio_upper"], 0.0)
        self.assertLess(d["widened_latent_vector_rotation_gain_upper"], 2.0)

    def test_scalar_identities_match_direct_rotation_geometry(self):
        for q in (1e-6, 0.1, 0.6, 1.0, 2.0, 4.0, 8.0):
            theta = 2.0 * math.atan(q / 2.0)
            # Unit vector perpendicular to the Cayley axis.
            y2 = (math.cos(theta) - 1.0) ** 2 + math.sin(theta) ** 2
            h_minus_y_2 = (math.cos(theta) - 1.0) ** 2 + (q - math.sin(theta)) ** 2
            self.assertAlmostEqual(h_minus_y_2 / y2, q*q/4.0, places=10)
            self.assertLessEqual(4.0 / (4.0 + q*q), 1.0)

    def test_annular_partition_is_outward_cover_and_monotone_safe(self):
        cells = self.d["annular_subdivision_cells"]
        self.assertEqual(len(cells), 64)
        self.assertEqual(cells[0]["q_interval"][0], 0.0)
        self.assertGreaterEqual(cells[-1]["q_interval"][1], self.d["widened_cayley_norm_upper"])
        previous_hi = 0.0
        previous_eta = 0.0
        for row in cells:
            lo, hi = row["q_interval"]
            self.assertLessEqual(lo, previous_hi)
            self.assertGreaterEqual(hi, previous_hi)
            self.assertGreater(row["exact_vector_residual_factor_lower"], 0.0)
            self.assertGreaterEqual(row["exact_eta_to_rotational_residual_information_ratio_upper"], previous_eta)
            previous_hi = hi
            previous_eta = row["exact_eta_to_rotational_residual_information_ratio_upper"]

    def test_eta_subdivision_does_not_promote_word_or_change_filter(self):
        d = self.d
        self.assertFalse(d["complete_word_covariance_reset_transport_closed_here"])
        self.assertFalse(d["global_packet_count_times_Lipschitz_defect_used"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertTrue(d["full_S_to_attitude_gain_retained"])


if __name__ == "__main__":
    unittest.main()
