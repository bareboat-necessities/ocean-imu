"""Sharp AW covariance ceiling under the isotropic stationary sync."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tools.stability.ou3_theorem.aw_covariance_ceiling import (
    SIGMA2_MAX, anisotropic_sync_counterexample, ceiling, certificate, diagonal,
    isotropic_sync_witness, prediction_upper, storage_radius)
from tools.stability.ou3_theorem.lin_path_certificate import small_x_source_defect
from tools.stability.ou3_theorem.matrix_certificates import add, matmul, transpose
from tools.stability.ou3_theorem.world_frame import quaternion_rotation


class AwCovarianceCeilingTest(unittest.TestCase):
    def test_committed_certificate_matches_exact_reproduction(self):
        path = ROOT/'reports/results/ou3_stability/aw-covariance-ceiling-certificate.json'
        record = json.loads(path.read_text())
        self.assertEqual(record, certificate())
        self.assertFalse(record['propagated_to_nuisance_upper_comparison'])
        self.assertFalse(record['uniform_AW_tracking_bound'])
        self.assertFalse(record['theorem_closed'])

    def test_prediction_keeps_the_ceiling_and_contracts_the_excess(self):
        eps = small_x_source_defect()[0]
        top = ceiling(eps)
        self.assertEqual(top, (1+eps)*16)
        for phi2 in (F(0), F(1, 7), F(999, 1000), F(1)):
            for m in (F(0), F(121, 25), top):
                self.assertLessEqual(prediction_upper(m, phi2, SIGMA2_MAX, eps), top)
            # Excess over any level s>=(1+eps)sigma^2 contracts by phi^2.
            for sigma2, level, m in ((F(1, 25), F(1, 20), F(3)), (F(9), (1+eps)*9, F(12))):
                excess = max(F(0), prediction_upper(m, phi2, sigma2, eps)-level)
                self.assertLessEqual(excess, phi2*max(F(0), m-level))
        with self.assertRaises(ValueError):
            prediction_upper(F(1), F(2), F(1), eps)

    def test_isotropic_sync_is_the_spectral_max(self):
        q = quaternion_rotation((3, -1, 4, 2))
        for beta, sigma2 in (([F(1), F(2), F(9)], F(4)), ([F(5), F(6), F(7)], F(1)),
                             ([F(0), F(0), F(0)], F(16))):
            p = matmul(q, matmul(diagonal(beta), transpose(q)))
            delta = matmul(q, matmul(diagonal([max(F(0), sigma2-b) for b in beta]), transpose(q)))
            expected = matmul(q, matmul(diagonal([max(b, sigma2) for b in beta]), transpose(q)))
            self.assertEqual(add(p, delta), expected)
        self.assertTrue(isotropic_sync_witness()['post_sync_floor_sigma2'])

    def test_anisotropic_sync_needs_the_isotropy_premise(self):
        witness = anisotropic_sync_counterexample()
        self.assertTrue(witness['isotropy_required'])
        self.assertGreater(F(witness['after_x_diagonal']), F(witness['lambda_max_operands_upper']))

    def test_storage_radius_is_a_valid_lower_bound(self):
        eps = small_x_source_defect()[0]
        threshold = F(112383, 100000)
        for sigma2 in (F(16), F(1)):
            radius = storage_radius(threshold, eps, sigma2)
            self.assertLessEqual(radius**2*(1+eps)*sigma2, threshold**2)
        record = certificate()
        self.assertGreater(F(record['corollary_A_storage_radius_lower_at_sigma_clamp']), F(28, 100))
        self.assertGreater(F(record['corollary_A_storage_radius_lower_at_unit_sigma']), 1)
        self.assertGreater(F(record['standard_deviation_improvement_factor_lower']), 38)
        with self.assertRaises(ValueError):
            storage_radius(threshold, eps, F(2))

    def test_shipping_source_premises(self):
        kalman = (ROOT/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h').read_text()
        fusion = (ROOT/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h').read_text()
        # Pending sync adds Pi_+(target-P_aw); the stationary std is diagonal.
        self.assertIn('for (int i = 0; i < 3; ++i) evals(i) = std::max(T(0), evals(i));', kalman)
        self.assertIn('Pext.template block<3,3>(OFF_AW, OFF_AW) += Delta;', kalman)
        self.assertIn('Sigma_aw_stat = s.array().square().matrix().asDiagonal();\n'
                      '        aw_process_correlated_ = false;', kalman)
        # Default profile: equal horizontal/vertical std, 4 m/s^2 clamp,
        # non-congruent periodic sync.
        self.assertIn('const Eigen::Vector3f aw_std(sH, sH, sZ);', fusion)
        self.assertIn('float S_factor_      = 1.0f;', fusion)
        self.assertIn('constexpr float MAX_SIGMA_A = 4.0f;', fusion)
        self.assertIn('bool congruent_aw_cov_sync_ = false;', fusion)


if __name__ == '__main__':
    unittest.main()
