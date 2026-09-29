"""Nominal-mean attitude columns and the literal AW correction loop."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tools.stability.ou3_theorem.aw_tracking import (  # noqa: E402
    G, attitude_gram, aw_increment_cap, aw_loop_identity, certificate, corollary_a_star,
    nominal_mean_representation, nonlocal_tracking_example, physical_transfer,
    world_innovation_factorization)
from tools.stability.ou3_theorem.aw_tracking_source_diagnostic import (  # noqa: E402
    CORR_A_PHYSICAL_THRESHOLD, CORR_A_STAR_THRESHOLD, PROFILES, TAP_ANCHOR, admissibility,
    verify_diagnostic)
from tools.stability.ou3_theorem.matrix_certificates import add, identity, is_psd  # noqa: E402

RESULTS = ROOT/'reports/results/ou3_stability'


class CorollaryAStarTest(unittest.TestCase):
    def test_committed_certificate_matches_exact_reproduction(self):
        record = json.loads((RESULTS/'aw-tracking-certificate.json').read_text())
        self.assertEqual(record, certificate())
        self.assertFalse(record['source_uniform_nominal_mean_bound'])
        self.assertFalse(record['uniform_AW_tracking_bound'])
        self.assertFalse(record['theorem_closed'])
        self.assertFalse(record['pointwise_tracking_premise_necessary'])

    def test_threshold_and_floor(self):
        self.assertEqual(corollary_a_star(0, 0)['normalized_gram_floor'], F(1, 50))
        self.assertEqual(corollary_a_star(0, 0)['nominal_transverse_mean_threshold'], G/5)
        self.assertFalse(corollary_a_star(G/5, G/5)['positive'])
        self.assertTrue(corollary_a_star(G/5-F(1, 10**6), G/5)['positive'])
        floors = [corollary_a_star(m, m)['normalized_gram_floor'] for m in (0, F(1, 2), 1, F(3, 2))]
        self.assertEqual(floors, sorted(floors, reverse=True))
        with self.assertRaises(ValueError):
            corollary_a_star(1, F(1, 2))

    def test_gram_dominates_the_mean_floor_on_supplied_rows(self):
        b = [F(7, 25), F(0), F(24, 25)]
        rows = [[F(3), F(-1), F(2)], [F(-2), F(4), F(-1)], [F(1, 2), F(-3), F(0)], [F(0), F(0), F(-7, 2)]]
        weights = [F(1, 8), F(3, 8), F(1, 4), F(1, 4)]
        mean = [sum(w*r[c] for w, r in zip(weights, rows)) for c in range(3)]
        cross = [mean[1]*b[2]-mean[2]*b[1], mean[2]*b[0]-mean[0]*b[2], mean[0]*b[1]-mean[1]*b[0]]
        m_perp2 = sum(x*x for x in cross)
        m2 = sum(x*x for x in mean)
        # Rational upper bounds of the square roots keep the floor conservative.
        m_perp, m = F(3, 2), F(3)
        self.assertLessEqual(m_perp2, m_perp**2)
        self.assertLessEqual(m2, m**2)
        floor = corollary_a_star(m_perp, m, F(7, 25))['normalized_gram_floor']
        self.assertGreater(floor, 0)
        self.assertTrue(is_psd(add(attitude_gram(weights, rows, b), identity(3), -floor)))
        self.assertTrue(nonlocal_tracking_example()['exact_gram_dominates_gamma_star'])

    def test_signed_transfer_equals_corollary_a(self):
        self.assertEqual(F(physical_transfer(16)['signed_mean_error_threshold_mps2']), F(112383, 100000))
        self.assertTrue(physical_transfer(64)['equals_pointwise_corollary_A_threshold'])


class LiteralLoopTest(unittest.TestCase):
    def test_nominal_mean_is_a_weighted_signed_correction_sum(self):
        ops = [('row', F(1, 3)), ('pred', F(9, 10)), ('corr', [F(1), F(-1), F(2)]),
               ('row', F(1, 3)), ('corr', [F(-1, 2), F(0), F(1, 3)]), ('pred', F(1, 2)),
               ('row', F(1, 3))]
        rep = nominal_mean_representation(ops, [F(2), F(0), F(-1)])
        self.assertTrue(all(0 <= w <= 1 for w in rep['correction_weights']))
        with self.assertRaises(ValueError):
            nominal_mean_representation([('pred', F(2))], [F(0)]*3)

    def test_exact_identities(self):
        self.assertEqual(aw_loop_identity()['exact_residual'], '0')
        self.assertEqual(world_innovation_factorization()['exact_residual'], '0')
        cap = aw_increment_cap()
        self.assertTrue(cap['prior_aw_dominates_increment'])
        self.assertGreater(F(cap['per_step_NIS_needed_to_follow_jerk_limit_at_floor']), 143)


class CarriedAuditTest(unittest.TestCase):
    def test_committed_audit_is_consistent_and_non_promoting(self):
        record = json.loads((RESULTS/'aw-tracking-source-feasibility.json').read_text())
        self.assertTrue(verify_diagnostic(record))
        self.assertTrue(record['pointwise_premise_refuted_on_admitted_history'])
        self.assertGreater(record['worst_pointwise_ratio'], 6)
        self.assertLess(record['worst_signed_mean_ratio'], 1)
        self.assertLess(record['worst_corollary_A_star_ratio'], 1)
        self.assertLess(record['worst_nominal_transverse_mean_mps2'], 0.4)
        self.assertLess(record['worst_nominal_force_L1'], 1.2)
        self.assertFalse(record['source_uniform_certificate'])
        self.assertTrue(record['G0_floor_below_literal_and_injection_free'])
        self.assertGreater(record['literal_over_injection_free_min'], 0.99)
        for row in record['profiles'].values():
            self.assertLess(row['A21_active_step'], 32000)
            array = row['literal_array']
            self.assertLess(array['A_tilde_minus_I_max'], 0.01)
            self.assertLessEqual(array['G0_floor_service_gap_1s'], array['literal_sigma_min'])
        self.assertEqual(CORR_A_STAR_THRESHOLD, G/5)
        self.assertEqual(CORR_A_PHYSICAL_THRESHOLD, F(112383, 100000))

    def test_every_profile_is_admitted_exactly(self):
        for profile in PROFILES:
            self.assertTrue(admissibility(*profile[1:])['admitted'], profile[0])
        # A plain 8.8 m/s^2 triangle with the onset ramp exceeds the envelope.
        self.assertFalse(admissibility('8.8', '2.75', '0', '1', '0', '1')['admitted'])

    def test_shipping_anchors(self):
        kalman = (ROOT/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h').read_text()
        fusion = (ROOT/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h').read_text()
        self.assertEqual(kalman.count(TAP_ANCHOR), 1)
        self.assertIn('const float sigma_floor = std::max(0.05f, band_noise_floor_sigma_());', fusion)
        self.assertIn('const Eigen::Vector3f aw_std(sH, sH, sZ);', fusion)
        common = (ROOT/'src/kalman_common/SeaStateAdaptationCommon.h').read_text()
        self.assertIn('return !(time_ - last_aw_cov_sync_sec_ <= adapt_every_secs_);', common)


if __name__ == '__main__':
    unittest.main()
