import copy
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / 'tools/stability')]

import numpy as np

from tools.stability.ou3_alt_contraction import common_storage_master as MASTER
from tools.stability.ou3_alt_contraction import finite_word_rho_diagnostic as RHO
from tools.stability.ou3_alt_contraction import proof_plan as PLAN

SHORT = 3


class WordCompositionTests(unittest.TestCase):
    def test_every_literal_event_kind_composes_in_both_modes(self):
        gauged = next(p for p in RHO.SOURCE_PHASES if p.gauged)
        for mode in RHO.MODES:
            with self.subTest(mode=mode):
                word = RHO.compose_word(gauged, RHO.GUARD_PROBE_SAMPLES, mode)
                self.assertEqual(word['A_joint24'].shape, (RHO.JOINT_DIM, RHO.JOINT_DIM))
                for kind in ('prediction', 'S_zero', 'accelerometer', 'magnetometer'):
                    self.assertIn(kind, word['literal_event_kinds'])
                self.assertFalse(word['replay_used'])

    def test_ungauged_word_contains_no_magnetic_event(self):
        ungauged = next(p for p in RHO.SOURCE_PHASES if not p.gauged)
        word = RHO.compose_word(ungauged, SHORT)
        self.assertNotIn('magnetometer', word['literal_event_kinds'])

    def test_H18_neutral_rows_carry_no_motion_column_so_the_projection_is_exact(self):
        for phase in RHO.SOURCE_PHASES:
            word = RHO.compose_word(phase, SHORT, 'H')
            certificate = RHO.block_triangular_certificate(word['A_joint24'])
            self.assertTrue(certificate['neutral_rows_have_no_motion_column'], phase.name)
            self.assertEqual(certificate['worst_neutral_to_motion_entry'], 0.0)

    def test_A21_releases_the_bias_so_the_motion_block_is_not_a_floor(self):
        """The released bias is corrected from motion, so ker(C) is not invariant."""
        word = RHO.compose_word(RHO.SOURCE_PHASES[0], SHORT, 'A')
        certificate = RHO.block_triangular_certificate(word['A_joint24'])
        self.assertFalse(certificate['projected_restriction_equals_motion_block'])
        self.assertGreater(certificate['worst_neutral_to_motion_entry'], 0.0)
        rho = RHO.word_rho(word['A_joint24'],
                           certificate['projected_restriction_equals_motion_block'])
        self.assertFalse(rho['motion_block_ratio_is_a_floor'])
        self.assertIsNone(rho['rho_floor'])

    def test_projected_finsler_restriction_equals_the_motion_block(self):
        """common_storage_master's reduction must agree with the motion block."""
        word = RHO.compose_word(RHO.SOURCE_PHASES[0], SHORT, 'H')
        A = word['A_joint24']
        rng = np.random.default_rng(20260917)
        root = rng.standard_normal((RHO.JOINT_DIM, RHO.JOINT_DIM))
        M = root @ root.T + RHO.JOINT_DIM * np.eye(RHO.JOINT_DIM)
        rho = 0.9
        restriction = MASTER.unsupplied_storage_form(A, M, rho)
        motion = A[:RHO.MOTION_DIM, :RHO.MOTION_DIM]
        M_motion = M[:RHO.MOTION_DIM, :RHO.MOTION_DIM]
        expected = motion.T @ M_motion @ motion - rho * M_motion
        self.assertTrue(np.allclose(restriction, (expected + expected.T) / 2, atol=1e-9))


class UngaugedObstructionTests(unittest.TestCase):
    def test_ungauged_yaw_pair_is_an_exactly_invariant_unipotent_block(self):
        """Checked at full joint24 width, so it holds in H18 and A21 alike."""
        phase = RHO.SOURCE_PHASES[0]
        for mode in RHO.MODES:
            with self.subTest(mode=mode):
                word = RHO.compose_word(phase, RHO.GUARD_PROBE_SAMPLES, mode)
                certificate = RHO.ungauged_yaw_certificate(
                    word['A_joint24'], word['word_horizon_s'])
                self.assertTrue(certificate['subspace_is_exactly_invariant'])
                self.assertEqual(certificate['worst_row_leak_outside_block'], 0.0)
                self.assertEqual(certificate['worst_column_leak_outside_block'], 0.0)
                self.assertTrue(certificate['unipotent_jordan_block_certified'])
                self.assertEqual(certificate['certified_spectral_radius'], 1.0)
                self.assertTrue(certificate['drift_matches_elapsed_horizon'])

    def test_bias_to_yaw_drift_equals_the_elapsed_word_horizon(self):
        phase = RHO.SOURCE_PHASES[0]
        for samples in (1, 2, 4):
            word = RHO.compose_word(phase, samples, 'H')
            certificate = RHO.ungauged_yaw_certificate(
                word['A_joint24'], word['word_horizon_s'])
            self.assertAlmostEqual(certificate['observed_bias_to_yaw_drift'],
                                   word['word_horizon_s'], places=6)

    def test_obstruction_is_independent_of_word_length(self):
        for mode in RHO.MODES:
            with self.subTest(mode=mode):
                induction = RHO.per_sample_unipotent_induction(mode=mode)
                self.assertTrue(induction['every_prefix_unipotent'])
                self.assertTrue(induction['both_literal_sample_shapes_exercised'])
                self.assertTrue(induction['length_independent'])

    def test_a_gauged_word_breaks_the_invariant_subspace(self):
        """Magnetic service is exactly what removes the obstruction."""
        gauged = next(p for p in RHO.SOURCE_PHASES if p.gauged)
        word = RHO.compose_word(gauged, RHO.GUARD_PROBE_SAMPLES, 'H')
        certificate = RHO.ungauged_yaw_certificate(word['A_joint24'], word['word_horizon_s'])
        self.assertFalse(certificate['subspace_is_exactly_invariant'])
        self.assertFalse(certificate['unipotent_jordan_block_certified'])


class ReportTests(unittest.TestCase):
    """One shared diagnose; composing a report is expensive and these read it."""

    @classmethod
    def setUpClass(cls):
        cls.built = RHO.diagnose(samples=RHO.GUARD_PROBE_SAMPLES,
                                 phases=(RHO.SOURCE_PHASES[0],), modes=('H', 'A'))

    def setUp(self):
        self.report = copy.deepcopy(self.built)

    def test_report_validates_and_reports_a_theorem_failure(self):
        self.assertEqual(RHO.validate(self.report), [])
        self.assertTrue(self.report['declared_joint24_contraction_falsified'])
        self.assertEqual(self.report['failure_classification'], 'theorem_failure')
        self.assertEqual(self.report['rho_floor_over_legal_words'], 1.0)
        self.assertEqual(self.report['rho_floor_source'], 'exact_unipotent_subspace_algebra')
        self.assertEqual(self.report['distance_to_rho_one'], 0.0)
        self.assertTrue(self.report['metric_independent'])

    def test_limiting_direction_is_the_ungauged_heading_pair(self):
        named = [x['coordinate'] for x in self.report['limiting_state_direction']]
        self.assertIn('theta_z', named)
        self.assertIn('bg_z', named)

    def test_measured_ratios_are_reported_beside_their_drift(self):
        self.assertIn('numerical_drift_proxy', self.report)
        self.assertLess(self.report['numerical_drift_proxy'], 1e-4)
        self.assertTrue(self.report['measured_ratios_are_meaningful_only_outside_drift'])

    def test_worst_ratio_is_reported_for_both_modes(self):
        worst = self.report['worst_motion_block_ratio_by_mode']
        self.assertEqual(set(worst), {'H', 'A'})
        for mode, value in worst.items():
            self.assertIsNotNone(value, mode)

    def test_both_modes_certify_the_obstruction(self):
        self.assertEqual(set(self.report['ungauged_unipotent_words']),
                         {'H:quiet_ungauged', 'A:quiet_ungauged'})

    def test_diagnostic_can_never_promote(self):
        for key in ('storage_search_allowed', 'ALT_STARTUP_PASS', 'ALT_LIVE_PASS',
                    'ALT_END_TO_END_PASS'):
            self.assertFalse(self.report[key])
        promoted = dict(self.report, storage_search_allowed=True)
        with self.assertRaises(RuntimeError):
            PLAN.assert_non_promoting_report(promoted)

    def test_declared_phase_is_non_promoting_and_not_a_promoting_phase(self):
        self.assertIn(RHO.PHASE, PLAN.NON_PROMOTING_PHASES)
        self.assertNotIn(RHO.PHASE, PLAN.PROMOTING_PHASES)
        self.assertNotIn(RHO.EVIDENCE_KIND, PLAN.DIAGNOSTIC_ONLY)

    def test_validate_rejects_a_falsification_without_length_independence(self):
        broken = dict(self.report, per_sample_unipotent_induction=None)
        self.assertIn('falsification was not shown independent of word length',
                      RHO.validate(broken))

    def test_validate_rejects_a_floor_claimed_without_an_invariant_kernel(self):
        broken = copy.deepcopy(self.report)
        row = broken['words']['H:quiet_ungauged']
        row['block_triangular_certificate']['projected_restriction_equals_motion_block'] = False
        self.assertIn('H:quiet_ungauged: H18 kernel is no longer invariant',
                      RHO.validate(broken))

    def test_validate_rejects_a_falsification_resting_on_a_measured_radius(self):
        broken = dict(self.report, rho_floor_source='measured_spectral_radius')
        self.assertIn('falsification rests on a measured radius rather than exact algebra',
                      RHO.validate(broken))


class SeedTests(unittest.TestCase):
    def test_live_entry_covariance_comes_from_the_shipping_seed_producer(self):
        P = RHO.live_entry_covariance()
        self.assertEqual(P.shape, (RHO.MOTION_DIM, RHO.MOTION_DIM))
        self.assertTrue(np.allclose(P, P.T))
        self.assertGreater(float(np.linalg.eigvalsh(P).min()), 0.0)
        # Ungauged yaw is not seeded as a small angle, and the translation block
        # keeps the shipping seed magnitudes.
        self.assertAlmostEqual(P[2, 2], 0.087 ** 2, places=12)
        self.assertAlmostEqual(P[9, 9], 400.0, places=9)
        self.assertAlmostEqual(P[12, 12], 2500.0, places=9)

    def test_A21_seed_adds_the_released_bias_block(self):
        P = RHO.live_entry_covariance('A')
        self.assertEqual(P.shape, (21, 21))
        self.assertGreater(float(np.linalg.eigvalsh(P).min()), 0.0)
        self.assertAlmostEqual(P[RHO.OFF_BA, RHO.OFF_BA], 1.6e-05, places=12)
        self.assertTrue(np.allclose(P[:RHO.MOTION_DIM, RHO.OFF_BA:], 0.0))


if __name__ == '__main__':
    unittest.main()
