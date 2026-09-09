"""Algebra of the per-entry-block every-prefix retention producer.

The producer itself consumes the attached physical capture; these tests pin the
parts that do not, so an indexing or bookkeeping regression is caught without
regenerating the word.
"""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np

PRODUCER = Path(__file__).resolve().parents[1] / "kalman_ou_iii"
sys.path.insert(0, str(PRODUCER))

import ou3_p4_entry_block_retention as R


def word(blocks):
    """One synthetic sample whose prefixes are the given (Phi, r) pairs."""
    return [{"prefixes": [(np.asarray(a, dtype=float), np.asarray(f, dtype=float), {})
                          for a, f in blocks]}]


class BlockRetentionAlgebraTests(unittest.TestCase):
    def setUp(self):
        self.radii = np.array([2*np.tan(np.pi/12), 0.01, 5.0, 20.0, 300.0, 2.941995, 0.4])

    def test_entry_and_group_orders_are_the_declared_ones(self):
        self.assertEqual(R.GROUPS, ("attitude", "gyro_bias", "velocity", "position",
                                    "integral_displacement", "latent_acceleration"))
        self.assertEqual(R.ENTRIES, R.GROUPS + ("accelerometer_bias",))
        self.assertEqual(len(R.ENTRIES), 7)
        self.assertEqual(R.DECLARED_CHART_CAYLEY_NORM_UPPER, 1.0)

    def test_single_ball_reaches_dominate_the_subadditive_total(self):
        rng = np.random.default_rng(20260909)
        blocks = [(rng.normal(size=(21, 21))/8 + np.eye(21), rng.normal(size=21)/50)
                  for _ in range(4)]
        out = R.block_retention(word(blocks), self.radii)
        self.assertEqual(out["completed_prefix_count"], 4)
        for g in range(len(R.GROUPS)):
            alone = out["single_ball_reach"][g].sum() + out["template_reach"][g]
            # Each maximum is taken over prefixes independently, so the sum of
            # the parts can only exceed the maximum of their sum.
            self.assertGreaterEqual(alone + 1e-12, out["subadditive_total"][g])
            self.assertGreaterEqual(out["subadditive_total"][g] + 1e-12,
                                    out["subadditive_total_without_independent_integral_ball"][g])

    def test_one_prefix_makes_the_decomposition_exact(self):
        rng = np.random.default_rng(7)
        phi = rng.normal(size=(21, 21))
        r = rng.normal(size=21)
        out = R.block_retention(word([(phi, r)]), self.radii)
        for g in range(len(R.GROUPS)):
            parts = out["single_ball_reach"][g].sum() + out["template_reach"][g]
            self.assertAlmostEqual(parts/out["subadditive_total"][g], 1.0, places=12)

    def test_identity_word_reaches_exactly_its_own_ball(self):
        out = R.block_retention(word([(np.eye(21), np.zeros(21))]), self.radii)
        alone = out["single_ball_reach"]
        for g in range(len(R.GROUPS)):
            self.assertAlmostEqual(alone[g][g], 1.0, places=12)
            for j in range(len(R.ENTRIES)):
                if j != g:
                    self.assertAlmostEqual(alone[g][j], 0.0, places=12)

    def test_summary_reports_the_chart_and_never_promotes(self):
        out = R.block_retention(word([(np.eye(21), np.zeros(21))]), self.radii)
        s = R.summarize(out, self.radii)
        self.assertTrue(s["accelerometer_bias_ball_limits_no_coordinate"])
        self.assertAlmostEqual(s["attitude_cayley_norm_reached_declared_box"],
                               float(self.radii[0]), places=12)
        self.assertTrue(s["declared_box_stays_inside_declared_chart"])
        self.assertFalse(s["entry_radii_reduced"])
        self.assertTrue(s["working_domain_is_larger_than_entry_set"])
        self.assertTrue(s["frozen_point_binary64_diagnostic_not_outward_certificate"])
        for k in ("P4_MOTION_PASS", "P4_PASS", "P5_MAY_START"):
            self.assertFalse(s[k])

    def test_chart_violation_is_reported_not_hidden(self):
        big = np.eye(21)
        big[0:3, 12:15] = np.eye(3)  # integral-displacement ball into attitude
        out = R.block_retention(word([(big, np.zeros(21))]), self.radii)
        s = R.summarize(out, self.radii)
        self.assertFalse(s["declared_box_stays_inside_declared_chart"])
        self.assertTrue(s["without_integral_ball_stays_inside_declared_chart"])
        self.assertEqual(s["worst_single_ball_per_coordinate"]["attitude"], "integral_displacement")

    def test_minimal_working_radius_never_shrinks_the_entry_set(self):
        rng = np.random.default_rng(11)
        blocks = [(rng.normal(size=(21, 21))/20 + np.eye(21), np.zeros(21)) for _ in range(3)]
        s = R.summarize(R.block_retention(word(blocks), self.radii), self.radii)
        for g in R.GROUPS:
            self.assertGreaterEqual(s["minimal_working_radius_inflation_declared_box"][g], 1.0)
            self.assertGreaterEqual(s["minimal_working_radius_inflation_without_integral_ball"][g], 1.0)


if __name__ == "__main__":
    unittest.main()
