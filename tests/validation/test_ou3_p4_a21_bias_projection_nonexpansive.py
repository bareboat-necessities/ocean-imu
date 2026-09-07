from __future__ import annotations

import math
import unittest

import ou3_p4_a21_bias_projection_nonexpansive as PROJ


class TestA21BiasProjectionNonexpansive(unittest.TestCase):
    def test_interior_point_is_unchanged(self):
        x = [0.1, -0.2, 0.05]
        self.assertEqual(PROJ.project_ball(x, 0.4), x)

    def test_exterior_point_projects_to_radius(self):
        y = PROJ.project_ball([1.0, 2.0, -2.0], 0.4)
        self.assertAlmostEqual(math.sqrt(sum(v * v for v in y)), 0.4, places=14)

    def test_distance_to_any_interior_truth_cannot_increase(self):
        truths = ([0.0, 0.0, 0.0], [0.35, 0.0, 0.0], [0.1, -0.2, 0.2])
        estimates = ([1.2, -0.4, 0.7], [-1.0, 0.2, 0.1], [0.1, 0.1, 0.1])
        for truth in truths:
            self.assertLessEqual(math.sqrt(sum(v * v for v in truth)), 0.35 + 1e-15)
            for estimate in estimates:
                self.assertGreaterEqual(PROJ.distance_nonincrease(estimate, truth, 0.4), -1e-15)

    def test_truth_outside_projection_ball_is_rejected(self):
        with self.assertRaises(ValueError):
            PROJ.distance_nonincrease([0.0, 0.0, 0.0], [0.41, 0.0, 0.0], 0.4)

    def test_repository_contract(self):
        d = PROJ.build()
        self.assertEqual(PROJ.validate(d), [])
        self.assertTrue(d["projection_bias_error_energy_nonincrease"])
        self.assertTrue(d["projection_can_be_omitted_from_adverse_P4_bias_budget"])
        self.assertFalse(d["projection_requires_correction_radius"])
        self.assertFalse(d["projection_requires_inactive_branch_assumption"])
        self.assertFalse(d["P4_promoted_here"])


if __name__ == "__main__":
    unittest.main()
