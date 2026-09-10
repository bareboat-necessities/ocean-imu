#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_p4_basin_frontier as FRONTIER


class BasinFrontierTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = FRONTIER.build()

    def test_validates_and_promotes_nothing(self):
        self.assertEqual(FRONTIER.validate(self.d), [])
        self.assertFalse(self.d["P4_promoted_here"])
        self.assertFalse(self.d["outward_nonlinear_certificate_closed_here"])
        self.assertFalse(self.d["entry_radii_reduced_for_proof_convenience"])
        self.assertTrue(self.d["basin_is_a_correlated_polytope_not_a_product_box"])
        self.assertTrue(self.d["frozen_map_binary64_diagnostic_not_outward_certificate"])

    def test_declared_box_overshoots_the_chart_in_H18(self):
        h = self.d["modes"]["H18"]
        self.assertFalse(h["declared_box_retains_chart"])
        self.assertGreater(h["declared_box_overshoot_factor"], 1.0)

    def test_max_uniform_scale_saturates_the_half_space(self):
        for mode, row in self.d["modes"].items():
            with self.subTest(mode=mode):
                spend = sum(row["max_uniform_scale"] * a
                            for a in row["half_space_coefficients"].values())
                self.assertLessEqual(spend, row["half_space_budget"] * (1 + 1e-9))
                self.assertGreaterEqual(spend, row["half_space_budget"] * (1 - 1e-6))

    def test_max_volume_point_beats_uniform_on_the_slack_axes(self):
        h = self.d["modes"]["H18"]
        uniform = h["max_uniform_scale"]
        mv = h["max_volume_scale"]
        # A coordinate whose single-ball reach is below average gets a larger
        # radius under the volume optimum than under a uniform inflation.
        self.assertGreater(mv["gyro_bias"], uniform)
        self.assertGreater(mv["latent_acceleration"], uniform)
        # The integral ball dominates the reach, so it is the axis the volume
        # optimum pays for; it still admits a positive radius.
        self.assertGreater(mv["integral_displacement"], 0.0)
        self.assertLess(mv["integral_displacement"], uniform)

    def test_max_volume_radii_are_physically_larger_than_the_certified_uniform_box(self):
        h = self.d["modes"]["H18"]
        radii = self.d["declared_entry_radii"]
        uniform_radii = h["max_uniform_radii"]
        mv_radii = h["max_volume_radii"]
        for axis in ("gyro_bias", "velocity", "latent_acceleration", "accelerometer_bias"):
            with self.subTest(axis=axis):
                self.assertGreater(mv_radii[axis], uniform_radii[axis])
                self.assertLessEqual(uniform_radii[axis], radii[axis])

    def test_integral_axis_uses_the_sharper_retained_joint_maximum(self):
        for mode, row in self.d["modes"].items():
            integral = row["per_axis_max_scale_with_others_declared"]["integral_displacement"]
            with self.subTest(mode=mode):
                self.assertTrue(integral["sharper_bound_used"])
                self.assertGreaterEqual(integral["sharper_max_scale_from_retained_joint_maximum"],
                                        integral["max_scale"])
        self.assertGreater(self.d["max_integral_radius_with_others_declared_m_s"], 0.0)
        self.assertLess(self.d["max_integral_radius_with_others_declared_m_s"],
                        self.d["declared_entry_radii"]["integral_displacement"])

    def test_binding_mode_is_reported_and_is_the_minimum(self):
        binding = self.d["binding_mode_for_uniform_scale"]
        for row in self.d["modes"].values():
            self.assertGreaterEqual(row["max_uniform_scale"],
                                    self.d["modes"][binding]["max_uniform_scale"])
        self.assertEqual(self.d["max_uniform_scale_over_modes"],
                         self.d["modes"][binding]["max_uniform_scale"])

    def test_rebuild_is_demanded_if_the_declared_box_ever_retains_the_chart(self):
        d = {k: v for k, v in self.d.items()}
        modes = {m: dict(row) for m, row in d["modes"].items()}
        modes["H18"]["declared_box_retains_chart"] = True
        d["modes"] = modes
        self.assertIn("H18 declared box now retains the chart; rebuild the frontier",
                      FRONTIER.validate(d))


if __name__ == "__main__":
    unittest.main()
