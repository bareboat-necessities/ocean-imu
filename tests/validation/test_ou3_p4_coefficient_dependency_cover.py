#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_brmm_tuner_scheduler_step as TUNER
import ou3_p4_coefficient_dependency_cover as COVER


class CoefficientDependencyCoverTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # A coarse cover is enough for every structural assertion; the joint
        # scale only tightens with more cells and stays an upper bound.
        cls.d = COVER.build(cells=48)

    def test_validates_and_promotes_nothing(self):
        self.assertEqual(COVER.validate(self.d), [])
        self.assertFalse(self.d["P4_promoted_here"])
        self.assertFalse(self.d["SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED"])
        self.assertFalse(self.d["independent_tau_sigma_TS_RS_rectangle_used"])
        self.assertFalse(self.d["reachable_f_sigma_set_certified_here"])
        self.assertFalse(self.d["composed_with_Riccati_Joseph_reset_graph_here"])

    def test_pseudo_cadence_is_the_clamped_affine_image_of_tau(self):
        # build() raises if the deployed transition disagrees; re-assert the
        # probe rows here so a silently weakened check is visible.
        c = TUNER.constants()
        for row in self.d["pinning_probe"]:
            tau = row["tau_s"]
            expect = TUNER.clamp_interval(TUNER.I(c.pseudo_ratio) * TUNER.I(tau),
                                          c.pseudo_min_s, c.pseudo_max_s)
            self.assertEqual([expect.lo, expect.hi], row["T_S_s"])

    def test_applied_anisotropic_R_S_is_a_single_scalar_ray(self):
        ray = self.d["applied_R_S_ray"]
        self.assertTrue(ray["one_scalar_generates_three_axes"])
        self.assertEqual(ray["x_factor"], ray["y_factor"])
        self.assertEqual(ray["z_factor"], 1.0)
        self.assertEqual(ray["probe_applied_z"], ray["probe_base"])

    def test_reachable_target_set_is_lower_dimensional_than_the_rectangle(self):
        interior = self.d["interior_invertibility"]
        self.assertTrue(interior["T_S_pins_tau_on_the_open_branch"])
        self.assertLess(interior["reachable_target_surface_dimension"],
                        interior["nominal_independent_rectangle_dimension"])
        self.assertTrue(interior["T_S_lower_rail_unreachable_for_targets"])

    def test_coefficient_box_is_forward_invariant_without_a_reachable_set_argument(self):
        box = self.d["forward_invariant_coefficient_box"]
        self.assertFalse(box["needs_reachable_set_argument"])
        tau_lo, tau_hi = box["tau_s"]
        self.assertGreater(tau_lo, 0.0)
        self.assertLess(tau_lo, tau_hi)
        ts_lo, ts_hi = box["T_S_s"]
        c = TUNER.constants()
        self.assertGreater(ts_lo, c.pseudo_min_s)
        self.assertLessEqual(ts_hi, c.pseudo_max_s)
        self.assertFalse(box["T_S_lower_rail_reachable"])

    def test_joint_cell_scale_is_sharper_than_the_independent_corner(self):
        cover = self.d["coefficient_cover"]
        self.assertLess(cover["joint_S_zero_residual_scale_upper"],
                        cover["independent_corner_residual_scale"])
        self.assertGreater(cover["independent_corner_over_approximation_factor"], 1.0)
        self.assertGreater(cover["independent_corner_energy_over_approximation_factor"],
                           cover["independent_corner_over_approximation_factor"])

    def test_every_adaptation_chain_is_a_strict_contraction_toward_its_target(self):
        rates = self.d["adaptation_rate_limits"]
        self.assertGreater(rates["steps_per_staged_commit_interval"], 1)
        chains = [k for k, v in rates.items() if isinstance(v, dict)]
        self.assertGreaterEqual(len(chains), 5)
        for name in chains:
            row = rates[name]
            with self.subTest(chain=name):
                self.assertTrue(row["strict_contraction_toward_target"])
                self.assertLess(row["per_step_alpha_upper"], 1.0)
                self.assertGreater(row["per_step_alpha_lower"], 0.0)
                self.assertLessEqual(row["target_gap_closed_per_commit_interval_upper"], 1.0)
        # The wave-period moments are the slowest chain, which is what makes the
        # coefficient trajectory rate limited rather than free.
        self.assertLess(rates["wave_period_moment_EWMA"]["per_step_alpha_upper"],
                        rates["candidate_tau_sigma_EMA"]["per_step_alpha_upper"])

    def test_missing_BRMM_S_m_primitive_is_reported_with_its_budget(self):
        cons = self.d["S_zero_correction_consequence"]
        self.assertIsNone(cons["contract_S_m_value"])
        self.assertFalse(cons["contract_S_m_instantiated"])
        self.assertTrue(cons["S_m_is_the_missing_source_primitive"])
        # The budget divides by the transverse attitude variance. That variance
        # was once taken as the prior-independent r/|f|^2, which is retracted;
        # under the conditional cap the magnitude-only route has no headroom, so
        # the honest budget is absent rather than 7.35 m*s. Either outcome is
        # admissible here, but the two must agree with each other.
        self.assertFalse(cons["transverse_attitude_variance_is_prior_independent"])
        self.assertGreater(cons["conditional_transverse_attitude_variance_rad2"],
                           cons["retracted_prior_independent_transverse_variance_rad2"])
        budget = cons["admissible_BRMM_S_m_upper_m_s"]
        self.assertEqual(budget is None, cons["no_admissible_S_m_on_the_magnitude_only_route"])
        self.assertEqual(budget is not None, cons["magnitude_only_route_has_headroom"])
        if budget is not None:
            self.assertGreater(budget, 0.0)
        else:
            self.assertLessEqual(cons["headroom_for_true_integral_displacement"], 0.0)
            self.assertIn("magnitude-only route has no headroom",
                          cons["fallback_if_qualified_S_m_exceeds_budget"])

    def test_tiling_is_a_gapless_closed_cover(self):
        tiles = COVER._geometric_tiling(0.5, 4.0, 7)
        self.assertEqual(len(tiles), 7)
        self.assertEqual(tiles[0][0], 0.5)
        self.assertEqual(tiles[-1][1], 4.0)
        for (_, b), (a2, _) in zip(tiles, tiles[1:]):
            self.assertEqual(a2, b)
        with self.assertRaises(ValueError):
            COVER._geometric_tiling(0.0, 1.0, 4)


if __name__ == "__main__":
    unittest.main()
